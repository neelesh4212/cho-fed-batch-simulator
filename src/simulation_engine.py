"""Continuous fed-batch balances with a sampled, zero-order-hold feed controller."""

from dataclasses import asdict, dataclass, replace
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
import json
import sqlite3
import uuid

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.special import expit

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODEL_VERSION = "2.0-ode"


@dataclass(frozen=True)
class Parameters:
    duration: float = 14.0  # days
    dt: float = 1.0  # days, controller and export interval; solver steps are smaller
    initial_volume: float = 1.0  # L
    initial_glucose: float = 8.0  # g/L
    initial_viable_cells: float = 1.0  # million cells/mL
    feed_glucose: float = 100.0  # g/L
    feed_rate_on: float = 0.1  # L/day
    glucose_target: float = 10.0  # g/L
    max_feed_volume: float = 1.5  # L, feeding cutoff
    reactor_capacity: float = 2.0  # L, physical capacity
    mu_max: float = 0.4  # 1/day
    Ks: float = 2.0  # g/L
    Kl: float = 4.0  # g/L
    Ka: float = 2.0  # g/L
    qG: float = 600.0  # pg/(cell.day), maximum glucose uptake
    qp: float = 30.0  # pg/(cell.day)
    qA: float = 40.0  # pg/(cell.day)
    qL_production: float = 250.0  # pg/(cell.day)
    qL_consumption: float = 50.0  # pg/(cell.day)
    mu_switch: float = 0.15  # 1/day
    switch_steepness: float = 20.0  # days, makes the sigmoid exponent dimensionless
    glucose_uptake_K: float = 0.1  # g/L, assumed uptake saturation constant
    lactate_uptake_K: float = 0.1  # g/L, assumed uptake saturation constant
    kd_min: float = 0.02  # 1/day
    kd_max: float = 0.20  # 1/day


def validate_parameters(p):
    values = asdict(p)
    if not all(np.isfinite(v) and v >= 0 for v in values.values()):
        raise ValueError("Parameters must be finite and nonnegative.")
    for key in ("duration", "dt", "initial_volume", "Ks", "Kl", "Ka",
                "glucose_uptake_K", "lactate_uptake_K"):
        if values[key] <= 0:
            raise ValueError(f"{key} must be positive.")
    if not p.initial_volume <= p.max_feed_volume <= p.reactor_capacity:
        raise ValueError("Require initial volume <= feeding cutoff <= capacity.")
    if p.kd_min > p.kd_max:
        raise ValueError("Require kd_min <= kd_max.")


def calculate_rates(state, p):
    glucose, product, lactate, ammonia, xv, xd, volume = np.maximum(state, 0)
    glucose_factor = glucose / (p.Ks + glucose)
    health = glucose_factor * p.Kl / (p.Kl + lactate) * p.Ka / (p.Ka + ammonia)
    mu = p.mu_max * health  # 1/day
    kd = p.kd_min + (p.kd_max - p.kd_min) * (1 - health)  # 1/day
    metabolic_state = expit(p.switch_steepness * (mu - p.mu_switch))
    availability = glucose / (p.glucose_uptake_K + glucose)
    # Production needs glucose; consumption cannot continue at zero lactate.
    ql = (metabolic_state * p.qL_production * availability
          - (1 - metabolic_state) * p.qL_consumption
          * lactate / (p.lactate_uptake_K + lactate))  # pg/(cell.day)
    return mu, kd, metabolic_state, ql, availability, health


def calculate_feed_rate(state, p, interval):
    glucose, *_, volume = state
    limit = min(p.max_feed_volume, p.reactor_capacity)
    if volume >= limit - 1e-10 or glucose >= p.glucose_target:
        return 0.0
    return min(p.feed_rate_on, max((limit - volume) / interval, 0.0))  # L/day


def derivatives(time, state, p, feed):
    glucose, product, lactate, ammonia, xv, xd, volume = np.maximum(state, 0)
    mu, kd, metabolic_state, ql, availability, health = calculate_rates(state, p)
    dilution = feed / volume  # 1/day; no separate dilution update is applied
    pg_to_concentration = xv * 1e9 * 1e-12  # cells/L * g/pg
    return np.array([
        dilution * (p.feed_glucose - glucose) - p.qG * availability * pg_to_concentration,
        p.qp * health * pg_to_concentration - dilution * product,
        ql * pg_to_concentration - dilution * lactate,
        p.qA * availability * pg_to_concentration - dilution * ammonia,
        (mu - kd - dilution) * xv,  # million cells/(mL.day)
        kd * xv - dilution * xd,  # million cells/(mL.day)
        feed,  # L/day
    ])


def simulate(p=Parameters(), scenario="Baseline", rtol=1e-7, max_step=0.05):
    """No file writes. Each exported feed is the command starting at that time."""
    validate_parameters(p)
    if max_step <= 0 or rtol <= 0:
        raise ValueError("Solver step and tolerance must be positive.")
    times = np.arange(0, p.duration, p.dt)
    times = np.append(times, p.duration)
    state = np.array([p.initial_glucose, 0., 0., 0., p.initial_viable_cells, 0., p.initial_volume])
    rows = []
    cumulative_feed = 0.0  # L, integrated applied flow, not a sum including final command
    for index, time in enumerate(times):
        interval = times[index + 1] - time if index + 1 < len(times) else p.dt
        feed = calculate_feed_rate(state, p, interval)
        mu, kd, metabolic, ql, _, _ = calculate_rates(state, p)
        glucose, product, lactate, ammonia, xv, xd, volume = state
        total = xv + xd
        rows.append(dict(Days=float(time), Glucose=glucose, Product=product,
                         Lactate=lactate, Ammonia=ammonia, Viable_Cells=xv,
                         Dead_Cells=xd, Total_Cells=total,
                         Viability=100 * xv / total if total > 0 else 0.,
                         Growth_Rate=mu, Death_Rate=kd, Metabolic_State=metabolic,
                         Dynamic_qL=ql, Reactor_Volume=volume, Feed_Rate=feed))
        if index + 1 == len(times):
            break
        solution = solve_ivp(derivatives, (time, times[index + 1]), state,
                             args=(p, feed), method="BDF", rtol=rtol,
                             atol=1e-10, max_step=min(max_step, interval))
        if not solution.success:
            raise RuntimeError(solution.message)
        state = solution.y[:, -1]
        if np.min(solution.y) < -1e-7:
            raise RuntimeError("Negative state exceeds numerical tolerance.")
        state = np.maximum(state, 0)  # roundoff only; depletion is handled by rates
        cumulative_feed += feed * interval
        if state[-1] > min(p.max_feed_volume, p.reactor_capacity) + 1e-8:
            raise RuntimeError("Volume limit exceeded.")
    results = pd.DataFrame(rows)
    last = results.iloc[-1]
    summary = dict(Simulation_Days=p.duration, Final_Product=last.Product,
                   Final_Lactate=last.Lactate, Final_Ammonia=last.Ammonia,
                   Final_Viability=last.Viability, Minimum_Glucose=results.Glucose.min(),
                   Maximum_Glucose=results.Glucose.max(), Final_Volume=last.Reactor_Volume,
                   Total_Feed=cumulative_feed, Scenario=scenario)
    return results, summary


def scenario_parameters():
    """Explicit new scenario definitions; historical parameter provenance is unavailable."""
    baseline = Parameters()
    overrides = {"Baseline": {}, "lower target glucose": {"glucose_target": 8.},
                 "higher target glucose": {"glucose_target": 12.},
                 "higher feed rate": {"feed_rate_on": .15},
                 "lower feed rate": {"feed_rate_on": .05},
                 "higher initial glucose": {"initial_glucose": 10.},
                 "lower feed glucose": {"feed_glucose": 60.},
                 "higher feed glucose": {"feed_glucose": 120.}}
    return {name: replace(baseline, **changes) for name, changes in overrides.items()}


def export_batch(scenarios, output_dir=None):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory = Path(output_dir or PROJECT_ROOT / "data" / "runs") / f"{stamp}_{uuid.uuid4().hex[:8]}"
    directory.mkdir(parents=True, exist_ok=False)
    all_results, summaries, metadata = [], [], []
    for simulation_id, (scenario, parameters) in enumerate(scenarios.items(), 1):
        results, summary = simulate(parameters, scenario)
        results.insert(0, "Simulation_ID", simulation_id)
        summary = {"Simulation_ID": simulation_id, **summary}
        all_results.append(results)
        summaries.append(summary)
        metadata.append(dict(Simulation_ID=simulation_id, Scenario=scenario,
                             Created_UTC=stamp, Model_Version=MODEL_VERSION,
                             Parameters_JSON=json.dumps(asdict(parameters), sort_keys=True)))
    results = pd.concat(all_results, ignore_index=True)
    summary = pd.DataFrame(summaries)
    with closing(sqlite3.connect(directory / "cho_digital_twin.db")) as connection:
        results.to_sql("simulation_results", connection, index=False)
        summary.to_sql("simulation_summary", connection, index=False)
        pd.DataFrame(metadata).to_sql("simulation_runs", connection, index=False)
        connection.commit()
    results.to_csv(directory / "simulation_results.csv", index=False)
    summary.to_csv(directory / "simulation_summary.csv", index=False)
    (directory / "manifest.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return directory, results
