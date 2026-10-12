"""In-memory adapter for the existing interval model; no exports or plotting."""
import numpy as np
import pandas as pd
from growth_simulator import (
    BASELINE_PARAMETERS, calculate_health, calculate_death_rate,
    calculate_metabolic_state, calculate_feed_rate, update_reactor_volume,
    update_metabolites, update_cells, classify_scenario,
)


def simulate(overrides=None):
    p = BASELINE_PARAMETERS.copy()
    if overrides:
        unknown = set(overrides) - set(p)
        if unknown:
            raise ValueError(f"Unknown parameters: {unknown}")
        p.update(overrides)
    if any(not np.isfinite(v) or v < 0 for v in p.values()):
        raise ValueError("Parameters must be finite and nonnegative.")
    for key in ("dt", "simulation_duration", "controller_interval", "initial_volume", "initial_cells", "Ks", "Kl", "Ka"):
        if p[key] <= 0:
            raise ValueError(f"{key} must be greater than zero.")
    if p["kd_min"] > p["kd_max"]:
        raise ValueError("Minimum death rate cannot exceed maximum death rate.")
    if p["initial_volume"] > p["max_feed_volume"]:
        raise ValueError("Initial volume cannot exceed the feeding cutoff.")
    duration, dt, interval = p["simulation_duration"], p["dt"], p["controller_interval"]
    if duration / dt + duration / interval > 10000:
        raise ValueError("Choose settings with fewer than 10,000 simulation steps.")
    commands = np.arange(0, duration, interval)
    days = np.unique(np.concatenate((np.arange(0, duration, dt), [duration], commands)))
    columns = ["Glucose", "Product", "Lactate", "Ammonia", "Viable_Cells", "Dead_Cells", "Reactor_Volume", "Growth_Rate", "Death_Rate", "Metabolic_State", "Dynamic_qL"]
    a = {name: np.empty(len(days)) for name in columns}
    for name, key in [("Glucose", "initial_glucose"), ("Product", "initial_product"),
                      ("Lactate", "initial_lactate"), ("Ammonia", "initial_ammonia"),
                      ("Viable_Cells", "initial_cells"), ("Dead_Cells", "initial_dead_cells"),
                      ("Reactor_Volume", "initial_volume")]:
        a[name][0] = p[key]
    feed = np.zeros(len(days))

    def rates(i):
        health = calculate_health(a["Glucose"][i], a["Lactate"][i], a["Ammonia"][i], p["Ks"], p["Kl"], p["Ka"])
        a["Growth_Rate"][i] = p["mu_max"] * health
        a["Death_Rate"][i] = calculate_death_rate(health, p["kd_min"], p["kd_max"])
        state = calculate_metabolic_state(a["Growth_Rate"][i], p["mu_switch"], p["switch_steepness"])
        a["Metabolic_State"][i] = state
        a["Dynamic_qL"][i] = state * p["qL_production"] - (1 - state) * p["qL_consumption"]
        return health

    health = rates(0)
    index, command = 0, 0.0
    for i in range(1, len(days)):
        h = days[i] - days[i - 1]  # days
        old_v = a["Reactor_Volume"][i - 1]
        if index < len(commands) and days[i - 1] >= commands[index]:
            command = calculate_feed_rate(old_v, p["max_feed_volume"], a["Glucose"][i - 1], p["glucose_target"], p["feed_rate_on"], min(interval, duration - days[i - 1]))
            index += 1
        feed[i - 1] = min(command, max((p["max_feed_volume"] - old_v) / h, 0))
        new_v = update_reactor_volume(old_v, feed[i - 1], h)
        a["Reactor_Volume"][i] = new_v
        metabolites = update_metabolites(
            a["Glucose"][i - 1], a["Product"][i - 1], a["Lactate"][i - 1], a["Ammonia"][i - 1],
            a["Viable_Cells"][i - 1], old_v, new_v, p["feed_glucose"], feed[i - 1],
            p["qG"], p["qp"], a["Dynamic_qL"][i - 1], p["qA"], health, h,
        )
        for name, value in zip(columns[:4], metabolites):
            a[name][i] = value
        a["Viable_Cells"][i], a["Dead_Cells"][i] = update_cells(
            a["Viable_Cells"][i - 1], old_v, new_v, a["Dead_Cells"][i - 1],
            a["Growth_Rate"][i - 1], a["Death_Rate"][i - 1], h,
        )
        health = rates(i)
    feed[-1] = np.nan
    result = pd.DataFrame({"Days": days, **a, "Feed_Rate": feed})
    result["Total_Cells"] = result.Viable_Cells + result.Dead_Cells
    result["Viability"] = result.Viable_Cells / result.Total_Cells * 100
    if not np.isfinite(result.drop(columns="Feed_Rate").to_numpy()).all():
        raise ValueError("These settings produced non-finite results.")
    summary = {
        "Scenario": classify_scenario(p), "Final_Product": result.Product.iloc[-1],
        "Final_Viability": result.Viability.iloc[-1], "Final_Lactate": result.Lactate.iloc[-1],
        "Final_Volume": result.Reactor_Volume.iloc[-1],
        "Total_Feed": float(np.sum(feed[:-1] * np.diff(days))),
    }
    return result, summary, p
