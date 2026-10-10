import numpy as np
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
from pathlib import Path
import json
from datetime import datetime, timezone


# =====================================================
# MODEL FUNCTIONS
# =====================================================

# CHANGED: fixed reference settings for automatic scenario classification.
BASELINE_PARAMETERS = {
    "simulation_duration": 14.0, "dt": 0.125, "controller_interval": 1.0,
    "initial_volume": 1.0, "initial_glucose": 8.0, "initial_cells": 1.0,
    "initial_product": 0.0, "initial_dead_cells": 0.0,
    "initial_lactate": 0.0, "initial_ammonia": 0.0,
    "feed_glucose": 100.0, "feed_rate_on": 0.1, "glucose_target": 10.0,
    "max_feed_volume": 1.5, "mu_max": 0.4, "Ks": 2.0, "Kl": 4.0,
    "Ka": 2.0, "qG": 600.0, "qp": 30.0, "qA": 40.0,
    "qL_production": 250.0, "qL_consumption": 50.0,
    "mu_switch": 0.15, "switch_steepness": 20.0,
    "kd_min": 0.02, "kd_max": 0.20,
}


def classify_scenario(parameters):
    labels = {
        "feed_glucose": "feed glucose", "feed_rate_on": "feed rate",
        "initial_glucose": "initial glucose", "glucose_target": "target glucose",
    }
    changes = []
    for name, baseline in BASELINE_PARAMETERS.items():
        value = parameters[name]
        if value != baseline:
            if name in labels:
                direction = "higher" if value > baseline else "lower"
                changes.append(f"{direction} {labels[name]}: {value:g}")
            else:
                changes.append(f"{name}={value:g}")
    return "; ".join(changes) if changes else "Baseline"

def calculate_health(glucose, lactate, ammonia, Ks, Kl, Ka):
    glucose_factor = glucose / (Ks + glucose)
    lactate_factor = Kl / (Kl + lactate)
    ammonia_factor = Ka / (Ka + ammonia)

    health_factor = glucose_factor * lactate_factor * ammonia_factor

    return health_factor


def calculate_death_rate(health_factor, kd_min, kd_max):
    kd = kd_min + (kd_max - kd_min) * (1 - health_factor)

    return kd


def calculate_metabolic_state(mu, mu_switch, switch_steepness):
    metabolic_state = 1 / (
        1 + np.exp(
            -switch_steepness * (mu - mu_switch)
        )
    )

    return metabolic_state


def update_reactor_volume(V_previous, feed_rate, dt):
    V_new = V_previous + feed_rate * dt

    return V_new


def calculate_growth(viable_cells, mu, dt):
    new_cells = viable_cells * (np.exp(mu * dt) - 1)

    return new_cells

def calculate_feed_rate(
        reactor_volume_previous, 
        max_feed_volume, 
        glucose_previous, 
        glucose_target, 
        feed_rate_on, 
        dt,
):
    if reactor_volume_previous >= max_feed_volume or glucose_previous >= glucose_target:
        feed_rate = 0.0

    else:
        feed_rate = feed_rate_on

    feed_rate = min(
        feed_rate, 
        max(((max_feed_volume - reactor_volume_previous) / dt), 0.0))

    return feed_rate

def update_metabolites(
    glucose_previous,
    product_previous,
    lactate_previous,
    ammonia_previous,
    viable_cells_previous,
    reactor_volume_previous,
    reactor_volume_new,
    feed_glucose,
    feed_rate,
    qG,
    qp,
    qL_dynamic,
    qA,
    health_factor,
    dt,
):
    cells_per_L = viable_cells_previous * 1e9  # million cells/mL to cells/L
    dynamic_qp = qp * health_factor  # pg/(cell.day)

    glucose_previous_mass = glucose_previous * reactor_volume_previous  # g
    glucose_added = feed_glucose * feed_rate * dt  # g
    glucose_demand = qG * cells_per_L * reactor_volume_previous * 1e-12 * dt  # g
    actual_glucose_consumed = min(glucose_previous_mass + glucose_added, glucose_demand)  # g

    glucose_new = (
        glucose_previous_mass + glucose_added - actual_glucose_consumed
    ) / reactor_volume_new  # g/L

    product_previous_mass = product_previous * reactor_volume_previous  # g
    product_formed = dynamic_qp * cells_per_L * 1e-12 * dt * reactor_volume_previous  # g
    product_new = (product_previous_mass + product_formed) / reactor_volume_new  # g/L

    lactate_previous_mass = lactate_previous * reactor_volume_previous  # g
    lactate_change = qL_dynamic * cells_per_L * 1e-12 * dt * reactor_volume_previous  # g
    lactate_new = (lactate_previous_mass + lactate_change) / reactor_volume_new  # g/L
    lactate_new = max(lactate_new, 0)  # Prevent negative lactate concentration

    ammonia_previous_mass = ammonia_previous * reactor_volume_previous  # g
    ammonia_formed = qA * cells_per_L * 1e-12 * dt * reactor_volume_previous  # g
    ammonia_new = (ammonia_previous_mass + ammonia_formed) / reactor_volume_new

    return glucose_new, product_new, lactate_new, ammonia_new


def update_cells(viable_cells_previous, reactor_volume_previous, reactor_volume_new, dead_cells_previous, mu, kd, dt):
    # CHANGED: solve growth and death together with rates held constant over dt.
    viable_cells_previous_mass = viable_cells_previous * 1e9 * reactor_volume_previous  # cells
    net_growth_rate = mu - kd  # 1/day
    viable_cells_new_mass = viable_cells_previous_mass * np.exp(net_growth_rate * dt)  # cells

    # CHANGED: integrate viable cell exposure to count deaths during the step.
    if abs(net_growth_rate * dt) < 1e-8:
        viable_cell_exposure = viable_cells_previous_mass * dt  # cell.day
    else:
        viable_cell_exposure = (
            viable_cells_previous_mass * np.expm1(net_growth_rate * dt) / net_growth_rate
        )  # cell.day
    dead_today_mass = kd * viable_cell_exposure  # cells

    # CHANGED: volume balance accounts for dilution once, after biological change.
    viable_cells_new = viable_cells_new_mass / reactor_volume_new / 1e9  # million cells/mL
    dead_cells_new = (
        dead_cells_previous * 1e9 * reactor_volume_previous + dead_today_mass
    ) / reactor_volume_new / 1e9  # million cells/mL

    return viable_cells_new, dead_cells_new

# =====================================================
# PLOTTING FUNCTIONS
# =====================================================

def plot_cell_population(days, viable_cells, dead_cells, total_cells):
    plt.figure(figsize=(8, 5))
    plt.plot(days, viable_cells, label="Viable Cells", color="blue")
    plt.plot(days, dead_cells, label="Dead Cells", color="red")
    plt.plot(days, total_cells, label="Total Cells", color="orange")
    plt.title("Cell Population Dynamics Over Time")
    plt.xlabel("Days")
    plt.ylabel("Million Cells/mL")
    plt.legend()
    plt.grid()
    plt.show()


def plot_substrate(days, glucose):
    plt.figure(figsize=(8, 5))
    plt.plot(days, glucose, label="Glucose", color="green")
    plt.title("Substrate Concentration Over Time")
    plt.xlabel("Days")
    plt.ylabel("Concentration (g/L)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_product(days, product):
    plt.figure(figsize=(8, 5))
    plt.plot(days, product, label="Product", color="purple")
    plt.title("Product Formation Over Time")
    plt.xlabel("Days")
    plt.ylabel("Concentration (g/L)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_viability(days, viability):
    plt.figure(figsize=(8, 5))
    plt.plot(days, viability, label="Viability", color="brown")
    plt.title("Cell Viability Over Time")
    plt.xlabel("Days")
    plt.ylabel("Viability (%)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_death_rate(days, kd):
    plt.figure(figsize=(8, 5))
    plt.plot(days, kd, label="Death Rate", color="black")
    plt.title("Cell Death Rate Over Time")
    plt.xlabel("Days")
    plt.ylabel("Death Rate (1/day)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_growth_rate(days, mu):
    plt.figure(figsize=(8, 5))
    plt.plot(days, mu, label="Growth Rate", color="darkgreen")
    plt.title("Cell Growth Rate Over Time")
    plt.xlabel("Days")
    plt.ylabel("Growth Rate (1/day)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_lactate(days, lactate):
    plt.figure(figsize=(8, 5))
    plt.plot(days, lactate, label="Lactate", color="cyan")
    plt.title("Lactate Concentration Over Time")
    plt.xlabel("Days")
    plt.ylabel("Lactate Concentration (g/L)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_ammonia(days, ammonia):
    plt.figure(figsize=(8, 5))
    plt.plot(days, ammonia, label="Ammonia", color="magenta")
    plt.title("Ammonia Production Over Time")
    plt.xlabel("Days")
    plt.ylabel("Ammonia Concentration (g/L)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_metabolic_state(days, metabolic_state):
    plt.figure(figsize=(8, 5))
    plt.plot(days, metabolic_state, label="Metabolic State", color="navy")
    plt.title("Metabolic State Over Time")
    plt.xlabel("Days")
    plt.ylabel("Metabolic State")
    plt.legend()
    plt.grid()
    plt.show()


def plot_dynamic_lactate_rate(days, dynamic_qL):
    plt.figure(figsize=(8, 5))
    plt.plot(days, dynamic_qL, label="Dynamic Lactate Rate", color="teal")
    plt.title("Dynamic Lactate Rate Over Time")
    plt.xlabel("Days")
    plt.ylabel("qL Dynamic pg/(cell.day)")
    plt.legend()
    plt.grid()
    plt.show()

def plot_reactor_volume(days, reactor_volume):
    plt.figure(figsize=(8, 5))
    plt.plot(days, reactor_volume, label="Reactor Volume", color="orange")
    plt.title("Reactor Volume Over Time")
    plt.xlabel("Days")
    plt.ylabel("Volume (L)")
    plt.legend()
    plt.grid()
    plt.show()

def plot_results(
    days,
    viable_cells,
    dead_cells,
    total_cells,
    glucose,
    product,
    viability,
    kd,
    mu,
    lactate,
    ammonia,
    metabolic_state,
    dynamic_qL, 
    reactor_volume
):
    plot_cell_population(days, viable_cells, dead_cells, total_cells)
    plot_substrate(days, glucose)
    plot_product(days, product)
    plot_viability(days, viability)
    plot_death_rate(days, kd)
    plot_growth_rate(days, mu)
    plot_lactate(days, lactate)
    plot_ammonia(days, ammonia)
    plot_metabolic_state(days, metabolic_state)
    plot_dynamic_lactate_rate(days, dynamic_qL)
    plot_reactor_volume(days, reactor_volume)

# =====================================================
# MAIN SIMULATION
# =====================================================

def main():

    # ----------------------------
    # Simulation settings
    # ----------------------------

    # CHANGED: timestamps follow dt and include the exact simulation endpoint.
    simulation_duration = 14.0  # days
    # CHANGED: label this run for the Power BI scenario slicer; update when parameters change.
    # Scenario is now generated automatically from the parameters below.
    dt = 0.125  # days, three-hour biological steps
    # CHANGED: feed decisions remain daily, independent of biological step size.
    controller_interval = 1.0  # days between feed-command calculations
    # FUTURE: online glucose sensing, measurement delays, and shorter control intervals.

    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be finite and greater than zero.")
    if not np.isfinite(simulation_duration) or simulation_duration <= 0:
        raise ValueError("simulation_duration must be finite and greater than zero.")
    if not np.isfinite(controller_interval) or controller_interval <= 0:
        raise ValueError("controller_interval must be finite and greater than zero.")

    days = np.arange(0.0, simulation_duration, dt)  # days
    days = np.append(days, simulation_duration)  # days
    # CHANGED: split biological intervals at every controller update time.
    controller_times = np.arange(0.0, simulation_duration, controller_interval)  # days
    days = np.unique(np.concatenate((days, controller_times)))  # days, sorted
    
    reactor_volume = np.empty(len(days))  # L
    feed_rate = np.zeros(len(days))  # L/day
    glucose = np.empty(len(days))  # g/L
    product = np.empty(len(days))  # g/L
    viable_cells = np.empty(len(days))  # million cells/mL
    dead_cells = np.empty(len(days))  # million cells/mL
    kd = np.empty(len(days))  # 1/day
    mu = np.empty(len(days))  # 1/day
    lactate = np.empty(len(days))  # g/L
    ammonia = np.empty(len(days))  # g/L
    metabolic_state = np.empty(len(days))  # dimensionless, 0-1
    dynamic_qL = np.empty(len(days))  # pg/(cell.day)

    # ----------------------------
    # Initial Conditions
    # ----------------------------

    reactor_volume[0] = 1.0  # L
    glucose[0] = 8.0  # g/L
    initial_cells = 1.0  # million cells/mL
    product[0] = 0.0  # g/L
    viable_cells[0] = initial_cells  # million cells/mL
    dead_cells[0] = 0.0  # million cells/mL
    lactate[0] = 0.0  # g/L
    ammonia[0] = 0.0  # g/L

    # ----------------------------
    # Feed Parameters
    # ----------------------------

    feed_glucose = 100.0  # g/L, glucose concentration in feed
    feed_rate_on = 0.1  # L/day
    glucose_target = 10.0  # g/L, target glucose concentration
    max_feed_volume = 1.5  # L, maximum reactor volume for feeding
    # reactor_capacity = 2.0  # L, maximum reactor volume

    # ----------------------------
    # Growth Parameters
    # ----------------------------

    mu_max = 0.4  # 1/day
    Ks = 2.0  # g/L, glucose saturation constant
    Kl = 4.0  # g/L, lactate inhibition constant
    Ka = 2.0  # g/L, ammonia inhibition constant

    # ----------------------------
    # Metabolic Parameters
    # ----------------------------

    qG = 600.0  # pg/(cell.day), glucose consumption rate
    qp = 30.0  # pg/(cell.day), product formation rate
    qA = 40.0  # pg/(cell.day), ammonia production rate

    # Lactate metabolism parameters
    qL_production = 250.0  # pg/(cell.day)
    qL_consumption = 50.0  # pg/(cell.day)

    mu_switch = 0.15  # 1/day
    switch_steepness = 20.0  # dimensionless, controls LP -> LC transition

    # ----------------------------
    # Death Parameters
    # ----------------------------

    kd_min = 0.02  # 1/day, death rate under healthier conditions
    kd_max = 0.20  # 1/day, death rate under stressed conditions

    # CHANGED: snapshot all settings before simulation changes the state arrays.
    settings = locals().copy()
    parameters = {
        name: float(settings[name])
        for name in BASELINE_PARAMETERS if name in settings
    }
    parameters.update({
        "initial_volume": float(reactor_volume[0]),
        "initial_glucose": float(glucose[0]),
        "initial_product": float(product[0]),
        "initial_dead_cells": float(dead_cells[0]),
        "initial_lactate": float(lactate[0]),
        "initial_ammonia": float(ammonia[0]),
    })
    scenario = classify_scenario(parameters)

    initial_health = calculate_health(
        glucose[0],
        lactate[0],
        ammonia[0],
        Ks,
        Kl,
        Ka,
    )

    mu[0] = mu_max * initial_health  # 1/day

    kd[0] = calculate_death_rate(
        initial_health,
        kd_min,
        kd_max,
    )

    metabolic_state[0] = calculate_metabolic_state(
        mu[0],
        mu_switch,
        switch_steepness,
    )

    dynamic_qL[0] = (
        metabolic_state[0] * qL_production
        - (1 - metabolic_state[0]) * qL_consumption
    )  # pg/(cell.day)

    # ----------------------------
    # Simulation Loop
    # ----------------------------

    # CHANGED: hold the feed command between scheduled controller updates.
    controller_index = 0
    feed_command = 0.0  # L/day

    for i in range(1, len(days)):

        # CHANGED: the last interval may be shorter than the configured dt.
        dt_step = days[i] - days[i-1]  # days

        if (
            controller_index < len(controller_times)
            and days[i-1] >= controller_times[controller_index]
        ):
            command_duration = min(
                controller_interval, simulation_duration - days[i-1]
            )  # days until the next controller decision or simulation endpoint
            feed_command = calculate_feed_rate(
                reactor_volume[i-1],
                max_feed_volume,
                glucose[i-1],
                glucose_target,
                feed_rate_on,
                command_duration,
            )  # L/day, held until the next daily decision
            controller_index += 1

        # CHANGED: volume protection still applies between controller decisions.
        max_allowed_feed_rate = max(
            (max_feed_volume - reactor_volume[i-1]) / dt_step, 0.0
        )  # L/day
        feed_rate[i-1] = min(feed_command, max_allowed_feed_rate)  # L/day, applied flow

        health_factor = calculate_health(
            glucose[i-1],
            lactate[i-1],
            ammonia[i-1],
            Ks,
            Kl,
            Ka,
        )

        # CHANGED: starting-condition rates advance the interval, not the plotted endpoint.
        mu_for_step = mu_max * health_factor  # 1/day

        current_metabolic_state = calculate_metabolic_state(
            mu_for_step,
            mu_switch,
            switch_steepness,
        )

        qL_dynamic = (
            current_metabolic_state * qL_production
            - (1 - current_metabolic_state) * qL_consumption
        )  # pg/(cell.day)

        kd_for_step = calculate_death_rate(
            health_factor,
            kd_min,
            kd_max,
        )

        reactor_volume[i] = update_reactor_volume(
            reactor_volume[i-1],
            feed_rate[i-1],
            dt_step,  # CHANGED: actual interval duration
        )

        glucose[i], product[i], lactate[i], ammonia[i] = update_metabolites(
            glucose[i-1],
            product[i-1],
            lactate[i-1],
            ammonia[i-1],
            viable_cells[i-1],
            reactor_volume[i-1],
            reactor_volume[i],
            feed_glucose,
            feed_rate[i-1],
            qG,
            qp,
            qL_dynamic,
            qA,
            health_factor,
            dt_step,  # CHANGED: actual interval duration
        )

        viable_cells[i], dead_cells[i] = update_cells(
            viable_cells[i-1],
            reactor_volume[i-1],
            reactor_volume[i],
            dead_cells[i-1],
            mu_for_step,
            kd_for_step,
            dt_step,  # CHANGED: actual interval duration
        )

        # CHANGED: store rates corresponding to the updated concentrations at day i.
        endpoint_health = calculate_health(
            glucose[i], lactate[i], ammonia[i], Ks, Kl, Ka,
        )
        mu[i] = mu_max * endpoint_health  # 1/day
        kd[i] = calculate_death_rate(endpoint_health, kd_min, kd_max)  # 1/day
        metabolic_state[i] = calculate_metabolic_state(
            mu[i], mu_switch, switch_steepness,
        )  # dimensionless, 0-1
        dynamic_qL[i] = (
            metabolic_state[i] * qL_production
            - (1 - metabolic_state[i]) * qL_consumption
        )  # pg/(cell.day)

    # ----------------------------
    # Derived Variables
    # ----------------------------

    # CHANGED: no feeding interval follows the final timestamp; export blank/SQL NULL.
    feed_rate[-1] = np.nan  # L/day, not applicable at the simulation endpoint

    total_cells = viable_cells + dead_cells  # million cells/mL
    viability = viable_cells / total_cells * 100  # %
    # CHANGED: integrate applied flow over each interval; exclude the final missing rate.
    total_feed = np.sum(feed_rate[:-1] * np.diff(days))  # L delivered

    # CHANGED: resolve output paths from the script, not the terminal directory.
    project_root = Path(__file__).resolve().parents[1]
    database_path = project_root / "data" / "database" / "cho_digital_twin.db"

    # CHANGED: create output folders if they do not already exist.
    database_path.parent.mkdir(parents=True, exist_ok=True)

    connection = sqlite3.connect(str(database_path))
    
    connection.execute("""
        CREATE TABLE IF NOT EXISTS simulation_results (
            Simulation_ID INTEGER,
            Days INTEGER,
            Glucose REAL,
            Product REAL,
            Lactate REAL,
            Ammonia REAL,
            Viable_Cells REAL,
            Dead_Cells REAL,
            Total_Cells REAL,
            Viability REAL,
            Growth_Rate REAL,
            Death_Rate REAL,
            Metabolic_State REAL,
            Dynamic_qL REAL,
            Reactor_Volume REAL,
            Feed_Rate REAL
        )
        """)
    
    next_sim_id = connection.execute("""
    SELECT COALESCE(MAX(Simulation_ID), 0) + 1 AS Next_Simulation_ID
    FROM simulation_results
    """)

    simulation_id = next_sim_id.fetchone()[0]

    # CHANGED: use the database simulation ID to give each run its own CSV filename.
    csv_path = project_root / "data" / "csv" / f"simulation_{simulation_id}.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    # =====================================================
    # CREATING DATAFRAMES
    # =====================================================
    
    results = pd.DataFrame({
        "Simulation_ID": simulation_id,
        "Days": days, 
        "Glucose (g/L)": glucose,
        "Product (g/L)": product,
        "Lactate (g/L)": lactate,
        "Ammonia (g/L)": ammonia,
        "Viable Cells (million cells/mL)": viable_cells,
        "Dead Cells (million cells/mL)": dead_cells,
        "Total Cells (million cells/mL)": total_cells,
        "Viability (%)": viability,
        "Growth Rate (1/day)": mu,
        "Death Rate (1/day)": kd,
        "Metabolic State": metabolic_state,
        "Dynamic qL (pg/(cell.day))": dynamic_qL,
        "Reactor Volume (L)": reactor_volume,
        "Feed Rate (L/day)": feed_rate
    })

    database_results = results.rename(columns={
        "Simulation_ID": "Simulation_ID",
        "Days": "Days",
        "Glucose (g/L)": "Glucose",
        "Product (g/L)": "Product",
        "Lactate (g/L)": "Lactate",
        "Ammonia (g/L)": "Ammonia",
        "Viable Cells (million cells/mL)": "Viable_Cells",
        "Dead Cells (million cells/mL)": "Dead_Cells",
        "Total Cells (million cells/mL)": "Total_Cells",
        "Viability (%)": "Viability",
        "Growth Rate (1/day)": "Growth_Rate",
        "Death Rate (1/day)": "Death_Rate",
        "Metabolic State": "Metabolic_State",
        "Dynamic qL (pg/(cell.day))": "Dynamic_qL",
        "Reactor Volume (L)": "Reactor_Volume",
        "Feed Rate (L/day)": "Feed_Rate"
    })

    database_results.to_sql(
        "simulation_results",
        connection,
        if_exists="append",
        index=False
    )

    # CHANGED: save one summary row with the same ID as this run's time-series results.
    simulation_summary = pd.DataFrame([{
        "Simulation_ID": simulation_id,
        "Simulation_Days": days[-1] - days[0],  # days
        "Final_Product": product[-1],  # g/L
        "Final_Lactate": lactate[-1],  # g/L
        "Final_Ammonia": ammonia[-1],  # g/L
        "Final_Viability": viability[-1],  # %
        "Minimum_Glucose": np.min(glucose),  # g/L, sampled timestamps
        "Maximum_Glucose": np.max(glucose),  # g/L, sampled timestamps
        "Final_Volume": reactor_volume[-1],  # L
        "Total_Feed": total_feed,  # L delivered
        "Scenario": scenario,
    }])
    simulation_summary.to_sql(
        "simulation_summary",
        connection,
        if_exists="append",
        index=False,
    )

    # CHANGED: retain reproducible settings without changing Power BI's existing tables.
    run_metadata = pd.DataFrame([{
        "Simulation_ID": simulation_id,
        "Scenario": scenario,
        "Created_UTC": datetime.now(timezone.utc).isoformat(),
        "Model_Version": "interval_exponential_v1",
        "Parameters_JSON": json.dumps(parameters, sort_keys=True),
    }])
    run_metadata.to_sql("simulation_metadata", connection, if_exists="append", index=False)

    connection.close()

    results.to_csv(
    csv_path,
    index=False
    )

    # ----------------------------
    # Plot
    # ----------------------------

    plot_results(
        days,
        viable_cells,
        dead_cells,
        total_cells,
        glucose,
        product,
        viability,
        kd,
        mu,
        lactate,
        ammonia,
        metabolic_state,
        dynamic_qL, 
        reactor_volume
    )


# =====================================================
# ENTRY POINT
# =====================================================

if __name__ == "__main__":
    main()
