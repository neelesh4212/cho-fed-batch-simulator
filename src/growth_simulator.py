import numpy as np
import matplotlib.pyplot as plt


# =====================================================
# MODEL FUNCTIONS
# =====================================================

def calculate_health(glucose, lactate, ammonia, Ks, Kl, Ka):
    glucose_factor = glucose / (Ks + glucose)
    lactate_factor = Kl / (Kl + lactate)
    ammonia_factor = Ka / (Ka + ammonia)

    health_factor = glucose_factor * lactate_factor * ammonia_factor

    return health_factor


def calculate_death_rate(health_factor, kd_min, kd_max):
    kd = kd_min + (kd_max - kd_min) * (1 - health_factor)

    return kd


def calculate_growth(viable_cells, mu, dt):
    new_cells = viable_cells * (np.exp(mu * dt) - 1)

    return new_cells


def update_metabolites(
    glucose_previous,
    product_previous,
    lactate_previous,
    ammonia_previous,
    viable_cells_previous,
    qG,
    qp,
    qL,
    qA,
    health_factor,
    dt,
):
    cells_per_L = viable_cells_previous * 1e9
    dynamic_qp = qp * health_factor

    glucose_new = glucose_previous - (qG * cells_per_L * 1e-12 * dt)
    product_new = product_previous + (dynamic_qp * cells_per_L * 1e-12 * dt)
    lactate_new = lactate_previous + (qL * cells_per_L * 1e-12 * dt)
    ammonia_new = ammonia_previous + (qA * cells_per_L * 1e-12 * dt)

    glucose_new = max(glucose_new, 0)

    return glucose_new, product_new, lactate_new, ammonia_new


def update_cells(viable_cells_previous, dead_cells_previous, mu, kd, dt):
    new_cells = calculate_growth(viable_cells_previous, mu, dt)
    dead_today = kd * viable_cells_previous * dt

    viable_cells_new = viable_cells_previous + new_cells - dead_today
    viable_cells_new = max(viable_cells_new, 0)

    dead_cells_new = dead_cells_previous + dead_today

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


def plot_substrate_and_product(days, glucose, product):
    plt.figure(figsize=(8, 5))
    plt.plot(days, glucose, label="Glucose", color="green")
    plt.plot(days, product, label="Product", color="purple")
    plt.title("Substrate and Product Over Time")
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


def plot_lactate(days, lactate):
    plt.figure(figsize=(8, 5))
    plt.plot(days, lactate, label="Lactate", color="cyan")
    plt.title("Lactate Production Over Time")
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


def plot_results(
    days,
    viable_cells,
    dead_cells,
    total_cells,
    glucose,
    product,
    viability,
    kd,
    lactate,
    ammonia,
):
    plot_cell_population(days, viable_cells, dead_cells, total_cells)
    plot_substrate_and_product(days, glucose, product)
    plot_viability(days, viability)
    plot_death_rate(days, kd)
    plot_lactate(days, lactate)
    plot_ammonia(days, ammonia)


# =====================================================
# MAIN SIMULATION
# =====================================================

def main():

    # ----------------------------
    # Simulation settings
    # ----------------------------

    days = np.arange(0, 15)

    glucose = np.empty(len(days))
    product = np.empty(len(days))
    viable_cells = np.empty(len(days))
    dead_cells = np.empty(len(days))
    kd = np.empty(len(days))
    lactate = np.empty(len(days))
    ammonia = np.empty(len(days))

    # ----------------------------
    # Initial Conditions
    # ----------------------------

    glucose[0] = 30.0                   # g/L
    initial_cells = 1.0                 # million cells/mL
    product[0] = 0.0                    # g/L
    viable_cells[0] = initial_cells     
    dead_cells[0] = 0.0                 
    lactate[0] = 0.0                    # g/L
    ammonia[0] = 0.0                    # g/L

    # ----------------------------
    # Growth Parameters
    # ----------------------------

    mu_max = 0.4
    Ks = 2.0                            # g/L, glucose saturation constant
    Kl = 4.0                            # g/L, lactate inhibition constant
    Ka = 2.0                            # g/L, ammonia inhibition constant

    # ----------------------------
    # Metabolic Parameters
    # ----------------------------

    qG = 600.0                          #picograms / (cell . day)
    qp = 30.0                           #picograms / (cell . day)
    qL = 250.0                          #picograms / (cell . day)
    qA = 40.0                           #picograms / (cell . day)

    dt = 1.0                            # time step in days

    # ----------------------------
    # Death Parameters
    # ----------------------------

    kd_min = 0.02                       # 1/day, death rate when glucose is high    
    kd_max = 0.20                       # 1/day, death rate when glucose is depleted

    kd[0] = kd_min                      # 1/day, initial death rate

    # ----------------------------
    # Simulation Loop
    # ----------------------------

    for i in range(1, len(days)):

        health_factor = calculate_health(
            glucose[i-1],
            lactate[i-1],
            ammonia[i-1],
            Ks,
            Kl,
            Ka,
        )

        mu = mu_max * health_factor

        kd[i] = calculate_death_rate(
            health_factor,
            kd_min,
            kd_max,
        )

        glucose[i], product[i], lactate[i], ammonia[i] = update_metabolites(
            glucose[i-1],
            product[i-1],
            lactate[i-1],
            ammonia[i-1],
            viable_cells[i-1],
            qG,
            qp,
            qL,
            qA,
            health_factor,
            dt,
        )

        viable_cells[i], dead_cells[i] = update_cells(
            viable_cells[i-1],
            dead_cells[i-1],
            mu,
            kd[i],
            dt,
        )

    # ----------------------------
    # Derived Variables
    # ----------------------------

    total_cells = viable_cells + dead_cells
    viability = viable_cells / total_cells * 100

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
        lactate,
        ammonia,
    )


# =====================================================
# ENTRY POINT
# =====================================================

if __name__ == "__main__":
    main()