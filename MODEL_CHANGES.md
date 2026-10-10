# ODE simulator release

## Run

From any directory, run the project Python interpreter with:

```powershell
python src/growth_simulator.py --all-scenarios --no-plots
python -m unittest discover -s tests -v
```

Install requirements.txt in a suitable Python environment. The simulation_engine
module separates parameters, pure simulation, and versioned exports. Existing
plot functions and legacy update helpers remain available in growth_simulator.py;
the command-line simulation uses the new engine, not those legacy helpers.

## Equations and assumptions

State order is glucose, product, lactate, ammonia (g/L), viable/dead cells
(million cells/mL), and reactor volume (L). Time is days. For cell-free feed and
no outlet, D = F/V (1/day):

- dXv/dt = (mu - kd - D) Xv
- dXd/dt = kd Xv - D Xd; dead-cell lysis is not included.
- dG/dt = D(Gfeed - G) - qG aG Xv * 0.001
- dP/dt = qp health Xv * 0.001 - D P
- dL/dt = qL_dynamic Xv * 0.001 - D L
- dA/dt = qA aG Xv * 0.001 - D A
- dV/dt = F

The factor 0.001 converts million cells/mL and pg/(cell.day) to g/(L.day).
The original glucose/lactate/ammonia health and death formulas are retained.
Sigmoid steepness has units days, not dimensionless.

New assumptions: aG = G/(0.1 + G) limits uptake, ammonia formation, and lactate
production as glucose depletes. Lactate consumption is multiplied by
L/(0.1 + L). Both 0.1 g/L saturation constants are placeholders requiring
experimental calibration. Product and growth already vanish with health at
zero glucose. Lactate consumption can still occur without glucose. No biomass
yield, elemental balance, glutamine limitation, or experimental validation is
claimed. This is an educational mechanistic simulator, not a validated twin.

The BDF solver reevaluates biology continuously, with maximum step 0.05 days.
The controller samples once per dt (default 1 day) and holds the capped command
until the next sample. A changing dt changes controller sampling as well as
export times, so it is not simply a numerical-accuracy setting. Solver tolerance
and max_step can be changed independently for numerical convergence checks.

Each row's rates are evaluated from that row's state. Feed_Rate is the outgoing
interval command. The final row contains the requested next command, not an
artificial zero; it is not applied beyond the simulation horizon. Total_Feed is
integrated over simulated intervals, excluding that last command. Min/max glucose
summaries refer to exported samples, not all internal solver evaluations.

## Scenario and output provenance

Eight named scenarios are explicitly defined in scenario_parameters(). Historical
settings were manually changed and not recorded; the new overrides must not be
presented as recovered historical settings. New values: targets 8/12 g/L, flow
0.05/0.15 L/day, initial glucose 10 g/L, feed glucose 60/120 g/L. Baseline retains
8 g/L initial glucose, 100 g/L feed, 0.1 L/day flow, and 10 g/L target.

Each batch creates data/runs/<UTC timestamp>_<unique suffix>/ containing a new
database, results CSV, summary CSV, and parameter manifest. Simulation_ID is local
to each database. There is exactly one summary row per scenario per batch.
Old database/CSV files and CHO.pbix are not overwritten.

The simulation_results and simulation_summary names and existing columns are
preserved. Fractional Days are supported. simulation_runs records model version,
parameters, scenario, and creation time.

For Power BI, change the existing database source to the chosen new database
(retain your existing connector), refresh both tables, and confirm the active
one-to-many relationship on Simulation_ID. Do not combine batches by Simulation_ID
alone. Keep the new manifest with the dashboard's source snapshot. Corrected
results will differ from the old daily mixed-update outputs.
