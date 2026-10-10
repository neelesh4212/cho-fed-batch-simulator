"""Isolated demonstration runs using syntax-tree overrides of main settings."""
import ast
import copy
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

SOURCE = Path(__file__).resolve().parents[2] / "src/growth_simulator.py"
SCENARIOS = [{}, {"feed_glucose": 150.0}, {"feed_glucose": 50.0},
             {"feed_rate_on": 0.15}, {"feed_rate_on": 0.05},
             {"initial_glucose": 12.0}, {"glucose_target": 12.0},
             {"glucose_target": 8.0}]


def main():
    source = SOURCE.read_text()
    original = ast.parse(source)
    output = Path(__file__).resolve().parent / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output.mkdir(parents=True, exist_ok=False)
    for overrides in SCENARIOS:
        tree = copy.deepcopy(original)
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
        found = set()
        for node in function.body:
            if not isinstance(node, ast.Assign):
                continue
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in overrides:
                node.value = ast.Constant(overrides[target.id])
                found.add(target.id)
            elif (isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name)
                  and target.value.id == "glucose" and isinstance(target.slice, ast.Constant)
                  and target.slice.value == 0 and "initial_glucose" in overrides):
                node.value = ast.Constant(overrides["initial_glucose"])
                found.add("initial_glucose")
        if found != set(overrides):
            raise ValueError(f"Unmatched settings: {set(overrides) - found}")
        # Redirect exports without changing the actual simulator file.
        namespace = {"__name__": "demo_model", "__file__": str(output / "src" / SOURCE.name)}
        exec(compile(ast.fix_missing_locations(tree), str(SOURCE), "exec"), namespace)
        namespace["plot_results"] = lambda *args: None
        namespace["main"]()
    database = output / "data/database/cho_digital_twin.db"
    with sqlite3.connect(database) as connection:
        summary = pd.read_sql_query("SELECT * FROM simulation_summary ORDER BY Simulation_ID", connection)
        assert connection.execute("SELECT COUNT(*), COUNT(DISTINCT Simulation_ID) FROM simulation_results").fetchone() == (904, 8)
        assert len(summary) == 8 and summary.Scenario.nunique() == 8
        assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert connection.execute("SELECT COUNT(*) FROM simulation_metadata").fetchone()[0] == 8
        for _, row in summary.iterrows():
            assert abs(row.Total_Feed - (row.Final_Volume - 1.0)) < 1e-10
        summary.to_csv(output / "simulation_summary.csv", index=False)
    (output / "manifest.json").write_text(json.dumps({
        "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "overrides": SCENARIOS,
        "description": "Illustrative comparisons, not experimental validation",
    }, indent=2))
    print("Output:", output)
    print(summary[["Scenario", "Final_Product", "Final_Viability", "Final_Volume"]].to_string(index=False))


if __name__ == "__main__":
    main()
