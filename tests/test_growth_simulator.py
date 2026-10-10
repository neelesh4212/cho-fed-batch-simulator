"""Tests for the active simulator; real project exports are never modified."""

import json
import math
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import growth_simulator as model


class ModelTests(unittest.TestCase):
    def test_volume_balance(self):
        self.assertAlmostEqual(model.update_reactor_volume(1.0, 0.1, 0.5), 1.05)

    def test_feed_switch_and_cap(self):
        # Arrange inputs, act by calling the function, assert the expected flow.
        cases = [(1.0, 8.0, 0.1), (1.0, 10.0, 0.0),
                 (1.5, 8.0, 0.0), (1.6, 8.0, 0.0), (1.46, 8.0, 0.08)]
        for volume, glucose, expected in cases:
            with self.subTest(volume=volume, glucose=glucose):
                flow = model.calculate_feed_rate(volume, 1.5, glucose, 10.0, 0.1, 0.5)
                self.assertAlmostEqual(flow, expected)

    def test_health_and_death_limits(self):
        self.assertAlmostEqual(model.calculate_health(8, 4, 2, 2, 4, 2), 0.2)
        self.assertEqual(model.calculate_health(0, 0, 0, 2, 4, 2), 0)
        self.assertAlmostEqual(model.calculate_death_rate(1, 0.02, 0.2), 0.02)
        self.assertAlmostEqual(model.calculate_death_rate(0, 0.02, 0.2), 0.2)

    def test_metabolic_switch(self):
        self.assertAlmostEqual(model.calculate_metabolic_state(0.15, 0.15, 20), 0.5)
        self.assertGreater(model.calculate_metabolic_state(0.3, 0.15, 20), 0.5)
        self.assertLess(model.calculate_metabolic_state(0.0, 0.15, 20), 0.5)

    def test_exponential_growth_without_death(self):
        viable, dead = model.update_cells(1, 1, 1, 0, 0.2, 0, 1)
        self.assertAlmostEqual(viable, math.exp(0.2))
        self.assertEqual(dead, 0)

    def test_death_transfers_cells_without_changing_total(self):
        viable, dead = model.update_cells(1, 1, 1, 0, 0, 0.2, 1)
        self.assertAlmostEqual(viable, math.exp(-0.2))
        self.assertAlmostEqual(dead, 1 - math.exp(-0.2))
        self.assertAlmostEqual(viable + dead, 1)

    def test_equal_growth_and_death(self):
        # Net growth is zero, but births and deaths still occur.
        viable, dead = model.update_cells(1, 1, 1, 0, 0.2, 0.2, 1)
        self.assertAlmostEqual(viable, 1)
        self.assertAlmostEqual(dead, 0.2)

    def test_cell_dilution_preserves_numbers(self):
        viable, dead = model.update_cells(1, 1, 2, 0.4, 0, 0, 1)
        self.assertAlmostEqual(viable, 0.5)
        self.assertAlmostEqual(dead, 0.2)

    def test_metabolite_mass_balances(self):
        glucose, product, lactate, ammonia = model.update_metabolites(
            8, 0.2, 1, 0.1, 1, 1, 1.1, 100, 0.1, 600, 30, 250, 40, 0.8, 1,
        )
        # Initial mass + feed/formation - consumption, divided by final volume.
        self.assertAlmostEqual(glucose, (8 + 10 - 0.6) / 1.1)
        self.assertAlmostEqual(product, (0.2 + 0.024) / 1.1)
        self.assertAlmostEqual(lactate, (1 + 0.25) / 1.1)
        self.assertAlmostEqual(ammonia, (0.1 + 0.04) / 1.1)

    def test_depletion_caps(self):
        glucose, _, lactate, _ = model.update_metabolites(
            0.01, 0, 0.01, 0, 1, 1, 1, 0, 0, 600, 30, -50, 40, 0, 1,
        )
        self.assertEqual(glucose, 0)
        self.assertEqual(lactate, 0)

    def test_scenario_labels(self):
        parameters = model.BASELINE_PARAMETERS.copy()
        self.assertEqual(model.classify_scenario(parameters), "Baseline")
        for flow in [0.15, 0.2]:
            parameters["feed_rate_on"] = flow
            self.assertEqual(model.classify_scenario(parameters), f"higher feed rate: {flow:g}")
        parameters["qG"] = 500
        self.assertEqual(model.classify_scenario(parameters), "higher feed rate: 0.2; qG=500")


class ExportTests(unittest.TestCase):
    def test_main_exports_and_daily_controller(self):
        # Retain the in-memory database for inspection after main closes it.
        class InspectableConnection(sqlite3.Connection):
            def close(self):
                pass

        connection = sqlite3.connect(":memory:", factory=InspectableConnection)
        try:
            with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as directory:
                fake_script = Path(directory) / "src" / "growth_simulator.py"
                with patch.object(model, "__file__", str(fake_script)), \
                     patch.object(model.sqlite3, "connect", return_value=connection), \
                     patch.object(model, "plot_results"), \
                     patch.object(model, "calculate_feed_rate", wraps=model.calculate_feed_rate) as feed:
                    model.main()
                self.assertEqual(feed.call_count, 14)
                self.assertTrue((Path(directory) / "data/csv/simulation_1.csv").is_file())
                self.assertEqual(connection.execute("SELECT COUNT(*) FROM simulation_results").fetchone()[0], 113)
                schema = dict((row[1], row[2]) for row in connection.execute("PRAGMA table_info(simulation_results)"))
                self.assertEqual(schema["Days"], "REAL")
                times = connection.execute("SELECT Days FROM simulation_results ORDER BY Days LIMIT 3").fetchall()
                self.assertEqual(times, [(0.0,), (0.125,), (0.25,)])
                last = connection.execute("SELECT Days, Feed_Rate FROM simulation_results ORDER BY Days DESC LIMIT 1").fetchone()
                self.assertEqual(last, (14.0, None))
                row = connection.execute("SELECT Simulation_ID, Scenario, Total_Feed, Final_Volume FROM simulation_summary").fetchone()
                self.assertEqual(row[:2], (1, "Baseline"))
                self.assertAlmostEqual(row[2], row[3] - 1)
                saved = connection.execute("SELECT Parameters_JSON FROM simulation_metadata").fetchone()[0]
                self.assertEqual(json.loads(saved), model.BASELINE_PARAMETERS)
                states = connection.execute("SELECT Glucose, Product, Lactate, Ammonia, Viable_Cells, Dead_Cells, Viability, Reactor_Volume FROM simulation_results").fetchall()
                values = np.array(states)
                self.assertTrue(np.isfinite(values).all())
                self.assertTrue((values >= 0).all())
                self.assertTrue((values[:, 6] <= 100).all())
                self.assertTrue((values[:, 7] <= 1.5 + 1e-12).all())
        finally:
            sqlite3.Connection.close(connection)


if __name__ == "__main__":
    unittest.main()
