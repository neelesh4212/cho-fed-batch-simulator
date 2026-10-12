import json
import sqlite3
import sys
import unittest
from contextlib import closing
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from web_simulation import simulate


class WebSimulationTests(unittest.TestCase):
    def test_matches_all_published_demo_runs(self):
        path = ROOT / "data/demo/20261010T184644691930Z/data/database/cho_digital_twin.db"
        with closing(sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)) as c:
            for run_id, saved in c.execute("SELECT Simulation_ID, Parameters_JSON FROM simulation_metadata"):
                with self.subTest(run_id=run_id):
                    data, summary, _ = simulate(json.loads(saved))
                    fields = list(data.columns)
                    expected = c.execute("SELECT " + ", ".join(fields) + " FROM simulation_results WHERE Simulation_ID=? ORDER BY Days", (run_id,)).fetchall()
                    expected = np.array([[np.nan if v is None else v for v in row] for row in expected])
                    np.testing.assert_allclose(data.to_numpy(), expected, rtol=1e-12, atol=1e-12, equal_nan=True)
                    self.assertAlmostEqual(summary["Total_Feed"], data.Reactor_Volume.iloc[-1] - 1)

    def test_invalid_settings(self):
        for changes in [{"dt":0}, {"initial_volume":2}, {"initial_cells":0}, {"qG":-1}, {"unknown":1}]:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                simulate(changes)

    def test_streamlit_run_and_submit(self):
        from streamlit.testing.v1 import AppTest
        app = AppTest.from_file(str(ROOT / "streamlit_app.py")).run(timeout=30)
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(app.metric[0].value, "0.38 g/L")
        next(x for x in app.number_input if x.label == "Feed glucose (g/L)").set_value(150.0)
        app.button[0].click().run(timeout=30)
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(app.metric[0].value, "0.42 g/L")


if __name__ == "__main__":
    unittest.main()
