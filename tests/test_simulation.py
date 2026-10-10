import sys
from contextlib import closing
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
import sqlite3
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from simulation_engine import Parameters, derivatives, simulate, export_batch, scenario_parameters


class SimulationTests(unittest.TestCase):
    def test_mass_balances(self):
        p = Parameters()
        y = np.array([8., .2, 1., .1, 2., .5, 1.2])
        flow = .1
        dy = derivatives(0., y, p, flow)
        mu = p.mu_max * 8 / 10 * 4 / 5 * 2 / 2.1
        self.assertAlmostEqual((dy[4] + dy[5]) * y[6] + (y[4] + y[5]) * flow,
                               mu * y[4] * y[6])
        uptake = p.qG * 8 / 8.1 * y[4] * .001
        self.assertAlmostEqual(dy[0] * y[6] + y[0] * flow,
                               flow * p.feed_glucose - uptake * y[6])
        health = mu / p.mu_max
        self.assertAlmostEqual(dy[1] * y[6] + y[1] * flow,
                               p.qp * health * y[4] * .001 * y[6])

    def test_depleted_substrates(self):
        dy = derivatives(0., [0., 0., 0., 0., 1., 0., 1.], Parameters(), 0.)
        np.testing.assert_allclose(dy[:4], 0.)
        self.assertLess(dy[4], 0.)
        result, _ = simulate(replace(Parameters(), initial_glucose=.001, feed_rate_on=0.))
        self.assertGreaterEqual(result.Glucose.min(), 0.)

    def test_pure_dilution(self):
        p = replace(Parameters(), mu_max=0., kd_min=0., kd_max=0.,
                    qG=0., qp=0., qA=0., qL_production=0., qL_consumption=0.)
        result, summary = simulate(p)
        np.testing.assert_allclose(result.Viable_Cells * result.Reactor_Volume, 1., atol=1e-7)
        self.assertAlmostEqual(summary["Total_Feed"], summary["Final_Volume"] - 1.)

    def test_scenarios_and_alignment(self):
        for name, p in scenario_parameters().items():
            result, summary = simulate(p, name)
            self.assertTrue(np.isfinite(result.to_numpy()).all())
            self.assertTrue((result.to_numpy() >= -1e-8).all())
            self.assertLessEqual(result.Reactor_Volume.max(), p.max_feed_volume + 1e-8)
            self.assertTrue(result.Viability.between(0, 100).all())
            expected_mu = p.mu_max * result.Glucose / (p.Ks + result.Glucose) * p.Kl / (p.Kl + result.Lactate) * p.Ka / (p.Ka + result.Ammonia)
            np.testing.assert_allclose(result.Growth_Rate, expected_mu)

    def test_grid_convergence_and_final_command(self):
        p = replace(Parameters(), duration=2.3, dt=.7, glucose_target=100., max_feed_volume=2.)
        coarse, _ = simulate(p)
        fine, _ = simulate(p, rtol=1e-9, max_step=.01)
        self.assertEqual(coarse.Days.iloc[-1], 2.3)
        self.assertGreater(coarse.Feed_Rate.iloc[-1], 0.)
        np.testing.assert_allclose(coarse.iloc[-1], fine.iloc[-1], rtol=2e-5, atol=1e-7)

    def test_export_schema(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            directory, result = export_batch({"Baseline": Parameters()}, tmp)
            with closing(sqlite3.connect(directory / "cho_digital_twin.db")) as c:
                self.assertEqual(c.execute("SELECT COUNT(*) FROM simulation_results").fetchone()[0], 15)
                self.assertEqual(c.execute("SELECT Scenario FROM simulation_summary").fetchone()[0], "Baseline")
                self.assertEqual(c.execute("SELECT COUNT(*) FROM simulation_runs").fetchone()[0], 1)

    def test_starvation_death_solution(self):
        p = replace(Parameters(), initial_glucose=0., feed_rate_on=0.)
        result, _ = simulate(p)
        np.testing.assert_allclose(result.Viable_Cells,
                                   np.exp(-p.kd_max * result.Days), rtol=2e-6)
        np.testing.assert_allclose(result.Total_Cells, 1., atol=1e-8)

    def test_partial_feed_cap(self):
        p = replace(Parameters(), initial_volume=1.46, glucose_target=100.)
        result, summary = simulate(p)
        self.assertAlmostEqual(result.Feed_Rate.iloc[0], .04)
        self.assertAlmostEqual(summary["Total_Feed"], .04)
        self.assertAlmostEqual(result.Reactor_Volume.iloc[-1], 1.5)


if __name__ == "__main__":
    unittest.main()
