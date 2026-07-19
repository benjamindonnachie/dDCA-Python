import unittest
import warnings

import numpy as np

from ddca import __version__, dDCA


class DDCAImplementationTests(unittest.TestCase):
    def test_version(self):
        self.assertEqual(__version__, "9.0.0")

    def test_first_antigen_matches_c_cell_selection(self):
        engine = dDCA(cells=3, antigenCount=3, maxMigration=10)
        engine.doAntigen(2)

        self.assertEqual(engine.DCs[0].antigen, [])
        self.assertEqual(engine.DCs[1].antigen, [2])
        self.assertEqual(engine.DCs[2].antigen, [])

    def test_single_cell_population(self):
        engine = dDCA(cells=1, antigenCount=1, maxMigration=4)
        engine.doAntigen(0)
        engine.doSignals(danger=3.0, safe=1.0)
        engine.results()

        np.testing.assert_array_equal(engine.m, [1])
        np.testing.assert_array_equal(engine.s, [0])
        np.testing.assert_allclose(engine.mcav, [1.0])
        np.testing.assert_allclose(engine.ka, [1.0])
        self.assertEqual(engine.Tk, 1.0)
        self.assertEqual(engine.Tmcav, 3.0)

    def test_reset_preserves_floating_migration_interval(self):
        engine = dDCA(cells=100, antigenCount=1, maxMigration=100)
        cell = engine.DCs[1]
        expected = 100 / 99

        self.assertAlmostEqual(cell.lifespan, expected)
        cell.lifespan = -1
        cell.reset()

        self.assertAlmostEqual(cell.lifespan, expected)
        self.assertNotEqual(cell.lifespan, 1)

    def test_safe_context_outputs(self):
        engine = dDCA(cells=1, antigenCount=1, maxMigration=2)
        engine.doAntigen(0)
        engine.doSignals(danger=1.0, safe=1.0)
        engine.results()

        np.testing.assert_array_equal(engine.m, [0])
        np.testing.assert_array_equal(engine.s, [1])
        np.testing.assert_allclose(engine.mcav, [0.0])
        np.testing.assert_allclose(engine.ka, [-1.0])
        self.assertEqual(engine.Tk, -1.0)
        self.assertEqual(engine.Tmcav, 1.0)

    def test_logging_sorts_and_clears_dynamic_antigen_list(self):
        engine = dDCA(cells=1, antigenCount=3, maxMigration=10)
        cell = engine.DCs[0]
        cell.antigen = [2, 0, 1]
        cell.k = 2.0

        engine.logAntigen(cell)

        np.testing.assert_array_equal(engine.m, [1, 1, 1])
        np.testing.assert_allclose(engine.k, [2.0, 2.0, 2.0])
        self.assertEqual(cell.antigen, [])

    def test_unobserved_antigen_outputs_nan(self):
        engine = dDCA(cells=1, antigenCount=2, maxMigration=2)
        engine.doAntigen(0)
        engine.doSignals(danger=1.0, safe=1.0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            engine.results()

        self.assertTrue(np.isnan(engine.mcav[1]))
        self.assertTrue(np.isnan(engine.ka[1]))


if __name__ == "__main__":
    unittest.main()
