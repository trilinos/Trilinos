import unittest

from pyrol.pyrol.Teuchos import ParameterList


class TestParameterList(unittest.TestCase):

    def test_supported_values_round_trip(self):
        p = ParameterList()
        p["bool"] = True
        p["int"] = 3
        p["float"] = 2.5
        p["str"] = "Trust Region"
        self.assertIs(p["bool"], True)
        self.assertEqual(p["int"], 3)
        self.assertEqual(p["float"], 2.5)
        self.assertEqual(p["str"], "Trust Region")

    def test_non_ascii_string(self):
        p = ParameterList()
        p["name"] = "café"
        self.assertEqual(p["name"], "café")

    def test_missing_entry_raises_key_error(self):
        p = ParameterList()
        with self.assertRaises(KeyError):
            p["missing"]
        with self.assertRaises(KeyError):
            p.get("missing")

    def test_unsupported_value_raises_type_error(self):
        p = ParameterList()
        for value in (None, [1, 2], object()):
            with self.subTest(value=value):
                with self.assertRaises(TypeError):
                    p["x"] = value
                with self.assertRaises(TypeError):
                    p.set("x", value)
        self.assertFalse(p.isParameter("x"))

    def test_sublists(self):
        p = ParameterList()
        p["Step"] = ParameterList()
        p.sublist("Step").set("Type", "Line Search")
        self.assertEqual(p["Step"]["Type"], "Line Search")
        with self.assertRaises(KeyError):
            p.sublist("missing")

    def test_contains_and_len(self):
        p = ParameterList()
        p["Step"] = ParameterList()
        p["Iteration Limit"] = 10
        self.assertIn("Step", p)
        self.assertIn("Iteration Limit", p)
        self.assertNotIn("missing", p)
        self.assertEqual(len(p), 2)


if __name__ == "__main__":
    unittest.main()
