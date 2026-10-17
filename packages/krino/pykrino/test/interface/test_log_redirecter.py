import unittest
import pykrino


class TestLogRedirecter(unittest.TestCase):
    """Unit tests for the pykrino LogRedirecter Python bindings."""

    def test_log_redirecter_import(self):
        """Test that LogRedirecter class can be imported"""
        self.assertTrue(hasattr(pykrino, 'LogRedirecter'))

    def test_log_redirecter_construction(self):
        """Test LogRedirecter construction"""
        redirecter = pykrino.LogRedirecter()
        self.assertIsNotNone(redirecter)
        self.assertIsInstance(redirecter, pykrino.LogRedirecter)

    def test_log_redirecter_clear(self):
        """Test LogRedirecter clear method"""
        redirecter = pykrino.LogRedirecter()
        # clear should not raise an exception
        redirecter.clear()

    def test_log_redirecter_get_log(self):
        """Test LogRedirecter get_log method"""
        redirecter = pykrino.LogRedirecter()
        # get_log should return a string
        log_output = redirecter.get_log()
        self.assertIsInstance(log_output, str)

    def test_log_redirecter_lifecycle(self):
        """Test LogRedirecter lifecycle - construction and destruction"""
        redirecter = pykrino.LogRedirecter()
        redirecter.clear()
        log_output = redirecter.get_log()
        self.assertIsInstance(log_output, str)
        # When redirecter goes out of scope, destructor should restore original buffer


if __name__ == "__main__":
    unittest.main()
