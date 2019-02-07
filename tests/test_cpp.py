import unittest
import subprocess
import os

class MainTest(unittest.TestCase):
    def test_cpp(self):
	print("\n\n Testing the C++ functions...")
	subprocess.check_call(os.path.join(os.path.dirname(os.path.relpath(__file__)),'bin','integer_polyomino_cpp'))
	print()

if __name__ == "__main__":
    unittest.main()
