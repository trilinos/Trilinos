import os
import sys


thisDir = os.path.dirname(os.path.abspath(__file__))
tribitsDir = os.path.abspath(os.path.join(thisDir, "..", "..", "..", "cmake", "tribits"))
print("tribitsDir = '"+tribitsDir+"'")
pythonUtilsDir = os.path.join(tribitsDir, "python_utils")
ciSupportDir = os.path.join(tribitsDir, "ci_support")

sys.path = [ciSupportDir] + [pythonUtilsDir] + sys.path
