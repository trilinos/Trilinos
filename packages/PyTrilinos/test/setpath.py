"""Set the python search path to include the library directory created by
distutils."""

# System includes
import commands
import os
import sys

# Get the path to the build directory
libDir = os.path.join("..", commands.getoutput("../pyLocate --build"))

# Insert the library directory name at the beginning of
# the python search path
if libDir:
    sys.path.insert(0,libDir)
