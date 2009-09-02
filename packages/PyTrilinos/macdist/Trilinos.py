# This script requires installed PyTrilinos (make install), so that
# all libraries are located in the same place.
# It should be used with non-MPI installations.
# The tarball is tight to specific major versions of MAC OS X.
#
# TODO:
# - support for NOX has to be done

# It seems that these modules have to be imported here, otherwise they
# will not be included in the tarball
import sys, os
import Numeric

# All the PyTrilinos import here
from PyTrilinos import Epetra, EpetraExt, Galeri, Amesos, IFPACK, AztecOO, ML

# Simple logic to run the input script.
args = sys.argv[1:]
if len(args) == 0:
  print "Usage: <path>/Trilinos script-name.py"
else:
  execfile(args[0])
  # FIXME: I am not sure about command line arguments
