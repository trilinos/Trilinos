# Importing mpi4py.MPI initializes MPI and registers MPI_Finalize to run at exit.  pykrino.finalize is registered
# after it, so it runs first (atexit is last in, first out), while MPI is still initialized.
import atexit as _atexit
from mpi4py import MPI as _MPI
from ._pykrino import *
from ._pykrino import finalize as _finalize
_atexit.register(_finalize)
from .compute_training_data import compute_training_data
from .mesh_adapt import adapt_mesh
from .mesh_adapt import parse_mesh_adapt_command_line
from .printP0 import printP0

__version__ = '0.0.1'
def version():
    return 'PyKrino version: ' + __version__

