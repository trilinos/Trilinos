# PyTrilinos2: Python 3 bindings for Trilinos
# Requires: Python 3, python3-dev, pybind11
# Set PYTHON_ROOT environment variable to Python installation path

set(PYTHON3_EXECUTABLE "$ENV{PYTHON_ROOT}/bin/python3" CACHE FILEPATH "")

SET_BOOL_CACHE_VAR(Trilinos_ENABLE_PyTrilinos2 ON)
SET_BOOL_CACHE_VAR(PyTrilinos2_ENABLE_TESTS ON)
SET_BOOL_CACHE_VAR(PyTrilinos2_ENABLE_EXAMPLES ON)

# Uncomment for verbose binder output during build
# SET_BOOL_CACHE_VAR(PyTrilinos2_BINDER_VERBOSE ON)
