# Debug build configuration
# Balances optimization with debug symbols for profiling

# RELWITHDEBINFO: optimized (-O2) with debug symbols (-g)
# Use DEBUG for maximum debug info (-O0, -g3) - slower builds
set(CMAKE_BUILD_TYPE "RELWITHDEBINFO" CACHE STRING "")

# Kokkos runtime bounds checking (catches array access violations)
SET_BOOL_CACHE_VAR(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK ON)

# Communication timers for profiling MPI operations
SET_BOOL_CACHE_VAR(HAVE_TEUCHOS_COMM_TIMERS ON)

# Enable timers in solver packages for performance analysis
SET_BOOL_CACHE_VAR(Belos_Tpetra_Timers ON)
SET_BOOL_CACHE_VAR(Anasazi_Tpetra_Timers ON)
