# Default Trilinos build configuration
# CPU-only build with MPI, shared libraries, debug-friendly settings

SET(SHARED_SNIPPETS_DIR "${CMAKE_CURRENT_LIST_DIR}/../shared")

INCLUDE("${SHARED_SNIPPETS_DIR}/utils.cmake")


##################################################
# Third-party libraries (TPLs)

# Ignore warnings from TPL headers (#1458)
SET_BOOL_CACHE_VAR(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS TRUE)

SET_BOOL_CACHE_VAR(TPL_ENABLE_MPI ON)

# Linear Algebra Libraries
SET_BOOL_CACHE_VAR(TPL_ENABLE_BLAS ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_LAPACK ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_Boost ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_BoostLib OFF)

# I/O Libraries
SET_BOOL_CACHE_VAR(TPL_ENABLE_Zlib ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_HDF5 ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_Netcdf ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_Pnetcdf ON)

# Testing Framework
SET_BOOL_CACHE_VAR(TPL_ENABLE_gtest ON)

# YAML file parser
SET_BOOL_CACHE_VAR(TPL_ENABLE_yamlcpp OFF)

# Disable unused TPLs
SET_BOOL_CACHE_VAR(TPL_ENABLE_BinUtils OFF)
SET_BOOL_CACHE_VAR(TPL_ENABLE_SuperLU OFF)
SET_BOOL_CACHE_VAR(TPL_ENABLE_Matio OFF)
SET_BOOL_CACHE_VAR(TPL_ENABLE_X11 OFF)
# ParMETIS: uncomment if needed for graph partitioning
# SET_BOOL_CACHE_VAR(TPL_ENABLE_ParMETIS ON)

# MPI: use compiler wrappers (mpicxx, mpicc, mpifort)
SET_BOOL_CACHE_VAR(MPI_USE_COMPILER_WRAPPERS ON)
INCLUDE("${SHARED_SNIPPETS_DIR}/mpi.cmake")


##################################################
# Trilinos package enable/disable configuration

# Enable tests and examples
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_EXAMPLES ON)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_TESTS ON)
# Instead of globally turning them on for all packages, the above variables
# can be set to OFF and individual packages can be turned ON, e.g.:
# SET_BOOL_CACHE_VAR(Ifpack2_ENABLE_TESTS ON)
# SET_BOOL_CACHE_VAR(Ifpack2_ENABLE_EXAMPLES ON)

# Kokkos profiling support and demangling for stack traces
SET_BOOL_CACHE_VAR(Teuchos_KOKKOS_PROFILING ON)
SET_BOOL_CACHE_VAR(Teuchos_ENABLE_GCC_DEMANGLE ON)
# Globally reduce output from unit tests (less noise in test runs)
SET_BOOL_CACHE_VAR(Teuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS ON)

# Linear Solvers: Belos (Iterative solvers: CG, GMRES, BiCGStab, etc.)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_Belos ON)

# Linear Solver Integration: Stratimikos (Thyra-based solver stacks)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_Stratimikos ON)

# Multi-Level Algebraic Multigrid: MueLu (algebraic multigrid preconditioner)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_MueLu ON)

# Block Preconditioning: Teko (block preconditioning framework)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_Teko ON)

# PDE Discretization: Panzer (finite element assembly framework)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_Panzer ON)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_PanzerMiniEM ON)

# Optimization: ROL (Rapid Optimization Library)
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_ROL OFF)

# TrilinosCouplings: example drivers
SET_BOOL_CACHE_VAR(Trilinos_ENABLE_TrilinosCouplings OFF)


##################################################
#

# Build shared libraries
SET_BOOL_CACHE_VAR(BUILD_SHARED_LIBS ON)

SET_CACHE_VAR(CMAKE_Fortran_FLAGS "-fPIE" STRING)

# Use ninja
INCLUDE("${SHARED_SNIPPETS_DIR}/ninja.cmake")

# Complex number support
# INCLUDE("${SHARED_SNIPPETS_DIR}/complex.cmake")

# Hide deprecated code
# INCLUDE("${SHARED_SNIPPETS_DIR}/deprecated.cmake")

# Python bindings
# INCLUDE("${SHARED_SNIPPETS_DIR}/pytrilinos2.cmake")

# Debug build settings (slower, more checks)
# INCLUDE("${SHARED_SNIPPETS_DIR}/debug.cmake")

# SET(Trilinos_WARNINGS_MODE "WARN" CACHE STRING "Trilinos warnings mode")

