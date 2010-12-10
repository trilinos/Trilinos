INCLUDE("../../../TrilinosVersion.cmake")

TRILINOS_DRIVER_SETUP()

#
# Zoltan only, global IDs are long long
#
TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_ZOLTAN_ONLY_64
    ctest_zoltan_only_ids64.cmake
    CTEST_INSTALLER_TYPE dev
    TIMEOUT_MINUTES 300
)

#
# Zoltan and packages that use Zoltan, global IDs are unsigned int
#
TRILINOS_DRIVER_ADD_DASHBOARD(
    MPI_ZOLTAN_PLUS_ZOLTAN_USERS
    ctest_zoltan_plus.cmake
    CTEST_INSTALLER_TYPE dev
    TIMEOUT_MINUTES 300
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
