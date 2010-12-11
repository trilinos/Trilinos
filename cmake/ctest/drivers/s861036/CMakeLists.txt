INCLUDE("../../../TrilinosVersion.cmake")

TRILINOS_DRIVER_SETUP()

#
# Zoltan only, with Parmetis and Scotch
#  Zoltan uses 64-bit global IDs (unless compiled on 32-bit machine)
#  Parmetis and Scotch use 32-bit IDs
#
TRILINOS_DRIVER_ADD_DASHBOARD(
    ZOLTAN_ONLY_64
    ctest_zoltan_only_ids64.cmake
    CTEST_INSTALLER_TYPE dev
    TIMEOUT_MINUTES 120
)

#
# Zoltan only, with Scotch
#  Zoltan uses 64-bit global IDs (unless compiled on 32-bit machine)
#  Scotch uses 64-bit IDs
#
TRILINOS_DRIVER_ADD_DASHBOARD(
    ZOLTAN_ONLY_64_SCOTCH_64
    ctest_zoltan_only_ids64_scotch_64.cmake
    CTEST_INSTALLER_TYPE dev
    TIMEOUT_MINUTES 120
)

#
# Zoltan and packages that use Zoltan, global IDs are unsigned int
#
TRILINOS_DRIVER_ADD_DASHBOARD(
    ZOLTAN_PLUS_USERS
    ctest_zoltan_plus.cmake
    CTEST_INSTALLER_TYPE dev
    TIMEOUT_MINUTES 240
)

TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
