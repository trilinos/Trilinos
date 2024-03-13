# if (TPL_ENABLE_CUDA AND NOT Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
#   if ("${${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho}" STREQUAL "")
#     message(
#       "-- " "NOTE: Setting ${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho=OFF by default since TPL_ENABLE_CUDA='${TPL_ENABLE_CUDA}' AND Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE='${Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE}'!\n"
#       "-- NOTE: To allow the enable of ShyLU_NodeTacho, please set Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON.")
#     set(${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho OFF)
#     # NOTE: Above we set the non-cache var
#     # ${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho so that each reconfigure will
#     # show this same note.
#   elseif (${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho)
#     message(FATAL_ERROR "ERROR: ${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho=ON but TPL_ENABLE_CUDA='${TPL_ENABLE_CUDA}' AND Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE='${Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE}' which is not allowed!")
#   endif()
# endif() 

SET(KDD_INT_INT OFF)     # Current default
SET(KDD_INT_UNSIGNED OFF)  # Current default
SET(KDD_INT_LONG OFF)      # Current default
SET(KDD_INT_LONG_LONG ON)  # Current default

IF (NOT ("${Tpetra_INST_INT_INT}" STREQUAL ""))
  SET(KDD_INT_INT ${Tpetra_INST_INT_INT})
  if (${KDD_INT_INT})
    SET(KDD_INT_LONG_LONG OFF)
  ENDIF()
ENDIF()

IF(NOT ("${Tpetra_INST_INT_UNSIGNED}" STREQUAL ""))
  SET(KDD_INT_UNSIGNED ${Tpetra_INST_INT_UNSIGNED})
  if (${KDD_INT_UNSIGNED})
    SET(KDD_INT_LONG_LONG OFF)
  ENDIF()
ENDIF()

IF(NOT ("${Tpetra_INST_INT_LONG}" STREQUAL ""))
  SET(KDD_INT_LONG ${Tpetra_INST_INT_LONG})
  if (${KDD_INT_LONG})
    SET(KDD_INT_LONG_LONG OFF)
  ENDIF()
ENDIF()

IF(NOT ("${Tpetra_INST_INT_LONG_LONG}" STREQUAL ""))
  SET(KDD_INT_LONG_LONG ${Tpetra_INST_INT_LONG_LONG})
ENDIF()

# Tpetra supports only one GO type at a time, and the default is long long.
# Epetra uses GO=int. To support both libraries, Xpetra requires they use
# the same GO. So if Tpetra's GO is not INT (either the default long long,
# or set explicitly to something else), turn off
# Xpetra_Epetra support. But if Xpetra_ENABLE_Epetra is explicitly on,
# throw an error.

# Note: if >1 Tpetra GO is explicitly enabled, the logic below could turn off
# Xpetra_ENABLE_Epetra even though Tpetra_INST_INT_INT. But multiple GOs are not
# allowed, so Tpetra configuration will error out before even getting to Xpetra.
IF(KDD_INT_LONG OR KDD_INT_LONG_LONG OR KDD_INT_UNSIGNED)
   # Tpetra will not using GO=int.
   # Is Xpetra_ENABLE_Tpetra explicitly OFF?
   SET(BMK_EXPLICIT_TPETRA_OFF OFF)
   IF(NOT "${Trilinos_ENABLE_Tpetra}" STREQUAL "" AND NOT ${Trilinos_ENABLE_Tpetra})
     SET(BMK_EXPLICIT_TPETRA_OFF ON)
   ENDIF()
   SET(BMK_EXPLICIT_XT_OFF OFF)
   IF(NOT "${Xpetra_ENABLE_Tpetra}" STREQUAL "" AND NOT ${Xpetra_ENABLE_Tpetra})
     SET(BMK_EXPLICIT_XT_OFF ON)
   ENDIF()
   # Assuming that Tpetra is always on by default (which it is, as a PT package)
   # If Xpetra is not enabled for any reason, nothing here will have an effect.
   #
   # Several cases to consider:
   #  -If Trilinos_ENABLE_Tpetra or Xpetra_ENABLE_Tpetra are explicitly OFF, nothing to do
   #  -If Xpetra_ENABLE_Epetra is explicitly set either way,
   #    let Xpetra error out (if ON) or be fine (if OFF) later.
   #  -Otherwise, turn off Xpetra_ENABLE_Epetra.
   IF("${Xpetra_ENABLE_Epetra}" STREQUAL "" AND NOT ${BMK_EXPLICIT_TPETRA_OFF} AND NOT ${BMK_EXPLICIT_XT_OFF})
     SET(Xpetra_ENABLE_Epetra OFF)
     SET(Xpetra_ENABLE_EpetraExt OFF)
     SET(MueLu_ENABLE_Epetra OFF)
   ENDIF()
ENDIF()

# Special logic to disable SEACAS subpackages depending on Fortran enabled or not

if (SEACAS_SOURCE_DIR)
  include("${SEACAS_SOURCE_DIR}/cmake/SeacasDisableSubpackagesDependingOnFortran.cmake")
  seacas_disable_subpackages_depending_on_fortran()
endif()

# Enable Build Status package and enable reporting the reporting test

set(bulidStats "${${PROJECT_NAME}_ENABLE_BUILD_STATS}")
if (
    bulidStats
    AND
    ("${${PROJECT_NAME}_ENABLE_TrilinosBuildStats}" STREQUAL "")
    )
  message("-- " "Setting ${PROJECT_NAME}_ENABLE_TrilinosBuildStats=ON"
    " by default because"
    " ${PROJECT_NAME}_ENABLE_BUILD_STATS=${bulidStats}"
    )
  set(${PROJECT_NAME}_ENABLE_TrilinosBuildStats ON)
endif()
