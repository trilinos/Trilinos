if (TPL_ENABLE_CUDA AND NOT Kokkos_ENABLE_Cuda_Relocatable_Device_Code)
  if ("${${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho}" STREQUAL "")
    message(
      "-- " "NOTE: Setting ${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho=OFF by default since TPL_ENABLE_CUDA='${TPL_ENABLE_CUDA}' AND Kokkos_ENABLE_Cuda_Relocatable_Device_Code='${Kokkos_ENABLE_Cuda_Relocatable_Device_Code}'!\n"
      "-- NOTE: To allow the enable of ShyLU_NodeTacho, please set Kokkos_ENABLE_Cuda_Relocatable_Device_Code=ON.")
    set(${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho OFF)
    # NOTE: Above we set the non-cache var
    # ${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho so that each reconfigure will
    # show this same note.
  elseif (${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho)
    message(FATAL_ERROR "ERROR: ${PROJECT_NAME}_ENABLE_ShyLU_NodeTacho=ON but TPL_ENABLE_CUDA='${TPL_ENABLE_CUDA}' AND Kokkos_ENABLE_Cuda_Relocatable_Device_Code='${Kokkos_ENABLE_Cuda_Relocatable_Device_Code}' which is not allowed!")
  endif()
endif() 

#########################################################################
# STKBalance does not work with GO=INT or GO=UNSIGNED
# Note and disable it if it was set as a dependence of STK
# Error out if it was explicitly requested by user

IF ((NOT ("${Tpetra_ENABLE_DEPRECATED_CODE}" STREQUAL "")) # Remove this test
    AND                                                    # when removing
    (NOT Tpetra_ENABLE_DEPRECATED_CODE))                   # Tpetra DEPRECATED
                                                           # code 
  SET(KDD_ENABLE_DEPRECATED OFF)  # Remove this line when DEPRECATED is removed
  SET(KDD_INT_INT OFF)     # Current default
ELSE() # Remove this "else" section when Tpetra DEPRECATED code is removed
  SET(KDD_ENABLE_DEPRECATED ON)   # Remove this line when DEPRECATED is removed
  SET(KDD_INT_INT ON)      # Deprecated default
ENDIF()

SET(KDD_INT_UNSIGNED OFF)  # Current default
SET(KDD_INT_LONG OFF)      # Current default
SET(KDD_INT_LONG_LONG ON)  # Current default

IF (NOT ("${Tpetra_INST_INT_INT}" STREQUAL ""))
  SET(KDD_INT_INT ${Tpetra_INST_INT_INT})
  if (${KDD_INT_INT}                      # Keep this test but
      AND (NOT ${KDD_ENABLE_DEPRECATED})) # remove this test when DEPRECATED
                                          # is removed
    SET(KDD_INT_LONG_LONG OFF)
  ENDIF()
ENDIF()

IF(NOT ("${Tpetra_INST_INT_UNSIGNED}" STREQUAL ""))
  SET(KDD_INT_UNSIGNED ${Tpetra_INST_INT_UNSIGNED})
  if (${KDD_INT_UNSIGNED}                 # Keep this test but
      AND (NOT ${KDD_ENABLE_DEPRECATED})) # remove this test when DEPRECATED
                                          # is removed
    SET(KDD_INT_LONG_LONG OFF)
  ENDIF()
ENDIF()

IF(NOT ("${Tpetra_INST_INT_LONG}" STREQUAL ""))
  SET(KDD_INT_LONG ${Tpetra_INST_INT_LONG})
  if (${KDD_INT_LONG}                     # Keep this test but
      AND (NOT ${KDD_ENABLE_DEPRECATED})) # remove this test when DEPRECATED
                                          # is removed
    SET(KDD_INT_LONG_LONG OFF)
  ENDIF()
ENDIF()

IF(NOT ("${Tpetra_INST_INT_LONG_LONG}" STREQUAL ""))
  SET(KDD_INT_LONG_LONG ${Tpetra_INST_INT_LONG_LONG})
ENDIF()

IF ((NOT ${KDD_INT_LONG}) AND (NOT ${KDD_INT_LONG_LONG}))
  IF ("${${PROJECT_NAME}_ENABLE_STKBalance}" STREQUAL "")
    # STKBalance may be enabled but only implicitly (as a dependence of STK);
    # give a message but turn off STKBalance support
    MESSAGE("NOTE:  int global indices are enabled in Trilinos. "
            "Because STKBalance requires long or long long "
            "global indices, STKBalance will be disabled.  "
            "To make this warning go away, do not request "
            "int global indices in Trilinos (that is,  do not "
            "set Tpetra_INST_INT_INT=ON or "
            "Tpetra_INST_INT_UNSIGNED=ON)." )
    SET(${PROJECT_NAME}_ENABLE_STKBalance OFF)
  ELSEIF (${PROJECT_NAME}_ENABLE_STKBalance)
    # STKBalance was explicitly enabled by the user, so error out
    MESSAGE(FATAL_ERROR 
            "STKBalance requires long or long long global indices, "
            "but Trilinos is using int indices "
            "(likely via Tpetra_INST_INT_INT or Tpetra_INST_INT_UNSIGNED).  "
            "Disable STKBalance or specify Tpetra_INST_INT_LONG_LONG.")
  ENDIF()
ENDIF()
  
