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