macro(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  assert_defined(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    message(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_WrapExternal=OFF"
      " because ${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES='${${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES}'!"
      "\n***\n"
      )
    set(${PROJECT_NAME}_ENABLE_WrapExternal OFF)
  endif()

  if ("${Python3_EXECUTABLE}" STREQUAL "")
    message(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_WrapExternal=OFF"
      " because Python3_EXECUTABLE=''!"
      "\n***\n"
      )
    set(${PROJECT_NAME}_ENABLE_WrapExternal OFF)
  endif()

  assert_defined(${PROJECT_NAME}_ENABLE_Fortran)

  if (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    message(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_MixedLang=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran='${${PROJECT_NAME}_ENABLE_Fortran}'!"
      "\n***\n"
      )
    set(${PROJECT_NAME}_ENABLE_MixedLang OFF)
  endif()

endmacro()
