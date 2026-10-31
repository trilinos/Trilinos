cmake_minimum_required(VERSION 3.23)

IF(COMMAND TRIBITS_PACKAGE_DECL)
    message("Building Krino as a Trilinos package")
    INCLUDE(TribitsAddOptionAndDefine)
    
    TRIBITS_PACKAGE(Krino)

    TRIBITS_ADD_DEBUG_OPTION()
    TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()
    
    TRIBITS_ADD_OPTION_AND_DEFINE(Krino_ENABLE_PyKrino
      ENABLE_PYBIND11_PYKRINO
      "Build Krino with PyKrino interface."
      OFF)

    if (${${PROJECT_NAME}_ENABLE_Krino})
      add_subdirectory(krino)
      add_subdirectory(delete_small_elements)
      IF( Krino_ENABLE_PyKrino ) 
        message("Enabling PyKrino")
        ADD_SUBDIRECTORY(pykrino)
      ENDIF()
    endif()

    TRIBITS_PACKAGE_POSTPROCESS()
ELSE()
    MESSAGE("*** Building Krino as a stand-alone cmake package. ***")
    project(krino LANGUAGES C CXX)

    include(SierraCmakeCode.cmake)
ENDIF()
