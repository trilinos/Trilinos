cmake_minimum_required(VERSION 3.23)

IF(COMMAND TRIBITS_PACKAGE_DECL)
    message("Building Krino as a Trilinos package")
    TRIBITS_PACKAGE(Krino)

    TRIBITS_ADD_DEBUG_OPTION()
    TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

    if (${${PROJECT_NAME}_ENABLE_Krino})
      add_subdirectory(krino)
      add_subdirectory(delete_small_elements)
    endif()

    TRIBITS_PACKAGE_POSTPROCESS()
ELSE()
    MESSAGE("*** Building Krino as a stand-alone cmake package. ***")

    include(SierraCmakeCode.cmake)
ENDIF()
