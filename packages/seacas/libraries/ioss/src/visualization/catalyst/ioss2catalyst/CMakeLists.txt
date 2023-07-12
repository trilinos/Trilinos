
SET(CATALYST_PLUGIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
SET(CATALYST_PLUGIN_BUILD_DIR ${PROJECT_BINARY_DIR})
CONFIGURE_FILE(CatalystPluginPaths.h.in CatalystPluginPaths.h @ONLY)

ADD_LIBRARY(catalystiossapp SHARED
    IossRegionReport.cxx
    IossApplication.cxx)

TARGET_INCLUDE_DIRECTORIES(catalystiossapp
    PUBLIC
        ${SEACAS_INCLUDE_DIRS}
        ${CMAKE_CURRENT_BINARY_DIR})

TARGET_LINK_LIBRARIES(catalystiossapp
    PUBLIC
        -L${SEACAS_INSTALL_DIR}/lib
        SEACASIoss::all_libs
        MPI::MPI_C
    PRIVATE
        ${Python3_LIBRARY})

ADD_EXECUTABLE(ioss2catalyst main.cxx)
TARGET_LINK_LIBRARIES(ioss2catalyst
    PRIVATE
        catalystiossapp
        MPI::MPI_C)

INSTALL(TARGETS ioss2catalyst RUNTIME DESTINATION bin)
INSTALL(TARGETS catalystiossapp LIBRARY DESTINATION lib)
