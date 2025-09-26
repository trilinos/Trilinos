cmake_minimum_required(VERSION 3.23)
project(krino LANGUAGES C CXX Fortran)
find_package(sierra-cmake-utils REQUIRED)
add_parser_commands(TARGET krino_commands  
	XML_FILES
		${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra/xml/Akri_Levelset.xml)
install(FILES ${CMAKE_BINARY_DIR}/krino_commands.xmldb DESTINATION xml)

add_library(krino_diagwriter)
FILE(GLOB krino_diagwriter_headers CONFIGURE_DEPENDS krino/diagwriter/*.hpp)
FILE(GLOB krino_diagwriter_sources CONFIGURE_DEPENDS krino/diagwriter/*.cpp)
target_sources(krino_diagwriter PRIVATE ${krino_diagwriter_sources})
find_package(stk REQUIRED)
target_link_libraries(krino_diagwriter PUBLIC stk::stk_util_diag)
target_include_directories(krino_diagwriter PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/diagwriter>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/diagwriter>)
target_sources(krino_diagwriter PUBLIC
    FILE_SET krino_diagwriter_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_diagwriter_headers})
target_compile_definitions(krino_diagwriter PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_diagwriter PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_diagwriter PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_diagwriter PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_diagwriter
    EXPORT krinoTargets
    FILE_SET krino_diagwriter_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_geometry)
FILE(GLOB krino_geometry_headers CONFIGURE_DEPENDS krino/geometry/*.hpp)
FILE(GLOB krino_geometry_sources CONFIGURE_DEPENDS krino/geometry/*.cpp)
target_sources(krino_geometry PRIVATE ${krino_geometry_sources})
find_package(stk REQUIRED)
find_package(Sacado REQUIRED)
target_link_libraries(krino_geometry PUBLIC stk::stk_math)
target_link_libraries(krino_geometry PUBLIC Sacado::all_libs)
target_link_libraries(krino_geometry PUBLIC stk::stk_util_parallel)
target_include_directories(krino_geometry PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/geometry>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/geometry>)
target_sources(krino_geometry PUBLIC
    FILE_SET krino_geometry_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_geometry_headers})
target_compile_definitions(krino_geometry PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_geometry PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_geometry PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_geometry PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_geometry
    EXPORT krinoTargets
    FILE_SET krino_geometry_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_quality_metric)
FILE(GLOB krino_quality_metric_headers CONFIGURE_DEPENDS krino/quality_metric/*.hpp)
FILE(GLOB krino_quality_metric_sources CONFIGURE_DEPENDS krino/quality_metric/*.cpp)
target_sources(krino_quality_metric PRIVATE ${krino_quality_metric_sources})
find_package(stk REQUIRED)
target_link_libraries(krino_quality_metric PUBLIC stk::stk_math)
target_link_libraries(krino_quality_metric PUBLIC stk::stk_topology)
target_include_directories(krino_quality_metric PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/quality_metric>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/quality_metric>)
target_sources(krino_quality_metric PUBLIC
    FILE_SET krino_quality_metric_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_quality_metric_headers})
target_compile_definitions(krino_quality_metric PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_quality_metric PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_quality_metric PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_quality_metric PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_quality_metric
    EXPORT krinoTargets
    FILE_SET krino_quality_metric_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_quality_metric_sens)
FILE(GLOB krino_quality_metric_sens_headers CONFIGURE_DEPENDS krino/quality_metric_sens/*.hpp)
FILE(GLOB krino_quality_metric_sens_sources CONFIGURE_DEPENDS krino/quality_metric_sens/*.cpp)
target_sources(krino_quality_metric_sens PRIVATE ${krino_quality_metric_sens_sources})
find_package(stk REQUIRED)
find_package(Sacado REQUIRED)
target_link_libraries(krino_quality_metric_sens PUBLIC MPI::MPI_C)
target_link_libraries(krino_quality_metric_sens PUBLIC stk::stk_math)
target_link_libraries(krino_quality_metric_sens PUBLIC Sacado::all_libs)
target_link_libraries(krino_quality_metric_sens PUBLIC krino_quality_metric)
target_include_directories(krino_quality_metric_sens PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/quality_metric_sens>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/quality_metric_sens>)
target_sources(krino_quality_metric_sens PUBLIC
    FILE_SET krino_quality_metric_sens_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_quality_metric_sens_headers})
target_compile_definitions(krino_quality_metric_sens PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_quality_metric_sens PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_quality_metric_sens PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_quality_metric_sens PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_quality_metric_sens
    EXPORT krinoTargets
    FILE_SET krino_quality_metric_sens_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_surface)
FILE(GLOB krino_surface_headers CONFIGURE_DEPENDS krino/surface/*.hpp)
FILE(GLOB krino_surface_sources CONFIGURE_DEPENDS krino/surface/*.cpp)
target_sources(krino_surface PRIVATE ${krino_surface_sources})
target_link_libraries(krino_surface PUBLIC krino_diagwriter)
target_link_libraries(krino_surface PUBLIC krino_geometry)
target_link_libraries(krino_surface PUBLIC krino_math_utils)
target_link_libraries(krino_surface PUBLIC krino_quality_metric)
find_package(stk REQUIRED)
target_link_libraries(krino_surface PUBLIC stk::stk_expreval)
target_include_directories(krino_surface PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/surface>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/surface>)
target_sources(krino_surface PUBLIC
    FILE_SET krino_surface_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_surface_headers})
target_compile_definitions(krino_surface PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_surface PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_surface PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_surface PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_surface
    EXPORT krinoTargets
    FILE_SET krino_surface_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_mesh_surface)
FILE(GLOB krino_mesh_surface_headers CONFIGURE_DEPENDS krino/mesh_surface/*.hpp)
FILE(GLOB krino_mesh_surface_sources CONFIGURE_DEPENDS krino/mesh_surface/*.cpp)
target_sources(krino_mesh_surface PRIVATE ${krino_mesh_surface_sources})
target_link_libraries(krino_mesh_surface PUBLIC krino_surface)
find_package(stk REQUIRED)
target_link_libraries(krino_mesh_surface PUBLIC stk::stk_mesh_base)
target_link_libraries(krino_mesh_surface PUBLIC stk::stk_tools_lib)
target_link_libraries(krino_mesh_surface PUBLIC stk::stk_topology)
target_include_directories(krino_mesh_surface PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_surface>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/mesh_surface>)
target_sources(krino_mesh_surface PUBLIC
    FILE_SET krino_mesh_surface_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_mesh_surface_headers})
target_compile_definitions(krino_mesh_surface PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_mesh_surface PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_mesh_surface PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_mesh_surface PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_mesh_surface
    EXPORT krinoTargets
    FILE_SET krino_mesh_surface_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_master_element)
FILE(GLOB krino_master_element_headers CONFIGURE_DEPENDS krino/master_element/*.hpp)
FILE(GLOB krino_master_element_sources CONFIGURE_DEPENDS krino/master_element/*.cpp)
target_sources(krino_master_element PRIVATE ${krino_master_element_sources})
find_package(Intrepid2 REQUIRED)
target_link_libraries(krino_master_element PUBLIC Intrepid2::all_libs)
find_package(stk REQUIRED)
target_link_libraries(krino_master_element PUBLIC stk::stk_mesh_base)
target_link_libraries(krino_master_element PUBLIC stk::stk_topology)
target_include_directories(krino_master_element PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/master_element>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/master_element>)
target_sources(krino_master_element PUBLIC
    FILE_SET krino_master_element_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_master_element_headers})
target_compile_definitions(krino_master_element PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_master_element PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_master_element PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_master_element PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_master_element
    EXPORT krinoTargets
    FILE_SET krino_master_element_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_math_utils)
FILE(GLOB krino_math_utils_headers CONFIGURE_DEPENDS krino/math_utils/*.hpp)
FILE(GLOB krino_math_utils_sources CONFIGURE_DEPENDS krino/math_utils/*.cpp)
target_sources(krino_math_utils PRIVATE ${krino_math_utils_sources})
target_link_libraries(krino_math_utils PUBLIC krino_diagwriter)
find_package(stk REQUIRED)
target_link_libraries(krino_math_utils PUBLIC stk::stk_math)
target_include_directories(krino_math_utils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/math_utils>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/math_utils>)
target_sources(krino_math_utils PUBLIC
    FILE_SET krino_math_utils_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_math_utils_headers})
target_compile_definitions(krino_math_utils PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_math_utils PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_math_utils PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_math_utils PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_math_utils
    EXPORT krinoTargets
    FILE_SET krino_math_utils_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_mesh_utils)
FILE(GLOB krino_mesh_utils_headers CONFIGURE_DEPENDS krino/mesh_utils/*.hpp)
FILE(GLOB krino_mesh_utils_sources CONFIGURE_DEPENDS krino/mesh_utils/*.cpp)
target_sources(krino_mesh_utils PRIVATE ${krino_mesh_utils_sources})
target_link_libraries(krino_mesh_utils PUBLIC krino_diagwriter)
find_package(stk REQUIRED)
target_link_libraries(krino_mesh_utils PUBLIC stk::stk_io)
target_link_libraries(krino_mesh_utils PUBLIC stk::stk_math)
target_link_libraries(krino_mesh_utils PUBLIC stk::stk_mesh_base)
target_link_libraries(krino_mesh_utils PUBLIC stk::stk_topology)
target_include_directories(krino_mesh_utils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_utils>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/mesh_utils>)
target_sources(krino_mesh_utils PUBLIC
    FILE_SET krino_mesh_utils_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_mesh_utils_headers})
target_compile_definitions(krino_mesh_utils PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_mesh_utils PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_mesh_utils PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_mesh_utils PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_mesh_utils
    EXPORT krinoTargets
    FILE_SET krino_mesh_utils_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_refinement)
FILE(GLOB krino_refinement_headers CONFIGURE_DEPENDS krino/refinement/*.hpp)
FILE(GLOB krino_refinement_sources CONFIGURE_DEPENDS krino/refinement/*.cpp)
target_sources(krino_refinement PRIVATE ${krino_refinement_sources})
target_link_libraries(krino_refinement PUBLIC krino_mesh_utils)
target_link_libraries(krino_refinement PUBLIC krino_quality_metric)
find_package(stk REQUIRED)
target_link_libraries(krino_refinement PUBLIC stk::stk_math)
target_link_libraries(krino_refinement PUBLIC stk::stk_mesh_base)
target_link_libraries(krino_refinement PUBLIC stk::stk_topology)
target_include_directories(krino_refinement PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/refinement>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/refinement>)
target_sources(krino_refinement PUBLIC
    FILE_SET krino_refinement_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_refinement_headers})
target_compile_definitions(krino_refinement PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_refinement PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_refinement PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_refinement PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_refinement
    EXPORT krinoTargets
    FILE_SET krino_refinement_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_refinement_rebalance)
FILE(GLOB krino_refinement_rebalance_headers CONFIGURE_DEPENDS krino/refinement_rebalance/*.hpp)
FILE(GLOB krino_refinement_rebalance_sources CONFIGURE_DEPENDS krino/refinement_rebalance/*.cpp)
target_sources(krino_refinement_rebalance PRIVATE ${krino_refinement_rebalance_sources})
target_link_libraries(krino_refinement_rebalance PUBLIC krino_refinement)
find_package(stk REQUIRED)
target_link_libraries(krino_refinement_rebalance PUBLIC stk::stk_balance_lib)
target_include_directories(krino_refinement_rebalance PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/refinement_rebalance>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/refinement_rebalance>)
target_sources(krino_refinement_rebalance PUBLIC
    FILE_SET krino_refinement_rebalance_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_refinement_rebalance_headers})
target_compile_definitions(krino_refinement_rebalance PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_refinement_rebalance PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_refinement_rebalance PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_refinement_rebalance PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_refinement_rebalance
    EXPORT krinoTargets
    FILE_SET krino_refinement_rebalance_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_lib)
FILE(GLOB krino_lib_headers CONFIGURE_DEPENDS krino/krino_lib/*.hpp)
FILE(GLOB krino_lib_sources CONFIGURE_DEPENDS krino/krino_lib/*.cpp)
target_sources(krino_lib PRIVATE ${krino_lib_sources})
find_package(MPI REQUIRED COMPONENTS C)
target_link_libraries(krino_lib PUBLIC MPI::MPI_C)
find_package(SEACAS REQUIRED COMPONENTS SEACASIoss)
target_link_libraries(krino_lib PUBLIC SEACASIoss::Ioss)
target_link_libraries(krino_lib PUBLIC krino_geometry)
target_link_libraries(krino_lib PUBLIC krino_master_element)
target_link_libraries(krino_lib PUBLIC krino_math_utils)
target_link_libraries(krino_lib PUBLIC krino_mesh_surface)
target_link_libraries(krino_lib PUBLIC krino_mesh_utils)
target_link_libraries(krino_lib PUBLIC krino_quality_metric)
target_link_libraries(krino_lib PUBLIC krino_quality_metric_sens)
target_link_libraries(krino_lib PUBLIC krino_refinement)
target_link_libraries(krino_lib PUBLIC krino_surface)
find_package(stk REQUIRED)
target_link_libraries(krino_lib PUBLIC stk::stk_emend)
target_link_libraries(krino_lib PUBLIC stk::stk_io)
target_link_libraries(krino_lib PUBLIC stk::stk_math)
target_link_libraries(krino_lib PUBLIC stk::stk_mesh_base)
target_link_libraries(krino_lib PUBLIC stk::stk_search)
target_link_libraries(krino_lib PUBLIC stk::stk_tools_lib)
target_link_libraries(krino_lib PUBLIC stk::stk_util_diag)
target_include_directories(krino_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/krino_lib>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/krino_lib>)
target_sources(krino_lib PUBLIC
    FILE_SET krino_lib_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_lib_headers})
target_compile_definitions(krino_lib PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_lib PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_lib PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_lib PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_lib
    EXPORT krinoTargets
    FILE_SET krino_lib_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_rebalance_utils)
FILE(GLOB krino_rebalance_utils_headers CONFIGURE_DEPENDS krino/rebalance_utils/*.hpp)
FILE(GLOB krino_rebalance_utils_sources CONFIGURE_DEPENDS krino/rebalance_utils/*.cpp)
target_sources(krino_rebalance_utils PRIVATE ${krino_rebalance_utils_sources})
target_link_libraries(krino_rebalance_utils PUBLIC krino_lib)
find_package(stk REQUIRED)
target_link_libraries(krino_rebalance_utils PUBLIC stk::stk_balance_lib)
target_link_libraries(krino_rebalance_utils PUBLIC stk::stk_mesh_base)
target_include_directories(krino_rebalance_utils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/rebalance_utils>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/rebalance_utils>)
target_sources(krino_rebalance_utils PUBLIC
    FILE_SET krino_rebalance_utils_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_rebalance_utils_headers})
target_compile_definitions(krino_rebalance_utils PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_rebalance_utils PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_rebalance_utils PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_rebalance_utils PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_rebalance_utils
    EXPORT krinoTargets
    FILE_SET krino_rebalance_utils_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_region)
FILE(GLOB krino_region_headers CONFIGURE_DEPENDS krino/region/*.hpp)
FILE(GLOB krino_region_sources CONFIGURE_DEPENDS krino/region/*.cpp)
target_sources(krino_region PRIVATE ${krino_region_sources})
target_link_libraries(krino_region PUBLIC krino_lib)
target_include_directories(krino_region PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/region>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/region>)
target_sources(krino_region PUBLIC
    FILE_SET krino_region_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_region_headers})
target_compile_definitions(krino_region PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_region PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_region PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_region PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_region
    EXPORT krinoTargets
    FILE_SET krino_region_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_sierra)
FILE(GLOB krino_sierra_headers CONFIGURE_DEPENDS krino_sierra/*.hpp)
FILE(GLOB krino_sierra_sources CONFIGURE_DEPENDS krino_sierra/*.cpp)
target_sources(krino_sierra PRIVATE ${krino_sierra_sources})
target_link_libraries(krino_sierra PUBLIC krino_lib)
find_package(sierra_common REQUIRED)
target_link_libraries(krino_sierra PUBLIC sierra_common::sierra)
target_link_libraries(krino_sierra PUBLIC sierra_common::sierra_util_sctl)
target_link_libraries(krino_sierra PUBLIC sierra_common::sierra_util_user_input_function_parser)
target_link_libraries(krino_sierra PUBLIC sierra_common::sierraparser)
target_include_directories(krino_sierra PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino_sierra>)
target_sources(krino_sierra PUBLIC
    FILE_SET krino_sierra_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_sierra_headers})
target_compile_definitions(krino_sierra PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_sierra PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_sierra PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_sierra PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_sierra
    EXPORT krinoTargets
    FILE_SET krino_sierra_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(krino_parser)
FILE(GLOB krino_parser_headers CONFIGURE_DEPENDS krino/parser/*.hpp)
FILE(GLOB krino_parser_sources CONFIGURE_DEPENDS krino/parser/*.cpp)
target_sources(krino_parser PRIVATE ${krino_parser_sources})
target_link_libraries(krino_parser PUBLIC krino_lib)
target_link_libraries(krino_parser PUBLIC krino_region)
find_package(yaml-cpp REQUIRED)
target_link_libraries(krino_parser PUBLIC yaml-cpp)
target_include_directories(krino_parser PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/parser>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino/parser>)
target_sources(krino_parser PUBLIC
    FILE_SET krino_parser_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${krino_parser_headers})
target_compile_definitions(krino_parser PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_parser PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(krino_parser PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(krino_parser PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS krino_parser
    EXPORT krinoTargets
    FILE_SET krino_parser_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_library(mesh_adapt_lib)
FILE(GLOB mesh_adapt_lib_headers CONFIGURE_DEPENDS krino_mesh_adapt/mesh_adapt_lib/*.hpp)
FILE(GLOB mesh_adapt_lib_sources CONFIGURE_DEPENDS krino_mesh_adapt/mesh_adapt_lib/*.cpp)
target_sources(mesh_adapt_lib PRIVATE ${mesh_adapt_lib_sources})
target_link_libraries(mesh_adapt_lib PUBLIC krino_refinement)
find_package(stk REQUIRED)
target_link_libraries(mesh_adapt_lib PUBLIC stk::stk_io)
target_link_libraries(mesh_adapt_lib PUBLIC stk::stk_mesh_base)
target_link_libraries(mesh_adapt_lib PUBLIC stk::stk_tools_lib)
target_link_libraries(mesh_adapt_lib PUBLIC stk::stk_util_diag)
target_include_directories(mesh_adapt_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino_mesh_adapt/mesh_adapt_lib>
    $<INSTALL_INTERFACE:include/krino>
    $<INSTALL_INTERFACE:include/krino/krino_mesh_adapt/mesh_adapt_lib>)
target_sources(mesh_adapt_lib PUBLIC
    FILE_SET mesh_adapt_lib_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${mesh_adapt_lib_headers})
target_compile_definitions(mesh_adapt_lib PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(mesh_adapt_lib PUBLIC Build64)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(mesh_adapt_lib PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>)
endif ()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(mesh_adapt_lib PUBLIC $<$<COMPILE_LANGUAGE:C>:-Wshadow>)
endif ()
install(
    TARGETS mesh_adapt_lib
    EXPORT krinoTargets
    FILE_SET mesh_adapt_lib_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

add_executable(krino)
target_sources(krino PRIVATE krino/Apps_krino.cpp)
target_link_libraries(krino PUBLIC krino_parser)
target_link_libraries(krino PUBLIC krino_region)
find_package(stk REQUIRED)
target_link_libraries(krino PUBLIC stk::stk_util_registry)
target_include_directories(krino PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/krino>)
target_compile_definitions(krino PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino PUBLIC Build64)
endif ()
install(TARGETS krino)

if (BUILD_TESTS)
        FILE(GLOB krino_unit_headers CONFIGURE_DEPENDS krino/unit_tests/*.hpp)
        FILE(GLOB krino_unit_sources CONFIGURE_DEPENDS krino/unit_tests/*.cpp)
	add_executable(krino_unit)
	target_sources(krino_unit PRIVATE ${krino_unit_sources})
	target_link_libraries(krino_unit PUBLIC krino_math_utils)
        target_link_libraries(krino_unit PUBLIC krino_quality_metric_sens)
	target_link_libraries(krino_unit PUBLIC krino_rebalance_utils)
	target_link_libraries(krino_unit PUBLIC krino_region)
	find_package(stk REQUIRED)
	target_link_libraries(krino_unit PUBLIC stk::stk_unit_test_utils)
	target_include_directories(krino_unit PUBLIC
	    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_tests>
	    $<INSTALL_INTERFACE:include/krino>
	    $<INSTALL_INTERFACE:include/krino/krino/unit_tests>)
	target_sources(krino_unit PUBLIC
	    FILE_SET krino_unit_headers
	    TYPE HEADERS
	    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
	    FILES ${krino_unit_headers})
	target_compile_definitions(krino_unit PUBLIC KRINO_BUILT_IN_SIERRA)
	if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
	    target_compile_definitions(krino_unit PUBLIC Build64)
	endif ()
	install(
	    TARGETS krino_unit
	    FILE_SET krino_unit_headers
	        DESTINATION include/krino
	        INCLUDES DESTINATION include/krino
	)
endif()

add_executable(delete_small_elements)
target_sources(delete_small_elements PRIVATE delete_small_elements/Akri_DeleteSmallElementsMain.cpp)
target_link_libraries(delete_small_elements PUBLIC krino_lib)
find_package(stk REQUIRED)
target_link_libraries(delete_small_elements PUBLIC stk::stk_util_command_line)
target_include_directories(delete_small_elements PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/krino>)
target_compile_definitions(delete_small_elements PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(delete_small_elements PUBLIC Build64)
endif ()
install(TARGETS delete_small_elements)

add_executable(krino_mesh_adapt)
target_sources(krino_mesh_adapt PRIVATE krino_mesh_adapt/KrinoMeshAdaptMain.cpp)
find_package(MPI REQUIRED COMPONENTS C)
target_link_libraries(krino_mesh_adapt PUBLIC MPI::MPI_C)
target_link_libraries(krino_mesh_adapt PUBLIC mesh_adapt_lib)
find_package(stk REQUIRED)
target_link_libraries(krino_mesh_adapt PUBLIC stk::stk_util_command_line)
target_include_directories(krino_mesh_adapt PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/krino>)
target_sources(krino_mesh_adapt PUBLIC
    FILE_SET krino_mesh_adapt_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(krino_mesh_adapt PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_mesh_adapt PUBLIC Build64)
endif ()
install(TARGETS krino_mesh_adapt)

install(
    EXPORT krinoTargets
    NAMESPACE krino::
    DESTINATION share/cmake/krino)

install(
    FILES cmake/krinoConfig.cmake
    DESTINATION share/cmake/krino
    )

