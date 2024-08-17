cmake_minimum_required(VERSION 3.23)
project(krino LANGUAGES C CXX Fortran)
cmake_path(GET CMAKE_CURRENT_LIST_DIR PARENT_PATH SIERRA_SOURCE_DIR)
list(PREPEND CMAKE_MODULE_PATH ${SIERRA_SOURCE_DIR}/modules)
include(${SIERRA_SOURCE_DIR}/modules/addParserCommands.cmake)
add_parser_commands(TARGET krino_commands  
	XML_FILES
		${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra/xml/Akri_Levelset.xml)
install(FILES ${CMAKE_BINARY_DIR}/krino_commands.xmldb DESTINATION xml)

add_library(krino_diagwriter)
target_sources(krino_diagwriter PRIVATE krino/diagwriter/Akri_DiagWriter.cpp)
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
    FILES krino/diagwriter/Akri_DiagWriter.hpp
        krino/diagwriter/Akri_DiagWriter_fwd.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_geometry PRIVATE krino/geometry/Akri_BoundingBox.cpp
    krino/geometry/Akri_BoundingBoxDistance.cpp
    krino/geometry/Akri_Plane_Intersections.cpp
    krino/geometry/Akri_WindingNumber.cpp)
find_package(stk REQUIRED)
target_link_libraries(krino_geometry PUBLIC stk::stk_math)
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
    FILES krino/geometry/Akri_BoundingBox.hpp
        krino/geometry/Akri_BoundingBoxDistance.hpp
        krino/geometry/Akri_Plane_Intersections.hpp
        krino/geometry/Akri_SearchTree.hpp
        krino/geometry/Akri_Segment.hpp
        krino/geometry/Akri_Triangle.hpp
        krino/geometry/Akri_WindingNumber.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
add_library(krino_surface)
target_sources(krino_surface PRIVATE krino/surface/Akri_AnalyticSurf.cpp
    krino/surface/Akri_Composite_Surface.cpp
    krino/surface/Akri_Facet.cpp
    krino/surface/Akri_FacetedSurfaceCalcs.cpp
    krino/surface/Akri_Faceted_Surface.cpp
    krino/surface/Akri_String_Function_Expression.cpp
    krino/surface/Akri_Surface.cpp
    krino/surface/Akri_SurfaceIntersectionFromSignedDistance.cpp
    krino/surface/Akri_Transformation.cpp)
target_link_libraries(krino_surface PUBLIC krino_diagwriter)
target_link_libraries(krino_surface PUBLIC krino_geometry)
target_link_libraries(krino_surface PUBLIC krino_math_utils)
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
    FILES krino/surface/Akri_AnalyticSurf.hpp
        krino/surface/Akri_Composite_Surface.hpp
        krino/surface/Akri_Facet.hpp
        krino/surface/Akri_FacetedSurfaceCalcs.hpp
        krino/surface/Akri_Faceted_Surface.hpp
        krino/surface/Akri_String_Function_Expression.hpp
        krino/surface/Akri_Surface.hpp
        krino/surface/Akri_SurfaceIntersectionFromSignedDistance.hpp
        krino/surface/Akri_Transformation.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_mesh_surface PRIVATE krino/mesh_surface/Akri_MeshSurface.cpp)
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
    FILES krino/mesh_surface/Akri_MeshSurface.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_master_element PRIVATE krino/master_element/Akri_MasterElementCalc.cpp
    krino/master_element/Akri_MasterElementHybrid.cpp
    krino/master_element/Akri_MasterElementIntrepid.cpp)
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
    FILES krino/master_element/Akri_MasterElement.hpp
        krino/master_element/Akri_MasterElementBasis.hpp
        krino/master_element/Akri_MasterElementCalc.hpp
        krino/master_element/Akri_MasterElementHybrid.hpp
        krino/master_element/Akri_MasterElementIntrepid.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_math_utils PRIVATE krino/math_utils/Akri_CramersRuleSolver.cpp
    krino/math_utils/Akri_CurvatureLeastSquares.cpp
    krino/math_utils/Akri_MathUtil.cpp)
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
    FILES krino/math_utils/Akri_CramersRuleSolver.hpp
        krino/math_utils/Akri_CurvatureLeastSquares.hpp
        krino/math_utils/Akri_MathUtil.hpp
        krino/math_utils/Akri_MortonIndex.hpp
        krino/math_utils/Akri_Sign.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_mesh_utils PRIVATE krino/mesh_utils/Akri_ChildNodeCreator.cpp
    krino/mesh_utils/Akri_Edge.cpp
    krino/mesh_utils/Akri_EntityIdPool.cpp
    krino/mesh_utils/Akri_FieldRef.cpp
    krino/mesh_utils/Akri_MeshHelpers.cpp
    krino/mesh_utils/Akri_ParallelErrorMessage.cpp
    krino/mesh_utils/Akri_QuadFace.cpp
    krino/mesh_utils/Akri_SideAttachedElements.cpp)
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
    FILES krino/mesh_utils/Akri_AllReduce.hpp
        krino/mesh_utils/Akri_ChildNodeCreator.hpp
        krino/mesh_utils/Akri_Edge.hpp
        krino/mesh_utils/Akri_EntityIdPool.hpp
        krino/mesh_utils/Akri_FieldRef.hpp
        krino/mesh_utils/Akri_MeshHelpers.hpp
        krino/mesh_utils/Akri_ParallelErrorMessage.hpp
        krino/mesh_utils/Akri_QuadFace.hpp
        krino/mesh_utils/Akri_ReportHandler.hpp
        krino/mesh_utils/Akri_SideAttachedElements.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
add_library(krino_quality_metric)
target_sources(krino_quality_metric PRIVATE krino/quality_metric/Akri_QualityMetric.cpp)
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
    FILES krino/quality_metric/Akri_QualityMetric.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
add_library(krino_refinement)
target_sources(krino_refinement PRIVATE krino/refinement/Akri_HexRefiner.cpp
    krino/refinement/Akri_MOAB_TetRefiner.cpp
    krino/refinement/Akri_NodeRefiner.cpp
    krino/refinement/Akri_QuadRefiner.cpp
    krino/refinement/Akri_Refinement.cpp
    krino/refinement/Akri_TransitionElementEdgeMarker.cpp
    krino/refinement/Akri_TriRefiner.cpp)
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
    FILES krino/refinement/Akri_HexRefiner.hpp
        krino/refinement/Akri_MOAB_TetRefiner.hpp
        krino/refinement/Akri_NodeRefiner.hpp
        krino/refinement/Akri_QuadRefiner.hpp
        krino/refinement/Akri_Refinement.hpp
        krino/refinement/Akri_RefinerUtils.hpp
        krino/refinement/Akri_TransitionElementEdgeMarker.hpp
        krino/refinement/Akri_TriRefiner.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_refinement_rebalance PRIVATE krino/refinement_rebalance/Akri_RefinementRebalance.cpp)
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
    FILES krino/refinement_rebalance/Akri_RefinementRebalance.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_lib PRIVATE krino/krino_lib/Akri_AdaptiveElementContour.cpp
    krino/krino_lib/Akri_AdaptivityHelpers.cpp
    krino/krino_lib/Akri_AnalyticSurfaceInterfaceGeometry.cpp
    krino/krino_lib/Akri_AuxMetaData.cpp
    krino/krino_lib/Akri_BoundingBoxMesh.cpp
    krino/krino_lib/Akri_BoundingSurface.cpp
    krino/krino_lib/Akri_CDFEM_Parent_Edge.cpp
    krino/krino_lib/Akri_CDFEM_Parent_Edges.cpp
    krino/krino_lib/Akri_CDFEM_Support.cpp
    krino/krino_lib/Akri_CDMesh.cpp
    krino/krino_lib/Akri_CDMesh_Debug.cpp
    krino/krino_lib/Akri_CDMesh_Refinement.cpp
    krino/krino_lib/Akri_CDMesh_Utils.cpp
    krino/krino_lib/Akri_ChildNodeStencil.cpp
    krino/krino_lib/Akri_Compute_Surface_Distance.cpp
    krino/krino_lib/Akri_ConformingPhaseParts.cpp
    krino/krino_lib/Akri_ContourElement.cpp
    krino/krino_lib/Akri_ContourSubElement.cpp
    krino/krino_lib/Akri_CreateInterfaceGeometry.cpp
    krino/krino_lib/Akri_Cutting_Surface.cpp
    krino/krino_lib/Akri_DecompositionHasChanged.cpp
    krino/krino_lib/Akri_DetermineElementSign.cpp
    krino/krino_lib/Akri_DistanceSweeper.cpp
    krino/krino_lib/Akri_Eikonal_Calc.cpp
    krino/krino_lib/Akri_Element.cpp
    krino/krino_lib/Akri_ElementCutterUtils.cpp
    krino/krino_lib/Akri_Element_Cutter.cpp
    krino/krino_lib/Akri_Element_Intersections.cpp
    krino/krino_lib/Akri_FastIterativeMethod.cpp
    krino/krino_lib/Akri_Fast_Marching.cpp
    krino/krino_lib/Akri_IC_Alg.cpp
    krino/krino_lib/Akri_IC_Calculator.cpp
    krino/krino_lib/Akri_IO_Helpers.cpp
    krino/krino_lib/Akri_Intersection_Points.cpp
    krino/krino_lib/Akri_LevelSet.cpp
    krino/krino_lib/Akri_LevelSetInterfaceGeometry.cpp
    krino/krino_lib/Akri_LevelSetPolicy.cpp
    krino/krino_lib/Akri_LevelSetShapeSensitivities.cpp
    krino/krino_lib/Akri_LevelSetSurfaceInterfaceGeometry.cpp
    krino/krino_lib/Akri_LowerEnvelope.cpp
    krino/krino_lib/Akri_MasterElementDeterminer.cpp
    krino/krino_lib/Akri_MeshClone.cpp
    krino/krino_lib/Akri_MeshDiagnostics.cpp
    krino/krino_lib/Akri_MeshFromFile.cpp
    krino/krino_lib/Akri_MeshInputOptions.cpp
    krino/krino_lib/Akri_NodalBoundingBox.cpp
    krino/krino_lib/Akri_NodalSurfaceDistance.cpp
    krino/krino_lib/Akri_NodeToCapturedDomains.cpp
    krino/krino_lib/Akri_OutputUtils.cpp
    krino/krino_lib/Akri_ParentsToChildMapper.cpp
    krino/krino_lib/Akri_PatchInterpolator.cpp
    krino/krino_lib/Akri_PhaseTag.cpp
    krino/krino_lib/Akri_Phase_Support.cpp
    krino/krino_lib/Akri_PostProcess.cpp
    krino/krino_lib/Akri_ProlongationData.cpp
    krino/krino_lib/Akri_Quality.cpp
    krino/krino_lib/Akri_RefineNearLevelSets.cpp
    krino/krino_lib/Akri_RefinementInterface.cpp
    krino/krino_lib/Akri_RefinementSupport.cpp
    krino/krino_lib/Akri_SemiLagrangian.cpp
    krino/krino_lib/Akri_SharpFeature.cpp
    krino/krino_lib/Akri_Snap.cpp
    krino/krino_lib/Akri_SnapInfo.cpp
    krino/krino_lib/Akri_SnapToNode.cpp
    krino/krino_lib/Akri_SubElement.cpp
    krino/krino_lib/Akri_SubElementChildNodeAncestry.cpp
    krino/krino_lib/Akri_SubElementNodeAncestry.cpp
    krino/krino_lib/Akri_Surface_Manager.cpp
    krino/krino_lib/Akri_VolumePreservingSnappingLimiter.cpp)
find_package(MPI REQUIRED COMPONENTS C Fortran)
target_link_libraries(krino_lib PUBLIC MPI::MPI_C)
target_link_libraries(krino_lib PUBLIC MPI::MPI_Fortran)
find_package(SEACAS REQUIRED COMPONENTS SEACASIoss)
target_link_libraries(krino_lib PUBLIC SEACASIoss::Ioss)
target_link_libraries(krino_lib PUBLIC krino_geometry)
target_link_libraries(krino_lib PUBLIC krino_master_element)
target_link_libraries(krino_lib PUBLIC krino_math_utils)
target_link_libraries(krino_lib PUBLIC krino_mesh_surface)
target_link_libraries(krino_lib PUBLIC krino_mesh_utils)
target_link_libraries(krino_lib PUBLIC krino_quality_metric)
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
    FILES krino/krino_lib/Akri_AdaptiveElementContour.hpp
        krino/krino_lib/Akri_AdaptivityHelpers.hpp
        krino/krino_lib/Akri_AnalyticSurfaceInterfaceGeometry.hpp
        krino/krino_lib/Akri_AuxMetaData.hpp
        krino/krino_lib/Akri_BoundingBoxMesh.hpp
        krino/krino_lib/Akri_BoundingSurface.hpp
        krino/krino_lib/Akri_CDFEM_Parent_Edge.hpp
        krino/krino_lib/Akri_CDFEM_Parent_Edges.hpp
        krino/krino_lib/Akri_CDFEM_Snapper.hpp
        krino/krino_lib/Akri_CDFEM_Support.hpp
        krino/krino_lib/Akri_CDMesh.hpp
        krino/krino_lib/Akri_CDMesh_Debug.hpp
        krino/krino_lib/Akri_CDMesh_Refinement.hpp
        krino/krino_lib/Akri_CDMesh_Utils.hpp
        krino/krino_lib/Akri_ChildNodeStencil.hpp
        krino/krino_lib/Akri_Compute_Surface_Distance.hpp
        krino/krino_lib/Akri_ConformingPhaseParts.hpp
        krino/krino_lib/Akri_ContourElement.hpp
        krino/krino_lib/Akri_ContourSubElement.hpp
        krino/krino_lib/Akri_CreateInterfaceGeometry.hpp
        krino/krino_lib/Akri_Cutting_Surface.hpp
        krino/krino_lib/Akri_DecompositionHasChanged.hpp
        krino/krino_lib/Akri_DetermineElementSign.hpp
        krino/krino_lib/Akri_DistanceSweeper.hpp
        krino/krino_lib/Akri_Eikonal_Calc.hpp
        krino/krino_lib/Akri_Element.hpp
        krino/krino_lib/Akri_ElementCutterUtils.hpp
        krino/krino_lib/Akri_Element_Cutter.hpp
        krino/krino_lib/Akri_Element_Intersections.hpp
        krino/krino_lib/Akri_FastIterativeMethod.hpp
        krino/krino_lib/Akri_Fast_Marching.hpp
        krino/krino_lib/Akri_IC_Alg.hpp
        krino/krino_lib/Akri_IC_Calculator.hpp
        krino/krino_lib/Akri_IO_Helpers.hpp
        krino/krino_lib/Akri_InterfaceGeometry.hpp
        krino/krino_lib/Akri_InterfaceID.hpp
        krino/krino_lib/Akri_Interface_Name_Generator.hpp
        krino/krino_lib/Akri_Intersection_Points.hpp
        krino/krino_lib/Akri_LevelSet.hpp
        krino/krino_lib/Akri_LevelSetInterfaceGeometry.hpp
        krino/krino_lib/Akri_LevelSetPolicy.hpp
        krino/krino_lib/Akri_LevelSetShapeSensitivities.hpp
        krino/krino_lib/Akri_LevelSetSurfaceInterfaceGeometry.hpp
        krino/krino_lib/Akri_LowerEnvelope.hpp
        krino/krino_lib/Akri_MasterElementDeterminer.hpp
        krino/krino_lib/Akri_MeshClone.hpp
        krino/krino_lib/Akri_MeshDiagnostics.hpp
        krino/krino_lib/Akri_MeshFromFile.hpp
        krino/krino_lib/Akri_MeshInputOptions.hpp
        krino/krino_lib/Akri_MeshInterface.hpp
        krino/krino_lib/Akri_NodalBoundingBox.hpp
        krino/krino_lib/Akri_NodalSurfaceDistance.hpp
        krino/krino_lib/Akri_NodeToCapturedDomains.hpp
        krino/krino_lib/Akri_OrderedIdPair.hpp
        krino/krino_lib/Akri_OutputUtils.hpp
        krino/krino_lib/Akri_ParentsToChildMapper.hpp
        krino/krino_lib/Akri_PatchInterpolator.hpp
        krino/krino_lib/Akri_PhaseTag.hpp
        krino/krino_lib/Akri_Phase_Support.hpp
        krino/krino_lib/Akri_PostProcess.hpp
        krino/krino_lib/Akri_ProlongationData.hpp
        krino/krino_lib/Akri_Quality.hpp
        krino/krino_lib/Akri_RefineNearLevelSets.hpp
        krino/krino_lib/Akri_RefinementInterface.hpp
        krino/krino_lib/Akri_RefinementSupport.hpp
        krino/krino_lib/Akri_SemiLagrangian.hpp
        krino/krino_lib/Akri_SharpFeature.hpp
        krino/krino_lib/Akri_Snap.hpp
        krino/krino_lib/Akri_SnapIndependentSetFinder.hpp
        krino/krino_lib/Akri_SnapInfo.hpp
        krino/krino_lib/Akri_SnapToNode.hpp
        krino/krino_lib/Akri_SubElement.hpp
        krino/krino_lib/Akri_SubElementChildNodeAncestry.hpp
        krino/krino_lib/Akri_SubElementNodeAncestry.hpp
        krino/krino_lib/Akri_Surface_Identifier.hpp
        krino/krino_lib/Akri_Surface_Manager.hpp
        krino/krino_lib/Akri_TopologyData.hpp
        krino/krino_lib/Akri_TypeDefs.hpp
        krino/krino_lib/Akri_Utility.hpp
        krino/krino_lib/Akri_VolumePreservingSnappingLimiter.hpp
        krino/krino_lib/Akri_config.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_rebalance_utils PRIVATE krino/rebalance_utils/Akri_RebalanceUtils.cpp
    krino/rebalance_utils/Akri_RebalanceUtils_Impl.cpp)
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
    FILES krino/rebalance_utils/Akri_RebalanceUtils.hpp
        krino/rebalance_utils/Akri_RebalanceUtils_Impl.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_region PRIVATE krino/region/Akri_Region.cpp
    krino/region/Akri_RegisterProduct.cpp
    krino/region/Akri_Simulation.cpp
    krino/region/Akri_Startup.cpp)
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
    FILES krino/region/Akri_Region.hpp
        krino/region/Akri_RegisterProduct.hpp
        krino/region/Akri_ResultsOutputOptions.hpp
        krino/region/Akri_Simulation.hpp
        krino/region/Akri_Startup.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_sierra PRIVATE krino_sierra/Akri_CDFEM_Options_SierraParser.cpp
    krino_sierra/Akri_Events.cpp
    krino_sierra/Akri_IC_SierraParser.cpp
    krino_sierra/Akri_LevelSet_Sctl.cpp
    krino_sierra/Akri_LevelSet_SierraParser.cpp
    krino_sierra/Akri_Motion_SierraParser.cpp
    krino_sierra/Akri_PerceptRefinementInterface.cpp
    krino_sierra/Akri_Phase_SierraParser.cpp
    krino_sierra/Akri_RegionInterface.cpp)
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
    FILES krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino_parser PRIVATE krino/parser/Akri_CDFEM_Options_Parser.cpp
    krino/parser/Akri_IC_Parser.cpp
    krino/parser/Akri_LevelSet_Parser.cpp
    krino/parser/Akri_MeshInput_Parser.cpp
    krino/parser/Akri_Parser.cpp
    krino/parser/Akri_Phase_Parser.cpp
    krino/parser/Akri_Region_Parser.cpp
    krino/parser/Akri_ResultsOutput_Parser.cpp
    krino/parser/Akri_Simulation_Parser.cpp
    krino/parser/Akri_Surface_Parser.cpp)
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
    FILES krino/parser/Akri_CDFEM_Options_Parser.hpp
        krino/parser/Akri_IC_Parser.hpp
        krino/parser/Akri_LevelSet_Parser.hpp
        krino/parser/Akri_MeshInput_Parser.hpp
        krino/parser/Akri_Parser.hpp
        krino/parser/Akri_Phase_Parser.hpp
        krino/parser/Akri_Region_Parser.hpp
        krino/parser/Akri_ResultsOutput_Parser.hpp
        krino/parser/Akri_Simulation_Parser.hpp
        krino/parser/Akri_Surface_Parser.hpp
        krino/parser/Akri_YAML.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(mesh_adapt_lib PRIVATE krino_mesh_adapt/mesh_adapt_lib/KrinoMeshAdapt.cpp
    krino_mesh_adapt/mesh_adapt_lib/KrinoMeshAdaptParser.cpp)
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
    FILES krino_mesh_adapt/mesh_adapt_lib/KrinoMeshAdapt.hpp
        krino_mesh_adapt/mesh_adapt_lib/KrinoMeshAdaptAlgorithmParameters.hpp
        krino_mesh_adapt/mesh_adapt_lib/KrinoMeshAdaptInputData.hpp
        krino_mesh_adapt/mesh_adapt_lib/KrinoMeshAdaptParser.hpp
        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(krino PUBLIC
    FILE_SET krino_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
target_compile_definitions(krino PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino PUBLIC Build64)
endif ()
install(
    TARGETS krino
    FILE_SET krino_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)
if (${SIERRA_DEVELOPER_BUILD})
	add_executable(krino_unit)
	target_sources(krino_unit PRIVATE krino/unit_tests/Akri_StkMeshBuilder.cpp
	    krino/unit_tests/Akri_UnitMathUtils.cpp
	    krino/unit_tests/Akri_UnitMeshUtils.cpp
	    krino/unit_tests/Akri_UnitTestUtils.cpp
	    krino/unit_tests/Akri_Unit_Analytic_CDMesh.cpp
	    krino/unit_tests/Akri_Unit_BoundingBoxDistance.cpp
	    krino/unit_tests/Akri_Unit_CDFEM_Parent_Edge.cpp
	    krino/unit_tests/Akri_Unit_CDMesh.cpp
	    krino/unit_tests/Akri_Unit_Constrained_Redistance.cpp
	    krino/unit_tests/Akri_Unit_ContourElement.cpp
	    krino/unit_tests/Akri_Unit_CurvatureLeastSquares.cpp
	    krino/unit_tests/Akri_Unit_DecomposeWithSensitivities.cpp
	    krino/unit_tests/Akri_Unit_Eikonal.cpp
	    krino/unit_tests/Akri_Unit_Element.cpp
	    krino/unit_tests/Akri_Unit_Element_Cutter.cpp
	    krino/unit_tests/Akri_Unit_Explicit_Hamilton_Jacobi.cpp
	    krino/unit_tests/Akri_Unit_FastMarching.cpp
	    krino/unit_tests/Akri_Unit_Geometry.cpp
	    krino/unit_tests/Akri_Unit_InterfaceGeometry.cpp
	    krino/unit_tests/Akri_Unit_LogRedirecter.cpp
	    krino/unit_tests/Akri_Unit_LowerEnvelope.cpp
	    krino/unit_tests/Akri_Unit_MeshHelpers.cpp
	    krino/unit_tests/Akri_Unit_MortonIndex.cpp
	    krino/unit_tests/Akri_Unit_OutputUtils.cpp
	    krino/unit_tests/Akri_Unit_ParallelErrorMessage.cpp
	    krino/unit_tests/Akri_Unit_Part_Decomposition_Fixture.cpp
	    krino/unit_tests/Akri_Unit_PatchInterpolator.cpp
	    krino/unit_tests/Akri_Unit_Phase_Support.cpp
	    krino/unit_tests/Akri_Unit_RebalanceUtils.cpp
	    krino/unit_tests/Akri_Unit_RefineInterval.cpp
	    krino/unit_tests/Akri_Unit_Refine_Beam.cpp
	    krino/unit_tests/Akri_Unit_Refine_CDMesh.cpp
	    krino/unit_tests/Akri_Unit_Refine_General.cpp
	    krino/unit_tests/Akri_Unit_Refine_Hex.cpp
	    krino/unit_tests/Akri_Unit_Refine_Quad.cpp
	    krino/unit_tests/Akri_Unit_Refine_Tet.cpp
	    krino/unit_tests/Akri_Unit_Refine_Tri.cpp
	    krino/unit_tests/Akri_Unit_SearchTree.cpp
	    krino/unit_tests/Akri_Unit_SemiLagrangian.cpp
	    krino/unit_tests/Akri_Unit_SideAttachedElements.cpp
	    krino/unit_tests/Akri_Unit_Single_Element_Fixtures.cpp
	    krino/unit_tests/Akri_Unit_Snap.cpp
	    krino/unit_tests/Akri_Unit_TriangleCalcs.cpp
	    krino/unit_tests/Akri_Unit_WindingNumber.cpp
	    krino/unit_tests/Akri_Unit_main.cpp)
	target_link_libraries(krino_unit PUBLIC krino_math_utils)
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
	    FILES krino/unit_tests/Akri_MeshSpecs.hpp
	        krino/unit_tests/Akri_StkMeshBuilder.hpp
	        krino/unit_tests/Akri_StkMeshFixture.hpp
	        krino/unit_tests/Akri_UnitMeshUtils.hpp
	        krino/unit_tests/Akri_UnitTestUtils.hpp
	        krino/unit_tests/Akri_Unit_BoundingBoxMesh.hpp
	        krino/unit_tests/Akri_Unit_DecompositionFixture.hpp
	        krino/unit_tests/Akri_Unit_InterfaceGeometry.hpp
	        krino/unit_tests/Akri_Unit_LogRedirecter.hpp
	        krino/unit_tests/Akri_Unit_MeshHelpers.hpp
	        krino/unit_tests/Akri_Unit_Part_Decomposition_Fixture.hpp
	        krino/unit_tests/Akri_Unit_RefinementFixture.hpp
	        krino/unit_tests/Akri_Unit_Single_Element_Fixtures.hpp
	        krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
	        krino_sierra/Akri_Events.hpp
	        krino_sierra/Akri_IC_SierraParser.hpp
	        krino_sierra/Akri_LevelSet_Sctl.hpp
	        krino_sierra/Akri_LevelSet_SierraParser.hpp
	        krino_sierra/Akri_Motion_SierraParser.hpp
	        krino_sierra/Akri_PerceptRefinementInterface.hpp
	        krino_sierra/Akri_Phase_SierraParser.hpp
	        krino_sierra/Akri_RegionInterface.hpp)
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
target_sources(delete_small_elements PUBLIC
    FILE_SET delete_small_elements_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
target_compile_definitions(delete_small_elements PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(delete_small_elements PUBLIC Build64)
endif ()
install(
    TARGETS delete_small_elements
    FILE_SET delete_small_elements_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)
add_executable(krino_mesh_adapt)
target_sources(krino_mesh_adapt PRIVATE krino_mesh_adapt/KrinoMeshAdaptMain.cpp)
find_package(MPI REQUIRED COMPONENTS C Fortran)
target_link_libraries(krino_mesh_adapt PUBLIC MPI::MPI_C)
target_link_libraries(krino_mesh_adapt PUBLIC MPI::MPI_Fortran)
target_link_libraries(krino_mesh_adapt PUBLIC mesh_adapt_lib)
find_package(stk REQUIRED)
target_link_libraries(krino_mesh_adapt PUBLIC stk::stk_util_command_line)
target_include_directories(krino_mesh_adapt PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/krino>)
target_sources(krino_mesh_adapt PUBLIC
    FILE_SET krino_mesh_adapt_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES krino_sierra/Akri_CDFEM_Options_SierraParser.hpp
        krino_sierra/Akri_Events.hpp
        krino_sierra/Akri_IC_SierraParser.hpp
        krino_sierra/Akri_LevelSet_Sctl.hpp
        krino_sierra/Akri_LevelSet_SierraParser.hpp
        krino_sierra/Akri_Motion_SierraParser.hpp
        krino_sierra/Akri_PerceptRefinementInterface.hpp
        krino_sierra/Akri_Phase_SierraParser.hpp
        krino_sierra/Akri_RegionInterface.hpp)
target_compile_definitions(krino_mesh_adapt PUBLIC KRINO_BUILT_IN_SIERRA)
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(krino_mesh_adapt PUBLIC Build64)
endif ()
install(
    TARGETS krino_mesh_adapt
    FILE_SET krino_mesh_adapt_headers
        DESTINATION include/krino
        INCLUDES DESTINATION include/krino
)

install(
    EXPORT krinoTargets
    NAMESPACE krino::
    DESTINATION share/cmake/krino)

install(
    FILES cmake/krinoConfig.cmake
    DESTINATION share/cmake/krino
    )

install(
	FILES ${SIERRA_SOURCE_DIR}/modules/createParserTarget.cmake
	DESTINATION share/cmake/krino
	)