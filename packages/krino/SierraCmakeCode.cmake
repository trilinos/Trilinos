cmake_minimum_required(VERSION 3.23)

project(krino LANGUAGES C CXX Fortran)

include(cmake/commonFunctions.cmake)
include(cmake/krinoFunctions.cmake)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

option(Krino_ENABLE_PIC
    "Build position-independent krino libraries (needed to link krino into \
shared modules such as the pykrino bindings)." OFF)
if (Krino_ENABLE_PIC)
    # Must be set before any target is declared.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------

find_package(sierra-cmake-utils REQUIRED)
find_package(stk REQUIRED)
find_package(Sacado REQUIRED)
find_package(MPI REQUIRED COMPONENTS C)
find_package(pocket_tensor REQUIRED)
find_package(Intrepid2 REQUIRED)
find_package(MiniTensor REQUIRED)
find_package(ROL REQUIRED)
find_package(SEACAS REQUIRED COMPONENTS SEACASIoss)
find_package(sierra_common REQUIRED)
find_package(yaml-cpp REQUIRED)

# ---------------------------------------------------------------------------
# Parser command database
# ---------------------------------------------------------------------------

add_parser_commands(
    TARGET krino_commands
    XML_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra/xml/Akri_Levelset.xml
)

install(
    FILES ${CMAKE_BINARY_DIR}/krino_commands.xmldb
    DESTINATION xml
)

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------

krino_add_interface_library(
    krino_config
    krino/config
)

krino_add_library(
    krino_diagwriter
    krino/diagwriter
    PUBLIC_stk::stk_util_diag
)

krino_add_library(
    krino_geometry
    krino/geometry
    PUBLIC_stk::stk_math
    PUBLIC_Sacado::all_libs
    PUBLIC_stk::stk_util_parallel
)

krino_add_library(
    krino_quality_metric
    krino/quality_metric
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_quality_metric_sens
    krino/quality_metric_sens
    PUBLIC_MPI::MPI_C
    PUBLIC_stk::stk_math
    PUBLIC_Sacado::all_libs
    PUBLIC_krino_quality_metric
)

krino_add_library(
    krino_math_utils
    krino/math_utils
    PUBLIC_krino_diagwriter
    PUBLIC_stk::stk_math
)

krino_add_library(
    krino_surface
    krino/surface
    PUBLIC_krino_config
    PUBLIC_krino_diagwriter
    PUBLIC_krino_geometry
    PUBLIC_krino_math_utils
    PUBLIC_krino_quality_metric
    PUBLIC_stk::stk_expreval
    PUBLIC_pocket_tensor::pocket_tensor
)

krino_add_library(
    krino_plato_workarounds
    krino/plato_workarounds
    PUBLIC_krino_interface_geometry
    PUBLIC_krino_surface
)

krino_add_library(
    krino_mesh_surface
    krino/mesh_surface
    PUBLIC_krino_surface
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_parametric_element
    krino/parametric_element
    PUBLIC_Intrepid2::all_libs
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_mesh_utils
    krino/mesh_utils
    PUBLIC_krino_diagwriter
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_foundation
    krino/foundation
    PUBLIC_krino_diagwriter
    PUBLIC_krino_math_utils
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_parametric_element
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_field_calcs
    krino/field_calcs
    PUBLIC_krino_foundation
    PUBLIC_krino_geometry
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_mesh_clone
    krino/mesh_clone
    PUBLIC_krino_diagwriter
    PUBLIC_krino_foundation
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_util_diag
    PUBLIC_stk::stk_util_parallel
)

krino_add_library(
    krino_eikonal
    krino/eikonal
    PUBLIC_krino_diagwriter
    PUBLIC_krino_geometry
    PUBLIC_krino_math_utils
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_util_diag
)

krino_add_library(
    krino_mesh_quality
    krino/mesh_quality
    PUBLIC_krino_quality_metric
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_metric
    krino/metric
    PUBLIC_MPI::MPI_C
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_util_util
    PUBLIC_MiniTensor::all_libs
)

krino_add_library(
    krino_smoothing
    krino/smoothing
    PUBLIC_krino_diagwriter
    PUBLIC_krino_math_utils
    PUBLIC_krino_metric
    PUBLIC_krino_mesh_quality
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_quality_metric_sens
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_MiniTensor::all_libs
)

krino_add_library(
    krino_refinement
    krino/refinement
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_quality_metric
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_refinement_rebalance
    krino/refinement_rebalance
    PUBLIC_krino_refinement
    PUBLIC_stk::stk_balance_lib
)

krino_add_library(
    krino_mesh_metric
    krino/mesh_metric
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_metric
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_refinement_drivers
    krino/refinement_drivers
    PUBLIC_krino_diagwriter
    PUBLIC_krino_foundation
    PUBLIC_krino_mesh_metric
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_metric
    PUBLIC_krino_refinement
    PUBLIC_krino_refinement_rebalance
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_topology
    PUBLIC_stk::stk_util_diag
    PUBLIC_stk::stk_util_env
    PUBLIC_stk::stk_util_parallel
)

krino_add_library(
    krino_contour
    krino/contour
    PUBLIC_krino_foundation
    PUBLIC_krino_diagwriter
    PUBLIC_krino_geometry
    PUBLIC_krino_math_utils
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_quality_metric
    PUBLIC_krino_refinement
    PUBLIC_krino_surface
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_support
    krino/support
    PUBLIC_krino_diagwriter
    PUBLIC_krino_foundation
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_phase
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_levelset
    krino/levelset
    PUBLIC_krino_support
    PUBLIC_krino_contour
    PUBLIC_krino_eikonal
    PUBLIC_krino_field_calcs
    PUBLIC_krino_foundation
    PUBLIC_krino_geometry
    PUBLIC_krino_math_utils
    PUBLIC_krino_mesh_surface
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_phase
    PUBLIC_krino_surface
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_search
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_util_diag
)

krino_add_library(
    krino_rol
    krino/rol
    PUBLIC_stk::stk_math
    PUBLIC_ROL::all_libs
    PUBLIC_krino_diagwriter
    PUBLIC_krino_math_utils
)

krino_add_library(
    krino_mesh_interface
    krino/mesh_interface
    PUBLIC_MPI::MPI_C
    PUBLIC_SEACASIoss::Ioss
    PUBLIC_krino_diagwriter
    PUBLIC_krino_foundation
    PUBLIC_krino_geometry
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_util_diag
)

krino_add_library(
    krino_phase
    krino/phase
    PUBLIC_krino_diagwriter
    PUBLIC_krino_foundation
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_topology
)

krino_add_library(
    krino_interface_geometry
    krino/interface_geometry
    PUBLIC_MPI::MPI_C
    PUBLIC_SEACASIoss::Ioss
    PUBLIC_krino_contour
    PUBLIC_krino_eikonal
    PUBLIC_krino_field_calcs
    PUBLIC_krino_foundation
    PUBLIC_krino_geometry
    PUBLIC_krino_levelset
    PUBLIC_krino_mesh_interface
    PUBLIC_krino_parametric_element
    PUBLIC_krino_math_utils
    PUBLIC_krino_mesh_surface
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_phase
    PUBLIC_krino_quality_metric
    PUBLIC_krino_support
    PUBLIC_krino_surface
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_search
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_util_diag
)

krino_add_library(
    krino_decomp
    krino/decomp
    PUBLIC_MPI::MPI_C
    PUBLIC_SEACASIoss::Ioss
    PUBLIC_krino_interface_geometry
    PUBLIC_krino_contour
    PUBLIC_krino_eikonal
    PUBLIC_krino_field_calcs
    PUBLIC_krino_foundation
    PUBLIC_krino_geometry
    PUBLIC_krino_levelset
    PUBLIC_krino_mesh_interface
    PUBLIC_krino_parametric_element
    PUBLIC_krino_math_utils
    PUBLIC_krino_mesh_quality
    PUBLIC_krino_mesh_surface
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_phase
    PUBLIC_krino_quality_metric
    PUBLIC_krino_refinement
    PUBLIC_krino_refinement_drivers
    PUBLIC_krino_refinement_rebalance
    PUBLIC_krino_support
    PUBLIC_krino_surface
    PUBLIC_stk::stk_emend
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_search
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_util_diag
)

krino_add_library(
    krino_post_process
    krino/post_process
    PUBLIC_krino_diagwriter
    PUBLIC_krino_foundation
    PUBLIC_krino_mesh_utils
    PUBLIC_krino_surface
    PUBLIC_stk::stk_math
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_levelset_shape_sens
    krino/levelset_shape_sens
    PUBLIC_krino_decomp
    PUBLIC_krino_interface_geometry
)

# krino_lib is now a source-free umbrella (INTERFACE) library.  Its former
# contents have been dispersed into decomp, interface_geometry, post_process,
# mesh_clone and levelset_shape_sens; it re-exports them so existing consumers
# (region, parser, rebalance_utils, sierra) keep resolving every symbol
# transitively through the single krino_lib dependency.  It also carries the
# higher-level libs (smoothing, rol, plato_workarounds, quality_metric_sens)
# that consumers reached through it.
krino_add_interface_library(
    krino_lib
    krino/krino_lib
    INTERFACE_krino_decomp
    INTERFACE_krino_interface_geometry
    INTERFACE_krino_post_process
    INTERFACE_krino_mesh_clone
    INTERFACE_krino_levelset_shape_sens
    INTERFACE_krino_plato_workarounds
    INTERFACE_krino_quality_metric_sens
    INTERFACE_krino_refinement_rebalance
    INTERFACE_krino_rol
    INTERFACE_krino_smoothing
)

krino_add_library(
    krino_rebalance_utils
    krino/rebalance_utils
    PUBLIC_krino_decomp
    PUBLIC_stk::stk_balance_lib
    PUBLIC_stk::stk_mesh_base
)

krino_add_library(
    krino_region
    krino/region
    PUBLIC_krino_decomp
    PUBLIC_krino_post_process
)

krino_add_library(
    krino_sierra
    krino_sierra
    PUBLIC_krino_decomp
    PUBLIC_krino_post_process
    PUBLIC_krino_region
    PUBLIC_krino_parser
    PUBLIC_krino_surface
    PUBLIC_krino_mesh_surface
    PUBLIC_sierra_common::sierra
    PUBLIC_sierra_common::sierra_util_sctl
    PUBLIC_sierra_common::sierra_util_user_input_function_parser
    PUBLIC_sierra_common::sierraparser
)

krino_add_library(
    krino_parser
    krino/parser
    PUBLIC_krino_region
    PUBLIC_krino_decomp
    PUBLIC_krino_post_process
    PUBLIC_krino_surface
    PUBLIC_krino_mesh_surface
    PUBLIC_yaml-cpp
)

krino_add_library(
    krino_unit_main
    krino/unit_main
    PUBLIC_SEACASIoss::Ioss
    PUBLIC_krino_diagwriter
    PUBLIC_stk::stk_unit_test_utils
)

krino_add_library(
    krino_unit_test_utils
    krino/unit_test_utils
    PUBLIC_krino_decomp
    PUBLIC_krino_refinement
    PUBLIC_krino_mesh_interface
    PUBLIC_krino_mesh_quality
    PUBLIC_krino_unit_main
)

krino_add_library(
    mesh_adapt_lib
    krino_mesh_adapt/mesh_adapt_lib
    PUBLIC_krino_refinement
    PUBLIC_stk::stk_io
    PUBLIC_stk::stk_mesh_base
    PUBLIC_stk::stk_tools_lib
    PUBLIC_stk::stk_util_diag
)

# ---------------------------------------------------------------------------
# Main executables
# ---------------------------------------------------------------------------

krino_add_executable(
    krino
    krino/Apps_krino.cpp
    PUBLIC_krino_parser
    PUBLIC_krino_region
    PUBLIC_stk::stk_util_registry
)

krino_add_executable(
    delete_small_elements
    delete_small_elements/Akri_DeleteSmallElementsMain.cpp
    PUBLIC_krino_mesh_utils
    PUBLIC_stk::stk_util_command_line
)

krino_add_executable(
    krino_mesh_adapt
    krino_mesh_adapt/KrinoMeshAdaptMain.cpp
    PUBLIC_MPI::MPI_C
    PUBLIC_mesh_adapt_lib
    PUBLIC_stk::stk_util_command_line
)

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

if (BUILD_TESTS)
    krino_add_executable(
        krino_decomp_utest
        "krino/decomp/unit_tests/*.cpp"
        PUBLIC_krino_decomp
        PUBLIC_krino_rebalance_utils
        PUBLIC_krino_mesh_clone
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/decomp/unit_tests
    )

    krino_add_executable(
        krino_levelset_utest
        "krino/levelset/unit_tests/*.cpp"
        PUBLIC_krino_levelset
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/levelset/unit_tests
    )

    krino_add_executable(
        krino_interface_geometry_utest
        "krino/interface_geometry/unit_tests/*.cpp"
        PUBLIC_krino_interface_geometry
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/interface_geometry/unit_tests
    )

    krino_add_executable(
        krino_unit_test_utils_utest
        "krino/unit_test_utils/unit_tests/*.cpp"
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils/unit_tests
    )

    krino_add_executable(
        krino_parser_utest
        "krino/parser/unit_tests/*.cpp"
        PUBLIC_krino_parser
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/parser/unit_tests
    )

    krino_add_executable(
        krino_math_utils_utest
        "krino/math_utils/unit_tests/*.cpp"
        PUBLIC_krino_math_utils
        PUBLIC_krino_rol
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/math_utils/unit_tests
    )

    krino_add_executable(
        krino_geometry_utest
        "krino/geometry/unit_tests/*.cpp"
        PUBLIC_krino_geometry
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/geometry/unit_tests
    )

    krino_add_executable(
        krino_eikonal_utest
        "krino/eikonal/unit_tests/*.cpp"
        PUBLIC_krino_eikonal
        PUBLIC_krino_post_process
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/eikonal/unit_tests
    )

    krino_add_executable(
        krino_contour_utest
        "krino/contour/unit_tests/*.cpp"
        PUBLIC_krino_contour
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/contour/unit_tests
    )

    krino_add_executable(
        krino_surface_utest
        "krino/surface/unit_tests/*.cpp"
        PUBLIC_krino_surface
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/surface/unit_tests
    )

    krino_add_executable(
        krino_mesh_surface_utest
        "krino/mesh_surface/unit_tests/*.cpp"
        PUBLIC_krino_mesh_surface
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_surface/unit_tests
    )

    krino_add_executable(
        krino_mesh_interface_utest
        "krino/mesh_interface/unit_tests/*.cpp"
        PUBLIC_krino_mesh_interface
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_interface/unit_tests
    )

    krino_add_executable(
        krino_levelset_shape_sens_utest
        "krino/levelset_shape_sens/unit_tests/*.cpp"
        PUBLIC_krino_levelset_shape_sens
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/levelset_shape_sens/unit_tests
    )

    krino_add_executable(
        krino_mesh_clone_utest
        "krino/mesh_clone/unit_tests/*.cpp"
        PUBLIC_krino_mesh_clone
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_clone/unit_tests
    )

    krino_add_executable(
        krino_mesh_utils_utest
        "krino/mesh_utils/unit_tests/*.cpp"
        PUBLIC_krino_mesh_utils
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_utils/unit_tests
    )

    krino_add_executable(
        krino_mesh_quality_utest
        "krino/mesh_quality/unit_tests/*.cpp"
        PUBLIC_krino_mesh_quality
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_quality/unit_tests
    )

    krino_add_executable(
        krino_parametric_element_utest
        "krino/parametric_element/unit_tests/*.cpp"
        PUBLIC_krino_parametric_element
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/parametric_element/unit_tests
    )

    krino_add_executable(
        krino_phase_utest
        "krino/phase/unit_tests/*.cpp"
        PUBLIC_krino_phase
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/phase/unit_tests
    )

    krino_add_executable(
        krino_quality_metric_utest
        "krino/quality_metric/unit_tests/*.cpp"
        PUBLIC_krino_quality_metric
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/quality_metric/unit_tests
    )

    krino_add_executable(
        krino_quality_metric_sens_utest
        "krino/quality_metric_sens/unit_tests/*.cpp"
        PUBLIC_krino_quality_metric_sens
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/quality_metric_sens/unit_tests
    )

    krino_add_executable(
        krino_rebalance_utils_utest
        "krino/rebalance_utils/unit_tests/*.cpp"
        PUBLIC_krino_rebalance_utils
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/rebalance_utils/unit_tests
    )

    krino_add_executable(
        krino_refinement_utest
        "krino/refinement/unit_tests/*.cpp"
        PUBLIC_krino_refinement
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/refinement/unit_tests
    )

    krino_add_executable(
        krino_refinement_drivers_utest
        "krino/refinement_drivers/unit_tests/*.cpp"
        PUBLIC_krino_refinement_drivers
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/refinement_drivers/unit_tests
    )

    krino_add_executable(
        krino_mesh_metric_utest
        "krino/mesh_metric/unit_tests/*.cpp"
        PUBLIC_krino_mesh_metric
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/mesh_metric/unit_tests
    )

    krino_add_executable(
        krino_metric_utest
        "krino/metric/unit_tests/*.cpp"
        PUBLIC_krino_metric
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/metric/unit_tests
    )

    krino_add_executable(
        krino_smoothing_utest
        "krino/smoothing/unit_tests/*.cpp"
        PUBLIC_krino_smoothing
        PUBLIC_krino_rol
        PUBLIC_krino_unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/smoothing/unit_tests
    )

    # The krino_sierra parser unit tests exercise the SIERRA-framework parser,
    # so they need the command database compiled in.  Build a unit-test command
    # database that concatenates the unit-test-only XML with the main krino
    # command database (mirrors the bake krino_unit_sierra_commands target), then
    # link the resulting object library into the test executable.
    add_parser_commands(
        TARGET krino_unit_sierra_commands
        XML_FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra/unit_tests/unit_test.xml
            ${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra/xml/Akri_Levelset.xml
    )

    krino_add_executable(
        krino_sierra_utest
        "krino_sierra/unit_tests/*.cpp"
        PUBLIC_krino_sierra
        PUBLIC_krino_unit_test_utils
        PUBLIC_krino_unit_sierra_commands.o
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino/unit_test_utils
        INCLUDE_${CMAKE_CURRENT_SOURCE_DIR}/krino_sierra/unit_tests
    )
endif()

# ---------------------------------------------------------------------------
# Package export
# ---------------------------------------------------------------------------

install(
    EXPORT krinoTargets
    NAMESPACE krino::
    DESTINATION share/cmake/krino
)

install(
    FILES cmake/krinoConfig.cmake
    DESTINATION share/cmake/krino
    )
