# Crated by mbenlioglu on Jul 12, 2020.
# -------------------------------------------------------------

find_package(SARMA QUIET)
get_target_property(TPL_SARMA_INCLUDE_DIRS SARMA::libsarma INTERFACE_INCLUDE_DIRECTORIES)
list(APPEND TPL_SARMA_INCLUDE_DIRS "${TPL_SARMA_INCLUDE_DIRS}/data_structures" "${TPL_SARMA_INCLUDE_DIRS}/algorithms" "${TPL_SARMA_INCLUDE_DIRS}/tools")

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES(SARMA
        REQUIRED_HEADERS "sarma.hpp;algorithm.hpp;mixed_integer_program.hpp;nicol1d.hpp;nicol2d.hpp;ordered_probe_a_load.hpp;patoh_cpm.hpp;probe_a_load.hpp;refine_a_cut.hpp;subgradient_method.hpp;uniform.hpp;csr_matrix.hpp;sparse_prefix_sum.hpp;progress_bar.hpp;timer.hpp;utils.hpp"
        MUST_FIND_ALL_HEADERS
        REQUIRED_LIBS_NAMES "mmio"
        )

SET(CMAKE_REQUIRED_INCLUDES ${TPL_SARMA_INCLUDE_DIRS})
SET(CMAKE_REQUIRED_LIBRARIES ${TPL_SARMA_LIBRARIES})
SET(CMAKE_REQUIRED_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
