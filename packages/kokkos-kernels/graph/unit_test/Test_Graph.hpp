// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef TEST_GRAPH_HPP
#define TEST_GRAPH_HPP

#include "Test_Graph_graph_color_deterministic.hpp"
#include "Test_Graph_graph_color_distance2.hpp"
#include "Test_Graph_graph_color.hpp"
#include "Test_Graph_mis2.hpp"
#if !defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#include "Test_Graph_coarsen.hpp"
#endif
#include "Test_Graph_rcm.hpp"
#include "Test_Graph_rcb.hpp"

#endif  // TEST_GRAPH_HPP
