//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
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
