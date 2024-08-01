// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#pragma once

template <int N, typename ExecSpace>
double time_dfad_hierarchical_team(int ncells, int num_basis, int num_points,
                                   int ndim, int ntrial, bool check);

template <int N, typename ExecSpace>
double time_dfad_hierarchical_team_scratch(int ncells, int num_basis,
                                           int num_points, int ndim, int ntrial,
                                           bool check);
