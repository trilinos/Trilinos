// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#pragma once

template <typename FadType, int N, typename ExecSpace>
double time_fad_hierarchical_flat(int ncells, int num_basis, int num_points,
                                  int ndim, int ntrial, bool check);

template <typename FadType, int N, typename ExecSpace>
double time_fad_hierarchical_team(int ncells, int num_basis, int num_points,
                                  int ndim, int ntrial, bool check);
