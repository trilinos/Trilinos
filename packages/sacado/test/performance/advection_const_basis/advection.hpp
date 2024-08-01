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
double time_fad_flat(int ncells, int num_basis, int num_points, int ndim,
                     int ntrial, bool check);

template <typename FadType, int N, typename ExecSpace>
double time_fad_scratch(int ncells, int num_basis, int num_points, int ndim,
                        int ntrial, bool check);

template <int N, typename ExecSpace>
double time_analytic_flat(int ncells, int num_basis, int num_points, int ndim,
                          int ntrial, bool check);

template <int N, typename ExecSpace>
double time_analytic_const(int ncells, int num_basis, int num_points, int ndim,
                           int ntrial, bool check);

template <int N, typename ExecSpace>
double time_analytic_team(int ncells, int num_basis, int num_points, int ndim,
                          int ntrial, bool check);
