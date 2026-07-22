// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This is a dummy utility to keep the build system happy for header only
// subpackages (these will be the only symbols in the library)

//this creates two local symbols without any "unused function" warnings
namespace _dummy
{
static void fastilu_dummy1();
static void fastilu_dummy2();
void fastilu_dummy1() {fastilu_dummy2();}
void fastilu_dummy2() {fastilu_dummy1();}
}

