// @HEADER
// *****************************************************************************
//                   Basker: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the Basker contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER


#define Int int
#define BASKER(name) basker_ ## name

#include "basker.c"

#undef Int
#undef BASKER
#define Int long
#define BASKER(name) basker_ ## name ## _l

#include "basker.c"
