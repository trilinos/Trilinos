// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* This file must be C-compatible. */

#ifndef token_valueH
#define token_valueH

#include "pamgen_code_types.h"

typedef union 
{
  char        cval;
  Real        fval;
  int         ival;
  char       *sval; 
}
Token_Value;

#endif
