// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* This file must be C-compatible. */

#ifndef token_enumH
#define token_enumH

typedef enum
{
  TK_EXIT       =   0,
  TK_LP         = '(',
  TK_RP         = ')',
  TK_PER        = '/',
  TK_CARET      = '^',
  TK_LT         = '<',
  TK_GT         = '>',
  TK_PLUS       = '+',
  TK_NONE,
  TK_IDENTIFIER = 257,
  TK_INTEGER,
  TK_REAL,
  TK_END,
  TK_STRING,
  TK_ERROR
}
Token_Type;

#endif
