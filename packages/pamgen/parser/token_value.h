/* $Id$ */

/* This file must be C-compatible. */

#ifndef token_valueH
#define token_valueH

#include "../asrc/code_types.h"

typedef union 
{
  char        cval;
  Real        fval;
  int         ival;
  char       *sval; 
}
Token_Value;

#endif

