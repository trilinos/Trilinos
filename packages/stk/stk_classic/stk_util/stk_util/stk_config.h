/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_config_h
#define stk_util_config_h

#ifdef STK_BUILT_IN_SIERRA
#define STK_HAS_MPI
#else
// This file gets created by cmake during a Trilinos build
// and will not be present in a sierra build using bjam or associated wrappers
#include <STKClassic_config.h>
#ifdef HAVE_MPI
#define STK_HAS_MPI
#endif
#endif

#define STK_PACKAGE stk
#define STK_HAS_SNL_EXODUSII

#endif /* stk_util_config_h */
