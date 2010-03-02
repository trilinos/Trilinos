/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_config_h
#define stk_util_config_h

#include <STK_config.h>

#define STK_PACKAGE stk
#define STK_VERSION 0.1a
#define STK_HAS_SNL_EXODUSII

#ifdef HAVE_MPI
#define STK_HAS_MPI
#endif

#endif /* stk_util_config_h */
