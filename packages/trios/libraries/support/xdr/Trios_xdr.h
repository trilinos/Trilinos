/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/*-------------------------------------------------------------------------*/
/**
 *   @file Trios_xdr.h
 *
 *   @brief Essential definitions and include files for code that uses XDR.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 */

#ifndef _TRIOS_XDR_H_
#define _TRIOS_XDR_H_

#include "Trios_config.h"
#include <stdint.h>
#include <sys/param.h>
#include <errno.h>

#ifdef __MTA__
#include <xdr.h>
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

/* Some systems do not have xdr functions for uint16_t, uint32_t, uint64_t.
 * We should have checked for these in the cmake config.  To fix the problem,
 * we redefine these function names.
 */
#ifdef HAVE_TRIOS_XDR_U_INT16_T
#define xdr_uint16_t xdr_u_int16_t
#endif
#ifdef HAVE_TRIOS_XDR_U_INT32_T
#define xdr_uint32_t xdr_u_int32_t
#endif
#ifdef HAVE_TRIOS_XDR_U_INT64_T
#define xdr_uint64_t xdr_u_int64_t
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)


/* On some systems (mac) xdr_sizeof is not in header file,
 * but it is in the library.
 */
#ifndef HAVE_TRIOS_XDR_SIZEOF
extern unsigned long xdr_sizeof (xdrproc_t func, void *data);
#endif


#else /* K&R C */
#endif



#ifdef __cplusplus
}
#endif

#endif

