/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#ifndef _PE_COMMON_H
#define _PE_COMMON_H

/*****************************************************************************/
/* This file contains defines that are common to both nem_spread, and        */
/* nem_join. Most of them were taken from rf_salsa.h which is unique         */
/* to nem_spread.                                                            */
/*****************************************************************************/

/*
 * Default value of the chunk size (in bytes) for use in I/O and message
 * passing
 */

#ifndef MAX_CHUNK_SIZE
#define MAX_CHUNK_SIZE 1073741824
/* Small message size for debugging purposes */
/* #define MAX_CHUNK_SIZE 16384 */
#endif

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

/* Redefine function "exit" for OSF and Puma so that */
/* the process is killed on all nodes when a single  */
/* node exits.  See routines osf_exit and smos_exit  */
/* in md_wrap_intel_c.c and md_wrap_smos_c.c,        */
/* respectively.                                     */
#ifdef PARA
#define exit(a) osf_exit(a, __FILE__, __LINE__)
#endif
#ifdef SMOS
#define exit(a) smos_exit(a, __FILE__, __LINE__)
#endif

# define PEX_MAX(x,y) (( x > y ) ? x : y)     /* max function */
# define PEX_MIN(x,y) (( x < y ) ? x : y)     /* min function */

#endif /* _PE_COMMON_H */
