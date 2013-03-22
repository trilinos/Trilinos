/*
 * Copyright (c) 2012 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
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
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

/*****************************************************************************
 *
 * exodusII_par.h - Exodus II parallel-aware API include file
 *
 *****************************************************************************/

#ifndef EXODUS_II_PAR_HDR
#define EXODUS_II_PAR_HDR

#if !defined(PARALLEL_NETCDF)

#include "exodusII.h"

/*
 * need following extern if this include file is used in a C++
 * program, to keep the C++ compiler from mangling the function names.
 */
#ifdef __cplusplus
extern "C" {
#endif

#define ex_open_par(path, mode, comp_ws, io_ws, version, comm, info) ex_open_par_int(path, mode, comp_ws, io_ws, version, comm, info, EX_API_VERS_NODOT)  
#define ex_create_par(path, mode, comp_ws, io_ws, comm, info) ex_create_par_int(path, mode, comp_ws, io_ws, comm, info, EX_API_VERS_NODOT)  

EXODUS_EXPORT int ex_open_par_int (const char  *path,
				   int    mode,
				   int   *comp_ws,
				   int   *io_ws,
				   float *version,
				   MPI_Comm comm,
				   MPI_Info info,
				   int my_version);

EXODUS_EXPORT int ex_create_par_int (const char *path, int cmode, int *comp_ws, int *io_ws,
				     MPI_Comm comm,
				     MPI_Info info,
				     int my_version);

#else
#error "Parallel-aware exodusII_par.h included in non-parallel context"
#endif

#ifdef __cplusplus
}                               /* close brackets on extern "C" declaration */
#endif

#endif

