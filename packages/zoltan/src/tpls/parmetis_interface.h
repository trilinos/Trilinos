/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#ifndef __PARMETIS_INTERFACE_H
#define __PARMETIS_INTERFACE_H

#include <limits.h>
#include "zoltan_comm.h"
#include "third_library_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Guess the version number of ParMetis if not defined */
/* PARMETIS_MAJOR_VERSION is only defined in version 3.0 and higher */
#if (!defined(PARMETIS_MAJOR_VERSION))
#define PARMETIS_MAJOR_VERSION 0
#define PARMETIS_MINOR_VERSION 0
#endif

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 * Both ParMetis 2.0 and 3.0 defs are included below.
 */
#define PMV3_OPTION_DBGLVL      1
#define PMV3_OPTION_SEED        2
#define PMV3_OPTION_IPART       3
#define PMV3_OPTION_PSR         3
#define PMV3_OPT_USE_OBJ_SIZE   9  /* Added by EB, not in ParMetis */
#define MAX_OPTIONS             10 /* Max number of options +1 */
/* Other ParMetis constants we may need */
#define GLOBAL_DBGLVL		0  /* Default debug level */
#define GLOBAL_SEED		15 /* Default random seed */
#define COUPLED                 1  /* Processors coupled to partitions? */
#define DISCOUPLED              2  /* Processors coupled to partitions? */


int Zoltan_ParMetis(
  ZZ *, float *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
