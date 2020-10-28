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


#ifndef __THIRD_LIBRARY_CONST_PRIDX_H
#define __THIRD_LIBRARY_CONST_PRIDX_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdint.h>
#include <sys/types.h>
#include <float.h>
#include "zoltan_util.h"

/*****************************************************************************/
/* Include appropriate files for TPLs */
#ifdef ZOLTAN_METIS
  #include "metis.h"
  #define __metis__ 1
#else
  #define __metis__ 0
#endif

#ifdef ZOLTAN_PARMETIS
  #include "parmetis.h"
  #define __parmetis__ 1
#else
  #define __parmetis__ 0
#endif

#ifdef ZOLTAN_PTSCOTCH
  #include "ptscotch.h"
  #define __ptscotch__ 1
#else
  #define __ptscotch__ 0
#endif

#ifdef ZOLTAN_SCOTCH
  #ifndef ZOLTAN_PTSCOTCH
    #include "scotch.h"
  #endif
  #define __scotch__ 1
#else
  #define __scotch__ 0
#endif

/****************************************************************************/
/*  TPL-Specific settings for data types                                    */
/*
 * "indextype" is the type used for global numbers and indices in the graph
 *  data structure.
 * "weighttype" is the type used for weights.
 * 
 * If there are no third party graph/ordering libraries, let indextype be
 * ZOLTAN_GNO_TYPE and let "weighttype" be float.
 *
 * If there is only one third party library used for graph algorithms,
 * define indextype and weighttype to match the types used by that library.
 *
 * If more than one library is linked in, arbitrarily choose one.
 * Check for compatibility between the libraries here; all should use the same
 * size integer for indices.
 *
 * At runtime, if
 * either the indextype or weighttype is not compatible with the graph library
 * API, return an error.
 */

#undef TPL_SCOTCH_DATATYPES
#undef TPL_METIS_DATATYPES
#undef TPL_ZOLTAN_DATATYPES

#define TPL_SCOTCH_DATATYPES   1
#define TPL_METIS_DATATYPES    2
#define TPL_ZOLTAN_DATATYPES   3

#undef TPL_USE_DATATYPE

/* Select the data types to use */
#if __parmetis__ + __metis__ + __ptscotch__ + __scotch__ == 0
  /* No graph TPLs used; use Zoltan values */
  #define TPL_USE_DATATYPE TPL_ZOLTAN_DATATYPES
#elif (__ptscotch__ + __scotch__ > 0) && (__parmetis__ + __metis__  == 0)
  /* Using only Scotch/PTScotch */
  #define TPL_USE_DATATYPE TPL_SCOTCH_DATATYPES
#elif (__parmetis__ + __metis__ > 0) && (__ptscotch__ + __scotch__ == 0)
  /* Using only METIS/ParMETIS */
  #define TPL_USE_DATATYPE TPL_METIS_DATATYPES

#else
  /* Using both METIS/ParMETIS and Scotch/PTScotch; let METIS datatypes rule */
  #define TPL_USE_DATATYPE TPL_METIS_DATATYPES
#endif

#if TPL_USE_DATATYPE == TPL_METIS_DATATYPES
  #if (PARMETIS_MAJOR_VERSION == 4 || METIS_VER_MAJOR == 5)
    #define TPL_IDX_SPEC "%"PRIDX
    #define TPL_WGT_SPEC "%"PRIDX
  #endif
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
