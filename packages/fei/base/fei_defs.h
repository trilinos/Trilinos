/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_defs_h_
#define _fei_defs_h_

/*
   In this file we set some #defines to use as parameters to
   some fei functions, and also some error-code returns.
   We also provide the typedef for 'GlobalID' which appears in
   many FEI function prototypes. Note that the default case is
   for GlobalID to simply be an int.
   This file is included by both C and C++ versions of the fei.
*/

#ifdef EIGHT_BYTE_GLOBAL_ID
    typedef long long   GlobalID;
    #define GlobalID_MAX LLONG_MAX
    #define GlobalID_MIN LLONG_MIN
#else
    typedef int GlobalID;
#endif


/* solveType (used in 'setSolveType'): */
#define FEI_SINGLE_SYSTEM     0
#define FEI_EIGEN_SOLVE       1
#define FEI_AGGREGATE_SUM     2
#define FEI_AGGREGATE_PRODUCT 3

/* IDType (used in coefficient-access functions) */
#define FEI_NODE 0
#define FEI_ELEMENT 1
#define FEI_ONLY_NODES 2
#define FEI_ONLY_ELEMENTS 3

/* elemFormat (used in 'sumInElem' and 'sumInElemMatrix'): */
#define FEI_DENSE_ROW      0
#define FEI_UPPER_SYMM_ROW 1
#define FEI_LOWER_SYMM_ROW 2
#define FEI_DENSE_COL      3
#define FEI_UPPER_SYMM_COL 4
#define FEI_LOWER_SYMM_COL 5
#define FEI_DIAGONAL       6
#define FEI_BLOCK_DIAGONAL_ROW 7
#define FEI_BLOCK_DIAGONAL_COL 8


/* interleaveStrategy (used in initElemBlock): */
#define FEI_NODE_MAJOR  0
#define FEI_FIELD_MAJOR 1


/* timingMode (used in cumulative_MPI_Wtimes): */
#define FEI_LOCAL_TIMES 0
#define FEI_MAX_TIMES   1
#define FEI_MIN_TIMES   2

/* FEI function return values */
#define FEI_SUCCESS         0
#define FEI_FATAL_ERROR    -1
#define FEI_ID_NOT_FOUND   -2

#endif

