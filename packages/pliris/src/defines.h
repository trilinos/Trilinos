/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef __DEFINESH__
#define __DEFINESH__


#undef CBLAS
#undef DEBUG
#define  OVERLAP
#undef  PRINT_STATUS
#undef   TIMING0 


#ifdef DREAL
#define MPI_DATA_TYPE MPI_DOUBLE
#define MPI_DATA_TYPE2 MPI_DOUBLE_INT
#endif
#ifdef SREAL
#define MPI_DATA_TYPE MPI_FLOAT
#define MPI_DATA_TYPE2 MPI_FLOAT_INT
#endif  

/*  ------------------ Define Data Types -------------------------- */
#ifdef ZCPLX
    typedef struct {
      double r;
      double i;
    } dcomplex ;
#define DATA_TYPE dcomplex
#define COMPLEX
#define MPI_DATA_TYPE MPI_DOUBLE
#define MPI_DATA_TYPE2 MPI_DOUBLE_INT
#endif

#ifdef SCPLX
    typedef struct {
      float r;
      float i;
    } scomplex ;
#define DATA_TYPE scomplex
#define COMPLEX
#define MPI_DATA_TYPE  MPI_FLOAT
#define MPI_DATA_TYPE2 MPI_FLOAT_INT
#endif

#ifdef DREAL

#define DATA_TYPE double
#undef COMPLEX
#endif
 
#ifdef SREAL

#define DATA_TYPE float
#undef COMPLEX 
#endif


/*  ------------------ Define Constants and Operations ------------ */

#ifdef ZCPLX 

#define CONST_ONE {1.0, 0.0}
#define CONST_MINUS_ONE {-1.0, 0.0}
#define CONST_ZERO {0.0, 0.0}
#define NEGATIVE(X,Y) (Y).r=-(X).r;(Y).i=-(X).i
#define ABS_VAL(X) ((X).r * (X).r + (X).i * (X).i)
#define INVERSE(X,W,Z) (Z).r=(X).r/(W);(Z).i=-(X).i/(W)
#define MULTIPLY(X,Y,Z) (Z).r=(X).r*(Y).r-(X).i*(Y).i;(Z).i=(X).r*(Y).i+(X).i*(Y).r
#define DIVIDE(X,Y,W,Z) (Z).r=((X).r*(Y).r+(X).i*(Y).i)/(W);(Z).i=((X).i*(Y).r-(X).r*(Y).i)/(W)

#endif

#ifdef  SCPLX

#define CONST_ONE {1.0, 0.0}
#define CONST_MINUS_ONE {-1.0, 0.0}
#define CONST_ZERO {0.0, 0.0}
#define NEGATIVE(X,Y) (Y).r=-(X).r;(Y).i=-(X).i
#define ABS_VAL(X) ((X).r * (X).r + (X).i * (X).i)
#define INVERSE(X,W,Z) (Z).r=(X).r/(W);(Z).i=-(X).i/(W)
#define MULTIPLY(X,Y,Z) (Z).r=(X).r*(Y).r-(X).i*(Y).i;(Z).i=(X).r*(Y).i+(X).i*(Y).r
#define DIVIDE(X,Y,W,Z) (Z).r=((X).r*(Y).r+(X).i*(Y).i)/(W);(Z).i=((X).i*(Y).r-(X).r*(Y).i)/(W)

#endif


#ifdef DREAL 

#define CONST_ONE 1.0
#define CONST_MINUS_ONE -1.0
#define CONST_ZERO 0.0
#define NEGATIVE(X,Y) (Y)=-(X)
#define ABS_VAL(X) (fabs(X))
#define INVERSE(X,W,Z) (Z)=CONST_ONE/(X)
#define MULTIPLY(X,Y,Z) (Z)=(X)*(Y)
#define DIVIDE(X,Y,W,Z) (Z)=(X)/(Y)

#endif

#ifdef SREAL 

#define CONST_ONE 1.0
#define CONST_MINUS_ONE -1.0
#define CONST_ZERO 0.0
#define NEGATIVE(X,Y) (Y)=-(X)
#define ABS_VAL(X) (fabs(X))
#define INVERSE(X,W,Z) (Z)=CONST_ONE/(X)
#define MULTIPLY(X,Y,Z) (Z)=(X)*(Y)
#define DIVIDE(X,Y,W,Z) (Z)=(X)/(Y)

#endif

/*  ------------------- Define BLAS prototype definitions ------------  */

#ifdef ZCPLX
#include "zblassp.h"
#define XGEMM_  zgemm_
#define XGEMMS_ zgemm_
#define XGEMM  zgemm
#define XLU_SOLVE_ zlusolve_
#define X_SOLVE_ zsolve_
#define X_FACTOR_ zfactor_
#define X_PERMUTE_ zpermute_
#endif

#ifdef SCPLX
#include "cblassp.h"
#define XGEMM_  cgemm_
#define XGEMM  cgemm
#define XLU_SOLVE_ clusolve_
#define X_SOLVE_ csolve_
#define X_FACTOR_ cfactor_
#define X_PERMUTE_ cpermute_
#endif

#ifdef DREAL
#include "dblassp.h"
#define XGEMM_  dgemm_
#define XGEMM  dgemm
#define XGEMMS_ dgemm_
#define XLU_SOLVE_ dlusolve_
#define X_SOLVE_ dsolve_
#define X_FACTOR_ dfactor_
#define X_PERMUTE_ dpermute_
#endif

#ifdef SREAL
#include "sblassp.h"
#define XGEMM_  sgemm_
#define XGEMM  sgemm
#define XGEMMS_  sgemm_
#define XLU_SOLVE_ slusolve_
#define X_SOLVE_ ssolve_
#define X_FACTOR_ sfactor_
#define X_PERMUTE_ spermute_
#endif

#endif
