/*
//@HEADER
// ************************************************************************
//
//                        Adelus v. 1.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// 3. Neither the name of NTESS nor the names of the contributors may be
// used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL NTESS OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Vinh Dang (vqdang@sandia.gov)
//                    Joseph Kotulski (jdkotul@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef __ADELUS_DEFINES_H__
#define __ADELUS_DEFINES_H__

#include "Adelus_config.h"

#undef CBLAS
#undef DEBUG
#define OVERLAP

//  ------------------ Define Data Types --------------------------
//  ------------------ Define Constants and Operations ------------
#ifdef SREAL

  #define ADELUS_DATA_TYPE float
  #undef ADELUS_COMPLEX
  #define ADELUS_MPI_DATA_TYPE MPI_FLOAT
  #define ADELUS_MPI_DATA_TYPE2 MPI_FLOAT_INT
  
  #define ADELUS_CONST_ONE 1.0
  #define ADELUS_CONST_MINUS_ONE -1.0
  #define ADELUS_CONST_ZERO 0.0
  #define ADELUS_NEGATIVE(X,Y) (Y)=-(X)
  #define ADELUS_ABS_VAL(X) (fabs(X))
  #define ADELUS_INVERSE(X,W,Z) (Z)=CONST_ONE/(X)
  #define ADELUS_MULTIPLY(X,Y,Z) (Z)=(X)*(Y)
  #define ADELUS_DIVIDE(X,Y,W,Z) (Z)=(X)/(Y)

#elif defined(DREAL)

  #define ADELUS_DATA_TYPE double
  #undef ADELUS_COMPLEX
  #define ADELUS_MPI_DATA_TYPE MPI_DOUBLE
  #define ADELUS_MPI_DATA_TYPE2 MPI_DOUBLE_INT
  
  #define ADELUS_CONST_ONE 1.0
  #define ADELUS_CONST_MINUS_ONE -1.0
  #define ADELUS_CONST_ZERO 0.0
  #define ADELUS_NEGATIVE(X,Y) (Y)=-(X)
  #define ADELUS_ABS_VAL(X) (fabs(X))
  #define ADELUS_INVERSE(X,W,Z) (Z)=CONST_ONE/(X)
  #define ADELUS_MULTIPLY(X,Y,Z) (Z)=(X)*(Y)
  #define ADELUS_DIVIDE(X,Y,W,Z) (Z)=(X)/(Y)

#elif defined(SCPLX)

  typedef struct {
    float r;
    float i;
  } scomplex ;
  #define ADELUS_DATA_TYPE scomplex
  #define ADELUS_COMPLEX
  #define ADELUS_MPI_DATA_TYPE MPI_COMPLEX
  #define ADELUS_MPI_DATA_TYPE2 MPI_FLOAT_INT
  
  #define ADELUS_CONST_ONE {1.0, 0.0}
  #define ADELUS_CONST_MINUS_ONE {-1.0, 0.0}
  #define ADELUS_CONST_ZERO {0.0, 0.0}
  #define ADELUS_NEGATIVE(X,Y) (Y).r=-(X).r;(Y).i=-(X).i
  #define ADELUS_ABS_VAL(X) ((X).r * (X).r + (X).i * (X).i)
  #define ADELUS_INVERSE(X,W,Z) (Z).r=(X).r/(W);(Z).i=-(X).i/(W)
  #define ADELUS_MULTIPLY(X,Y,Z) (Z).r=(X).r*(Y).r-(X).i*(Y).i;(Z).i=(X).r*(Y).i+(X).i*(Y).r
  #define ADELUS_DIVIDE(X,Y,W,Z) (Z).r=((X).r*(Y).r+(X).i*(Y).i)/(W);(Z).i=((X).i*(Y).r-(X).r*(Y).i)/(W)

#else//ZCPLX

  typedef struct {
    double r;
    double i;
  } dcomplex ;
  #define ADELUS_DATA_TYPE dcomplex
  #define ADELUS_COMPLEX
  #define ADELUS_MPI_DATA_TYPE MPI_DOUBLE_COMPLEX
  #define ADELUS_MPI_DATA_TYPE2 MPI_DOUBLE_INT
  
  #define ADELUS_CONST_ONE {1.0, 0.0}
  #define ADELUS_CONST_MINUS_ONE {-1.0, 0.0}
  #define ADELUS_CONST_ZERO {0.0, 0.0}
  #define ADELUS_NEGATIVE(X,Y) (Y).r=-(X).r;(Y).i=-(X).i
  #define ADELUS_ABS_VAL(X) ((X).r * (X).r + (X).i * (X).i)
  #define ADELUS_INVERSE(X,W,Z) (Z).r=(X).r/(W);(Z).i=-(X).i/(W)
  #define ADELUS_MULTIPLY(X,Y,Z) (Z).r=(X).r*(Y).r-(X).i*(Y).i;(Z).i=(X).r*(Y).i+(X).i*(Y).r
  #define ADELUS_DIVIDE(X,Y,W,Z) (Z).r=((X).r*(Y).r+(X).i*(Y).i)/(W);(Z).i=((X).i*(Y).r-(X).r*(Y).i)/(W)
#endif

#endif
