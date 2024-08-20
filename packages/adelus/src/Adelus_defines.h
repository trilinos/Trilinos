/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef __ADELUS_DEFINES_H__
#define __ADELUS_DEFINES_H__

#include "Adelus_config.h"

#undef CBLAS
#undef DEBUG
#define OVERLAP

//#define ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST //NOTE: for perf comparison only

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
