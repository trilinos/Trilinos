// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef _TEUCHOS_SCALARTRAITS_CUDA_HPP_
#define _TEUCHOS_SCALARTRAITS_CUDA_HPP_

/*! \file Teuchos_ScalarTraitsCUDA.hpp
    \brief Defines basic traits for the scalar field type, appropriate for compilation under the NVIDIA CUDA C compiler.
*/
 
#include "Teuchos_ScalarTraitsDecl.hpp"

namespace Teuchos {

template<>
struct ScalarTraits<int>
{
  typedef int magnitudeType;
  typedef int halfPrecision;
  typedef int doublePrecision;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline __device__ __host__ magnitudeType magnitude(int a) { return (int)fabsf((float)a); }
  static inline __device__ __host__ int zero()  { return 0; }
  static inline __device__ __host__ int one()   { return 1; }
  static inline __device__ __host__ int conjugate(int x) { return x; }
  static inline __device__ __host__ int real(int x) { return x; }
  static inline __device__ __host__ int imag(int)   { return 0; }
  static inline __device__ __host__ bool isnaninf(int) { return false; }
  static inline __device__ __host__ int squareroot(int x) { return (int)sqrtf((float)x); }          // perhaps this cast should be replaced by an explicit call like __float2int_rn
  static inline __device__ __host__ int pow(int x, int y) { return (int)powf((float)x,(float)y); }  // perhaps this cast should be replaced by an explicit call like __float2int_rn
};

#ifdef HAVE_KOKKOS_CUDA_FLOAT
template<>
struct ScalarTraits<float>
{
  typedef float magnitudeType;
  typedef float halfPrecision; // should become IEEE754-2008 binary16 or fp16 later, perhaps specified at configure according to architectural support
#ifdef HAVE_KOKKOS_CUDA_DOUBLE
  typedef double doublePrecision;
#else
  typedef float  doublePrecision;
#endif
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline __device__ __host__ magnitudeType magnitude(float a) { return fabsf(a); }
  static inline __device__ __host__ float zero()  { return(0.0f); }
  static inline __device__ __host__ float one()   { return(1.0f); }    
  static inline __device__ __host__ float conjugate(float x)   { return(x); }    
  static inline __device__ __host__ float real(float x) { return x; }
  static inline __device__ __host__ float imag(float)   { return zero(); }
  static inline __device__ __host__ bool  isnaninf(float x) { return isnan(x) || isinf(x); }
  static inline __device__ __host__ float squareroot(float x) { return sqrtf(x); }
  static inline __device__ __host__ float pow(float x, float y) { return powf(x,y); }
};
#endif // HAVE_KOKKOS_CUDA_FLOAT

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
template<>
struct ScalarTraits<double>
{
  typedef double magnitudeType;
#ifdef HAVE_KOKKOS_CUDA_FLOAT
  typedef float  halfPrecision;
#else
  typedef double halfPrecision;
#endif
  typedef double doublePrecision;
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline __device__ __host__ magnitudeType magnitude(double a) { return abs(a); }
  static inline __device__ __host__ double zero()  { return(0.0); }
  static inline __device__ __host__ double one()   { return(1.0); }    
  static inline __device__ __host__ double conjugate(double x)   { return(x); }    
  static inline __device__ __host__ double real(double x) { return x; }
  static inline __device__ __host__ double imag(double)   { return zero(); }
  static inline __device__ __host__ bool  isnaninf(double x) { return isnan(x) || isinf(x); }
  static inline __device__ __host__ double squareroot(double x) { return sqrt(x); }
  static inline __device__ __host__ double pow(double x, double y) { return pow(x,y); }
};
#endif // HAVE_KOKKOS_CUDA_DOUBLE

} // Teuchos namespace

#endif // _TEUCHOS_SCALARTRAITS_CUDA_HPP_
