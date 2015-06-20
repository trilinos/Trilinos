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
#include "Teuchos_ScalarTraits.hpp"
#include <cfloat>
namespace Teuchos {

#ifdef __CUDA_ARCH__
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

  // Dummy operations, need to exist for parsing when compiling everything with NVCC
  static inline __device__ __host__ int random() { return 9; }
  static inline __device__ __host__ void seedrandom(unsigned int ) {}
};

template<>
struct ScalarTraits<unsigned int>
{
  typedef unsigned int magnitudeType;
  typedef unsigned int halfPrecision;
  typedef unsigned int doublePrecision;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline __device__ __host__ magnitudeType magnitude(unsigned int a) { return (unsigned int)fabsf((float)a); }
  static inline __device__ __host__ unsigned int zero()  { return 0; }
  static inline __device__ __host__ unsigned int one()   { return 1; }
  static inline __device__ __host__ unsigned int conjugate(unsigned int x) { return x; }
  static inline __device__ __host__ unsigned int real(unsigned int x) { return x; }
  static inline __device__ __host__ unsigned int imag(unsigned int)   { return 0; }
  static inline __device__ __host__ bool isnaninf(unsigned int) { return false; }
  static inline __device__ __host__ unsigned int squareroot(unsigned int x) { return (unsigned int)sqrtf((float)x); }          // perhaps this cast should be replaced by an explicit call like __float2int_rn
  static inline __device__ __host__ unsigned int pow(unsigned int x, unsigned int y) { return (unsigned int)powf((float)x,(float)y); }  // perhaps this cast should be replaced by an explicit call like __float2int_rn

  // Dummy operations, need to exist for parsing when compiling everything with NVCC
  static inline __device__ __host__ unsigned int random() { return 9; }
  static inline __device__ __host__ void seedrandom(unsigned int ) {}
};

template<>
struct ScalarTraits<long int>
{
  typedef long int magnitudeType;
  typedef long int halfPrecision;
  typedef long int doublePrecision;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline __device__ __host__ magnitudeType magnitude(long int a) { return (long int)fabsf((float)a); }
  static inline __device__ __host__ long int zero()  { return 0; }
  static inline __device__ __host__ long int one()   { return 1; }
  static inline __device__ __host__ long int conjugate(long int x) { return x; }
  static inline __device__ __host__ long int real(long int x) { return x; }
  static inline __device__ __host__ long int imag(long int) { return 0; }
  static inline __device__ __host__ bool isnaninf(int) { return false; }
  static inline __device__ __host__ long int squareroot(long int x) { return (long int)sqrtf((float)x); }          // perhaps this cast should be replaced by an explicit call like __float2int_rn
  static inline __device__ __host__ long int pow(long int x, long int y) { return (long int)powf((float)x,(float)y); }  // perhaps this cast should be replaced by an explicit call like __float2int_rn

  // Dummy operations, need to exist for parsing when compiling everything with NVCC
  static inline __device__ __host__ long int random() { return 9; }
  static inline __device__ __host__ void seedrandom(unsigned int ) {}
};

template<>
struct ScalarTraits<long unsigned int>
{
  typedef long unsigned int magnitudeType;
  typedef long unsigned int halfPrecision;
  typedef long unsigned int doublePrecision;
  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline __device__ __host__ magnitudeType magnitude(long unsigned int a) { return (long unsigned int)fabsf((float)a); }
  static inline __device__ __host__ long unsigned int zero()  { return 0; }
  static inline __device__ __host__ long unsigned int one()   { return 1; }
  static inline __device__ __host__ long unsigned int conjugate(long unsigned int x) { return x; }
  static inline __device__ __host__ long unsigned int real(long unsigned int x) { return x; }
  static inline __device__ __host__ long unsigned int imag(long unsigned int) { return 0; }
  static inline __device__ __host__ bool isnaninf(int) { return false; }
  static inline __device__ __host__ long unsigned int squareroot(long unsigned int x) { return (long unsigned int)sqrtf((float)x); }          // perhaps this cast should be replaced by an explicit call like __float2int_rn
  static inline __device__ __host__ long unsigned int pow(long unsigned int x, long unsigned int y) { return (long unsigned int)powf((float)x,(float)y); }  // perhaps this cast should be replaced by an explicit call like __float2int_rn

  // Dummy operations, need to exist for parsing when compiling everything with NVCC
  static inline __device__ __host__ long unsigned int random() { return 9; }
  static inline __device__ __host__ void seedrandom(unsigned int ) {}
};

template<>
struct ScalarTraits<float>
{
  typedef float magnitudeType;
  typedef float halfPrecision; // should become IEEE754-2008 binary16 or fp16 later, perhaps specified at configure according to architectural support
  typedef double doublePrecision;
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
  static inline __device__ __host__ float eps() { return FLT_EPSILON; }
  static inline __device__ __host__ float t() { return FLT_MANT_DIG; }
  static inline __device__ __host__ float base() { return FLT_RADIX; }
  static inline __device__ __host__ float log10(float x ) { return ::log10f(x); }

  static inline __device__ __host__ float prec()  { return eps()*base(); }
  static inline __device__ __host__ float rnd()   { return 1.0f; }
  static inline __device__ __host__ float sfmin() { return FLT_MIN; }
  static inline __device__ __host__ float emin()  { return FLT_MIN_EXP; }
  static inline __device__ __host__ float rmin()  { return FLT_MIN; }
  static inline __device__ __host__ float emax()  { return FLT_MAX_EXP; }
  static inline __device__ __host__ float rmax()  { return FLT_MAX; }
  static inline __device__ __host__ float nan()   { return 0.0f/std::sin(0.0f); }
  static inline __device__ __host__ const char* name() { return "float"; }

  // Dummy operations, need to exist for parsing when compiling everything with NVCC
  static inline __device__ __host__ float random() { return 9.0f; }
  static inline __device__ __host__ void seedrandom(unsigned int ) {}
};

template<>
struct ScalarTraits<double>
{
  typedef double magnitudeType;
  typedef float  halfPrecision;
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
  static inline __device__ __host__ double eps() { return DBL_EPSILON; }
  static inline __device__ __host__ double t() { return DBL_MANT_DIG; }
  static inline __device__ __host__ double base() { return FLT_RADIX; }
  static inline __device__ __host__ double log10(double x ) { return ::log10(x); }

  static inline __device__ __host__ double prec()  { return eps()*base(); }
  static inline __device__ __host__ double rnd()   { return 1.0; }
  static inline __device__ __host__ double sfmin() { return DBL_MIN; }
  static inline __device__ __host__ double emin()  { return DBL_MIN_EXP; }
  static inline __device__ __host__ double rmin()  { return DBL_MIN; }
  static inline __device__ __host__ double emax()  { return DBL_MAX_EXP; }
  static inline __device__ __host__ double rmax()  { return DBL_MAX; }
  static inline __device__ __host__ double nan()   { return 0.0/std::sin(0.0); }
  static inline __device__ __host__ const char* name() { return "double"; }

  // Dummy operations, need to exist for parsing when compiling everything with NVCC
  static inline __device__ __host__ double random() { return 9.0; }
  static inline __device__ __host__ void seedrandom(unsigned int ) {}
};
#endif // __CUDA_ARCH__

} // Teuchos namespace

#endif // _TEUCHOS_SCALARTRAITS_CUDA_HPP_
