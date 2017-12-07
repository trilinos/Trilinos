// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_NO_SIMD_H
#define STK_NO_SIMD_H

#include <stdio.h>
#include <Kokkos_Macros.hpp>

namespace stk {
namespace simd {

typedef double Double;
typedef bool Bool;
typedef float Float;
typedef bool Boolf;

constexpr int ndoubles = 1; 
constexpr int nfloats = 1; 

KOKKOS_INLINE_FUNCTION Double load_aligned(const double* x) {
  return *x;
}

KOKKOS_INLINE_FUNCTION Double load(const double* x) {
  return *x;
}
    
KOKKOS_INLINE_FUNCTION Double load(const double* x, const int offset) {
  return *x;
}

KOKKOS_INLINE_FUNCTION void store_aligned(double* x, const Double& z) {
  *x = z;
}

KOKKOS_INLINE_FUNCTION void store(double* x, const Double& z) {
  *x = z;
}
  
KOKKOS_INLINE_FUNCTION void store(double* x, const Double& z, const int offset) {
  *x = z;
}

KOKKOS_INLINE_FUNCTION Float load_aligned(const float* x) {
  return *x;
}

KOKKOS_INLINE_FUNCTION Float load(const float* x) {
  return *x;
}
    
KOKKOS_INLINE_FUNCTION Float load(const float* x, const int offset) {
  return *x;
}

KOKKOS_INLINE_FUNCTION void store_aligned(float* x, const Float& z) {
  *x = z;
}

KOKKOS_INLINE_FUNCTION void store(float* x, const Float& z) {
  *x = z;
}
  
KOKKOS_INLINE_FUNCTION void store(float* x, const Float& z, const int offset) {
  *x = z;
}

} // namespace simd
} // namespace stk

#endif // STK_NO_SIMD_H

