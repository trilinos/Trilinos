#ifndef KOKKOS_FUNCTIONAL_HPP
#define KOKKOS_FUNCTIONAL_HPP

#include <Kokkos_Macros.hpp>
#include <stdint.h>

namespace Kokkos {

namespace Impl {

// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
KOKKOS_FORCEINLINE_FUNCTION
uint32_t getblock32 ( const uint8_t * p, int i )
{
// used to avoid aliasing error which could cause errors with
// forced inlining
  return    ((uint32_t)p[i*4+0])
          | ((uint32_t)p[i*4+1] << 8)
          | ((uint32_t)p[i*4+2] << 16)
          | ((uint32_t)p[i*4+3] << 24);
}

KOKKOS_FORCEINLINE_FUNCTION
uint32_t rotl32 ( uint32_t x, int8_t r )
{ return (x << r) | (x >> (32 - r)); }

KOKKOS_FORCEINLINE_FUNCTION
uint32_t fmix32 ( uint32_t h )
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

KOKKOS_INLINE_FUNCTION
uint32_t MurmurHash3_x86_32 ( const void * key, int len, uint32_t seed )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;

  uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  for(int i=0; i<nblocks; ++i)
  {
    uint32_t k1 = getblock32(data,i);

    k1 *= c1;
    k1 = rotl32(k1,15);
    k1 *= c2;

    h1 ^= k1;
    h1 = rotl32(h1,13);
    h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

  uint32_t k1 = 0;

  switch(len & 3)
  {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;

  h1 = fmix32(h1);

  return h1;
}

} // namespace Impl


// These should work for most types

template <typename T>
struct hash
{
  typedef T argument_type;
  typedef T first_argument_type;
  typedef uint32_t second_argument_type;
  typedef uint32_t result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const & t) const
  { return Impl::MurmurHash3_x86_32( &t, sizeof(T), 0); }

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const & t, uint32_t seed) const
  { return Impl::MurmurHash3_x86_32( &t, sizeof(T), seed); }
};



template <typename T>
struct equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a == b; }
};

template <typename T>
struct not_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a != b; }
};


template <typename T>
struct greater
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a > b; }
};


template <typename T>
struct less
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a < b; }
};

template <typename T>
struct greater_equal
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a >= b; }
};


template <typename T>
struct less_equal
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a <= b; }
};

} // namespace Kokkos


#endif //KOKKOS_FUNCTIONAL_HPP


