//------------------------------------------------------------------------*/
//   Copyright (c) 2013, Sandia Corporation.
//   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//   the U.S. Governement retains certain rights in this software.
//   
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//   
//       * Redistributions of source code must retain the above copyright
//         notice, this list of conditions and the following disclaimer.
//   
//       * Redistributions in binary form must reproduce the above
//         copyright notice, this list of conditions and the following
//         disclaimer in the documentation and/or other materials provided
//         with the distribution.
//   
//       * Neither the name of Sandia Corporation nor the names of its
//         contributors may be used to endorse or promote products derived
//         from this software without specific prior written permission.
//   
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//                 

// Portions of this source code file are:
// Copyright (c) 2010, Victor J. Duvanenko.
// All rights reserved.
// Used with permission.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * The name of Victor J. Duvanenko may not be used to endorse or
//     promote products derived from this software without specific prior
//     written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_UTIL_UTIL_RADIXSORT_HPP
#define STK_UTIL_UTIL_RADIXSORT_HPP

#include <stdexcept>

//
// Internal implementation
//
namespace { // anonymous namespace

// Swap that does not check for self-assignment.
template< class T >
inline void swap_impl( T& a, T& b )
{
  T tmp = a;
  a     = b;
  b     = tmp;
}

// insertion sort similar to STL with no self-assignment
template <class T>
inline void insertion_sort_impl( T* a, size_t a_size )
{
  for ( size_t i = 1; i < a_size; ++i )
  {
    if ( a[ i ] < a[ i - 1 ] )    // no need to do (j > 0) compare for the first iteration
    {
      T currentElement = a[ i ];
      a[ i ] = a[ i - 1 ];
      size_t j;
      for ( j = i - 1; j > 0 && currentElement < a[ j - 1 ]; --j )
      {
        a[ j ] = a[ j - 1 ];
      }
      a[ j ] = currentElement;  // always necessary work/write
    }
    // Perform no work at all if the first comparison fails - i.e. never assign an element to itself!
  }
}

// Recursive implementation of radix sort for unsigned integer types only.
//
//    PowerOfTwoRadix - must be a power of 2, 256 is suggested for current hardware as of 03/20/2011.
//    Log2ofPowerOfTwoRadix - for example log( 256 ) = 8
//    Threshold - length below which array sections are sorted with an insertion sort
//    a - pointer to start of an array of integer type
//    last - one less than the length of the a array
//    bitMask - controls how many and which bits we process at a time
//    shiftRightAmount - sizeof(T) * 8 - Log2ofPowerOfTwoRadix
template< class T, unsigned long PowerOfTwoRadix, unsigned long Log2ofPowerOfTwoRadix, long Threshold >
inline void radix_sort_unsigned_impl( T* a, long last, T bitMask, unsigned long shiftRightAmount )
{
  unsigned long count[ PowerOfTwoRadix ];
  for( unsigned long i = 0; i < PowerOfTwoRadix; ++i )
    count[ i ] = 0;
  for ( long _current = 0; _current <= last; ++_current )
    count[( a[ _current ] & bitMask ) >> shiftRightAmount]++;  // Scan the array and count the number of times each value appears

  long startOfBin[ PowerOfTwoRadix + 1 ], endOfBin[ PowerOfTwoRadix ], nextBin = 1;
  startOfBin[ 0 ] = endOfBin[ 0 ] = 0;    startOfBin[ PowerOfTwoRadix ] = -1;     // sentinel
  for( unsigned long i = 1; i < PowerOfTwoRadix; ++i )
    startOfBin[ i ] = endOfBin[ i ] = startOfBin[ i - 1 ] + count[ i - 1 ];

  for ( long _current = 0; _current <= last; )
  {
    unsigned long digit;
    T _current_element = a[ _current ]; // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
    while( endOfBin[ digit = (( _current_element & bitMask ) >> shiftRightAmount )] != _current )  swap_impl( _current_element, a[ endOfBin[ digit ]++ ] );
    a[ _current ] = _current_element;

    endOfBin[ digit ]++;
    while( endOfBin[ nextBin - 1 ] == startOfBin[ nextBin ] )  ++nextBin; // skip over empty and full bins, when the end of the current bin reaches the start of the next bin
    _current = endOfBin[ nextBin - 1 ];
  }
  bitMask >>= Log2ofPowerOfTwoRadix;
  if ( bitMask != 0 )           // end recursion when all the bits have been processes
  {
    if ( shiftRightAmount >= Log2ofPowerOfTwoRadix )  shiftRightAmount -= Log2ofPowerOfTwoRadix;
    else                        shiftRightAmount  = 0;

    for( unsigned long i = 0; i < PowerOfTwoRadix; ++i )
    {
      long numberOfElements = endOfBin[ i ] - startOfBin[ i ];
      if ( numberOfElements >= Threshold )    // endOfBin actually points to one beyond the bin
        radix_sort_unsigned_impl< T, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( &a[ startOfBin[ i ]], numberOfElements - 1, bitMask, shiftRightAmount );
      else if ( numberOfElements >= 2 )
        insertion_sort_impl( &a[ startOfBin[ i ]], numberOfElements );
    }
  }
}

} // end anonymous namespace

//
// Public API
//
namespace stk {
namespace util {

/**
 * Sort a contiguous array of unsigned 8-, 16-, 32- or 64-bit integers in place,
 * using an unstable MSD radix sort algorithm.
 * <p>
 * This MSD radix sort algorithm is O(N k/d), where N is the number of items to be sorted,
 * k, the size of each key, and d, the digit size used by the implementation. Additional
 * memory used is (k/d 2^d) (k/d recursion levels, 2^d for count array).
 * </p>
 * <p>
 * For comparison, the std::sort algorithm is O(N logN), and uses (log N) additional memory.
 * The MSD radix sort on unsigned integers typically outperforms the std::sort by using less
 * than half the CPU time for 64-bit integers.
 * </p>
 * <p>
 * TODO: Also consider the threaded TBB implementation below of the MSD radix sort.
 * </p>
 *
 * @param a       pointer to the start of a contiguous array of unsigned or signed type
 * @param a_size  length of the array (must be less than MAX_LONG)
 */
template< class T >
void radix_sort_unsigned( T* a, size_t a_size )
{
  const long Threshold                      =  25;  // smaller array sections sorted by insertion sort
  const unsigned long PowerOfTwoRadix       = 256;  // Radix - must be a power of 2
  const unsigned long Log2ofPowerOfTwoRadix =   8;  // log( 256 ) = 8

  if ( a_size >= static_cast<size_t>(Threshold) ) {
    if (a_size > static_cast<size_t>(std::numeric_limits<long>::max())) {
      std::ostringstream msg ;
      msg << "stk::utility::radix_sort() exceeded allowable array size (";
      msg << a_size << " < " << std::numeric_limits<long>::max() << ")";
      throw std::runtime_error( msg.str() );
    }
    unsigned long shiftRightAmount = sizeof(T) * 8 - Log2ofPowerOfTwoRadix;
    T bitMask = static_cast<T>( (static_cast<T>( PowerOfTwoRadix - 1 )) << shiftRightAmount );  // bitMask controls/selects how many and which bits we process at a time
    radix_sort_unsigned_impl< T, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( a, a_size - 1, bitMask, shiftRightAmount );
  }
  else
    insertion_sort_impl( a, a_size );
}

} // namespace util
} // namespace stk

#endif /* STK_UTIL_UTIL_RADIXSORT_HPP */
