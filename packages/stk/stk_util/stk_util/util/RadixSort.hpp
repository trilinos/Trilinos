/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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

namespace stk {
namespace util {

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
void insertion_sort_impl( T* a, size_t a_size )
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
  for( unsigned long i = 0; i < PowerOfTwoRadix; ++i )     count[ i ] = 0;
  for ( long _current = 0; _current <= last; ++_current )  count[ (unsigned long)(( a[ _current ] & bitMask ) >> shiftRightAmount ) ]++;  // Scan the array and count the number of times each value appears

  long startOfBin[ PowerOfTwoRadix + 1 ], endOfBin[ PowerOfTwoRadix ], nextBin = 1;
  startOfBin[ 0 ] = endOfBin[ 0 ] = 0;    startOfBin[ PowerOfTwoRadix ] = -1;     // sentinel
  for( unsigned long i = 1; i < PowerOfTwoRadix; ++i )
    startOfBin[ i ] = endOfBin[ i ] = startOfBin[ i - 1 ] + count[ i - 1 ];

  for ( long _current = 0; _current <= last; )
  {
    unsigned long digit;
    T _current_element = a[ _current ]; // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
    while( endOfBin[ digit = (unsigned long)(( _current_element & bitMask ) >> shiftRightAmount )] != _current )  swap_impl( _current_element, a[ endOfBin[ digit ]++ ] );
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

/**
 * Sort a contiguous array of unsigned 8-, 16-, 32- or 64-bit integers in place,
 * using the radix sort algorithm.
 *
 * Performance is O(N). Runtime is typically 2-5 times faster than the comparison sort used
 * by the std::sort().
 *
 * Internal storage required is 3*(M+1)*sizeof(long), where M is PowerOfTwoRadix.
 *
 * @param a       pointer to the start of a contiguous array of unsigned or signed type
 * @param a_size  length of the array (must be less than MAX_LONG)
 */
template< class T >
inline void radix_sort_unsigned( T* a, size_t a_size )
{
  const long Threshold                      =  25;  // smaller array sections sorted by insertion sort
  const unsigned long PowerOfTwoRadix       = 256;  // Radix - must be a power of 2
  const unsigned long Log2ofPowerOfTwoRadix =   8;  // log( 256 ) = 8

  if ( a_size >= (size_t)Threshold ) {
    if (a_size > (size_t)std::numeric_limits<long>::max()) {
      std::ostringstream msg ;
      msg << "stk::utility::radix_sort() exceeded allowable array size (";
      msg << a_size << " < " << std::numeric_limits<long>::max() << ")";
      throw std::runtime_error( msg.str() );
    }
    unsigned long shiftRightAmount = sizeof(T) * 8 - Log2ofPowerOfTwoRadix;
    T bitMask = (T)( ((T)( PowerOfTwoRadix - 1 )) << shiftRightAmount );  // bitMask controls/selects how many and which bits we process at a time
    radix_sort_unsigned_impl< T, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( a, a_size - 1, bitMask, shiftRightAmount );
  }
  else
    insertion_sort_impl( a, a_size );
}

//--------------------------------------
// PARALLEL RADIX SORT using TBB library
//--------------------------------------

/* TODO: this code works with TBB, it just needs minor modifications and some renaming
template< unsigned long PowerOfTwoRadix, unsigned long Log2ofPowerOfTwoRadix, long Threshold, class T >
class RadixInPlaceOperation_3
{
  T* my_a;                    // a local copy to the input array to provide a pointer to each parallel task
  long*  my_startOfBin;
  long*  my_endOfBin;
  T  my_bitMask;                  // a local copy of the bitMask
  unsigned long my_shiftRightAmount;
  static const unsigned long Threshold_P = 10000;   // threshold when to switch between parallel and non-parallel implementations
public:
  static const unsigned long numberOfBins = PowerOfTwoRadix;
  unsigned long my_count[ numberOfBins ];       // the count for this task

  RadixInPlaceOperation_3(  T a[], long* startOfBin, long* endOfBin, T bitMask, unsigned long shiftRightAmount ) :
                my_a( a ), my_startOfBin( startOfBin ), my_endOfBin( endOfBin ),
                my_bitMask( bitMask ), my_shiftRightAmount( shiftRightAmount )
  {}
  void operator()( const blocked_range< long >& r ) const
  {
    for( long i = r.begin(); i != r.end(); ++i )
    {
      long numOfElements = my_endOfBin[ i ] - my_startOfBin[ i ];
      if ( numOfElements >= Threshold_P )
        _RadixSort_Unsigned_PowerOf2Radix_3< PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( &my_a[ my_startOfBin[ i ]], numOfElements - 1, my_bitMask, my_shiftRightAmount );
      else if ( numOfElements >= Threshold )
        _RadixSort_Unsigned_PowerOf2Radix_1< PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( &my_a[ my_startOfBin[ i ]], numOfElements - 1, my_bitMask, my_shiftRightAmount );
      else if ( numOfElements >= 2 )
        insertion_sort_impl( &my_a[ my_startOfBin[ i ]], numOfElements );
    }
  }
};

template< unsigned long PowerOfTwoRadix, unsigned long Log2ofPowerOfTwoRadix, long Threshold, class T >
inline void _RadixSort_Unsigned_PowerOf2Radix_3( T* a, long last, T bitMask, unsigned long shiftRightAmount )
{
  const unsigned long numberOfBins = PowerOfTwoRadix;

  CountingType_1< PowerOfTwoRadix, T >  count( a, bitMask, shiftRightAmount );  // contains the count array, which is initialized to all zeros
  // Scan the array and count the number of times each value appears
  parallel_reduce( blocked_range< unsigned long >( 0, last + 1, 1000 ), count );

  long startOfBin[ numberOfBins ], endOfBin[ numberOfBins ], nextBin;
  startOfBin[ 0 ] = endOfBin[ 0 ] = nextBin = 0;
  for( unsigned long i = 1; i < numberOfBins; i++ )
    startOfBin[ i ] = endOfBin[ i ] = startOfBin[ i - 1 ] + count.my_count[ i - 1 ];

  for ( long _current = 0; _current <= last; )
  {
    unsigned long digit;
    T tmp = a[ _current ];  // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
    while ( true ) {
      digit = (unsigned long)(( tmp & bitMask ) >> shiftRightAmount );  // extract the digit we are sorting based on
      if ( endOfBin[ digit ] == _current )
        break;
      swap_impl( tmp, a[ endOfBin[ digit ] ] );
      endOfBin[ digit ]++;
    }
    a[ _current ] = tmp;

    endOfBin[ digit ]++;          // leave the element at its location and grow the bin
    _current++;               // advance the current pointer to the next element
    while( _current >= startOfBin[ nextBin ] && nextBin < numberOfBins )
      nextBin++;
    while( endOfBin[ nextBin - 1 ] == startOfBin[ nextBin ] && nextBin < numberOfBins )
      nextBin++;
    if ( _current < endOfBin[ nextBin - 1 ] )
       _current = endOfBin[ nextBin - 1 ];
  }
  bitMask >>= Log2ofPowerOfTwoRadix;
  if ( bitMask != 0 )           // end recursion when all the bits have been processes
  {
    if ( shiftRightAmount >= Log2ofPowerOfTwoRadix )  shiftRightAmount -= Log2ofPowerOfTwoRadix;
    else                        shiftRightAmount  = 0;

    RadixInPlaceOperation_3< PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, T >  radixOp( a, startOfBin, endOfBin, bitMask, shiftRightAmount );
    parallel_for( blocked_range< long >( 0, numberOfBins ), radixOp );
  }
}

template< unsigned long PowerOfTwoRadix, unsigned long Log2ofPowerOfTwoRadix, long Threshold, class T >
inline void _RadixSort_Unsigned_PowerOf2Radix_1( T* a, long last, T bitMask, unsigned long shiftRightAmount )
{
  const unsigned long numberOfBins = PowerOfTwoRadix;
  unsigned long count[ numberOfBins ];
// Counting occurrence of digits within the array (related to Counting Sort)
  for( unsigned long i = 0; i < numberOfBins; i++ )
    count[ i ] = 0;

  for ( long _current = 0; _current <= last; _current++ ) // Scan the array and count the number of times each value appears
  {
    unsigned long digit = (unsigned long)(( a[ _current ] & bitMask ) >> shiftRightAmount );  // extract the digit we are sorting based on
    count[ digit ]++;
  }
// Moving array elements into their bins
  long startOfBin[ numberOfBins ], endOfBin[ numberOfBins ], nextBin;
  startOfBin[ 0 ] = endOfBin[ 0 ] = nextBin = 0;
  for( unsigned long i = 1; i < numberOfBins; i++ )
    startOfBin[ i ] = endOfBin[ i ] = startOfBin[ i - 1 ] + count[ i - 1 ];

  for ( long _current = 0; _current <= last; )
  {
    unsigned long digit;
    T tmp = a[ _current ];  // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
    while ( true ) {
      digit = (unsigned long)(( tmp & bitMask ) >> shiftRightAmount );  // extract the digit we are sorting based on
      if ( endOfBin[ digit ] == _current )
        break;
      swap_impl( tmp, a[ endOfBin[ digit ] ] );
      endOfBin[ digit ]++;
    }
    a[ _current ] = tmp;

    endOfBin[ digit ]++;          // leave the element at its location and grow the bin
    _current++;               // advance the current pointer to the next element
    while( _current >= startOfBin[ nextBin ] && nextBin < numberOfBins )
      nextBin++;
    while( endOfBin[ nextBin - 1 ] == startOfBin[ nextBin ] && nextBin < numberOfBins )
      nextBin++;
    if ( _current < endOfBin[ nextBin - 1 ] )
       _current = endOfBin[ nextBin - 1 ];
  }
// Recursion for each bin
  bitMask >>= Log2ofPowerOfTwoRadix;
  if ( bitMask != 0 )           // end recursion when all the bits have been processes
  {
    if ( shiftRightAmount >= Log2ofPowerOfTwoRadix )  shiftRightAmount -= Log2ofPowerOfTwoRadix;
    else                        shiftRightAmount  = 0;

    for( unsigned long i = 0; i < numberOfBins; i++ )
    {
      long numberOfElements = endOfBin[ i ] - startOfBin[ i ];
      if ( numberOfElements >= Threshold ) {    // endOfBin actually points to one beyond the bin
        _RadixSort_Unsigned_PowerOf2Radix_1< T, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( &a[ startOfBin[ i ]], numberOfElements - 1, bitMask, shiftRightAmount );
    }
    else if ( numberOfElements >= 2 ) {
    insertion_sort_impl( &a[ startOfBin[ i ]], numberOfElements );
    }
    }
  }
}

// Top-level template function
template< class T >
inline void RadixSortParallel( T* a, unsigned long a_size )
{
  if ( a_size < 2 ) return;

  const unsigned long Threshold             =  32;  // Threshold of when to switch to using Insertion Sort
  const unsigned long PowerOfTwoRadix       = 256;  // Radix - must be a power of 2
  const unsigned long Log2ofPowerOfTwoRadix =   8;
  unsigned long shiftRightAmount = sizeof( T ) * 8 - Log2ofPowerOfTwoRadix;   // Create bit-mask and shift right amount
  T bitMask = (T)( ((T)( PowerOfTwoRadix - 1 )) << shiftRightAmount );  // bitMask controls/selects how many and which bits we process at a time

  if ( a_size >= Threshold )
    _RadixSort_Unsigned_PowerOf2Radix_3< PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold >( a, a_size - 1, bitMask, shiftRightAmount );
  else
    insertion_sort_impl( a, a_size );
}

*/

} // namespace util
} // namespace stk


#endif /* STK_UTIL_UTIL_RADIXSORT_HPP */
