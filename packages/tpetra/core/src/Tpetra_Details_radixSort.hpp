// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RADIXSORT_HPP
#define TPETRA_DETAILS_RADIXSORT_HPP

#include "TpetraCore_config.h"
#include "Kokkos_Macros.hpp"

namespace Tpetra
{
namespace Details
{

/// \brief Radix sort the input array \c keys, and permute values identically to the keys.
/// 
/// Radix sort may be significantly faster (60%) than Details::shellsort but only works for integers
///
/// \pre KeyType is an integer (can be signed or unsigned)
///
/// \param keys [in/out] Input array of keys to sort.
/// \param keysAux [in] Scratch space (double buffer) for keys (must be allocated to same size as keys)
/// \param values [in/out] Input array of values to permute (must have same number of elements as keys)
/// \param valuesAux [in] Scratch space (double buffer) for values (must be allocated to same size as values)
/// \param n [in] Length of all 4 input arrays keys, keysAux, values and valuesAux.
template<typename KeyType, typename ValueType, typename IndexType>
KOKKOS_INLINE_FUNCTION void
radixSortKeysAndValues(KeyType* keys, KeyType* keysAux, ValueType* values, ValueType* valuesAux, IndexType n, IndexType upperBound)
{
  if(n <= 1)
    return;
  KeyType mask = 0xF;
  bool inAux = false;
  //maskPos counts the low bit index of mask (0, 4, 8, ...)
  IndexType maskPos = 0;
  //Count number of bits required to sort (8 * sizeof(KeyType) - lzcnt(maxKey - minKey))
  IndexType sortBits = 0;
  while(upperBound)
  {
    upperBound >>= 1;
    sortBits++;
  }
  for(IndexType s = 0; s < (sortBits + 3) / 4; s++)
  {
    //Count the number of elements in each bucket
    IndexType count[16] = {0};
    IndexType offset[17] = {0};
    if(!inAux)
    {
      for(IndexType i = 0; i < n; i++)
      {
        count[(keys[i] & mask) >> maskPos]++;
      }
    }
    else
    {
      for(IndexType i = 0; i < n; i++)
      {
        count[(keysAux[i] & mask) >> maskPos]++;
      }
    }
    //get offset as the prefix sum for count
    for(IndexType i = 0; i < 16; i++)
    {
      offset[i + 1] = offset[i] + count[i];
    }
    //now for each element in [lo, hi), move it to its offset in the other buffer
    //this branch should be ok because whichBuf is the same on all threads
    if(!inAux)
    {
      //copy from *Over to *Aux
      for(IndexType i = 0; i < n; i++)
      {
        IndexType bucket = (keys[i] & mask) >> maskPos;
        keysAux[offset[bucket + 1] - count[bucket]] = keys[i];
        valuesAux[offset[bucket + 1] - count[bucket]] = values[i];
        count[bucket]--;
      }
    }
    else
    {
      //copy from *Aux to *Over
      for(IndexType i = 0; i < n; i++)
      {
        IndexType bucket = (keysAux[i] & mask) >> maskPos;
        keys[offset[bucket + 1] - count[bucket]] = keysAux[i];
        values[offset[bucket + 1] - count[bucket]] = valuesAux[i];
        count[bucket]--;
      }
    }
    inAux = !inAux;
    mask = mask << 4;
    maskPos += 4;
  }
  if(inAux)
  {
    //need to deep copy from aux arrays to main
    for(IndexType i = 0; i < n; i++)
    {
      keys[i] = keysAux[i];
      values[i] = valuesAux[i];
    }
  }
}

} //namespace Details
} //namespace Tpetra

#endif  //include guard

