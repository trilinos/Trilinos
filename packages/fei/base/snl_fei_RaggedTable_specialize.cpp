/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <fei_macros.hpp>
#include <snl_fei_RaggedTable_specialize.hpp>

namespace snl_fei {

/** specialization for MapContig<fei::ctg_set<int>*> */
RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::RaggedTable(int firstKey, int lastKey)
  : map_(firstKey, lastKey),
    poolAllocatorSet_(),
    dummy()
{
  int len = lastKey-firstKey+1;
  if (len > 0) {
    map_type::value_type val;
    for(int i=0; i<len; ++i) {
      val.first = firstKey+i;
      row_type* row = poolAllocatorSet_.allocate(1);
      poolAllocatorSet_.construct(row,dummy);
      val.second = row;
      map_.insert(val);
    }
  }
}

RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::RaggedTable(const RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >& src)
 : map_(src.map_),
   poolAllocatorSet_()
{
}

void RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::addIndices(int row,
                             int numIndices,
                             const int* indices)
{
  iterator m_end = map_.end();
  iterator m_iter = map_.lower_bound(row);

  map_type::mapped_type mapped_indices = (*m_iter).second;

  if (mapped_indices == NULL) {
    throw std::runtime_error("RaggedTable<MapContig>, NULL row.");
  }

  for(int i=0; i<numIndices; ++i) {
    mapped_indices->insert2(indices[i]);
  }
}

void
RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::addIndices(int numRows,
                             const int* rows,
                             int numIndices,
                             const int* indices)
{
  iterator m_end = map_.end();
  map_type::mapped_type mapped_indices = NULL;

  for(int i=0; i<numRows; ++i) {
    int row = rows[i];
    iterator m_iter = map_.lower_bound(row);

    mapped_indices = (*m_iter).second;
    if (mapped_indices == NULL) {
      throw std::runtime_error("RaggedTable<MapContig>, NULL row.");
    }

    for(int j=0; j<numIndices; ++j) {
      mapped_indices->insert2(indices[j]);
    }
  }
}

void RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::addDiagonals(int numIndices,
                                                         const int* indices)
{
  for(int i=0; i<numIndices; ++i) {
    int ind = indices[i];
    addIndices(ind, 1, &ind);
  }
}

}//namespace snl_fei

