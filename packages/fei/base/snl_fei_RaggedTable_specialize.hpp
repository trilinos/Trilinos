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

#ifndef _snl_fei_RaggedTable_specialize_hpp_
#define _snl_fei_RaggedTable_specialize_hpp_


#include <fei_macros.hpp>
#include <snl_fei_RaggedTable.hpp>
#include <fei_ctg_set.hpp>
#include <snl_fei_MapContig.hpp>

namespace snl_fei {

/** specialization for MapContig<fei::ctg_set<int> > */
template<>
class RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >
  : public fei::IndexTable {
public:
  RaggedTable(int firstKey, int lastKey);

  RaggedTable(const RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >& src);

  typedef MapContig<fei::ctg_set<int>*> map_type;
  typedef fei::ctg_set<int> row_type;
  typedef map_type::iterator iterator;

  virtual ~RaggedTable()
  {
    iterator it = begin();
    iterator it_end = end();
    for(; it!=it_end; ++it) {
      poolAllocatorSet_.destroy( (*it).second );
      poolAllocatorSet_.deallocate( (*it).second, 1 );
    }
  }

  void addDiagonals(int numIndices,
                    const int* indices);

  void addIndices(int row,
                  int numIndices,
                  const int* indices);

  void addIndices(int numRows,
                  const int* rows,
                  int numIndices,
                  const int* indices);

  map_type& getMap() { return( map_ ); }

  row_type* getRow(int row)
  {
    iterator m_end = map_.end();
    iterator m_iter = map_.find(row);
    return( m_end == m_iter ? NULL : (*m_iter).second );
  }

  iterator begin() { return( map_.begin() ); }

  iterator end() { return( map_.end() ); }

 private:
  map_type map_;
  fei_Pool_alloc<row_type> poolAllocatorSet_;
  row_type dummy;
};//RaggedTable<MapContig<fei::ctg_set<int>*> >

}//namespace snl_fei
#endif

