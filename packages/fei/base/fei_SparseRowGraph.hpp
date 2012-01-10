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



#ifndef _fei_SparseRowGraph_hpp_
#define _fei_SparseRowGraph_hpp_

#include <fei_macros.hpp>
#include <vector>
#include <algorithm>

namespace fei {
  /** Extremely simple data container that represents a sparse row-oriented
      matrix-graph. Purely serial. If it is used to store the local portion of
      a distributed matrix-graph, the calling code is responsible for all
      knowledge related to parallelism.
  */
  class SparseRowGraph {
  public:
    /** Default constructor */
    SparseRowGraph()
      : rowNumbers(), rowOffsets(), packedColumnIndices(), blockEntries(false)
      {}

    /** Copy constructor */
    SparseRowGraph(const SparseRowGraph& src)
      : rowNumbers(src.rowNumbers), rowOffsets(src.rowOffsets),
      packedColumnIndices(src.packedColumnIndices), blockEntries(src.blockEntries)
      {}

    /** Destructor */
    virtual ~SparseRowGraph() {}

    /** comparison operator */
    bool operator==(const fei::SparseRowGraph& othergraph) const;

    /** not-equal operator */
    bool operator!=(const fei::SparseRowGraph& othergraph) const;

    /** Local row-numbers. */
    std::vector<int> rowNumbers;

    /** The starting offset of each row in the packedColumnIndices
        vector. The i-th row corresponds to positions rowOffsets[i] through
        rowOffsets[i+1]-1.  Note that rowOffsets should have length
        rowNumbers.size()+1, and
        rowOffsets[rowNumbers.size()] == packedColumnIndices.size().
    */
    std::vector<int> rowOffsets;

    /** Contiguous array of column-indices for all local rows.
       See the comments for the 'rowOffsets' attribute, for information
       about accessing column-indices for a particular row, etc.*/
    std::vector<int> packedColumnIndices;

    /** whether this graph represents a block-entry matrix. */
    bool blockEntries;
  };//class SparseRowGraph

inline bool SparseRowGraph::operator==(const fei::SparseRowGraph& othergraph) const
{
  if (rowNumbers != othergraph.rowNumbers) return(false);
  if (rowOffsets != othergraph.rowOffsets) return(false);
  if (packedColumnIndices != othergraph.packedColumnIndices) return(false);
  return(true);
}

inline bool SparseRowGraph::operator!=(const fei::SparseRowGraph& othergraph) const
{
  return( !(*this == othergraph) );
}

/** Given a row-number and a SparseRowGraph object, return the offset at which
    that row's column-indices start in the SparseRowGraph object's
    packedColumnIndices vector.

    If the given row-number is not found in the SparseRowGraph object's vector
    of row-numbers, return -1.
*/
inline
int find_row_start(int row, const SparseRowGraph& srg)
{
  std::vector<int>::const_iterator rowNumbers_iter =
    std::lower_bound(srg.rowNumbers.begin(), srg.rowNumbers.end(), row);
  if (rowNumbers_iter == srg.rowNumbers.end() || *rowNumbers_iter != row) {
    return -1;
  }

  size_t offset = rowNumbers_iter - srg.rowNumbers.begin();
  return srg.rowOffsets[offset];
}

}//namespace fei

#endif

