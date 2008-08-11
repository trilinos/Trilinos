
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_SparseRowGraph_hpp_
#define _fei_SparseRowGraph_hpp_

#include "fei_macros.hpp"
#include <vector>

namespace fei {
  /** Extremely simple data container (basically a struct) that represents
      the local portion of a sparse row-oriented matrix-graph (may or may not
      include rows that are shared but not owned locally).
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
    bool operator==(const fei::SparseRowGraph& othergraph);

    /** not-equal operator */
    bool operator!=(const fei::SparseRowGraph& othergraph);

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

inline bool SparseRowGraph::operator==(const fei::SparseRowGraph& othergraph)
{
  if (rowNumbers != othergraph.rowNumbers) return(false);
  if (rowOffsets != othergraph.rowOffsets) return(false);
  if (packedColumnIndices != othergraph.packedColumnIndices) return(false);
  return(true);
}

inline bool SparseRowGraph::operator!=(const fei::SparseRowGraph& othergraph)
{
  return( !(*this == othergraph) );
}

}//namespace fei

#endif
