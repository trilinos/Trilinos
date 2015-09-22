
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_IndexTable_hpp_
#define _fei_IndexTable_hpp_

#include <fei_macros.hpp>

namespace fei {
  /** Abstract interface for adding index mappings to a table of indices,
      such as an algebraic matrix-graph.
  */
  class IndexTable {
  public:
    /** Constructor */
    IndexTable() {}

    /** Destructor */
    virtual ~IndexTable() {}

    /** Input function to add diagonals to the index table.
     */
    virtual void addDiagonals(int numIndices,
			      const int* indices) = 0;

    /** Input function 'addIndices' specifies the row of the table to be
	operated on, and a list of indices to be added to that row.
    */
    virtual void addIndices(int row,
			    int numIndices,
			    const int* indices) = 0;

    /** Input function for adding a list of indices to multiple rows.
     */
    virtual void addIndices(int numRows,
			    const int* rows,
			    int numIndices,
			    const int* indices) = 0;
  };//class IndexTable

}//namespace fei

#endif
