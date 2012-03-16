/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_FECRSGRAPH_H
#define EPETRA_FECRSGRAPH_H

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"

  // TODO this file needs to be changed for long long

#include <map>

/**
  Epetra Finite-Element CrsGraph. This class provides the ability to insert
  indices into a matrix-graph, where the indices represent dense submatrices
  such as element-stiffnesses that might arise from a finite-element
  application.

  In a parallel setting, indices may be submitted on the local processor
  for rows that do not reside in the local portion of the row-map. After
  all indices have been submitted, the GlobalAssemble method gathers all
  non-local graph rows to the appropriate 'owning' processors (an owning
  processor is a processor which has the row in its row-map).
 */
class EPETRA_LIB_DLL_EXPORT Epetra_FECrsGraph : public Epetra_CrsGraph {
  friend class Epetra_FECrsMatrix;

  public:

  /** Constructor */
  Epetra_FECrsGraph(Epetra_DataAccess CV,
		    const Epetra_BlockMap& RowMap,
		    int* NumIndicesPerRow,
		    bool ignoreNonLocalEntries=false,
        bool buildNonlocalGraph=false);

  /** Constructor */
  Epetra_FECrsGraph(Epetra_DataAccess CV,
		    const Epetra_BlockMap& RowMap,
		    int NumIndicesPerRow,
		    bool ignoreNonLocalEntries=false,
        bool buildNonlocalGraph=false);

  /** Constructor */
  Epetra_FECrsGraph(Epetra_DataAccess CV,
		    const Epetra_BlockMap& RowMap, 
		    const Epetra_BlockMap& ColMap,
		    int* NumIndicesPerRow,
		    bool ignoreNonLocalEntries=false,
        bool buildNonlocalGraph=false);

  /** Constructor */
  Epetra_FECrsGraph(Epetra_DataAccess CV,
		    const Epetra_BlockMap& RowMap, 
		    const Epetra_BlockMap& ColMap,
		    int NumIndicesPerRow,
		    bool ignoreNonLocalEntries=false,
        bool buildNonlocalGraph=false);

  /** Constructor */
  Epetra_FECrsGraph(const Epetra_FECrsGraph& Graph);

  /** Destructor */
  virtual ~Epetra_FECrsGraph();

  //Let the compiler know we intend to overload the base-class function
  //InsertGlobalIndices rather than hide it.
  using Epetra_CrsGraph::InsertGlobalIndices;

  /** Insert a rectangular, dense 'submatrix' of entries (matrix nonzero
      positions) into the graph.

    @param numRows Number of rows in the submatrix.
    @param rows List of row-numbers for the submatrix.
    @param numCols Number of columns in the submatrix.
    @param cols List of column-indices that will be used for each row in
        the 'rows' list.
  */
  int InsertGlobalIndices(int numRows, const int* rows,
			  int numCols, const int* cols);

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this matrix at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.

      ***NOTE***: When GlobalAssemble() calls FillComplete(), it passes the
      arguments 'DomainMap()' and 'RangeMap()', which are the map attributes
      held by the base-class CrsMatrix and its graph. If a rectangular matrix
      is being assembled, the domain-map and range-map must be specified by
      calling the other overloading of this method. Otherwise, GlobalAssemble()
      has no way of knowing what these maps should really be.


      @param callFillComplete option argument, defaults to true.
        Determines whether GlobalAssemble() internally calls the
        FillComplete() method on this matrix.

      @return error-code 0 if successful, non-zero if some error occurs
   */
  int GlobalAssemble(bool callFillComplete=true);

  /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this matrix at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.

      ***NOTE***: When GlobalAssemble() (the other overloading of this method)
      calls FillComplete(), it passes the arguments 'DomainMap()' and
      'RangeMap()', which are the map attributes already held by the base-class
      CrsMatrix and its graph. If a rectangular matrix is being assembled, the
      domain-map and range-map must be specified. Otherwise, GlobalAssemble()
      has no way of knowing what these maps should really be.


      @param domain_map user-supplied domain map for this matrix

      @param range_map user-supplied range map for this matrix

      @param callFillComplete option argument, defaults to true.
        Determines whether GlobalAssemble() internally calls the
        FillComplete() method on this matrix.

      @return error-code 0 if successful, non-zero if some error occurs
   */
  int GlobalAssemble(const Epetra_Map& domain_map,
                     const Epetra_Map& range_map,
                     bool callFillComplete=true);

  bool UseNonlocalGraph () const {return buildNonlocalGraph_; };

 private:     
  void DeleteMemory();
  int InsertNonlocalRow(int row, int offset);
  int InputNonlocalIndices(int row,
			   int numCols,
			   const int* cols);
  int InputNonlocalIndex(int rowoffset,
			 int col);

  long long myFirstRow_;
  int myNumRows_;
  bool ignoreNonLocalEntries_;

  /**
   * This STL map holds all non-local data in format of Entries in the
   * individual rows together with the row number.
   */
//TODO FIXME std::map<int,Epetra_CrsGraphData::EntriesInOneRow> nonlocalRowData_;

  /**
   * A CrsGraph holding non-local data in case the respective flag is set in
   * the constructor.
   */
  Epetra_CrsGraph* nonlocalGraph_;
  bool buildNonlocalGraph_;

  Epetra_FECrsGraph & operator=(const Epetra_FECrsGraph& Graph);
     

};//class Epetra_FECrsGraph

#endif
