
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_FECrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"
#include "Epetra_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_HashTable.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_OffsetIndex.h"

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     int* NumIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, RowMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     int NumIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, RowMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     const Epetra_BlockMap& ColMap,
				     int* NumIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, RowMap, ColMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     const Epetra_BlockMap& ColMap,
				     int NumIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, RowMap, ColMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::~Epetra_FECrsGraph()
{
  DeleteMemory();
}

//----------------------------------------------------------------------------
void Epetra_FECrsGraph::DeleteMemory()
{
  if (nonlocalGraph_ != 0)
    delete nonlocalGraph_;
  // nothing else to do here, since the STL map has an appropriate
  // destructor
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const int* rows,
					   int numCols, const int* cols)
{
  int returncode = 0;
  int err = 0;

  Epetra_CrsGraph* thisgraph = static_cast<Epetra_CrsGraph*>(this);

  for(int i=0; i<numRows; ++i) {
    const int LID = thisgraph->LRID(rows[i]);
    if (LID > -1) {
      thisgraph->SetIndicesAreGlobal(true);
      err = thisgraph->InsertIndicesIntoSorted(LID, numCols,
          const_cast<int*>(cols));
    }
    else {
      nonlocalRowData_[rows[i]].AddEntries(numCols,cols);
    }

    if (err < 0) return (err);
    if (err > 0) returncode = err;
  }

  return(returncode);
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::GlobalAssemble(bool callFillComplete)
{
  return GlobalAssemble (static_cast<Epetra_Map&>(this->CrsGraphData_->RowMap_), 
      static_cast<Epetra_Map&>(this->CrsGraphData_->RowMap_),
      callFillComplete);
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::GlobalAssemble(const Epetra_Map& domain_map,
                                      const Epetra_Map& range_map,
                                      bool callFillComplete)
{
  if (Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    if (callFillComplete) {
      EPETRA_CHK_ERR( FillComplete(domain_map, range_map) );
    }
    return(0);
  }

  //In this method we need to gather all the non-local (overlapping) data
  //that's been input on each processor, into the
  //non-overlapping distribution defined by the map that 'this' graph was
  //constructed with.

  // First build a map that describes our nonlocal data.
  // We'll use the arbitrary distribution constructor of Map.
  // Start by extracting the column numbers from the STL map.

  const int numRows = nonlocalRowData_.size();
  int * presentRowIndices = new int[numRows];
  std::map<int,Epetra_CrsGraphData::EntriesInOneRow>::iterator nonlocalRows 
    = nonlocalRowData_.begin();
  for (int i=0 ; nonlocalRows != nonlocalRowData_.end(); ++nonlocalRows, ++i)
    presentRowIndices[i] = nonlocalRows->first;

  Epetra_Map* sourceMap = new Epetra_Map(-1, nonlocalRowData_.size(),
          presentRowIndices,
                                         Map().IndexBase(), Map().Comm());

  //If sourceMap has global size 0, then no nonlocal data exists and we can
  //skip most of this function.
  if (sourceMap->NumGlobalElements() < 1) {
    if (callFillComplete) {
      EPETRA_CHK_ERR( FillComplete(domain_map, range_map) );
    }
    delete [] presentRowIndices;
    delete sourceMap;
    return(0);
  }

  //We also need to build a column-map, containing the columns in our
  //nonlocal data. To do that, create a list of all column-indices that
  //occur in our nonlocal rows. This is most easily done using the
  //EntriesInOneRow struct, since that is sorted.
  Epetra_CrsGraphData::EntriesInOneRow allColumns;
  for (nonlocalRows = nonlocalRowData_.begin(); 
       nonlocalRows != nonlocalRowData_.end(); ++nonlocalRows)
    allColumns.AddEntries(nonlocalRows->second.entries_.size(),
       &nonlocalRows->second.entries_[0]);

  Epetra_Map* colMap = new Epetra_Map(-1, allColumns.entries_.size(),
             &allColumns.entries_[0],
                                      Map().IndexBase(), Map().Comm());

  //now we need to create a graph with sourceMap and colMap, and fill it with
  //our nonlocal data so we can then export it to the correct owning processors

  int * rowLengths = new int[numRows];
  {
    int i = 0;
    for (nonlocalRows = nonlocalRowData_.begin(); 
  nonlocalRows != nonlocalRowData_.end() ; ++nonlocalRows, ++i)
      rowLengths[i] = nonlocalRows->second.entries_.size();
  }

  Epetra_CrsGraph* tempGrph = NULL;
  if (buildNonlocalGraph_) {
    nonlocalGraph_ = new Epetra_CrsGraph(Copy, *sourceMap, *colMap, rowLengths);
    tempGrph = nonlocalGraph_;
  }
  else
    tempGrph = new Epetra_CrsGraph(Copy, *sourceMap, *colMap, rowLengths);

  //Next we need to make sure the 'indices-are-global' attribute of tempGrph
  //is set to true, in case this processor doesn't end up calling the
  //InsertGlobalIndices method...

  tempGrph->SetIndicesAreGlobal(true);

  for (nonlocalRows = nonlocalRowData_.begin(); 
       nonlocalRows != nonlocalRowData_.end(); ++nonlocalRows)
    EPETRA_CHK_ERR( tempGrph->InsertGlobalIndices(nonlocalRows->first,
             nonlocalRows->second.entries_.size(),
             &nonlocalRows->second.entries_[0]) );


  //Now we need to call FillComplete on our temp graph. We need to
  //pass a DomainMap and RangeMap.

  EPETRA_CHK_ERR(tempGrph->FillComplete(domain_map, range_map));

  if (buildNonlocalGraph_)
    tempGrph->OptimizeStorage();

  Epetra_Export* exporter = new Epetra_Export(*sourceMap, RowMap());

  EPETRA_CHK_ERR(Export(*tempGrph, *exporter, Add));

  if(callFillComplete) {
    EPETRA_CHK_ERR(FillComplete(domain_map, range_map));
  }

  //now reset the values in our nonlocal data
  for (nonlocalRows = nonlocalRowData_.begin(); 
       nonlocalRows != nonlocalRowData_.end(); ++nonlocalRows)
    nonlocalRows->second.entries_.clear();
  nonlocalRowData_.clear();

  delete [] rowLengths;
  delete [] presentRowIndices;
  delete exporter;
  if (!buildNonlocalGraph_)
    delete tempGrph;
  delete sourceMap;
  delete colMap;

  return(0);
}

