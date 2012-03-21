
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
				     const Epetra_BlockMap& rowMap,
				     int* numIndicesPerRow,
				     bool ignoreNonLocalEntries,
				     bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, rowMap, numIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& rowMap,
				     int numIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, rowMap, numIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& rowMap,
				     const Epetra_BlockMap& colMap,
				     int* numIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, rowMap, colMap, numIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& rowMap,
				     const Epetra_BlockMap& colMap,
				     int numIndicesPerRow,
				     bool ignoreNonLocalEntries,
             bool buildNonlocalGraph)
  : Epetra_CrsGraph(CV, rowMap, colMap, numIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    nonlocalGraph_ (NULL),
    buildNonlocalGraph_ (buildNonlocalGraph)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
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
template<typename int_type>
int Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const int_type* rows,
					   int numCols, const int_type* cols)
{
  int returncode = 0;
  int err = 0;

  Epetra_CrsGraph* thisgraph = static_cast<Epetra_CrsGraph*>(this);

  for(int i=0; i<numRows; ++i) {
    const int LID = thisgraph->LRID(rows[i]);
    if (LID > -1) {
      thisgraph->SetIndicesAreGlobal(true);
      err = thisgraph->InsertIndicesIntoSorted(LID, numCols,
          const_cast<int_type*>(cols));
    }
    else {
       nonlocalRowData<int_type>()[rows[i]].AddEntries(numCols,cols);
    }

    if (err < 0) return (err);
    if (err > 0) returncode = err;
  }

  return(returncode);
}

int Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const int* rows, int numCols, const int* cols)
{
  if(RowMap().GlobalIndicesInt())
	return InsertGlobalIndices<int>(numRows, rows, numCols, cols);
  else
	throw ReportError("Epetra_FECrsGraph::InsertGlobalIndices int version called for a matrix that is not int.", -1);
}

int Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const long long* rows, int numCols, const long long* cols)
{
  if(RowMap().GlobalIndicesLongLong())
	return InsertGlobalIndices<long long>(numRows, rows, numCols, cols);
  else
	throw ReportError("Epetra_FECrsGraph::InsertGlobalIndices long long version called for a matrix that is not long long.", -1);
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::GlobalAssemble(bool callFillComplete)
{
  return GlobalAssemble (static_cast<Epetra_Map&>(this->CrsGraphData_->RowMap_), 
      static_cast<Epetra_Map&>(this->CrsGraphData_->RowMap_),
      callFillComplete);
}

//----------------------------------------------------------------------------
template<typename int_type>
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

  std::map<int_type,Epetra_CrsGraphData::EntriesInOneRow<int_type> >& nonlocalRowData_var = nonlocalRowData<int_type>();

  const int numRows = (int) nonlocalRowData_var.size();
  int * presentRowIndices = new int[numRows];
  typename std::map<int_type,Epetra_CrsGraphData::EntriesInOneRow<int_type> >::iterator nonlocalRows 
    = nonlocalRowData<int_type>().begin();
  for (int i=0 ; nonlocalRows != nonlocalRowData_var.end(); ++nonlocalRows, ++i)
    presentRowIndices[i] = (int) nonlocalRows->first;

  Epetra_Map* sourceMap = new Epetra_Map((int_type) -1, (int) nonlocalRowData<int_type>().size(),
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
  Epetra_CrsGraphData::EntriesInOneRow<int_type> allColumns;
  for (nonlocalRows = nonlocalRowData_var.begin(); 
       nonlocalRows != nonlocalRowData_var.end(); ++nonlocalRows)
    allColumns.AddEntries((int) nonlocalRows->second.entries_.size(),
       &nonlocalRows->second.entries_[0]);

  Epetra_Map* colMap = new Epetra_Map((int_type) -1, (int) allColumns.entries_.size(),
             &allColumns.entries_[0],
                                      Map().IndexBase(), Map().Comm());

  //now we need to create a graph with sourceMap and colMap, and fill it with
  //our nonlocal data so we can then export it to the correct owning processors

  int * rowLengths = new int[numRows];
  {
    int i = 0;
    for (nonlocalRows = nonlocalRowData_var.begin(); 
  nonlocalRows != nonlocalRowData_var.end() ; ++nonlocalRows, ++i)
      rowLengths[i] = (int) nonlocalRows->second.entries_.size();
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

  for (nonlocalRows = nonlocalRowData_var.begin(); 
       nonlocalRows != nonlocalRowData_var.end(); ++nonlocalRows)
    EPETRA_CHK_ERR( tempGrph->InsertGlobalIndices(nonlocalRows->first,
             (int) nonlocalRows->second.entries_.size(),
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
  for (nonlocalRows = nonlocalRowData_var.begin(); 
       nonlocalRows != nonlocalRowData_var.end(); ++nonlocalRows)
    nonlocalRows->second.entries_.clear();
  nonlocalRowData_var.clear();

  delete [] rowLengths;
  delete [] presentRowIndices;
  delete exporter;
  if (!buildNonlocalGraph_)
    delete tempGrph;
  delete sourceMap;
  delete colMap;

  return(0);
}

int Epetra_FECrsGraph::GlobalAssemble(const Epetra_Map& domain_map,
                                      const Epetra_Map& range_map,
                                      bool callFillComplete)
{
  if(!domain_map.GlobalIndicesTypeMatch(range_map))
     throw ReportError("Epetra_FECrsGraph::GlobalAssemble: cannot be called with different indices types for domainMap and rangeMap", -1);

  if(!RowMap().GlobalIndicesTypeMatch(domain_map))
    throw ReportError("Epetra_FECrsGraph::GlobalAssemble: cannot be called with different indices types for row map and incoming rangeMap", -1);

  if(RowMap().GlobalIndicesInt())
	  return GlobalAssemble<int>(domain_map, range_map, callFillComplete);
  else if(RowMap().GlobalIndicesLongLong())
	  return GlobalAssemble<long long>(domain_map, range_map, callFillComplete);

  throw ReportError("Epetra_FECrsGraph::GlobalAssemble: cannot determine global index type", -1);
}
