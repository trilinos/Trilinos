// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __READMATRIXFROMBINARYFILE_HPP
#define __READMATRIXFROMBINARYFILE_HPP

/////////////////////////////////////////////////////////////////////////////
// These utilities read a sparse matrix from a binary file that was 
// created with the corresponding "largestComponent2Binary.cpp" code.
// Specific structure may be assumed.
// Note: This is research code. We do not guarantee it works in all cases.
/////////////////////////////////////////////////////////////////////////////
#include <climits> 
#include "Tpetra_Details_makeColMap.hpp"

template <typename global_ordinal_type, typename local_ordinal_type, typename scalar_type, typename map_type>
void
distribute (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
	    Teuchos::ArrayRCP<size_t>& myRowPtr,
	    Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
	    Teuchos::ArrayRCP<scalar_type>& myValues,
	    const Teuchos::RCP<const map_type>& pRowMap,
	    global_ordinal_type *rowPtr,
	    global_ordinal_type *colInd,
	    const bool debug=false)
{
  
  int maxNumEnt = INT_MAX/sizeof(global_ordinal_type);

  Teuchos::RCP<const Teuchos::Comm<int>> pComm = pRowMap->getComm ();
  const int numProcs = pComm->getSize ();
  const int myRank = pComm->getRank ();
  const int rootRank = 0;

  Teuchos::ArrayView<const global_ordinal_type> myRows = pRowMap->getLocalElementList();
  const size_t myNumRows = myRows.size();
  
  myNumEntriesPerRow = Teuchos::ArrayRCP<size_t> (myNumRows);

  if (myRank != rootRank) {

    Teuchos::send (*pComm, myNumRows, rootRank);
    
    if (myNumRows != 0) {

      Teuchos::send (*pComm, static_cast<int> (myNumRows), myRows.getRawPtr(), rootRank);      
      Teuchos::receive (*pComm, rootRank, static_cast<int> (myNumRows), myNumEntriesPerRow.getRawPtr());

      const global_ordinal_type myNumEntries =
	std::accumulate (myNumEntriesPerRow.begin(),
			 myNumEntriesPerRow.end(), 0);

      myColInd = Teuchos::ArrayRCP<global_ordinal_type> (myNumEntries);
      myValues = Teuchos::ArrayRCP<scalar_type> (myNumEntries, 1.0);
      if (myNumEntries > 0) {

        if(myNumEntries < maxNumEnt)
          Teuchos::receive (*pComm, rootRank, static_cast<int> (myNumEntries), &myColInd[0]); 
        else {
          int nchunks = myNumEntries/maxNumEnt;
          if(myNumEntries % maxNumEnt != 0)
            nchunks ++;
          for(int i = 0; i < nchunks-1; i++) {
            Teuchos::receive (*pComm, rootRank, maxNumEnt, &myColInd[maxNumEnt*i]);
            std::cout << "Chunk " << i << " received by myRank "<< myRank << ", size: " << maxNumEnt << "\n";
          }
          int lastsize = (int)(myNumEntries - (nchunks-1)*maxNumEnt);
          Teuchos::receive (*pComm, rootRank, lastsize, &myColInd[maxNumEnt*(nchunks-1)]);
          std::cout << "Chunk " << nchunks-1 << " received by myRank " << myRank << ", size: " << lastsize << "\n";
        }

      }

    } // If I own at least one row

  } // If I am not the root processor
  else { // I _am_ the root processor
    if (debug) {
      std::cout << "-- Proc 0: Copying my data from global arrays" << std::endl;
    }

    for (size_t k = 0; k < myNumRows; ++k) {
      const global_ordinal_type myCurRow = myRows[k];
      const global_ordinal_type numEntriesInThisRow = rowPtr[myCurRow+1] - rowPtr[myCurRow];
      myNumEntriesPerRow[k] = numEntriesInThisRow;
      
    }

    size_t myNumEntries = std::accumulate (myNumEntriesPerRow.begin(),
					   myNumEntriesPerRow.end(), 0);
    if (debug) {
      std::cout << "-- Proc 0: I own " << myNumRows << " rows and "
    		<< myNumEntries << " entries" << std::endl;
    }
    myColInd = Teuchos::ArrayRCP<global_ordinal_type> (myNumEntries);
    myValues = Teuchos::ArrayRCP<scalar_type> (myNumEntries, 1.0);

    global_ordinal_type myCurPos = 0;
    for (size_t k = 0; k < myNumRows; ++k) {
      const global_ordinal_type curNumEntries = myNumEntriesPerRow[k];
      const global_ordinal_type myRow = myRows[k];
      global_ordinal_type curPos = rowPtr[myRow];

      if (curNumEntries > 0) {
      	for(global_ordinal_type ii = 0; ii < curNumEntries; ++ii) {
      	  myColInd[myCurPos++] = colInd[curPos++];
      	}
      }
    }
    
    for (int p = 1; p < numProcs; ++p) {
      if (debug) {
	std::cout << "-- Proc 0: Processing proc " << p << std::endl;
      }

      size_t theirNumRows = 0;
      Teuchos::receive (*pComm, p, &theirNumRows);
      if (debug) {
	std::cout << "-- Proc 0: Proc " << p << " owns "
		  << theirNumRows << " rows" << std::endl;
      }
     
      if (theirNumRows != 0) {
	Teuchos::ArrayRCP<global_ordinal_type> theirRows(theirNumRows);
	Teuchos::receive (*pComm, p, Teuchos::as<int> (theirNumRows), theirRows.getRawPtr ());
	
	Teuchos::ArrayRCP<size_t> theirNumEntriesPerRow = Teuchos::ArrayRCP<size_t> (theirNumRows);
	for (size_t k = 0; k < theirNumRows; ++k) {
	  theirNumEntriesPerRow[k] = rowPtr[theirRows[k]+1] - rowPtr[theirRows[k]];
	}

	Teuchos::send (*pComm, static_cast<int> (theirNumRows), theirNumEntriesPerRow.getRawPtr(), p);

	const global_ordinal_type theirNumEntries =
	  std::accumulate (theirNumEntriesPerRow.begin(),
			   theirNumEntriesPerRow.end(), 0);

	if (debug) {
	  std::cout << "-- Proc 0: Proc " << p << " owns "
		    << theirNumEntries << " entries" << std::endl;
	}

	if (theirNumEntries == 0) {
	  continue;
	}

	Teuchos::ArrayRCP<global_ordinal_type> theirColInd (theirNumEntries);

	global_ordinal_type theirCurPos = 0;
	for (size_t k = 0; k < theirNumRows; k++) {
	  const global_ordinal_type curNumEntries = theirNumEntriesPerRow[k];
	  const global_ordinal_type theirRow = theirRows[k];
	  global_ordinal_type curPos = rowPtr[theirRow];

	  if (curNumEntries > 0) {

	    for(global_ordinal_type ii = 0; ii < curNumEntries; ++ii) {
	      theirColInd[theirCurPos++] = colInd[curPos++];
	    }
	    
	  }
	}

	if(theirNumEntries < maxNumEnt)
	  Teuchos::send (*pComm, static_cast<int> (theirNumEntries), &theirColInd[0], p);
	else {
	  int nchunks = theirNumEntries/maxNumEnt;
	  if(theirNumEntries % maxNumEnt != 0)
	    nchunks ++;
	  for(int i = 0; i < nchunks-1; i++) {
	    Teuchos::send (*pComm, maxNumEnt, &theirColInd[maxNumEnt*i], p);
	    std::cout << "Chunk " << i << " sent to Rank "<< p << ", size: " << maxNumEnt << "\n";
	  }
	  int lastsize = (int)(theirNumEntries - (nchunks-1)*maxNumEnt);
	  Teuchos::send (*pComm, lastsize, &theirColInd[maxNumEnt*(nchunks-1)], p);
	  std::cout << "Chunk " << nchunks-1 << " sent to Rank "<< p << ", size: " << lastsize << "\n";
	}
	  
	if (debug) {
	  std::cout << "-- Proc 0: Finished with proc " << p << std::endl;
	}
	
      } // If proc p owns at least one row
    } // For each proc p not the root proc 0
  } // If I'm (not) the root proc 0

  if (debug && myRank == 0) {
    std::cout << "-- Proc 0: About to fill in myRowPtr" << std::endl;
  }

  myRowPtr = Teuchos::ArrayRCP<size_t> (myNumRows+1);
  myRowPtr[0] = 0;
  for (size_t k = 1; k < myNumRows+1; ++k) {
    myRowPtr[k] = myRowPtr[k-1] + myNumEntriesPerRow[k-1];
  }
  if (debug && myRank == 0) {
    std::cout << "-- Proc 0: Done with distribute" << std::endl;
  }
}

template <typename crs_matrix_type>
Teuchos::RCP<crs_matrix_type>
readBinaryFile(std::string filename, const Teuchos::RCP<const Teuchos::Comm<int>> pComm, bool callFillComplete=true, bool debug=false)
{
  typedef typename crs_matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename crs_matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename crs_matrix_type::scalar_type scalar_type;
  typedef typename crs_matrix_type::node_type node_type;

  typedef typename Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  const int myRank = pComm->getRank();
  const int rootRank = 0;

  if (debug && myRank == rootRank) {
    std::cout << "Binary CRS reader: readSparse:" << std::endl
  	      << "-- Reading started" << std::endl;
  }

  Teuchos::RCP<std::ifstream> in;
  global_ordinal_type globalNumRows;
  global_ordinal_type globalNumNonzeros;
  
  if (myRank == rootRank) {
    
    // Open the file
    in = Teuchos::RCP<std::ifstream>(new std::ifstream(filename, std::ios::in | std::ios::binary));

    // Read number of vertices and number of edges
    in->read((char *)&globalNumRows, sizeof(global_ordinal_type));
    in->read((char *)&globalNumNonzeros, sizeof(global_ordinal_type));

    TEUCHOS_TEST_FOR_EXCEPTION(globalNumRows <= 0 || globalNumNonzeros <= 0, std::invalid_argument,
    			       "Global number of rows or nonzeros have nonpositive value." << globalNumRows << " " << globalNumNonzeros << " " << sizeof(global_ordinal_type) );
  }
  
  broadcast (*pComm, rootRank, 1, &globalNumRows);
  broadcast (*pComm, rootRank, 1, &globalNumNonzeros);

  global_ordinal_type *rowPtr = 0;
  global_ordinal_type *colInd = 0;
  
  if (myRank == rootRank) {

    rowPtr = new global_ordinal_type[globalNumRows+1];
    colInd = new global_ordinal_type[globalNumNonzeros];

    in->read((char*)rowPtr, sizeof(global_ordinal_type)*(globalNumRows+1));
    in->read((char*)colInd, sizeof(global_ordinal_type)*(globalNumNonzeros));
  }
  
  Teuchos::RCP<const map_type> pRowMap = Teuchos::rcp (new map_type (static_cast<Tpetra::global_size_t> (globalNumRows),
  								     static_cast<global_ordinal_type> (0),
  								     pComm, Tpetra::GloballyDistributed));
  
  Teuchos::RCP<const map_type> pRangeMap = Teuchos::rcp (new map_type (static_cast<Tpetra::global_size_t> (globalNumRows),
  								       static_cast<global_ordinal_type> (0),
  								       pComm, Tpetra::GloballyDistributed));
  
  Teuchos::RCP<const map_type> pDomainMap = pRangeMap;

  Teuchos::ArrayView<const global_ordinal_type> myRows = pRowMap->getLocalElementList ();
  const size_t myNumRows = myRows.size ();

  Teuchos::ArrayRCP<size_t> myNumEntriesPerRow(myNumRows);
  Teuchos::ArrayRCP<size_t> myRowPtr;
  Teuchos::ArrayRCP<global_ordinal_type> myColInd;
  Teuchos::ArrayRCP<scalar_type> myValues;
   
  distribute<global_ordinal_type, local_ordinal_type, scalar_type, map_type>(myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap, rowPtr, colInd, debug);
  pComm->barrier();
  
  if (debug && myRank == rootRank) {
    std::cout << "-- Inserting matrix entries on each processor";
    if (callFillComplete) {
      std::cout << " and calling fillComplete()";
    }
    std::cout << std::endl;
  }
  
  Teuchos::RCP<crs_matrix_type> pMatrix = Teuchos::rcp (new crs_matrix_type (pRowMap, myNumEntriesPerRow()));

  const global_ordinal_type indexBase = pRowMap->getIndexBase ();
  for (size_t i = 0; i < myNumRows; ++i) {
    const size_t myCurPos = myRowPtr[i];
    const local_ordinal_type curNumEntries = myNumEntriesPerRow[i];
    Teuchos::ArrayView<global_ordinal_type> curColInd = myColInd.view (myCurPos, curNumEntries);
    Teuchos::ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

    for (size_t k = 0; k < curNumEntries; ++k) {
      curColInd[k] += indexBase;
    }
    
    if (curNumEntries > 0) {
      pMatrix->insertGlobalValues (myRows[i], curColInd, curValues);
    }
  }
  pComm->barrier();
  if (debug && myRank == rootRank) {
    std::cout << "-- Done with inserting." << std::endl;
  }

  myNumEntriesPerRow = Teuchos::null;
  myRowPtr = Teuchos::null;
  myColInd = Teuchos::null;
  myValues = Teuchos::null;
  
  if (callFillComplete) {
    pMatrix->fillComplete (pDomainMap, pRangeMap);
  }
  pComm->barrier();
  if (debug && myRank == rootRank) {
    std::cout << "-- Done with fill complete." << std::endl;
  }
  

  if(myRank == rootRank) {
    delete [] rowPtr;
    delete [] colInd;
  }
    
  return pMatrix;
}

template <typename crs_matrix_type>
Teuchos::RCP<crs_matrix_type>
readBinaryFileFast(std::string filename, const Teuchos::RCP<const Teuchos::Comm<int>> pComm, bool callFillComplete=true, bool debug=false)
{
  typedef typename crs_matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename crs_matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename crs_matrix_type::scalar_type scalar_type;
  typedef typename crs_matrix_type::node_type node_type;
  typedef typename crs_matrix_type::crs_graph_type crs_graph_type;

  typedef typename Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  const int myRank = pComm->getRank();
  const int rootRank = 0;

  if (debug && myRank == rootRank) {
    std::cout << "Binary CRS reader: readSparse:" << std::endl
  	      << "-- Reading started" << std::endl;
  }

  Teuchos::RCP<std::ifstream> in;
  global_ordinal_type globalNumRows;
  global_ordinal_type globalNumNonzeros;
  
  if (myRank == rootRank) {
    
    // Open the file
    in = Teuchos::RCP<std::ifstream>(new std::ifstream(filename, std::ios::in | std::ios::binary));

    // Read number of vertices and number of edges
    in->read((char *)&globalNumRows, sizeof(global_ordinal_type));
    in->read((char *)&globalNumNonzeros, sizeof(global_ordinal_type));

    TEUCHOS_TEST_FOR_EXCEPTION(globalNumRows <= 0 || globalNumNonzeros <= 0, std::invalid_argument,
    			       "Global number of rows or nonzeros have nonpositive value." << globalNumRows << " " << globalNumNonzeros << " " << sizeof(global_ordinal_type) );
  }
  
  broadcast (*pComm, rootRank, 1, &globalNumRows);
  broadcast (*pComm, rootRank, 1, &globalNumNonzeros);

  global_ordinal_type *rowPtr = 0;
  global_ordinal_type *colInd = 0;
  
  if (myRank == rootRank) {

    rowPtr = new global_ordinal_type[globalNumRows+1];
    colInd = new global_ordinal_type[globalNumNonzeros];

    in->read((char*)rowPtr, sizeof(global_ordinal_type)*(globalNumRows+1));
    in->read((char*)colInd, sizeof(global_ordinal_type)*(globalNumNonzeros));

  }

  
  Teuchos::RCP<const map_type> pRowMap = Teuchos::rcp (new map_type (static_cast<Tpetra::global_size_t> (globalNumRows),
  								     static_cast<global_ordinal_type> (0),
  								     pComm, Tpetra::GloballyDistributed));
  
  Teuchos::RCP<const map_type> pRangeMap = Teuchos::rcp (new map_type (static_cast<Tpetra::global_size_t> (globalNumRows),
  								       static_cast<global_ordinal_type> (0),
  								       pComm, Tpetra::GloballyDistributed));
  
  Teuchos::RCP<const map_type> pDomainMap = pRangeMap;

  Teuchos::ArrayView<const global_ordinal_type> myRows = pRowMap->getLocalElementList ();
  const size_t myNumRows = myRows.size ();

  Teuchos::ArrayRCP<size_t> myNumEntriesPerRow(myNumRows);
  Teuchos::ArrayRCP<size_t> myRowPtr;
  Teuchos::ArrayRCP<global_ordinal_type> myColInd;
  Teuchos::ArrayRCP<scalar_type> myValues;
 
   
  distribute<global_ordinal_type, local_ordinal_type, scalar_type, map_type>(myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap, rowPtr, colInd, debug);
  pComm->barrier();
  
  if (debug && myRank == rootRank) {
    std::cout << "-- Inserting matrix entries on each processor";
    if (callFillComplete) {
      std::cout << " and calling fillComplete()";
    }
    std::cout << std::endl;
  }
  
  // get the colIds
  std::vector<bool> mark(globalNumRows, false);
  size_t myNumEntries = myRowPtr[myNumRows]; 
  for(size_t i = 0; i < myNumEntries; i++)
    mark[myColInd[i]] = true;

  local_ordinal_type myNumCols = 0;
  for(global_ordinal_type i = 0; i < globalNumRows; i++)
    if(mark[i] == true)
      myNumCols++;

  Kokkos::View<global_ordinal_type*, typename node_type::memory_space> myColGIDs("myColGIDs", myNumCols);
  auto myColGIDs_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), myColGIDs);

  myNumCols = 0;
  for(global_ordinal_type i = 0; i < globalNumRows; i++)
    if(mark[i] == true) {
      myColGIDs_host(myNumCols)= i;
      myNumCols++;
    };

  Kokkos::deep_copy(myColGIDs, myColGIDs_host);

  Teuchos::RCP<const map_type> pColumnMap; 
  Tpetra::Details::makeColMap(pColumnMap, pDomainMap, myColGIDs);

  std::vector<local_ordinal_type> map(globalNumRows);
  for(global_ordinal_type i = 0; i < globalNumRows; i++) {
    if(mark[i] == true)
      map[i] = pColumnMap->getLocalElement(i);
  }

  Teuchos::ArrayRCP<local_ordinal_type> myLclColInd(myNumEntries);
  for(size_t i = 0; i < myNumEntries; i++)
    myLclColInd[i] = map[myColInd[i]];

  local_ordinal_type *cur = myLclColInd.getRawPtr();
  
  for(size_t i = 0; i < myNumRows; i++) {
    size_t start = myRowPtr[i];
    size_t end = myRowPtr[i+1];

    std::sort(&cur[start], &cur[end]);
  }
  
  Teuchos::RCP<crs_graph_type> graph(new crs_graph_type(pRowMap, pColumnMap, myRowPtr, myLclColInd));
  graph->fillComplete(pDomainMap, pRangeMap);

  Kokkos::View<scalar_type*, typename node_type::memory_space> values("values", myNumEntries);
  Kokkos::deep_copy(values, 1.0);
  
  Teuchos::RCP<crs_matrix_type> pMatrix (new crs_matrix_type(graph, values));
  pMatrix->fillComplete(pDomainMap, pRangeMap);
  
  pComm->barrier();
  if (debug && myRank == rootRank) {
    std::cout << "-- Done with fill complete." << std::endl;
  }
  
  if(myRank == rootRank) {
    delete [] rowPtr;
    delete [] colInd;
  }
    
  return pMatrix;
}

template <typename crs_matrix_type>
Teuchos::RCP<crs_matrix_type>
readMatrixFromBinaryFile(std::string filename, const Teuchos::RCP<const Teuchos::Comm<int>> pComm, bool binary=true, bool debug=false)
{
  return readBinaryFileFast<crs_matrix_type>(filename, pComm, true, debug);
}

#endif

