#ifndef __READMATRIXFROMFILE_HPP
#define __READMATRIXFROMFILE_HPP

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
  
  Teuchos::RCP<const Teuchos::Comm<int>> pComm = pRowMap->getComm ();
  const int numProcs = pComm->getSize ();
  const int myRank = pComm->getRank ();
  const int rootRank = 0;

  Teuchos::ArrayView<const global_ordinal_type> myRows = pRowMap->getNodeElementList();
  const size_t myNumRows = myRows.size();
  
  myNumEntriesPerRow = Teuchos::ArrayRCP<size_t> (myNumRows);

  if (myRank != rootRank) {

    Teuchos::send (*pComm, myNumRows, rootRank);
    
    if (myNumRows != 0) {

      Teuchos::send (*pComm, static_cast<int> (myNumRows), myRows.getRawPtr(), rootRank);      
      Teuchos::receive (*pComm, rootRank, static_cast<int> (myNumRows), myNumEntriesPerRow.getRawPtr());

      const local_ordinal_type myNumEntries =
	std::accumulate (myNumEntriesPerRow.begin(),
			 myNumEntriesPerRow.end(), 0);

      myColInd = Teuchos::ArrayRCP<global_ordinal_type> (myNumEntries);
      myValues = Teuchos::ArrayRCP<scalar_type> (myNumEntries, 1.0);
      if (myNumEntries > 0) {
	receive (*pComm, rootRank, static_cast<int> (myNumEntries), myColInd.getRawPtr());
	// receive (*pComm, rootRank,
	// 	 static_cast<int> (myNumEntries),
	// 	 myValues.getRawPtr());
      }
      

    } // If I own at least one row

  } // If I am not the root processor
  else { // I _am_ the root processor
    if (debug) {
      std::cout << "-- Proc 0: Copying my data from global arrays" << std::endl;
    }

    for (size_t k = 0; k < myNumRows; ++k) {
      const global_ordinal_type myCurRow = myRows[k];
      const local_ordinal_type numEntriesInThisRow = rowPtr[myCurRow+1] - rowPtr[myCurRow];
      myNumEntriesPerRow[k] = numEntriesInThisRow;
      
    }

    const local_ordinal_type myNumEntries =
      std::accumulate (myNumEntriesPerRow.begin(),
    		       myNumEntriesPerRow.end(), 0);
    if (debug) {
      std::cout << "-- Proc 0: I own " << myNumRows << " rows and "
    		<< myNumEntries << " entries" << std::endl;
    }
    myColInd = Teuchos::ArrayRCP<global_ordinal_type> (myNumEntries);
    myValues = Teuchos::ArrayRCP<scalar_type> (myNumEntries);

    local_ordinal_type myCurPos = 0;
    for (size_t k = 0; k < myNumRows; ++k) {
      const local_ordinal_type curNumEntries = myNumEntriesPerRow[k];
      const global_ordinal_type myRow = myRows[k];
      global_ordinal_type curPos = rowPtr[myRow];

      if (curNumEntries > 0) {
      	for(size_t ii = 0; ii < curNumEntries; ++ii) {
      	  myColInd[myCurPos++] = colInd[curPos++];
      	}
      }
    }


    myValues = Teuchos::ArrayRCP<scalar_type>(myNumEntries, 1.0);
    
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

	const local_ordinal_type theirNumEntries =
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
	Teuchos::ArrayRCP<scalar_type> theirValues (theirNumEntries);

	local_ordinal_type theirCurPos = 0;
	for (size_t k = 0; k < theirNumRows; k++) {
	  const local_ordinal_type curNumEntries = theirNumEntriesPerRow[k];
	  const global_ordinal_type theirRow = theirRows[k];
	  local_ordinal_type curPos = rowPtr[theirRow];

	  if (curNumEntries > 0) {

	    for(size_t ii = 0; ii < curNumEntries; ++ii) {
	      theirColInd[theirCurPos++] = colInd[curPos++];
	    }
	    
	    //create theirValues if you are going to send

	  }
	}

	Teuchos::send (*pComm, static_cast<int> (theirNumEntries), theirColInd.getRawPtr(), p);
	// send (*pComm, static_cast<int> (theirNumEntries),
	//       theirValues.getRawPtr(), p);

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

  //std::cout << globalNumRows << " " << globalNumNonzeros << std::endl;

  global_ordinal_type *rowPtr;
  global_ordinal_type *colInd;

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

  Teuchos::ArrayRCP<size_t> myNumEntriesPerRow;
  Teuchos::ArrayRCP<size_t> myRowPtr;
  Teuchos::ArrayRCP<global_ordinal_type> myColInd;
  Teuchos::ArrayRCP<scalar_type> myValues;
 
  distribute<global_ordinal_type, local_ordinal_type, scalar_type, map_type>(myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap, rowPtr, colInd, debug);




  if (debug && myRank == rootRank) {
    std::cout << "-- Inserting matrix entries on each processor";
    if (callFillComplete) {
      std::cout << " and calling fillComplete()";
    }
    std::cout << std::endl;
  }
  
  Teuchos::RCP<crs_matrix_type> pMatrix = Teuchos::rcp (new crs_matrix_type (pRowMap, myNumEntriesPerRow, Tpetra::DynamicProfile));

  Teuchos::ArrayView<const global_ordinal_type> myRows = pRowMap->getNodeElementList ();
  const size_t myNumRows = myRows.size ();

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

  myNumEntriesPerRow = Teuchos::null;
  myRowPtr = Teuchos::null;
  myColInd = Teuchos::null;
  myValues = Teuchos::null;
  
  if (callFillComplete) {
    pMatrix->fillComplete (pDomainMap, pRangeMap);
  }
  
  
  if(myRank == rootRank) {
    delete [] rowPtr;
    delete [] colInd;
  }
    
  return pMatrix;

}


template <typename crs_matrix_type>
Teuchos::RCP<crs_matrix_type>
readMatrixFromFile(std::string filename, const Teuchos::RCP<const Teuchos::Comm<int>> pComm, bool binary=true, bool debug=false)
{
  if(binary)
    return readBinaryFile<crs_matrix_type>(filename, pComm, true, debug);
  else
    return Tpetra::MatrixMarket::Reader<crs_matrix_type>::readSparseFile(filename, pComm);

}



#endif

