#ifndef TPETRA_MATRIX_IO_DEF
#define TPETRA_MATRIX_IO_DEF

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include <iostream>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
Tpetra::Utils::generateMatrix(const Teuchos::RCP<Teuchos::ParameterList> &plist,
                              const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                              const Teuchos::RCP<Node> &node,
                              Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPTION( plist == Teuchos::null, std::runtime_error,
      "Tpetra::Utils::generateMatrix(): ParameterList is null.");
  TEST_FOR_EXCEPTION( Teuchos::isParameterType<std::string>(*plist,"mat_type") == false, std::runtime_error,
      "Tpetra::Utils::generateMatrix(): ParameterList did not contain string parameter ""mat_type"".");
  std::string mat_type = plist->get<std::string>("mat_type");
  if (mat_type == "Lap3D") {
    // 3D Laplacian, grid is a cube with dimension gridSize x gridSize x gridSize
    const int gridSize = plist->get<int>("gridSize",100);
    const GlobalOrdinal gS2 = (GlobalOrdinal)gridSize*(GlobalOrdinal)gridSize;
    const GlobalOrdinal numRows = gS2*(GlobalOrdinal)gridSize;
    Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap;
    rowMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((global_size_t)numRows,(GlobalOrdinal)0,comm,GloballyDistributed,node));
    // create with DynamicProfile, so that the fillComplete(DoOptimizeStorage) can do first-touch reallocation 
    A = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap,7,Tpetra::DynamicProfile));
    // fill matrix, one row at a time
    Teuchos::Array<GlobalOrdinal> neighbors;
    Teuchos::Array<Scalar> values(7, -ST::one());
    values[0] = (Scalar)6;
    for (GlobalOrdinal r = rowMap->getMinGlobalIndex(); r <= rowMap->getMaxGlobalIndex(); ++r) {
      neighbors.clear();
      neighbors.push_back(r); // add diagonal
      GlobalOrdinal ixy, iz, ix, iy;  // (x,y,z) coords and index in xy plane
      ixy = r%gS2;
      iz = (r - ixy)/gS2;
      ix = ixy%gridSize;
      iy = (ixy - ix)/gridSize;
      //
      if ( ix != 0 )          neighbors.push_back( r-1 );
      if ( ix != gridSize-1 ) neighbors.push_back( r+1 );
      if ( iy != 0 )          neighbors.push_back( r-gridSize );
      if ( iy != gridSize-1 ) neighbors.push_back( r+gridSize );
      if ( iz != 0 )          neighbors.push_back( r-gS2 );
      if ( iz != gridSize-1 ) neighbors.push_back( r+gS2 );
      A->insertGlobalValues( r, neighbors(), values(0,neighbors.size()) );
    }
    A->fillComplete(DoOptimizeStorage);
  }
  else {
    TEST_FOR_EXCEPTION( true, std::runtime_error, 
        "Tpetra::Utils::generateMatrix(): ParameterList specified unsupported ""mat_type"".");
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
Tpetra::Utils::readHBMatrix(const std::string &filename, 
                             const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                             const Teuchos::RCP<Node> &node,
                             Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A,
                             Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap)
{
  const int myRank = comm->getRank();
  int numRows,numCols,numNZ;
  Teuchos::ArrayRCP<Scalar> svals;
  Teuchos::ArrayRCP<GlobalOrdinal> colinds;
  Teuchos::ArrayRCP<int>           rowptrs;
  Teuchos::ArrayRCP<size_t>        nnzPerRow;
  Teuchos::ArrayRCP<char>          type;
  int fail = 0;
  if (myRank == 0) {
    bool isSymmetric=false;
    Teuchos::ArrayRCP<double> dvals;
    Teuchos::ArrayRCP<int> colptrs, rowinds;
    std::string type;
    Tpetra::Utils::readHBMatDouble(filename,numRows,numCols,numNZ,type,colptrs,rowinds,dvals);
    TEST_FOR_EXCEPT(type.size() != 3);
    if (type[0] != 'R' && type[0] != 'r') {
      // only real matrices right now
      fail = 1;
    }
    if (fail == 0 && numNZ > 0) {
      if (type[1] == 'S' || type[1] == 's') {
        isSymmetric = true;
      }
      else {
        isSymmetric = false;
      }
    }
    if (fail == 0 && numNZ > 0) {
      // find num non-zero per row
      nnzPerRow = Teuchos::arcp<size_t>(numRows);
      std::fill(nnzPerRow.begin(), nnzPerRow.end(), 0);
      for (Teuchos::ArrayRCP<int>::const_iterator ri=rowinds.begin(); ri != rowinds.end(); ++ri) {
        // count each row index towards its row
        ++nnzPerRow[*ri-1];
      }
      if (isSymmetric) {
        // count each column toward the corresponding row as well
        for (int c=0; c < numCols; ++c) {
          // the diagonal was already counted; neglect it, if it exists
          for (int i=colptrs[c]-1; i != colptrs[c+1]-1; ++i) {
            if (rowinds[i] != c+1) {
              ++nnzPerRow[c];
              ++numNZ;
            }
          }
        }
      }
      // allocate/set new matrix data
      svals = Teuchos::arcp<Scalar>(numNZ);
      colinds = Teuchos::arcp<GlobalOrdinal>(numNZ);
      rowptrs = Teuchos::arcp<int>(numRows+1);
      rowptrs[0] = 0;
#ifdef HAVE_TPETRA_DEBUG
      Teuchos::ArrayRCP<size_t> nnzPerRow_debug(nnzPerRow.size());
      std::copy(nnzPerRow.begin(), nnzPerRow.end(), nnzPerRow_debug.begin());
#endif
      for (int j=1; j <= numRows; ++j) {
        rowptrs[j] = rowptrs[j-1] + nnzPerRow[j-1];
        nnzPerRow[j-1] = 0;
      }
      // translate from column-oriented to row-oriented
      for (int col=0; col<numCols; ++col) {
        for (int i=colptrs[col]-1; i != colptrs[col+1]-1; ++i) {
          const int row = rowinds[i]-1;
          // add entry to (row,col), with value dvals[i]
          const size_t entry = rowptrs[row] + nnzPerRow[row];
          svals[entry] = Teuchos::as<Scalar>(dvals[i]);
          colinds[entry] = Teuchos::as<GlobalOrdinal>(col);
          ++nnzPerRow[row];
          if (isSymmetric && row != col) {
            // add entry to (col,row), with value dvals[i]
            const size_t symentry = rowptrs[col] + nnzPerRow[col];
            svals[symentry] = Teuchos::as<Scalar>(dvals[i]);
            colinds[symentry] = Teuchos::as<GlobalOrdinal>(row);
            ++nnzPerRow[col];
          }
        }
      }
#ifdef HAVE_TPETRA_DEBUG
      {
        bool isequal = true;
        Teuchos::ArrayRCP<size_t>::const_iterator it1, it2;
        for (it1 = nnzPerRow.begin(), it2 = nnzPerRow_debug.begin(); it1 != nnzPerRow.end(); ++it1, ++it2) {
          if (*it1 != *it2) {
            isequal = false; 
            break;
          }
        }
        TEST_FOR_EXCEPTION(!isequal || nnzPerRow.size() != nnzPerRow_debug.size(), std::logic_error,
            "Tpetra::Utils::readHBMatrix(): Logic error.");
      }
#endif
    }
    // std::cout << "Matrix " << filename << " of type " << type << ": " << numRows << " by " << numCols << ", " << numNZ << " nonzeros" << std::endl;
  }
  // check for read errors
  broadcast(*comm,0,&fail);
  TEST_FOR_EXCEPTION(fail == 1, std::runtime_error, "Tpetra::Utils::readHBMatrix() can only read Real matrices.");
  // distribute global matrix info
  broadcast(*comm,0,&numRows);
  broadcast(*comm,0,&numCols);
  // create map with uniform partitioning
  if (rowMap == Teuchos::null) {
    rowMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((global_size_t)numRows,(GlobalOrdinal)0,comm,GloballyDistributed,node));
  }
  else {
    TEST_FOR_EXCEPTION( rowMap->getGlobalNumElements() != (global_size_t)numRows, std::runtime_error,
        "Tpetra::Utils::readHBMatrix(): specified map has incorrect number of elements.");
    TEST_FOR_EXCEPTION( rowMap->isDistributed() == false && comm->getSize() > 1, std::runtime_error,
        "Tpetra::Utils::readHBMatrix(): specified map is not distributed.");
  }
  Teuchos::ArrayRCP<size_t> myNNZ;
  if (rowMap->getNodeNumElements()) {
    myNNZ = Teuchos::arcp<size_t>(rowMap->getNodeNumElements());
  }
  if (myRank == 0) {
    size_t numRowsAlreadyDistributed = rowMap->getNodeNumElements();
    std::copy(nnzPerRow.begin(), nnzPerRow.begin()+numRowsAlreadyDistributed,myNNZ);
    for (int p=1; p < Teuchos::size(*comm); ++p) {
      size_t numRowsForP;
      Teuchos::receive(*comm,p,&numRowsForP);
      if (numRowsForP) {
        Teuchos::send<int,size_t>(*comm,numRowsForP,nnzPerRow(numRowsAlreadyDistributed,numRowsForP).getRawPtr(),p);
        numRowsAlreadyDistributed += numRowsForP;
      }
    }
  }
  else {
    const size_t numMyRows = rowMap->getNodeNumElements();
    Teuchos::send(*comm,numMyRows,0);
    if (numMyRows) {
      Teuchos::receive<int,size_t>(*comm,0,numMyRows,myNNZ(0,numMyRows).getRawPtr());
    }
  }
  nnzPerRow = Teuchos::null;
  // create column map
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domMap;
  if (numRows == numCols) {
    domMap = rowMap;
  }
  else {
    domMap = createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numCols,comm,node);
  }
  // create with DynamicProfile, so that the fillComplete(DoOptimizeStorage) can do first-touch reallocation 
  A = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap,myNNZ,Tpetra::DynamicProfile));
  // free this locally, A will keep it allocated as long as it is needed by A (up until allocation of nonzeros)
  myNNZ = Teuchos::null;
  if (myRank == 0 && numNZ > 0) {
    for (int r=0; r < numRows; ++r) {
      const size_t nnz = rowptrs[r+1] - rowptrs[r];
      if (nnz > 0) {
        Teuchos::ArrayView<const GlobalOrdinal> inds = colinds(rowptrs[r],nnz);
        Teuchos::ArrayView<const        Scalar> vals = svals(  rowptrs[r],nnz);
        A->insertGlobalValues(r, inds, vals);
      }
    }
  }
  // don't need these anymore
  colinds = Teuchos::null;
  svals   = Teuchos::null;
  rowptrs = Teuchos::null;
  A->fillComplete(domMap,rowMap,Tpetra::DoOptimizeStorage);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
Tpetra::Utils::readMatrixMarketMatrix(const std::string &filename,
  Teuchos::RCP<const Teuchos::Comm<int> > comm,
  Teuchos::RCP<Node> node,
 Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >& A,
 bool transpose,
 Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap,
 bool verbose,
 std::ostream* outstream
 )
{
  const int chunk_read = 500000;  //  Modify this variable to change the
                                  //  size of the chunks read from the file.
  const int headerlineLength = 257;
  const int lineLength = 81;
  const int tokenLength = 35;
  char line[lineLength];
  char token1[tokenLength];
  char token2[tokenLength];
  char token3[tokenLength];
  char token4[tokenLength];
  char token5[tokenLength];
  int M, N, NZ;      // Matrix dimensions
  int me = comm->getRank();

  Teuchos::Time timer("Mtx reader timer", true);

  // Make sure domain and range maps are either both defined or undefined 
  /*if ((domainMap!=0 && rangeMap==0) || (domainMap==0 && rangeMap!=0)) {
    EPETRA_CHK_ERR(-3);
  }*/

  // check maps to see if domain and range are 1-to-1

  /*if (domainMap!=0) {
    if (!domainMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
    if (!rangeMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
  }
  else {
    // If domain and range are not specified, 
    // rowMap becomes both and must be 1-to-1 if specified
    if (rowMap!=0) {
      if (!rowMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
    }
  }*/
      
  FILE * handle = 0;
  if (me == 0) {
    if (verbose) *outstream << 
      "Reading MatrixMarket file " << filename << std::endl;
    handle = fopen(filename.c_str(),"r");  // Open file
	TEST_FOR_EXCEPTION(
    handle == 0, 
    std::runtime_error, 
    "Couldn't open file: " << filename <<std::endl);


    // Check first line, which should be 
    // %%MatrixMarket matrix coordinate real general
	try{
	  TEST_FOR_EXCEPTION(fgets(line, headerlineLength, handle)==0, std::runtime_error,
	  "Bad file header in: " << filename << std::endl);
      TEST_FOR_EXCEPTION(sscanf(line, "%s %s %s %s %s", token1,token2,token3,token4,token5)==0,
	  std::runtime_error,
	  "Bad file header in: " << filename << std::endl);
      TEST_FOR_EXCEPTION(strcmp(token1, "%%MatrixMarket") ||
        strcmp(token2, "matrix") ||
        strcmp(token3, "coordinate") ||
        strcmp(token4, "real") ||
        strcmp(token5, "general"),
	  std::runtime_error,
	  "Bad file header in: " << filename << std::endl);
	}
	catch(std::runtime_error& e){
      if (handle!=0) fclose(handle); 
	  throw std::runtime_error(e.what());
	}
    // Next, strip off header lines (which start with "%")
    do {
	  try{
        TEST_FOR_EXCEPTION(fgets(line, headerlineLength, handle)==0,
	      std::runtime_error,
	      "Bad file header in: " << filename << std::endl);
	  }
	  catch(std::runtime_error& e){
        if (handle!=0) fclose(handle); 
	    throw std::runtime_error(e.what());
	  }
    } while (line[0] == '%');

    // Next get problem dimensions: M, N, NZ
	try{
      TEST_FOR_EXCEPTION(sscanf(line, "%d %d %d", &M, &N, &NZ)==0,
	    std::runtime_error,
	    "Couldn't get problem dimensions in file: " << filename << std::endl);
	}
	catch(std::runtime_error& e){
      if (handle!=0) fclose(handle); 
	  throw std::runtime_error(e.what());
	}
  }
  Teuchos::broadcast(*comm,  0, 1, &M);
  Teuchos::broadcast(*comm,  0, 1, &N);
  Teuchos::broadcast(*comm,  0, 1, &NZ);

  // Now create matrix using user maps if provided.


  // Now read in chunks of triplets and broadcast them to other processors.
  char *buffer = new char[chunk_read*lineLength];
  int nchunk; 
  int nmillion = 0;
  int nread = 0;
  int rlen;

  // Storage for this processor's nonzeros.
  const int localblock = 100000;
  size_t localsize = NZ / comm->getSize() + localblock;

  Teuchos::ArrayRCP<GlobalOrdinal> iv(localsize);
  Teuchos::ArrayRCP<GlobalOrdinal> jv(localsize);
  Teuchos::ArrayRCP<Scalar> vv(localsize);
  size_t lnz = Teuchos::OrdinalTraits<size_t>::zero();   //  Number of non-zeros on this processor.

  Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap1 = Teuchos::null;
  bool allocatedHere=false;
  if (!rowMap.is_null()) 
    rowMap1 = Teuchos::rcp_const_cast<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >(rowMap);
  else {
	  rowMap1 = Teuchos::rcp(
      new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(
        (global_size_t)M, 
        Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), 
        comm));
    allocatedHere = true;
  }
  GlobalOrdinal ioffset = rowMap1->getIndexBase()-Teuchos::OrdinalTraits<GlobalOrdinal>::one();
  int joffset = ioffset; 

  int rowmajor = 1;  // Assume non-zeros are listed in row-major order, 
  int prevrow = -1;  // but test to detect otherwise.  If non-zeros
                     // are row-major, we can avoid the sort.


  while (nread < NZ) {
    if (NZ-nread > chunk_read) nchunk = chunk_read;
    else nchunk = NZ - nread;

    if (me == 0) {
      char *eof;
      rlen = 0;
      for (int i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[rlen],lineLength,handle);
        if (eof == NULL) {
          fprintf(stderr, "%s", "Unexpected end of matrix file.");
		  throw std::runtime_error("Unexpected end of matrix file.");
        }
        rlen += strlen(&buffer[rlen]);
      }
      buffer[rlen++]='\n';
    }
	  Teuchos::broadcast(*comm, 0, 1, &rlen);
      comm->broadcast(0, rlen, buffer);

    buffer[rlen++] = '\0';
    nread += nchunk;

    // Process the received data, saving non-zeros belonging on this proc.
    char *lineptr = buffer;
    for (rlen = 0; rlen < nchunk; rlen++) {
      char *next = strchr(lineptr,'\n');
      GlobalOrdinal I = atoi(strtok(lineptr," \t\n")) + ioffset;
      GlobalOrdinal J = atoi(strtok(NULL," \t\n")) + joffset;
      Scalar V = atof(strtok(NULL," \t\n"));
      lineptr = next + 1;
      /*if (transpose) {
        // swap I and J indices.
        int tmp = I;
        I = J;
        J = tmp;
      }*/
      if (rowMap1->isNodeGlobalElement(I)) {
        //  This processor keeps this non-zero.
        if (lnz >= localsize) {  
          // Need more memory to store nonzeros.
          localsize += localblock;
          iv.resize(localsize);
          jv.resize(localsize);
          vv.resize(localsize);
        }
        iv[lnz] = I;
        jv[lnz] = J;
        vv[lnz] = V;
        lnz++;
        if (I < prevrow) rowmajor = 0;
        prevrow = I;
      }
    }
    // Status check
    if (nread / 1000000 > nmillion) {
      nmillion++;
      if (verbose && me == 0) *outstream << nmillion << "M ";
    }
  }

  delete [] buffer;

  if (!rowmajor) {
    // Sort into row-major order (by iv) so can insert entire rows at once.
    // Reorder jv and vv to parallel iv.
    if (verbose && me == 0) *outstream << std::endl << "   Sorting local nonzeros" << std::endl;
	sort3(iv.begin(), iv.end(), jv.begin(), vv.begin());
  }

  //  Compute number of entries per local row for use in constructor;
  //  saves reallocs in FillComplete.

  //  Now construct the matrix.
  //  Count number of entries in each row so can use StaticProfile=2.
  if (verbose && me == 0) *outstream << std::endl << "   Constructing the matrix" << std::endl;
  size_t numRows = rowMap1->getNodeNumElements();
  Teuchos::ArrayRCP<size_t> numNonzerosPerRow(numRows);
  for (size_t i = Teuchos::OrdinalTraits<size_t>::zero(); i < lnz; i++) 
    numNonzerosPerRow[rowMap1->getLocalElement(iv[i])]++;

  A = Teuchos::rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap1, numNonzerosPerRow));
  //A->SetTracebackMode(1);

  // Rows are inserted in ascending global number, and the insertion uses numNonzerosPerRow.
  // Therefore numNonzerosPerRow must be permuted, as it was created in ascending local order.
  //Teuchos::ArrayView<GlobalOrdinal> gidList = rowMap1->getNodeElementList(); 
  //sort2(gidList.begin(), gidList.end(), numNonzerosPerRow.begin());

  //  Insert the global values into the matrix row-by-row.
  if (verbose && me == 0) *outstream << "   Inserting global values" << std::endl;
  for (size_t sum = Teuchos::OrdinalTraits<size_t>::zero(), i = Teuchos::OrdinalTraits<size_t>::zero(); i < numRows; i++) {
    if (numNonzerosPerRow[i]) {
       A->insertGlobalValues(iv[sum],  
                             jv(sum,numNonzerosPerRow[i]),
                             vv(sum, numNonzerosPerRow[i]));
      sum += numNonzerosPerRow[i];
    }
  }

  if (verbose && me == 0) *outstream << "   Completing matrix fill" << std::endl;
  /*if (rangeMap != 0 && domainMap != 0) {
    EPETRA_CHK_ERR(A->FillComplete(*domainMap, *rangeMap));
  }
  else if (M!=N) {
    Epetra_Map newDomainMap(N,rowMap1->IndexBase(), comm);
    EPETRA_CHK_ERR(A->FillComplete(newDomainMap, *rowMap1));
  }
  else {
    EPETRA_CHK_ERR(A->FillComplete());
  }*/

  A->fillComplete();
  
  if (handle!=0) fclose(handle);
  timer.stop();
  double dt = timer.totalElapsedTime();
  if (verbose && me == 0) *outstream << "File Read time (secs):  " << dt << std::endl;
  //return(0);

}


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Utils namespace!
//

#define TPETRA_MATRIXIO_INSTANT(SCALAR,LO,GO,NODE)                                                                                                      \
  template                                                                                                                                              \
  void                                                                                                                                                  \
  readHBMatrix<SCALAR,LO,GO,NODE,Kokkos::DefaultKernels<SCALAR,LO,NODE>::SparseOps>(                                                                    \
          const std::string &, const Teuchos::RCP<const Teuchos::Comm<int> > &, const Teuchos::RCP<NODE > &,                                            \
          Teuchos::RCP< CrsMatrix<SCALAR,LO,GO,NODE,Kokkos::DefaultKernels<SCALAR,LO,NODE>::SparseOps > > &,                                            \
          Teuchos::RCP< const Tpetra::Map<LO,GO,NODE> >);                                                                                               \
                                                                                                                                                        \
  template                                                                                                                                              \
  void                                                                                                                                                  \
  generateMatrix<SCALAR,LO,GO,NODE,Kokkos::DefaultKernels<SCALAR,LO,NODE>::SparseOps>(                                                                  \
          const Teuchos::RCP<Teuchos::ParameterList> &plist, const Teuchos::RCP<const Teuchos::Comm<int> > &, const Teuchos::RCP<NODE > &,              \
          Teuchos::RCP< CrsMatrix<SCALAR,LO,GO,NODE,Kokkos::DefaultKernels<SCALAR,LO,NODE>::SparseOps > > &);

#endif
