// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIX_IO_DEF
#define TPETRA_MATRIX_IO_DEF

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include <iostream>

namespace Tpetra {
namespace Utils {


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
readHBMatrix (const std::string &filename,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &A,
              Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap,
              const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  const int myRank = comm->getRank();
  int numRows,numCols,numNZ;
  Teuchos::ArrayRCP<Scalar> svals;
  Teuchos::ArrayRCP<GlobalOrdinal> colinds;
  Teuchos::ArrayRCP<int>           rowptrs;
  Teuchos::ArrayRCP<size_t>        nnzPerRow;
  int fail = 0;
  if (myRank == 0) {
    bool isSymmetric=false;
    Teuchos::ArrayRCP<double> dvals;
    Teuchos::ArrayRCP<int> colptrs, rowinds;
    std::string type;
    Tpetra::Utils::readHBMatDouble(filename,numRows,numCols,numNZ,type,colptrs,rowinds,dvals);
    TEUCHOS_TEST_FOR_EXCEPT(type.size() != 3);
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
        typename Teuchos::ArrayRCP<size_t>::const_iterator it1, it2;
        for (it1 = nnzPerRow.begin(), it2 = nnzPerRow_debug.begin(); it1 != nnzPerRow.end(); ++it1, ++it2) {
          if (*it1 != *it2) {
            isequal = false;
            break;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(!isequal || nnzPerRow.size() != nnzPerRow_debug.size(), std::logic_error,
            "Tpetra::Utils::readHBMatrix(): Logic error.");
      }
#endif
    }
    // std::cout << "Matrix " << filename << " of type " << type << ": " << numRows << " by " << numCols << ", " << numNZ << " nonzeros" << std::endl;
  }
  // check for read errors
  broadcast(*comm,0,&fail);
  TEUCHOS_TEST_FOR_EXCEPTION(fail == 1, std::runtime_error, "Tpetra::Utils::readHBMatrix() can only read Real matrices.");
  // distribute global matrix info
  broadcast(*comm,0,&numRows);
  broadcast(*comm,0,&numCols);
  // create map with uniform partitioning
  if (rowMap == Teuchos::null) {
    rowMap = Teuchos::rcp (new map_type (static_cast<global_size_t> (numRows),
                                         static_cast<GlobalOrdinal> (0),
                                         comm, GloballyDistributed));
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap->getGlobalNumElements() != (global_size_t)numRows, std::runtime_error,
        "Tpetra::Utils::readHBMatrix(): specified map has incorrect number of elements.");
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap->isDistributed() == false && comm->getSize() > 1, std::runtime_error,
        "Tpetra::Utils::readHBMatrix(): specified map is not distributed.");
  }
  Teuchos::Array<size_t> myNNZ (rowMap->getLocalNumElements ());
  if (myRank == 0) {
    LocalOrdinal numRowsAlreadyDistributed = rowMap->getLocalNumElements();
    std::copy(nnzPerRow.begin(), nnzPerRow.begin()+numRowsAlreadyDistributed, myNNZ.begin());
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
    const size_t numMyRows = rowMap->getLocalNumElements();
    Teuchos::send(*comm,numMyRows,0);
    if (numMyRows) {
      Teuchos::receive<int,size_t>(*comm,0,numMyRows,myNNZ(0,numMyRows).getRawPtr());
    }
  }
  nnzPerRow = Teuchos::null;
  // create column map
  Teuchos::RCP<const map_type> domMap;
  if (numRows == numCols) {
    domMap = rowMap;
  }
  else {
    domMap = createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numCols,comm);
  }
  A = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap, myNNZ ()));
  // free this locally, A will keep it allocated as long as it is needed by A (up until allocation of nonzeros)
  {
    // Classic idiom for freeing an std::vector; resize doesn't
    // promise deallocation.
    Teuchos::Array<size_t> empty;
    std::swap (myNNZ, empty);
  }
  if (myRank == 0 && numNZ > 0) {
    for (int r=0; r < numRows; ++r) {
      const LocalOrdinal nnz = rowptrs[r+1] - rowptrs[r];
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
  A->fillComplete(domMap,rowMap,params);
}


} // namespace Utils
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Utils namespace!
//


#define TPETRA_MATRIXIO_INSTANT(SCALAR,LO,GO,NODE) \
  template void \
  readHBMatrix< SCALAR, LO, GO, NODE > (const std::string&, \
                                        const Teuchos::RCP<const Teuchos::Comm<int> > &, \
                                        Teuchos::RCP<CrsMatrix< SCALAR, LO, GO, NODE > >&, \
                                        Teuchos::RCP<const Tpetra::Map< LO, GO, NODE> >, \
                                        const Teuchos::RCP<Teuchos::ParameterList>& ); 



#endif
