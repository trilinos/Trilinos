/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

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
  using Teuchos::as;
  TEUCHOS_TEST_FOR_EXCEPTION( plist == Teuchos::null, std::runtime_error,
      "Tpetra::Utils::generateMatrix(): ParameterList is null.");
  TEUCHOS_TEST_FOR_EXCEPTION( Teuchos::isParameterType<std::string>(*plist,"mat_type") == false, std::runtime_error,
      "Tpetra::Utils::generateMatrix(): ParameterList did not contain string parameter ""mat_type"".");
  std::string mat_type = plist->get<std::string>("mat_type");
  if (mat_type == "Lap3D") {
    // 3D Laplacian, grid is a cube with dimension gridSize x gridSize x gridSize
    const GlobalOrdinal gridSize = as<GlobalOrdinal>(plist->get<int>("gridSize",100));
    const GlobalOrdinal gS2 = gridSize*gridSize;
    const GlobalOrdinal numRows = gS2*gridSize;
    Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap;
    rowMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(as<global_size_t>(numRows),as<GlobalOrdinal>(0),comm,GloballyDistributed,node));
    A = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap,7,Tpetra::StaticProfile));
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
    A->fillComplete();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, 
        "Tpetra::Utils::generateMatrix(): ParameterList specified unsupported ""mat_type"".");
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
Tpetra::Utils::readHBMatrix(const std::string &filename, 
                             const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                             const Teuchos::RCP<Node> &node,
                             Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A,
                             Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap, 
                             const Teuchos::RCP<ParameterList> &params)
{
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
    rowMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((global_size_t)numRows,(GlobalOrdinal)0,comm,GloballyDistributed,node));
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap->getGlobalNumElements() != (global_size_t)numRows, std::runtime_error,
        "Tpetra::Utils::readHBMatrix(): specified map has incorrect number of elements.");
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap->isDistributed() == false && comm->getSize() > 1, std::runtime_error,
        "Tpetra::Utils::readHBMatrix(): specified map is not distributed.");
  }
  Teuchos::ArrayRCP<size_t> myNNZ;
  if (rowMap->getNodeNumElements()) {
    myNNZ = Teuchos::arcp<size_t>(rowMap->getNodeNumElements());
  }
  if (myRank == 0) {
    LocalOrdinal numRowsAlreadyDistributed = rowMap->getNodeNumElements();
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
  A = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap,myNNZ,Tpetra::StaticProfile));
  // free this locally, A will keep it allocated as long as it is needed by A (up until allocation of nonzeros)
  myNNZ = Teuchos::null;
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
          Teuchos::RCP< const Tpetra::Map<LO,GO,NODE> >,                                                                                                \
          const Teuchos::RCP< Teuchos::ParameterList > & );                                                                                             \
                                                                                                                                                        \
  template                                                                                                                                              \
  void                                                                                                                                                  \
  generateMatrix<SCALAR,LO,GO,NODE,Kokkos::DefaultKernels<SCALAR,LO,NODE>::SparseOps>(                                                                  \
          const Teuchos::RCP<Teuchos::ParameterList> &plist, const Teuchos::RCP<const Teuchos::Comm<int> > &, const Teuchos::RCP<NODE > &,              \
          Teuchos::RCP< CrsMatrix<SCALAR,LO,GO,NODE,Kokkos::DefaultKernels<SCALAR,LO,NODE>::SparseOps > > &);

#endif
