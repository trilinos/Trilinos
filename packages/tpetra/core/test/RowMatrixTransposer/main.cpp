// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include <Teuchos_FancyOStream.hpp>
#include <cmath>

template<class CrsMatrix_t>
typename CrsMatrix_t::scalar_type getNorm(CrsMatrix_t& matrix){
  typedef typename CrsMatrix_t::scalar_type Scalar;
  Scalar mySum = 0;

  for(int i =0; ((size_t)i)<matrix.getLocalNumRows(); ++i){
    size_t numRowEnts = matrix.getNumEntriesInLocalRow(i);
    typename CrsMatrix_t::local_inds_host_view_type indsView;
    typename CrsMatrix_t::values_host_view_type valsView;
    matrix.getLocalRowView(i, indsView, valsView);
    for(size_t j=0; ((size_t)j)<numRowEnts; ++j){
      mySum += valsView[j]*valsView[j];
    }
  }
  Scalar totalSum = 0;
  Teuchos::reduceAll(*(matrix.getComm()), Teuchos::REDUCE_SUM, 1, &mySum, &totalSum);
  return sqrt(totalSum);
}



int
main (int argc, char* argv[])
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::tuple;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::CrsMatrix<> crs_matrix_type;
  typedef typename crs_matrix_type::scalar_type Scalar;
  typedef Tpetra::Map<> map_type;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<Teuchos::FancyOStream> out =
    Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  //out->setOutputToRootOnly(comm->getRank());

  size_t myRank = comm->getRank();
  size_t numProc = comm->getSize();
  bool verbose = (myRank==0);

  std::cout << *comm;

  const global_size_t numGlobalElements = 4;
  if (numGlobalElements < numProc) {
    if (verbose) {
      std::cout << "numGlobalBlocks = " << numGlobalElements
                << " cannot be less than the number of processors = " << numProc << std::endl;
    }
    return -1;
  }

  // Construct a Map that puts approximately the same number of equations on each processor.

  RCP<const map_type> map =
    Tpetra::createUniformContigMap<LO, GO> (numGlobalElements, comm);

  // Get update list and number of local equations from newly created map.

  const size_t numMyElements = map->getLocalNumElements();

  Teuchos::ArrayView<const GO> myGlobalElements = map->getLocalElementList();

  // Create an OTeger vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor

  Teuchos::ArrayRCP<size_t> NumNz = Teuchos::arcp<size_t>(numMyElements);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (size_t i=0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0 ||
        static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements - 1) {
      // boundary
      NumNz[i] = 2;
    }
    else {
      NumNz[i] = 3;
    }
  }

  // Create a Tpetra::CrsMatrix using the Map, with a static allocation dictated by NumNz
  crs_matrix_type A (map, NumNz ());
  crs_matrix_type AT(map, NumNz ());
  RCP< crs_matrix_type > TestMatrix = Teuchos::null;

  // We are done with NumNZ
  NumNz = Teuchos::null;

  // Add  rows one-at-a-time
  // Off diagonal values will always be -1
  const Scalar two    = static_cast<Scalar>( 2.0);
  const Scalar negOne = static_cast<Scalar>(-1.0);
  const Scalar three  = static_cast<Scalar>(3.0);
  for (size_t i = 0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0) {
      A.insertGlobalValues (myGlobalElements[i],
                            tuple<GO> (myGlobalElements[i],
                                       myGlobalElements[i]+1),
                            tuple<Scalar> (two, negOne));
    }
    else if (static_cast<global_size_t> (myGlobalElements[i]) ==
             numGlobalElements - 1) {
      A.insertGlobalValues (myGlobalElements[i],
                            tuple<GO> (myGlobalElements[i]-1,
                                       myGlobalElements[i]),
                            tuple<Scalar> (negOne, two));
    }
    else {
      A.insertGlobalValues (myGlobalElements[i],
                            tuple<GO> (myGlobalElements[i]-1,
                                       myGlobalElements[i],
                                       myGlobalElements[i]+1),
                            tuple<Scalar> (three, two, negOne));
    }
  }

  for (size_t i=0; i<numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      AT.insertGlobalValues (myGlobalElements[i],
                             tuple<GO> (myGlobalElements[i],
                                        myGlobalElements[i] + 1),
                             tuple<Scalar> (two, three));
    }
    else if (static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
      AT.insertGlobalValues (myGlobalElements[i],
                             tuple<GO> (myGlobalElements[i] - 1,
                                        myGlobalElements[i]),
                             tuple<Scalar> (negOne, two));
    }
    else if(static_cast<global_size_t>(myGlobalElements[i])==1){
      AT.insertGlobalValues( myGlobalElements[i],
                             tuple<GO> (myGlobalElements[i]-1,
                                        myGlobalElements[i],
                                        myGlobalElements[i]+1),
                             tuple<Scalar> (negOne, two, three));
    }
    else if(static_cast<global_size_t>(myGlobalElements[i])==2){
      AT.insertGlobalValues( myGlobalElements[i],
                             tuple<GO> (myGlobalElements[i]-1,
                                        myGlobalElements[i],
                                        myGlobalElements[i]+1),
                             tuple<Scalar> (negOne, two, negOne));
    }
  }

  // Finish up
  A.fillComplete();
  AT.fillComplete();

  Tpetra::RowMatrixTransposer<Scalar, LO, GO> transposer (Teuchos::rcpFromRef (A));
  TestMatrix = transposer.createTranspose(); //, TestMatrix/*, tMap*/);

  RCP<crs_matrix_type > diffMatrix =
    Tpetra::createCrsMatrix<Scalar, LO, GO> (TestMatrix->getRowMap (),
                                             AT.getGlobalMaxNumRowEntries());

  // Apparently there is a problem with ADD because while these two matrices are
  // identical when I add them together I don't get 0 like I should. In fact
  // I just get a matrix that has the exact same entries and sparsity structure.
  // I'll have to come back to this later. But RowMatrixTransposer is working right.
  // And all my other tests are telling me ADD works right too.
  // KLN 06/14/2011
  Tpetra::MatrixMatrix::Add(AT,false,static_cast<Scalar>(-1.0),
                            *TestMatrix,false,static_cast<Scalar>(1.0),
                            diffMatrix);
  diffMatrix->fillComplete();
  //diffMatrix->describe(*out, Teuchos::VERB_EXTREME);
  Scalar diffNorm = getNorm(*diffMatrix);
  Scalar realNorm = getNorm(AT);
  Scalar epsilon = diffNorm/realNorm;
  if(epsilon > 1e-10){
    *out << "The calculated A transpose and the real one don't match!" << std::endl;
    *out << "Diff Norm: " << diffNorm << std::endl;
    *out << "Real Norm: " << realNorm << std::endl;
    *out << "Epsilon: " << epsilon << std::endl;
    return 1;
  }

  return 0;
}




