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

#include "Tpetra_Core.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include <Teuchos_FancyOStream.hpp>
#include <cmath>

template<class CrsMatrix_t> double getNorm(CrsMatrix_t& matrix){
  double mySum = 0;

  typedef typename CrsMatrix_t::local_ordinal_type LO;

  Teuchos::Array<LO> inds(matrix.getNodeMaxNumRowEntries());
  Teuchos::Array<double> vals(matrix.getNodeMaxNumRowEntries());
  for(int i =0; ((size_t)i)<matrix.getNodeNumRows(); ++i){
    size_t numRowEnts = matrix.getNumEntriesInLocalRow(i);
    Teuchos::ArrayView<const LO> indsView = inds();
    Teuchos::ArrayView<const double> valsView = vals();
    matrix.getLocalRowView(i, indsView, valsView);
    for(size_t j=0; ((size_t)j)<numRowEnts; ++j){
      mySum += valsView[j]*valsView[j];
    }
  }
  double totalSum = 0;
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
  typedef double Scalar;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> crs_matrix_type;
  typedef Tpetra::Map<LO, GO> map_type;

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

  const size_t numMyElements = map->getNodeNumElements();

  Teuchos::ArrayView<const GO> myGlobalElements = map->getNodeElementList();

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

  // Create a Tpetra::Matrix using the Map, with a static allocation dictated by NumNz
  crs_matrix_type A (map, NumNz, Tpetra::StaticProfile);
  crs_matrix_type AT(map, NumNz, Tpetra::StaticProfile);
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
    Tpetra::createCrsMatrix<Scalar, LO, GO> (TestMatrix->getRowMap ());

  // Apparently there is a problem with ADD because while these two matrices are
  // identical when I add them together I don't get 0 like I should. In fact
  // I just get a matrix that has the exact same entries and sparsity structure.
  // I'll have to come back to this later. But RowMatrixTransposer is working right.
  // And all my other tests are telling me ADD works right too.
  // KLN 06/14/2011
  Tpetra::MatrixMatrix::Add(AT,false,-1.0,*TestMatrix,false, 1.0,diffMatrix);
  diffMatrix->fillComplete();
  //diffMatrix->describe(*out, Teuchos::VERB_EXTREME);
  double diffNorm = getNorm(*diffMatrix);
  double realNorm = getNorm(AT);
  double epsilon = diffNorm/realNorm;
  if(epsilon > 1e-10){
    *out << "The calculated A transpose and the real one don't match!" << std::endl;
    *out << "Diff Norm: " << diffNorm << std::endl;
    *out << "Real Norm: " << realNorm << std::endl;
    *out << "Epsilon: " << epsilon << std::endl;
    return 1;
  }

  return 0;
}




