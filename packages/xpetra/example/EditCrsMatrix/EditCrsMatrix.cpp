// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Xpetra_Parameters.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVector.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>

//#include <Xpetra_EditCrsMatrix.hpp>

typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;

int main(int argc, char *argv[]) {
  GlobalOrdinal numGlobalElements = 10; // problem size

  using Teuchos::RCP;
  using Teuchos::rcp;
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
 
  //
  // Parse command line arguments
  //

  Xpetra::UnderlyingLib lib;

  {
    Teuchos::CommandLineProcessor clp(false); // Note: use --help to list available options.
    Xpetra::Parameters xpetraParameters(clp);
    switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }
    
    lib = xpetraParameters.GetLib();
  }

  //
  // Build matrix
  //

  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal>::createUniformContigMap(lib, numGlobalElements, comm);

  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getNodeElementList();

  RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> > A =  Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (size_t i = 0; i < numMyElements; i++) {
     if (myGlobalElements[i] == 0) {
       A->insertGlobalValues(myGlobalElements[i], 
                             Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] +1), 
                             Teuchos::tuple<Scalar> (2.0, -1.0));
     }
     else if (myGlobalElements[i] == numGlobalElements - 1) {
       A->insertGlobalValues(myGlobalElements[i], 
                             Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]), 
                             Teuchos::tuple<Scalar> (-1.0, 2.0));
     }
     else {
       A->insertGlobalValues(myGlobalElements[i], 
                             Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1), 
                             Teuchos::tuple<Scalar> (-1.0, 2.0, -1.0));
     }
  }

#define TEST_REPLACEGLOBALVALUES
#ifdef  TEST_REPLACEGLOBALVALUES
  {
    Teuchos::ArrayView<const GlobalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;

    A->getGlobalRowView(0, indices, values);

    const size_t nnz = indices.size();
    Teuchos::Array<Scalar> newValues(nnz, 0.0);

    for (size_t j = 0; j < nnz; j++) {
      newValues[j] = 42;
    }

    A->replaceGlobalValues(0, indices, newValues);
  }
#endif // TEST_REPLACEGLOBALVALUES

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  //
  // Edit matrix
  //

  const LocalOrdinal nodeNumRows = A->getNodeNumRows();

  A->resumeFill();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

#define TEST_REPLACELOCALVALUES
#ifdef  TEST_REPLACELOCALVALUES
  for (LocalOrdinal i = 0; i < nodeNumRows; i++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;

    A->getLocalRowView(i, indices, values);

    const size_t nnz = indices.size();
    Teuchos::Array<Scalar> newValues(nnz, 0.0);

    for (size_t j = 0; j < nnz; j++) {
      if (indices[j] == i) /* diagonal term */
        newValues[j] = values[j] * 2;
      else
        newValues[j] = values[j];
    }

    A->replaceLocalValues(i, indices, newValues);
  }
#endif // TEST_REPLACELOCALVALUES

#ifdef TEST2_REPLACELOCALVALUES
  for (LocalOrdinal i = 0; i < nodeNumRows; i++) {
    Teuchos::Array<LocalOrdinal> newIndices(1, i);
    Teuchos::Array<Scalar> newValues(1, 4.0);

    A->replaceLocalValues(i, newIndices, newValues);
  }
#endif // TEST2_REPLACELOCALVALUES

  //
  //
  //
 
  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  //
  //
  //

  A->describe(*out,Teuchos::VERB_EXTREME);
  
  return EXIT_SUCCESS;
}

//TODO: transform this example to a unit test
