// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  typedef Tpetra::CrsMatrix<>::scalar_type Scalar;
  typedef Tpetra::CrsMatrix<Scalar> crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  bool success = false;
  bool verbose = false;
  try {
    const GO numGlobalElements          = 256;  // problem size
    const GO indexBase                  = 0;
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<const map_type> map             = rcp(new map_type(numGlobalElements, indexBase, comm));

    crs_matrix_type A(map, 3);

    Teuchos::ArrayView<const GO> myGlobalElements = map->getLocalElementList();
    const size_t numMyElements                    = map->getLocalNumElements();
    for (size_t i = 0; i < numMyElements; ++i) {
      if (myGlobalElements[i] == 0) {
        A.insertGlobalValues(myGlobalElements[i],
                             tuple<GO>(myGlobalElements[i], myGlobalElements[i] + 1),
                             tuple<Scalar>(2.0, -1.0));
      } else if (myGlobalElements[i] == numGlobalElements - 1) {
        A.insertGlobalValues(myGlobalElements[i],
                             tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i]),
                             tuple<Scalar>(-1.0, 2.0));
      } else {
        A.insertGlobalValues(myGlobalElements[i],
                             tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i], myGlobalElements[i] + 1),
                             tuple<Scalar>(-1.0, 2.0, -1.0));
      }
    }

    A.fillComplete();
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
