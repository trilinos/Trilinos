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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>

int main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  typedef double Scalar;
  typedef Tpetra::CrsMatrix<Scalar> crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  bool success = false;
  bool verbose = false;
  try {
    const GO numGlobalElements = 256; // problem size
    const GO indexBase = 0;
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
    RCP<const map_type> map = rcp (new map_type (numGlobalElements, indexBase, comm));

    crs_matrix_type A (map, 3);

    Teuchos::ArrayView<const GO> myGlobalElements = map->getNodeElementList();
    const size_t numMyElements = map->getNodeNumElements ();
    for (size_t i = 0; i < numMyElements; ++i) {
      if (myGlobalElements[i] == 0) {
        A.insertGlobalValues (myGlobalElements[i],
                              tuple<GO> (myGlobalElements[i], myGlobalElements[i] + 1),
                              tuple<Scalar> (2.0, -1.0));
      }
      else if (myGlobalElements[i] == numGlobalElements - 1) {
        A.insertGlobalValues (myGlobalElements[i],
                              tuple<GO> (myGlobalElements[i] - 1, myGlobalElements[i]),
                              tuple<Scalar> (-1.0, 2.0));
      }
      else {
        A.insertGlobalValues (myGlobalElements[i],
                              tuple<GO> (myGlobalElements[i] - 1, myGlobalElements[i], myGlobalElements[i] + 1),
                              tuple<Scalar> (-1.0, 2.0, -1.0));
      }
    }

    A.fillComplete ();
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
