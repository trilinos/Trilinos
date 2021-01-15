// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_Utilities.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TestFactory, CreatePoisson1DMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);

    out << "version: " << MueLu::Version() << std::endl;

    const GO numNodes = 29;
    RCP<Matrix> matrix = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(numNodes);
    TEST_ASSERT(!matrix.is_null());

    TEST_EQUALITY_CONST(matrix->getGlobalNumRows(), numNodes);

    const size_t numEntriesPerRow = 3;
    const global_size_t expectedGlobalNumEntries = numEntriesPerRow * numNodes - 2;
    TEST_EQUALITY_CONST(matrix->getGlobalNumEntries(), expectedGlobalNumEntries);

    // Check every single matrix entry
    for (size_t lRow = 0; lRow < matrix->getNodeNumRows(); ++lRow)
    {
      ArrayView<const LocalOrdinal> cols;
      ArrayView<const Scalar> vals;
      matrix->getLocalRowView(1, cols, vals);

      TEST_EQUALITY_CONST(cols.size(), numEntriesPerRow);
      TEST_EQUALITY_CONST(vals.size(), numEntriesPerRow);
      TEST_EQUALITY_CONST(vals[0], Teuchos::as<Scalar>(-1.0));
      TEST_EQUALITY_CONST(vals[1], Teuchos::as<Scalar>(2.0));
      TEST_EQUALITY_CONST(vals[2], Teuchos::as<Scalar>(-1.0));
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TestFactory, CreatePoisson2DMeshtyingMatrix_Obtain2x2SaddlePointSystem, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra);
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);

    out << "version: " << MueLu::Version() << std::endl;

    const Scalar SC_zero = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar SC_one = Teuchos::ScalarTraits<Scalar>::one();

    const GO numNodesX = 3;
    const GO numNodesY = 3;
    RCP<const BlockedCrsMatrix> blockedMatrix = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildPoissonSaddlePointMatrix(numNodesX, numNodesY);
    TEST_ASSERT(!blockedMatrix.is_null())

    const size_t numSubDomains = 2;
    const GO numNodesPrimal = numSubDomains * (numNodesX * numNodesY);
    const GO numNodesDual = numNodesX; // Interface is oriented along the global x-axis

    TEST_EQUALITY_CONST(blockedMatrix->Rows(), numSubDomains);
    TEST_EQUALITY_CONST(blockedMatrix->Cols(), numSubDomains);

    TEST_EQUALITY_CONST(blockedMatrix->getMatrix(0, 0)->getGlobalNumRows(), numNodesPrimal);
    TEST_EQUALITY_CONST(blockedMatrix->getMatrix(0, 1)->getGlobalNumRows(), numNodesPrimal);
    TEST_EQUALITY_CONST(blockedMatrix->getMatrix(1, 0)->getGlobalNumRows(), numNodesDual);

    // Saddle-point system needs to have a zero block
    TEST_ASSERT(blockedMatrix->getMatrix(1, 1).is_null());

    TEST_EQUALITY_CONST(blockedMatrix->getGlobalNumRows(), numNodesPrimal + numNodesDual);

    // Kinematic coupling: rows have entries [1, -1] or [-1, 1]
    RCP<const Matrix> A10 = blockedMatrix->getMatrix(1, 0);
    TEST_ASSERT(!A10.is_null());
    for (LocalOrdinal lDualRow = 0; lDualRow < Teuchos::as<LocalOrdinal>(A10->getNodeNumRows()); ++lDualRow)
    {
      ArrayView<const LocalOrdinal> cols;
      ArrayView<const Scalar> vals;
      A10->getLocalRowView(lDualRow, cols, vals);

      TEST_EQUALITY_CONST(cols.size(), 2);
      TEST_EQUALITY_CONST(vals.size(), 2);
      TEST_INEQUALITY_CONST(vals[0], SC_zero);
      TEST_INEQUALITY_CONST(vals[1], SC_zero);
      TEST_EQUALITY_CONST(vals[0] * vals[0], SC_one);
      TEST_EQUALITY_CONST(vals[1] * vals[1], SC_one);
      TEST_EQUALITY_CONST(vals[0] + vals[1], SC_zero);
    }

    // Coupling of fluxes: rows have entries [1] or [-1]
    RCP<const Matrix> A01 = blockedMatrix->getMatrix(0, 1);
    TEST_ASSERT(!A01.is_null());
    for (LocalOrdinal lPrimalRow = 0; lPrimalRow < Teuchos::as<LocalOrdinal>(A01->getNodeNumRows()); ++lPrimalRow)
    {
      ArrayView<const LocalOrdinal> cols;
      ArrayView<const Scalar> vals;
      A01->getLocalRowView(lPrimalRow, cols, vals);

      if (cols.size() > 0)
      {
        TEST_EQUALITY_CONST(cols.size(), 1);
        TEST_EQUALITY_CONST(vals.size(), 1);
        TEST_EQUALITY_CONST(vals[0] * vals[0], SC_one);
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TestFactory, CreatePoisson2DMeshtyingMatrix_Obtain2x2SaddlePointSystemAndInterfaceDofMap, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra);
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);

    out << "version: " << MueLu::Version() << std::endl;

    const Scalar SC_zero = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar SC_one = Teuchos::ScalarTraits<Scalar>::one();

    const GO numNodesX = 3;
    const GO numNodesY = 3;
    RCP<BlockedCrsMatrix> blockedMatrix = Teuchos::null;
    RCP<const Map> interfaceDofMap = Teuchos::null;
    TestHelpers::TestFactory<SC, LO, GO, NO>::BuildPoissonSaddlePointMatrix(blockedMatrix, interfaceDofMap, numNodesX, numNodesY);
    TEST_ASSERT(!blockedMatrix.is_null());
    TEST_ASSERT(!interfaceDofMap.is_null());

    const size_t numSubDomains = 2;
    const GO numNodesPrimal = numSubDomains * (numNodesX * numNodesY);
    const GO numNodesDual = numNodesX; // Interface is oriented along the global x-axis

    // Tests for block matrix
    {
      TEST_EQUALITY_CONST(blockedMatrix->Rows(), numSubDomains);
      TEST_EQUALITY_CONST(blockedMatrix->Cols(), numSubDomains);

      TEST_EQUALITY_CONST(blockedMatrix->getMatrix(0, 0)->getGlobalNumRows(), numNodesPrimal);
      TEST_EQUALITY_CONST(blockedMatrix->getMatrix(0, 1)->getGlobalNumRows(), numNodesPrimal);
      TEST_EQUALITY_CONST(blockedMatrix->getMatrix(1, 0)->getGlobalNumRows(), numNodesDual);

      // Saddle-point system needs to have a zero block
      TEST_ASSERT(blockedMatrix->getMatrix(1, 1).is_null());

      TEST_EQUALITY_CONST(blockedMatrix->getGlobalNumRows(), numNodesPrimal + numNodesDual);

      // Kinematic coupling: rows have entries [1, -1] or [-1, 1]
      RCP<const Matrix> A10 = blockedMatrix->getMatrix(1, 0);
      TEST_ASSERT(!A10.is_null());
      for (LocalOrdinal lDualRow = 0; lDualRow < static_cast<LocalOrdinal>(A10->getNodeNumRows()); ++lDualRow)
      {
        ArrayView<const LocalOrdinal> cols;
        ArrayView<const Scalar> vals;
        A10->getLocalRowView(lDualRow, cols, vals);

        TEST_EQUALITY_CONST(cols.size(), 2);
        TEST_EQUALITY_CONST(vals.size(), 2);
        TEST_INEQUALITY_CONST(vals[0], SC_zero);
        TEST_INEQUALITY_CONST(vals[1], SC_zero);
        TEST_EQUALITY_CONST(vals[0] * vals[0], SC_one);
        TEST_EQUALITY_CONST(vals[1] * vals[1], SC_one);
        TEST_EQUALITY_CONST(vals[0] + vals[1], SC_zero);
      }

      // Coupling of fluxes: rows have entries [1] or [-1]
      RCP<const Matrix> A01 = blockedMatrix->getMatrix(0, 1);
      TEST_ASSERT(!A01.is_null());
      for (LocalOrdinal lPrimalRow = 0; lPrimalRow < static_cast<LocalOrdinal>(A01->getNodeNumRows()); ++lPrimalRow)
      {
        ArrayView<const LocalOrdinal> cols;
        ArrayView<const Scalar> vals;
        A01->getLocalRowView(lPrimalRow, cols, vals);

        if (cols.size() > 0)
        {
          TEST_EQUALITY_CONST(cols.size(), 1);
          TEST_EQUALITY_CONST(vals.size(), 1);
          TEST_EQUALITY_CONST(vals[0] * vals[0], SC_one);
        }
      }
    }

    // Tests for interface dof map
    {
      TEST_EQUALITY_CONST(interfaceDofMap->getGlobalNumElements(), numNodesX);

      const auto comm = interfaceDofMap->getComm();
      if (comm->getSize() == 1)
      {
        Teuchos::Array<GlobalOrdinal> interfaceGIDs;
        interfaceGIDs.push_back(6);
        interfaceGIDs.push_back(7);
        interfaceGIDs.push_back(8);
        TEST_COMPARE_ARRAYS(interfaceDofMap->getNodeElementList(), interfaceGIDs);
      }
      else if (comm->getSize() == 4)
      {
        if (comm->getRank() == 0 || comm->getRank() == 1)
        {
          TEST_EQUALITY_CONST(interfaceDofMap->getNodeNumElements(), 0);
        }
        else if (comm->getRank() == 2)
        {
          TEST_EQUALITY_CONST(interfaceDofMap->getNodeNumElements(), 2);

          Teuchos::Array<GlobalOrdinal> myInterfaceGIDs;
          myInterfaceGIDs.push_back(6);
          myInterfaceGIDs.push_back(7);
          TEST_COMPARE_ARRAYS(interfaceDofMap->getNodeElementList(), myInterfaceGIDs);
        }
        else if (comm->getRank() == 3)
        {
          TEST_EQUALITY_CONST(interfaceDofMap->getNodeNumElements(), 1);

          Teuchos::Array<GlobalOrdinal> myInterfaceGIDs;
          myInterfaceGIDs.push_back(8);
          TEST_COMPARE_ARRAYS(interfaceDofMap->getNodeElementList(), myInterfaceGIDs);
        }
      }
    }
  }

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TestFactory, CreatePoisson1DMatrix, SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TestFactory, CreatePoisson2DMeshtyingMatrix_Obtain2x2SaddlePointSystem, SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TestFactory, CreatePoisson2DMeshtyingMatrix_Obtain2x2SaddlePointSystemAndInterfaceDofMap, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}