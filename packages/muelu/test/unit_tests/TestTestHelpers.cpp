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
    const GO expectedGlobalNumEntries = numEntriesPerRow * numNodes - 2;
    TEST_EQUALITY_CONST(matrix->getGlobalNumEntries(), expectedGlobalNumEntries);

    // Check every single matrix entry
    for (LocalOrdinal lRow = 0; lRow < matrix->getNodeNumRows(); ++lRow)
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
    RCP<const Matrix> matrix = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildPoissonSaddlePointMatrix(numNodesX, numNodesY);
    TEST_ASSERT(!matrix.is_null());

    const size_t numSubDomains = 2;
    const GO numNodesPrimal = numSubDomains * (numNodesX * numNodesY);
    const GO numNodesDual = numNodesX; // Interface is oriented along the global x-axis

    RCP<const BlockedCrsMatrix> blockedMatrix = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(matrix);
    TEST_ASSERT(!blockedMatrix.is_null())
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
    for (LocalOrdinal lDualRow = 0; lDualRow < A10->getNodeNumRows(); ++lDualRow)
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
    for (LocalOrdinal lPrimalRow = 0; lPrimalRow < A01->getNodeNumRows(); ++lPrimalRow)
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

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TestFactory, CreatePoisson1DMatrix, SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TestFactory, CreatePoisson2DMeshtyingMatrix_Obtain2x2SaddlePointSystem, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}
