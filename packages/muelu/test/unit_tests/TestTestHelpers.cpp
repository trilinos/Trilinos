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

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

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

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TestFactory, CreatePoisson1DMatrix, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}


