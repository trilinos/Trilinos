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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "Galeri_XpetraUtils.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_DistanceLaplacianFactory.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(DistanceLaplacian, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<DistanceLaplacianFactory> dlFact = rcp(new DistanceLaplacianFactory);
    TEST_EQUALITY(dlFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(DistanceLaplacian, Build)
  {

    out << "version: " << MueLu::Version() << std::endl;
    out << "Test DistanceLaplacianFactory with 2D Poisson" << std::endl;

    Level fineLevel, coarseLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);

    int nx=5; int ny=5;
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build2DPoisson(nx,ny);
    Teuchos::ParameterList list;
    list.set("nx",nx); list.set("ny",ny);
    RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D",A->getRowMap(),list);
    A->SetFixedBlockSize(1);
    fineLevel.Set("A",A);
    fineLevel.Set("Coordinates",coordinates);
    RCP<DistanceLaplacianFactory> dlFact = rcp(new DistanceLaplacianFactory);
    fineLevel.Request("Distance Laplacian",dlFact.get());
    dlFact->Build(fineLevel);
    RCP<Matrix> DL = fineLevel.Get<RCP<Matrix> >("Distance Laplacian", dlFact.get());

    TEST_EQUALITY(A->getGlobalNumEntries()==DL->getGlobalNumEntries(), true);
    TEST_EQUALITY(A->getGlobalNumRows()==DL->getGlobalNumRows(), true);
    TEST_EQUALITY(A->getGlobalNumCols()==DL->getGlobalNumCols(), true);

    RCP<MultiVector> onesVector = MultiVectorFactory::Build(A->getDomainMap(),1,false);
    RCP<MultiVector> result = MultiVectorFactory::Build(A->getRangeMap(),1,false);
    onesVector->putScalar(Teuchos::ScalarTraits<SC>::one());
    DL->apply(*onesVector,*result,Teuchos::NO_TRANS,1.0,0.0);
    Teuchos::Array<ST::magnitudeType> norms(1);
    result->norm2(norms);
    TEST_EQUALITY(norms[0]<1e-12, true);

  } //MakeTentative  Lapack QR

} // namespace MueLuTests
