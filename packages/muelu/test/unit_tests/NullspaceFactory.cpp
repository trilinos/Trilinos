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
#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_NullspaceFactory.hpp"
//#include "MueLu_FactoryManager.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(NullspaceFactory, 3DElasticity, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  GO nx = 2, ny = 3, nz = 4, mx = comm->getSize(), my = 1, mz = 1;
  ;
  LO blkSize = 3;
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", ny);
  matrixList.set("nz", nz);
  matrixList.set("mx", mx);
  matrixList.set("my", my);
  matrixList.set("mz", mz);
  matrixList.set("matrixType", "Elasticity3D");
  RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::BuildMatrix(matrixList, TestHelpers::Parameters::getLib());
  Op->SetFixedBlockSize(blkSize);

  fineLevel.Set("A", Op);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("nz", nz);
  RCP<Matrix> fake                       = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx * ny * nz);  // need scalar PDE to produce proper coordinate array
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", fake->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<NullspaceFactory> nspace = rcp(new NullspaceFactory());
  fineLevel.Request("Nullspace", &(*nspace));
  Teuchos::ParameterList nspaceList = *(nspace->GetValidParameterList());
  nspaceList.set("nullspace: calculate rotations", true);
  nspace->SetParameterList(nspaceList);
  fineLevel.Request(*nspace);
  nspace->Build(fineLevel);
  RCP<MultiVector> nullSpace = fineLevel.Get<RCP<MultiVector> >("Nullspace", nspace.get());

  ArrayRCP<real_type> xvals          = coordinates->getDataNonConst(0);
  ArrayRCP<real_type> yvals          = coordinates->getDataNonConst(1);
  ArrayRCP<real_type> zvals          = coordinates->getDataNonConst(2);
  size_t nullDim                     = 6;
  RCP<MultiVector> handCompNullSpace = MultiVectorFactory::Build(Op->getRowMap(), nullDim);
  handCompNullSpace->putScalar(0.0);
  {
    ArrayRCP<Scalar> nsValues;
    nsValues = handCompNullSpace->getDataNonConst(0);
    for (int j = 0; j < nsValues.size(); j += 3) nsValues[j] = 1.;
    nsValues = handCompNullSpace->getDataNonConst(1);
    for (int j = 1; j < nsValues.size(); j += 3) nsValues[j] = 1.;
    nsValues = handCompNullSpace->getDataNonConst(2);
    for (int j = 2; j < nsValues.size(); j += 3) nsValues[j] = 1.;
    nsValues = handCompNullSpace->getDataNonConst(3);
    for (int j = 0; j < nsValues.size(); j += 3) nsValues[j] = -(yvals[j / 3] - .5);
    for (int j = 1; j < nsValues.size(); j += 3) nsValues[j] = (xvals[j / 3] - .5);
    nsValues = handCompNullSpace->getDataNonConst(4);
    for (int j = 1; j < nsValues.size(); j += 3) nsValues[j] = -(zvals[j / 3] - .5);
    for (int j = 2; j < nsValues.size(); j += 3) nsValues[j] = (yvals[j / 3] - .5);
    nsValues = handCompNullSpace->getDataNonConst(5);
    for (int j = 0; j < nsValues.size(); j += 3) nsValues[j] = -(zvals[j / 3] - .5);
    for (int j = 2; j < nsValues.size(); j += 3) nsValues[j] = (xvals[j / 3] - .5);
  }

  RCP<MultiVector> diff = MultiVectorFactory::Build(Op->getRowMap(), 6);
  diff->putScalar(0.0);
  // diff = fineNS + (-1.0)*(handCompNullSpace) + 0*diff
  diff->update(1.0, *nullSpace, -1.0, *handCompNullSpace, 0.0);

  Teuchos::Array<typename TST::magnitudeType> norms(6);
  diff->norm2(norms);
  for (LocalOrdinal i = 0; i < 6; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_EQUALITY(norms[i] < 100 * TMT::eps(), true);
  }

}  // NullSpace test

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NullspaceFactory, 3DElasticity, Scalar, LO, GO, Node)
#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
