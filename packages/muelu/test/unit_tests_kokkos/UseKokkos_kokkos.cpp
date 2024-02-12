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

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_Hierarchy.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UseKokkos_kokkos, FactoryManager, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using real_type             = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  const bool useKokkos = false;

  // Build the problem
  RCP<Matrix> A                          = MueLuTests::TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2001);
  RCP<const Map> map                     = A->getMap();
  RCP<RealValuedMultiVector> coordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(map, 1);
  RCP<MultiVector> nullspace             = MultiVectorFactory::Build(map, 1);
  nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());

  // Setting paramters for the hierarchy
  Teuchos::ParameterList paramList;
  paramList.set("verbosity", "test");
  paramList.set("use kokkos refactor", useKokkos);
  RCP<HierarchyManager> mueluFactory = rcp(new ParameterListInterpreter(paramList));
  mueluFactory->CheckConfig();
  const RCP<const FactoryManagerBase> manager = mueluFactory->GetFactoryManager(0);
  RCP<const FactoryManager> factoryManager    = Teuchos::rcp_dynamic_cast<const FactoryManager>(manager, true);

  TEST_EQUALITY(factoryManager->GetKokkosRefactor() == useKokkos, true);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UseKokkos_kokkos, FactoryManager, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
