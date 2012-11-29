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

#include "MueLu_TentativePFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include "MueLu_TwoKeyMap.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(TwoKeyMap, ptrTests)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory> sapFactory = rcp(new TentativePFactory);
    TEST_EQUALITY(sapFactory != Teuchos::null, true);

    RCP<TentativePFactory> sapFactory2 = rcp(new TentativePFactory);
    TEST_EQUALITY(sapFactory2 != Teuchos::null, true);

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // build Matrix
    Teuchos::ParameterList params;
    const RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), 20, 0, comm);
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap> > Pr = Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap>("Laplace1D", map, params);
    RCP<Matrix> Op = Pr->BuildMatrix();

    // a TwoKeyMap
    //TODO    RCP<MueLu::UTILS::TwoKeyMap> exh = Teuchos::rcp(new MueLu::UTILS::TwoKeyMap());

    //TODO
//     exh->Set<RCP<Matrix> >("op", sapFactory.get(), Op);
//     RCP<Matrix> test = Teuchos::null;
//     test = exh->Get<RCP<Matrix> >("op", sapFactory.get());
//     TEST_EQUALITY_CONST( test, Op );

//     exh->Set("op", sapFactory.get(), 22);
//     int test2 = exh->Get<int> ("op", sapFactory.get());
//     TEST_EQUALITY_CONST( test2, 22 );

//     exh->Set("op", sapFactory2.get(), Op);
//     RCP<Matrix> test3 = exh->Get<RCP<Matrix> > ("op", sapFactory2.get());
//     TEST_EQUALITY_CONST( test3, Op );
//     TEST_EQUALITY_CONST( exh->GetType("op", sapFactory.get()), "int" );

//     exh->Remove("op", sapFactory2.get());
//     TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory.get()),  true );
//     TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory2.get()), false );

//     exh->Remove("op", sapFactory.get());
//     TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory.get()),  false );
//     TEST_EQUALITY_CONST( exh->IsKey("op", sapFactory2.get()), false );

  }

}//namespace MueLuTests

