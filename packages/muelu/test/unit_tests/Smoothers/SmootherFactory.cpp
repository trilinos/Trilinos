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
#include <MueLu_Version.hpp>
#include <MueLu_TestHelpers.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_FakeSmootherPrototype.hpp>

namespace MueLuTests {

namespace SmootherFactory {

using Teuchos::RCP;
using Teuchos::rcp;
using namespace MueLu::Exceptions;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class SmootherFactoryHelpers {
#include <MueLu_UseShortNames.hpp>

 public:
  ////////////////////////////////////////////////////////////////////////////////
  // Helper functions: class invariant
  ////////////////////////////////////////////////////////////////////////////////

  // Check the consistency of the internal state of a SmootherFactory object
  static void check(RCP<const SmootherFactory> smootherFactory, RCP<const SmootherPrototype> preSmooProto, RCP<const SmootherPrototype> postSmooProto, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> getPreSmooProto, getPostSmooProto;
    smootherFactory->GetSmootherPrototypes(getPreSmooProto, getPostSmooProto);

    TEST_EQUALITY(preSmooProto, getPreSmooProto);
    TEST_EQUALITY(postSmooProto, getPostSmooProto);

    // Check if internal prototypes are still prototypes (IsSetup() == false)
    if (preSmooProto != Teuchos::null) TEST_EQUALITY(preSmooProto->IsSetup(), false);
    if (postSmooProto != Teuchos::null) TEST_EQUALITY(postSmooProto->IsSetup(), false);
  }

  // Similar to check(), but with INEQUALITY.
  // Do not check if internal prototypes verifies IsSetup() == false
  static void ineqcheck(RCP<const SmootherFactory> smootherFactory, RCP<const SmootherBase> preSmooProto, RCP<const SmootherBase> postSmooProto, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> getPreSmooProto, getPostSmooProto;
    smootherFactory->GetSmootherPrototypes(getPreSmooProto, getPostSmooProto);
    TEST_INEQUALITY(preSmooProto, getPreSmooProto);
    TEST_INEQUALITY(postSmooProto, getPostSmooProto);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Helper functions: apply a test to a collection of test cases
  ////////////////////////////////////////////////////////////////////////////////

  // Apply a test to a list of cases (one argument)
  static void testCollection(void (*func)(RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> smooProto = rcp(new FakeSmootherPrototype());

    // tests
    func(smooProto, out, success);
    func(Teuchos::null, out, success);
  }

  // Apply a test to a list of invalid cases (one argument)
  static void testInvalidCollection(void (*func)(RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> invalidSmooProto = rcp(new FakeSmootherPrototype());
    Level level;
    invalidSmooProto->Setup(level);

    // tests
    func(invalidSmooProto, out, success);
  }

  // Apply a test to a list of cases (two argument)
  static void testCollection(
      void (*func)(RCP<SmootherPrototype>, RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &),
      Teuchos::FancyOStream &out,
      bool &success) {
    RCP<SmootherPrototype> smooProtoA = rcp(new FakeSmootherPrototype());
    RCP<SmootherPrototype> smooProtoB = rcp(new FakeSmootherPrototype());

    // tests
    func(Teuchos::null, Teuchos::null, out, success);
    func(smooProtoA, smooProtoB, out, success);
    func(smooProtoA, Teuchos::null, out, success);
    func(Teuchos::null, smooProtoB, out, success);
    func(smooProtoA, smooProtoA, out, success);
  }

  // Apply a test to a list of invalid cases (two argument)
  static void testInvalidCollection(void (*func)(RCP<SmootherPrototype>, RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> smooProtoA       = rcp(new FakeSmootherPrototype());
    RCP<SmootherPrototype> smooProtoB       = rcp(new FakeSmootherPrototype());
    RCP<SmootherPrototype> invalidSmooProto = rcp(new FakeSmootherPrototype());
    Level level;
    invalidSmooProto->Setup(level);

    // tests
    func(invalidSmooProto, smooProtoB, out, success);
    func(smooProtoA, invalidSmooProto, out, success);
    func(invalidSmooProto, Teuchos::null, out, success);
    func(Teuchos::null, invalidSmooProto, out, success);
    func(invalidSmooProto, invalidSmooProto, out, success);
  }

  // Build an object and test internal state (one argument)
  static void testConstructor1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream &out, bool &success) {
    check(rcp(new SmootherFactory(smooProto)), smooProto, smooProto, out, success);  // One argument == same pre and post smoother prototype
  }

  // Test with invalid input argument
  static void testInvalidConstructor1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream &out, bool &success) {
    TEST_THROW(rcp(new SmootherFactory(smooProto)), RuntimeError);
  }

  // Build an object and test internal state (two arguments)
  static void testConstructor2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    check(rcp(new SmootherFactory(smooProtoA, smooProtoB)), smooProtoA, smooProtoB, out, success);
  }

  // Test with invalid input arguments
  static void testInvalidConstructor2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    TEST_THROW(rcp(new SmootherFactory(smooProtoA, smooProtoB)), RuntimeError);
  }

  static void testSetSmootherPrototypes1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> smooProto0 = rcp(new FakeSmootherPrototype());
    RCP<SmootherFactory> smooFact     = rcp(new SmootherFactory(smooProto0));

    ineqcheck(smooFact, smooProto, smooProto, out, success);
    smooFact->SetSmootherPrototypes(smooProto);
    check(smooFact, smooProto, smooProto, out, success);  // One argument == same pre and post smoother prototype
  }

  static void testInvalidSetSmootherPrototypes1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream &out, bool &success) {
    SmootherFactory smooFact(Teuchos::null);
    TEST_THROW(smooFact.SetSmootherPrototypes(smooProto), RuntimeError);
  }

  static void testSetSmootherPrototypes2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherPrototype> smooProto0 = rcp(new FakeSmootherPrototype());
    RCP<SmootherFactory> smooFact     = rcp(new SmootherFactory(smooProto0));

    ineqcheck(smooFact, smooProtoA, smooProtoB, out, success);
    smooFact->SetSmootherPrototypes(smooProtoA, smooProtoB);
    check(smooFact, smooProtoA, smooProtoB, out, success);
  }

  static void testInvalidSetSmootherPrototypes2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    SmootherFactory smooFact(Teuchos::null);
    TEST_THROW(smooFact.SetSmootherPrototypes(smooProtoA, smooProtoB), RuntimeError);
  }

  // test if a smoother created by Build() is correct (check if it corresponds to the prototype)
  static void testBuildCheckOutput(RCP<SmootherFactory> smooFact, Level &level, RCP<SmootherPrototype> smooProto, const std::string &tag, Teuchos::FancyOStream &out, bool &success) {
    if (smooProto == Teuchos::null) {
      TEST_EQUALITY(level.IsAvailable(tag, smooFact.get()), false);
    } else {
      RCP<SmootherBase> smoother;
      TEST_NOTHROW(smoother = level.Get<RCP<SmootherBase> >(tag, smooFact.get()));

      TEST_INEQUALITY(smoother, Teuchos::null);
      TEST_INEQUALITY(smoother, smooProto);

      if (smooProto != Teuchos::null) {
        RCP<FakeSmootherPrototype> smootherF;

        // ouput test: smoothers same derived class as prototypes
        TEST_NOTHROW(smootherF = rcp_dynamic_cast<FakeSmootherPrototype>(smoother, true));

        if (smootherF != Teuchos::null) {
          // output test: smoother parameters == prototype parameters
          RCP<FakeSmootherPrototype> smooProtoF = rcp_dynamic_cast<FakeSmootherPrototype>(smooProto, true);
          TEST_EQUALITY(smootherF->GetParam(), smooProtoF->GetParam());

          // output test: smoothers are ready to be apply
          TEST_EQUALITY(smootherF->IsSetup(), true);

          // setup done only once.
          TEST_EQUALITY(smootherF->GetNumOfSetupCall(), 1);
          TEST_EQUALITY(smootherF->GetNumOfSetup(), 1);
        }
      }
    }
  }

  // test outpout of Build() and BuildSmoother()
  static void testBuildCheck(RCP<SmootherFactory> smooFact, Level &level, RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, MueLu::PreOrPost preOrPost, Teuchos::FancyOStream &out, bool &success) {
    // invariant: smoother prototypes kept unchanged
    check(smooFact, smooProtoA, smooProtoB, out, success);

    // output test
    if (preOrPost == MueLu::BOTH) {
      testBuildCheckOutput(smooFact, level, smooProtoA, "PreSmoother", out, success);
      testBuildCheckOutput(smooFact, level, smooProtoB, "PostSmoother", out, success);

      // ReUse: if pre and post prototype are the same, then pre smoother == post smoother
      // otherwise, they are different (have not been tested by previous tests)
      RCP<SmootherBase> smooA, smooB;
      if (smooProtoA != Teuchos::null) {
        smooA = level.Get<RCP<SmootherBase> >("PreSmoother", smooFact.get());
      }
      if (smooProtoB != Teuchos::null) {
        smooB = level.Get<RCP<SmootherBase> >("PostSmoother", smooFact.get());
      }
      if (smooProtoA == smooProtoB) {
        TEST_EQUALITY(smooA, smooB);
      } else {
        TEST_INEQUALITY(smooA, smooB);
      }

    } else if (preOrPost == MueLu::PRE) {
      testBuildCheckOutput(smooFact, level, smooProtoA, "PreSmoother", out, success);
      TEST_EQUALITY(level.IsAvailable("PostSmoother", smooFact.get()), false);

    } else if (preOrPost == MueLu::POST) {
      TEST_EQUALITY(level.IsAvailable("PreSmoother", smooFact.get()), false);
      testBuildCheckOutput(smooFact, level, smooProtoB, "PostSmoother", out, success);

    } else {
      TEST_EQUALITY(true, false);
    }
  }

  //
  static void testBuild(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smooProtoA, smooProtoB));

    Level level;  // level.SetupPhase(true);
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);

    level.Request("PreSmoother", smooFact.get());
    level.Request("PostSmoother", smooFact.get());

    smooFact->Build(level);

    testBuildCheck(smooFact, level, smooProtoA, smooProtoB, MueLu::BOTH, out, success);
  }

  static void testBuildSmootherPreOrPost(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, MueLu::PreOrPost preOrPost, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smooProtoA, smooProtoB));

    Level level;  // level.SetupPhase(true);
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);

    level.Request("PreSmoother", smooFact.get());
    level.Request("PostSmoother", smooFact.get());

    smooFact->BuildSmoother(level, preOrPost);

    testBuildCheck(smooFact, level, smooProtoA, smooProtoB, preOrPost, out, success);
  }

  static void testBuildSmoother(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    testBuildSmootherPreOrPost(smooProtoA, smooProtoB, MueLu::PRE, out, success);
    testBuildSmootherPreOrPost(smooProtoA, smooProtoB, MueLu::POST, out, success);
    testBuildSmootherPreOrPost(smooProtoA, smooProtoB, MueLu::BOTH, out, success);
  }

  static void testBuildSmootherDefaultArg(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream &out, bool &success) {
    RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smooProtoA, smooProtoB));

    Level level;  // level.SetupPhase(true);
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);

    level.Request("PreSmoother", smooFact.get());
    level.Request("PostSmoother", smooFact.get());

    smooFact->BuildSmoother(level);

    testBuildCheck(smooFact, level, smooProtoA, smooProtoB, MueLu::BOTH, out, success);
  }

};  // class SmootherFactoryHelpers

////////////////////////////////////////////////////////////////////////////////
// Test: Constructor_OneArg
////////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, Constructor_OneArg, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SFH::testCollection(&SFH::testConstructor1, out, success);                // TEST: Valid input parameter
  SFH::testInvalidCollection(&SFH::testInvalidConstructor1, out, success);  // TEST: Invalid input parameter
}

////////////////////////////////////////////////////////////////////////////////
// Test: Constructor_TwoArg
////////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, Constructor_TwoArg, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SFH::testCollection(&SFH::testConstructor2, out, success);                // TEST: Valid input parameter
  SFH::testInvalidCollection(&SFH::testInvalidConstructor2, out, success);  // TEST: Valid input parameter
}

////////////////////////////////////////////////////////////////////////////////
// Test: SetSmootherPrototypes_OneArg
////////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, SetSmootherPrototypes_OneArg, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SFH::testCollection(&SFH::testSetSmootherPrototypes1, out, success);                // TEST: Valid input parameter
  SFH::testInvalidCollection(&SFH::testInvalidSetSmootherPrototypes1, out, success);  // TEST: Invalid input parameter
}

////////////////////////////////////////////////////////////////////////////////
// Test: SetSmootherPrototypes_TwoArg
////////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, SetSmootherPrototypes_TwoArg, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SFH::testCollection(&SFH::testSetSmootherPrototypes2, out, success);                // TEST: Valid input parameter
  SFH::testInvalidCollection(&SFH::testInvalidSetSmootherPrototypes2, out, success);  // TEST: Invalid input parameter
}

////////////////////////////////////////////////////////////////////////////////
// Test: GetSmootherPrototypes
////////////////////////////////////////////////////////////////////////////////

/*
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, GetSmootherPrototypes, Scalar, LocalOrdinal, GlobalOrdinal, Node)
{
  // Already heavily tested by previous tests.
}
*/

////////////////////////////////////////////////////////////////////////////////
// Test: DeclareInput
////////////////////////////////////////////////////////////////////////////////

/*
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, DeclareInput, Scalar, LocalOrdinal, GlobalOrdinal, Node)
{
  // TODO:
  // Create a level, call DeclareInput() and Build() and then check that every counter == 0 at the end
  // Same test have to be done for each factory...
}
*/

////////////////////////////////////////////////////////////////////////////////
// Test: Build
////////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SFH::testCollection(&SFH::testBuild, out, success);
}

////////////////////////////////////////////////////////////////////////////////
// Test: BuildSmoother
////////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, BuildSmoother, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SFH::testCollection(&SFH::testBuildSmootherDefaultArg, out, success);  // TEST: default arg
  SFH::testCollection(&SFH::testBuildSmoother, out, success);            // TEST: PRE, POST and BOTH
}

////////////////////////////////////////////////////////////////////////////////
// Test: Compilation
////////////////////////////////////////////////////////////////////////////////

// Make sure that Constructors and SetSmootherPrototypes() are not using references as input arguments
// (to allow Teuchos::null as input arguments)
// TODO: const ref?
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SmootherFactory, Compilation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef SmootherFactoryHelpers<SC, LO, GO, NO> SFH;
  SmootherFactory smooFact1(Teuchos::null);
  SmootherFactory smooFact2(Teuchos::null, Teuchos::null);
  smooFact1.SetSmootherPrototypes(Teuchos::null);
  smooFact1.SetSmootherPrototypes(Teuchos::null, Teuchos::null);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, Constructor_OneArg, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, Constructor_TwoArg, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, SetSmootherPrototypes_OneArg, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, SetSmootherPrototypes_TwoArg, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, Build, SC, LO, GO, NO)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, BuildSmoother, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SmootherFactory, Compilation, SC, LO, GO, NO)
#include <MueLu_ETI_4arg.hpp>

}  // namespace SmootherFactory

}  // namespace MueLuTests

// TODO: advance reuse test
// TODO: remove some useless RCP in helper function prototypes?
