// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_VariableContainer.hpp"

namespace MueLuTests {

using MueLu::KeepType;
using MueLu::VariableContainer;

TEUCHOS_UNIT_TEST(VariableContainer, Constructor) {
  VariableContainer vc;

  KeepType keepFlag = vc.GetKeepFlag();
  TEST_EQUALITY(keepFlag, 0);
}

TEUCHOS_UNIT_TEST(VariableContainer, AddKeepFlag_RemoveKeepFlag) {
  VariableContainer vc;

  // Test: default input parameter is UserData
  vc.AddKeepFlag();
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData);
  vc.RemoveKeepFlag();
  TEST_EQUALITY(vc.GetKeepFlag(), 0);

  // Test: add and remove one flag
  vc.AddKeepFlag(MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::Keep);
  vc.RemoveKeepFlag(MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), 0);

  // Test: add and remove are using bitwise operations (| ^ instead of + -)
  // Test of add
  vc.AddKeepFlag(MueLu::Keep);  // 2 add
  vc.AddKeepFlag(MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::Keep);  // failed if AddKeepFlag uses + instead of |
  vc.RemoveKeepFlag(MueLu::Keep);                // 1 remove
  TEST_EQUALITY(vc.GetKeepFlag(), 0);
  // Test of remove (removing something twice or something that do not exist)
  vc.AddKeepFlag(MueLu::Keep);  // 1 add
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::Keep);
  vc.RemoveKeepFlag(MueLu::Keep);  // 2 remove
  vc.RemoveKeepFlag(MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), 0);

  // Test: internal flag can be a combinaison of flags.
  vc.AddKeepFlag(MueLu::UserData);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData);
  vc.AddKeepFlag(MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData | MueLu::Keep);
  vc.AddKeepFlag(MueLu::Final);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData | MueLu::Keep | MueLu::Final);
  vc.RemoveKeepFlag(MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData | MueLu::Final);
  vc.RemoveKeepFlag(MueLu::UserData);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::Final);
  vc.RemoveKeepFlag(MueLu::Final);
  TEST_EQUALITY(vc.GetKeepFlag(), 0);

  // Test: input of add and remove can be a combinaison of flags
  vc.AddKeepFlag(MueLu::UserData | MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData | MueLu::Keep);
  vc.RemoveKeepFlag(MueLu::UserData | MueLu::Keep);
  TEST_EQUALITY(vc.GetKeepFlag(), 0);

  // Test: MueLu::All == MueLu::UserData | MueLu::Keep | MueLu::Final
  TEST_EQUALITY(MueLu::UserData | MueLu::Keep | MueLu::Final, MueLu::All);

  // Test: a very common use case: set one flag, remove all flags
  // (this is what is actually done by Level::Keep() followed by Level::Delete())
  vc.AddKeepFlag(MueLu::UserData);
  TEST_EQUALITY(vc.GetKeepFlag(), MueLu::UserData);
  vc.RemoveKeepFlag(MueLu::All);
  TEST_EQUALITY(vc.GetKeepFlag(), 0);

  // TODO: Test: add remove check if input argument is a valid flag
  //     TEST_THROW(vc.AddKeepFlag(42))
  //     TEST_THROW(vc.DeleteKeepFlag(42))
}

}  // namespace MueLuTests
