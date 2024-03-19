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
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TestingHelpers.hpp>

//#include "MueLu_Hierarchy.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include <unistd.h>

/**********************************************************************************/

namespace MueLu {
/*!
  @class MiniFactory classes.
  @brief Mini factories for demonstration purposes only
*/
class MiniFactory : public SingleLevelFactoryBase {
 public:
  MiniFactory(std::string factoryName, std::string dependsOn) {
    factoryName_ = factoryName;
    dependsOn_   = dependsOn;
  }
  virtual ~MiniFactory() {}
  RCP<const ParameterList> GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set<std::string>("internal-data", "", "internal data for factory");
    validParamList->set<bool>("CallDependencyBuild", true, "If true, the dependencies are built");
    validParamList->set<RCP<const FactoryBase> >(dependsOn_, Teuchos::null, "Generating factory of data this factory depends on");
    return validParamList;
  };
  void DeclareInput(Level &Level) const {
    if (dependsOn_ != "")
      Input(Level, dependsOn_);
  };
  void Build(Level &Level) const {
    std::cout << "MiniFactory" << factoryName_ << "->Build" << std::endl;
    const ParameterList &pL = GetParameterList();
    std::string data        = pL.get<std::string>("internal-data");
    std::cout << "Data of MiniFactory" << factoryName_ << " = " << data << std::endl;

    // call build for dependencies:
    bool bBuild = pL.get<bool>("CallDependencyBuild");
    if (bBuild) {
      std::string dep_data = Get<std::string>(Level, dependsOn_);
      std::cout << "Content of variable " << dependsOn_ << " = " << dep_data << std::endl;
    }

    Set(Level, factoryName_, std::string("Data produced by MiniFactory" + factoryName_));
  };

 private:
  std::string factoryName_;
  std::string dependsOn_;
};  // class MiniFactory

/*!
  @class MiniSwitchFactory classes.
  @brief Mini factories for demonstration purposes only
*/
class MiniSwitchFactory : public SingleLevelFactoryBase {
 public:
  MiniSwitchFactory(std::string factoryName, std::string dependsOnA, std::string dependsOnB) {
    factoryName_ = factoryName;
    dependsOnA_  = dependsOnA;
    dependsOnB_  = dependsOnB;
  }
  virtual ~MiniSwitchFactory() {}
  RCP<const ParameterList> GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set<std::string>("internal-data", "", "internal data for factory");
    validParamList->set<bool>("UseA", true, "If true, the dependencies for A are built");
    validParamList->set<bool>("UseB", true, "If true, the dependencies for B are built");
    validParamList->set<RCP<const FactoryBase> >(dependsOnA_, Teuchos::null, "Generating factory of data this factory depends on");
    validParamList->set<RCP<const FactoryBase> >(dependsOnB_, Teuchos::null, "Generating factory of data this factory depends on");
    return validParamList;
  };
  void DeclareInput(Level &Level) const {
    if (dependsOnA_ != "")
      Input(Level, dependsOnA_);
    if (dependsOnB_ != "")
      Input(Level, dependsOnB_);
  };
  void Build(Level &Level) const {
    std::cout << "MiniFactory" << factoryName_ << "->Build" << std::endl;
    const ParameterList &pL = GetParameterList();
    std::string data        = pL.get<std::string>("internal-data");
    std::cout << "Data of MiniFactory" << factoryName_ << " = " << data << std::endl;

    // call build for dependencies:
    bool bBuild = pL.get<bool>("UseA");
    if (bBuild) {
      std::string dep_data = Get<std::string>(Level, dependsOnA_);
      std::cout << "Content of variable " << dependsOnA_ << " = " << dep_data << std::endl;
    }
    bBuild = pL.get<bool>("UseB");
    if (bBuild) {
      std::string dep_data = Get<std::string>(Level, dependsOnB_);
      std::cout << "Content of variable " << dependsOnB_ << " = " << dep_data << std::endl;
    }

    Set(Level, factoryName_, std::string("Data produced by MiniFactory" + factoryName_));
  };

 private:
  std::string factoryName_;
  std::string dependsOnA_;
  std::string dependsOnB_;
};  // class MiniSwitchFactory

}  // namespace MueLu

/**********************************************************************************/

int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success  = true;
  bool bSuccess = true;  // collect success results
  bool verbose  = true;
  try {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<MueLu::Level> finest = Teuchos::rcp(new MueLu::Level());

    // B depends on A
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      finest->Request("B", facB.get());
      facB->Build(*finest);
      // finest->print(std::cout,MueLu::Debug);
      finest->Release(*facB);
      // finest->print(std::cout,MueLu::Debug);
    }

    // B depends on A but does not use A
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));
      facB->SetFactory("A", facA);

      finest->Request("B", facB.get());
      facB->Build(*finest);
      // finest->print(std::cout,MueLu::Debug);
      finest->Release(*facB);
      // finest->print(std::cout,MueLu::Debug);
    }

    // C depends on B which depends on A
    //              B uses A
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));
      facC->SetFactory("B", facB);

      finest->Request("C", facC.get());
      facC->Build(*finest);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facC);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
    }

    // C depends on B which depends on A
    // C uses       B which uses A
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facC->SetFactory("B", facB);

      finest->Request("C", facC.get());
      facC->Build(*finest);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facC);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      // finest->print(std::cout,MueLu::Debug);
    }

    // C depends on B which depends on A
    // C uses       B
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facC->SetFactory("B", facB);

      finest->Request("C", facC.get());
      facC->Build(*finest);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facC);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
    }

    // D depends on A and C; C depends on B which depends on A
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facC->SetFactory("B", facB);

      Teuchos::RCP<MueLu::MiniSwitchFactory> facD = Teuchos::rcp(new MueLu::MiniSwitchFactory("D", "A", "C"));
      facD->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for D")));
      facD->SetParameter("UseA", Teuchos::ParameterEntry(true));
      facD->SetParameter("UseB", Teuchos::ParameterEntry(true));  // B means only second factory
      facD->SetFactory("A", facA);
      facD->SetFactory("C", facC);

      finest->Request("D", facD.get());
      facD->Build(*finest);
      // finest->print(std::cout,MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facD);
      // finest->print(std::cout,MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
    }

    // D depends on B and C; C depends on B which depends on A
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facC->SetFactory("B", facB);

      Teuchos::RCP<MueLu::MiniSwitchFactory> facD = Teuchos::rcp(new MueLu::MiniSwitchFactory("D", "B", "C"));
      facD->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for D")));
      facD->SetParameter("UseA", Teuchos::ParameterEntry(true));
      facD->SetParameter("UseB", Teuchos::ParameterEntry(true));  // B means only second factory
      facD->SetFactory("B", facB);
      facD->SetFactory("C", facC);

      finest->Request("D", facD.get());
      facD->Build(*finest);
      // finest->print(std::cout,MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facD);
      // finest->print(std::cout,MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
    }

    // D depends on B and C; C depends on B which depends on A
    // D only uses B (not C)
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facC->SetFactory("B", facB);

      Teuchos::RCP<MueLu::MiniSwitchFactory> facD = Teuchos::rcp(new MueLu::MiniSwitchFactory("D", "B", "C"));
      facD->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for D")));
      facD->SetParameter("UseA", Teuchos::ParameterEntry(true));
      facD->SetParameter("UseB", Teuchos::ParameterEntry(false));  // B means only second factory
      facD->SetFactory("B", facB);
      facD->SetFactory("C", facC);

      finest->Request("D", facD.get());
      facD->Build(*finest);
      // finest->print(std::cout,MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facD);
      // finest->print(std::cout,MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
    }

    // D depends on B and C; C depends on B which depends on A
    // D only uses C (not B)
    {
      finest = Teuchos::null;
      finest = Teuchos::rcp(new MueLu::Level());

      Teuchos::RCP<MueLu::MiniFactory> facA = Teuchos::rcp(new MueLu::MiniFactory("A", ""));
      facA->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for A")));
      facA->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(false));

      Teuchos::RCP<MueLu::MiniFactory> facB = Teuchos::rcp(new MueLu::MiniFactory("B", "A"));
      facB->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for B")));
      facB->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facB->SetFactory("A", facA);

      Teuchos::RCP<MueLu::MiniFactory> facC = Teuchos::rcp(new MueLu::MiniFactory("C", "B"));
      facC->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for C")));
      facC->SetParameter("CallDependencyBuild", Teuchos::ParameterEntry(true));
      facC->SetFactory("B", facB);

      Teuchos::RCP<MueLu::MiniSwitchFactory> facD = Teuchos::rcp(new MueLu::MiniSwitchFactory("D", "B", "C"));
      facD->SetParameter("internal-data", Teuchos::ParameterEntry(std::string("data for D")));
      facD->SetParameter("UseA", Teuchos::ParameterEntry(false));
      facD->SetParameter("UseB", Teuchos::ParameterEntry(true));  // B means only second factory
      facD->SetFactory("B", facB);
      facD->SetFactory("C", facC);

      finest->Request("D", facD.get());
      facD->Build(*finest);
      finest->print(std::cout, MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("C", facC.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      finest->Release(*facD);
      finest->print(std::cout, MueLu::Debug);
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("C", facC.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("B", facB.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsAvailable("A", facA.get()), false, std::cout, success);
      if (!success) bSuccess = false;
      TEUCHOS_TEST_EQUALITY(finest->IsRequested("D", facD.get()), true, std::cout, success);
      if (!success) bSuccess = false;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, bSuccess);

  return (bSuccess ? EXIT_SUCCESS : EXIT_FAILURE);
}
