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
#include <Teuchos_TabularOutputter.hpp>

#include "MueLu_Needs.hpp"

#include "MueLu_NoFactory.hpp"

namespace MueLu {

  Needs::Needs() { }

  Needs::~Needs() { }

  void Needs::Request(const std::string & ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    if(!dataTable_.IsKey(factory, ename)) {
      Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
      dataTable_.Set(factory, ename, newVar);
    }

    Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    var->Request(requestedBy);
  } //Request

  void Needs::Release(const std::string & ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");

    Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    var->Release(requestedBy);

    if (var->IsRequested() == false && var->GetKeepFlag() == 0) {
      var = Teuchos::null; // free data
      dataTable_.Remove(factory, ename);
    }
  } //Release

  void Needs::AddKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep) {
    if (!dataTable_.IsKey(factory, ename)) {
      // If the entry does not exist, create it to store the keep flag
      Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
      dataTable_.Set(factory, ename, newVar);
    }

    // Set the flag
    Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    var->AddKeepFlag(keep);
  }

  void Needs::RemoveKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep) {
    // No entry = nothing to do
    if (!dataTable_.IsKey(factory,ename)) return;

    // Remove the flag
    Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    var->RemoveKeepFlag(keep);

    // Remove data if no keep flag left and counter == 0
    if ((var->IsRequested() == false) && (var->GetKeepFlag() == 0)) {
      var = Teuchos::null; // free data
      dataTable_.Remove(factory, ename);
    }

  }

  KeepType Needs::GetKeepFlag(const std::string & ename, const FactoryBase* factory) const {
    if(!dataTable_.IsKey(factory,ename)) return false;
    //TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    return var->GetKeepFlag();
  }

  bool Needs::IsKey(const std::string & ename, const FactoryBase* factory) const {
    return dataTable_.IsKey(factory, ename);
  }

  bool Needs::IsAvailable(const std::string & ename, const FactoryBase* factory) const {
    //TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsAvailable(): " + ename + " not found. Do a request first.");
    if(!dataTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    return var->IsAvailable();
  }

  bool Needs::IsRequested(const std::string & ename, const FactoryBase* factory) const {
    //TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): " + ename + " not found. Do a request first.");
    if(!dataTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->GetKeepFlag() == 0, Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->IsRequested();
  }

  bool Needs::IsRequestedBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const {
    //TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsRequestedBy(): " + ename + " not found. Do a request first.");
    if(!dataTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->GetKeepFlag() == 0, Exceptions::RuntimeError, "MueLu::Needs::IsRequestedBy(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->IsRequested(requestedBy);
  }

  bool Needs::IsRequestedFactory(const FactoryBase* factory) {
    if(dataTable_.IsKey(factory) == false) return false;

    std::vector<std::string> ekeys = RequestedKeys(factory);
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); ++it) {
      if(IsRequested(*it,factory)) return true;
    }
    return false;
  }

  int Needs::CountRequestedFactory(const FactoryBase* factory) {
    int cnt = 0;
    if(dataTable_.IsKey(factory) == false) return cnt;
    std::vector<std::string> ekeys = RequestedKeys(factory);
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); ++it) {
      cnt += NumRequests(*it, factory);
    }
    return cnt;
  }

  bool Needs::IsAvailableFactory(const FactoryBase* factory) {
    if(dataTable_.IsKey(factory) == false) return false;
    std::vector<std::string> ekeys = RequestedKeys(factory);
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); ++it) {
      if(IsAvailable(*it, factory))
        return true;
    }
    return false;
  }

  int Needs::NumRequests(const std::string & ename, const FactoryBase* factory) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): " + ename + " not found. Do a request first."); //TODO: fix message when factory==NoFactory() (no request needed)?
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->GetKeepFlag() == 0, Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->NumAllRequests();
  }

  int Needs::NumRequestsBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::NumRequestsBy(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->GetKeepFlag() == 0, Exceptions::RuntimeError, "MueLu::Needs::NumRequestsBy(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->NumRequests(requestedBy);
  }

  std::vector<std::string> Needs::RequestedKeys(const FactoryBase* factory) const {
    return dataTable_.GetKey2List(factory);
  }

  std::vector<const FactoryBase*> Needs::RequestedFactories() const {
    return dataTable_.GetKeyList();
  }

  std::string Needs::GetType(const std::string & ename, const FactoryBase* factory) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory, ename), Exceptions::RuntimeError, "MueLu::Needs::Get(): " + ename + " not found in dataTable_");
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    return var->GetData().getAny(true).typeName();
  }

  void Needs::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    Teuchos::TabularOutputter outputter(out);

    outputter.pushFieldSpec("name",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 32);
    outputter.pushFieldSpec("gen. factory addr.", Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
    outputter.pushFieldSpec("req",                Teuchos::TabularOutputter::INT,    Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 3);
    outputter.pushFieldSpec("keep",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 5);
    outputter.pushFieldSpec("type",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 10);
    outputter.pushFieldSpec("data",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 20);
    outputter.outputHeader();

    std::vector<const MueLu::FactoryBase*> ehandles = RequestedFactories();
    for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
      std::vector<std::string> enames = RequestedKeys(*kt);
      for (std::vector<std::string>::iterator it = enames.begin(); it != enames.end(); ++it) {
        outputter.outputField(*it);                    // variable name
        outputter.outputField(*kt);                    // factory ptr
        int reqcount = NumRequests(*it, *kt);          // request counter
        outputter.outputField(reqcount);
        if (GetKeepFlag(*it, *kt) != 0) outputter.outputField("true");
        else outputter.outputField("false");
        // variable type
        std::string strType = GetType(*it, *kt);
        if (strType.find("Xpetra::Matrix") != std::string::npos) {
          outputter.outputField("Matrix" );
          outputter.outputField("");
        } else if (strType.find("Xpetra::MultiVector") != std::string::npos) {
          outputter.outputField("Vector");
          outputter.outputField("");
        } else if (strType == "int") {
          outputter.outputField(strType);
          int data = Get<int>(*it, *kt);
          outputter.outputField(data);
        } else if (strType == "double") {
          outputter.outputField(strType);
          double data = Get<double>(*it, *kt);
          outputter.outputField(data);
        } else if (strType == "string") {
          outputter.outputField(strType);
          std::string data = Get<std::string>(*it, *kt);
          outputter.outputField(data);
        } else {
          outputter.outputField(strType);
          outputter.outputField("unknown");
        }

        outputter.nextRow();
      }
    }
  }

} //namespace MueLu
