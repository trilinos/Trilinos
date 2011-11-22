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

    if (var->IsRequested() == false && var->IsKept() == false ) {
      var = Teuchos::null; // free data
      dataTable_.Remove(factory, ename);
    }
  } //Release

  void Needs::Keep(const std::string & ename, const FactoryBase* factory, bool keep ) {
    if (keep) {
      if(!dataTable_.IsKey(factory, ename)) {
        Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
        dataTable_.Set(factory, ename, newVar);
      }

      Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
      var->Keep();
    } else {
      if(!dataTable_.IsKey(factory,ename)) return;
      Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
      var->Keep(false);
      if (var->IsRequested() == false) {
        var = Teuchos::null; // free data
        dataTable_.Remove(factory, ename);
      }
    }
  }

  bool Needs::IsKept(const std::string & ename, const FactoryBase* factory) const {
    if(!dataTable_.IsKey(factory,ename)) return false;
    //TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    return var->IsKept();
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
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->IsRequested();
  }

  bool Needs::IsRequestedBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const {
    //TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsRequestedBy(): " + ename + " not found. Do a request first.");
    if(!dataTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::IsRequestedBy(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->IsRequested(requestedBy);
  }

  bool Needs::IsRequestedFactory(const FactoryBase* factory) {
    if(dataTable_.IsKey(factory) == false) return false;

    std::vector<std::string> ekeys = RequestedKeys(factory);
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      if(IsRequested(*it,factory)) return true;
    }
    return false;
  }

  int Needs::CountRequestedFactory(const FactoryBase* factory) {
    int cnt = 0;
    if(dataTable_.IsKey(factory) == false) return cnt;
    std::vector<std::string> ekeys = RequestedKeys(factory);
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      cnt += NumRequests(*it, factory);
    }
    return cnt;
  }

  bool Needs::IsAvailableFactory(const FactoryBase* factory) {
    if(dataTable_.IsKey(factory) == false) return false;
    std::vector<std::string> ekeys = RequestedKeys(factory);
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      if(IsAvailable(*it, factory))
        return true;
    }
    return false;
  }

  int Needs::NumRequests(const std::string & ename, const FactoryBase* factory) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->NumAllRequests();
  }

  int Needs::NumRequestsBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::NumRequestsBy(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::NumRequestsBy(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
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
    outputter.pushFieldSpec("keep",                Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 5);
    outputter.pushFieldSpec("type",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 10);
    outputter.pushFieldSpec("data",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 20);
    outputter.outputHeader();

    std::vector<const MueLu::FactoryBase*> ehandles = RequestedFactories();
    for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
      std::vector<std::string> enames = RequestedKeys(*kt);
      for (std::vector<std::string>::iterator it = enames.begin(); it != enames.end(); it++) {
        outputter.outputField(*it);                    // variable name
        outputter.outputField(*kt);                    // factory ptr          
        int reqcount = NumRequests(*it, *kt);          // request counter
        outputter.outputField(reqcount);
        if (IsKept(*it, *kt)) outputter.outputField("true");
        else outputter.outputField("false");
        // variable type
        std::string strType = GetType(*it, *kt);
        if (strType.find("Xpetra::Operator") != std::string::npos) {
          outputter.outputField("Operator" );
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
