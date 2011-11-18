#include <Teuchos_TabularOutputter.hpp>

#include "MueLu_Needs.hpp"

#include "MueLu_NoFactory.hpp"

namespace MueLu {
  
  Needs::Needs() { }

  Needs::~Needs() { }

  void Needs::Request(const std::string & ename, const FactoryBase* factory) {
#if OLD
    // If it's the first request for 'ename', create a new entry in the map
    if (!countTable_.IsKey(ename, factory))
      countTable_.Set(ename, factory, 0);

    // Increment counter
    int currentCount = countTable_.Get(ename, factory);
    countTable_.Set(ename, factory, ++currentCount);
#else
    if(!variableTable_.IsKey(factory, ename)) {
      Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
      variableTable_.Set(factory, ename, newVar);
    }

    Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    var->Request(MueLu::NoFactory::get());
#endif
  } //Request

  void Needs::Request(const std::string & ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
#if OLD
    // Increment counter (called after Request(..,..) has been called, so counter is already incremented
    if (!requestTable_.IsKey(factory, ename)) {
      std::map<const FactoryBase*, int> newrequests;
      newrequests[requestedBy] = 0;
      requestTable_.Set(factory, ename, newrequests);
    }

    std::map<const FactoryBase*,int>& requests = requestTable_.Get(factory,ename);
    int currentCount = requests[requestedBy];
    requests[requestedBy] = ++currentCount;
#else
    if(!variableTable_.IsKey(factory, ename)) {
      Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
      variableTable_.Set(factory, ename, newVar);
    }

    Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    var->Request(requestedBy);
#endif
  } //Request

  void Needs::Release(const std::string & ename, const FactoryBase* factory) {
#if OLD
    // We can only release data if the key 'ename' exists in countTable (and dataTable)
    TEUCHOS_TEST_FOR_EXCEPTION(!countTable_.IsKey(ename, factory), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");

    // Decrement reference counter
    int currentCount = countTable_.Get(ename, factory);
    countTable_.Set(ename, factory, --currentCount);

    // Desallocation if counter becomes zero and data shall not be kept
    if (currentCount == 0 && IsKept(ename, factory) == false) {
      countTable_.Remove(ename, factory);   // Desallocation of the counter
      if (dataTable_.IsKey(ename, factory))
        dataTable_.Remove(ename, factory);  // Desallocation of the data
    }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");

    Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    var->Release(MueLu::NoFactory::get());

    if (var->IsRequested() == false && var->IsKept() == false ) {
      var = Teuchos::null; // free data
      variableTable_.Remove(factory, ename);
    }
#endif
  }

  void Needs::Release(const std::string & ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
#if OLD
    // Decrement reference counter
    if (!requestTable_.IsKey(factory, ename)) {
      std::cout << "STRANGE: " << ename << " gen by " << factory << " is not in the request table any more!" << std::endl;
      std::cout << "how can this be released?" << std::endl;
      return;
    }

    std::map<const FactoryBase*,int>& requests = requestTable_.Get(factory,ename);
    int currentCount = requests[requestedBy];
    requests[requestedBy] = --currentCount;

    // Desallocation for requestTable
    if (currentCount == 0 && IsKept(ename, factory) == false) {
      int nElementErased = requests.erase(requestedBy);
      TEUCHOS_TEST_FOR_EXCEPTION(nElementErased != 1, Exceptions::RuntimeError, "cannot remove data");

      if (requests.size() == 0)
        requestTable_.Remove(factory,ename);
    }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");

    Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    var->Release(requestedBy);

    if (var->IsRequested() == false && var->IsKept() == false ) {
      var = Teuchos::null; // free data
      variableTable_.Remove(factory, ename);
    }
#endif
  } //Release

  void Needs::Keep(const std::string & ename, const FactoryBase* factory, bool keep ) {
#if 0
    //Request(ename, factory);
    if (keep) {
      // if not listed in countTable_ yet, add it to count Table
      if (!countTable_.IsKey(ename, factory))
        countTable_.Set(ename, factory, 0);

      keepTable_.Set(ename, factory, true);
    }
    else {
      if(keepTable_.IsKey(ename, factory) == false) return; // variable is not kept
      keepTable_.Set(ename, factory, false);

      // Desallocation if request counter is zero
      int currentCount = countTable_.Get(ename, factory);
      if (currentCount == 0) {
        countTable_.Remove(ename, factory);   // Desallocation of the counter
        if (dataTable_.IsKey(ename, factory))
          dataTable_.Remove(ename, factory);  // Desallocation of the data
      }
    }
#else
    if (keep) {
      if(!variableTable_.IsKey(factory, ename)) {
        Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
        variableTable_.Set(factory, ename, newVar);
      }

      Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
      var->Keep();
    } else {
      if(!variableTable_.IsKey(factory,ename)) return;
      Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
      var->Keep(false);
      if (var->IsRequested() == false) {
        var = Teuchos::null; // free data
        variableTable_.Remove(factory, ename);
      }
    }
#endif
  }

  bool Needs::IsKept(const std::string & ename, const FactoryBase* factory) const {
#if OLD
    if (keepTable_.IsKey(ename, factory))
      return keepTable_.Get(ename,factory);
    return false;
#else
    if(!variableTable_.IsKey(factory,ename)) return false;
    //TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::Release(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    return var->IsKept();
#endif
  }

  bool Needs::IsAvailable(const std::string & ename, const FactoryBase* factory) const {
#if OLD
    return dataTable_.IsKey(ename, factory);
#else
    //TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsAvailable(): " + ename + " not found. Do a request first.");
    if(!variableTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    return var->IsAvailable();
#endif
  }

  bool Needs::IsRequested(const std::string & ename, const FactoryBase* factory) const {
#if OLD
    TEUCHOS_TEST_FOR_EXCEPTION(countTable_.IsKey(ename, factory) && (countTable_.Get(ename, factory) == 0) && (IsKept(ename, factory)==false), Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    if (countTable_.IsKey(ename, factory) == false) return false;
    bool bRequested = NumRequests(ename, factory) > 0 ? true : false;
    return bRequested;

    //TEUCHOS_TEST_FOR_EXCEPTION(countTable_.IsKey(ename, factory) && (countTable_.Get(ename, factory) == 0), Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    //return countTable_.IsKey(ename, factory);
#else
    //TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): " + ename + " not found. Do a request first.");
    if(!variableTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->IsRequested();
#endif
  }

  bool Needs::IsRequestedBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const {
#if OLD
    //if(!IsRequested(ename,factory)) return false;
    if (requestTable_.IsKey(factory,ename) == false) return false;
    std::map<const FactoryBase*,int> requests = requestTable_.Get(factory,ename);

    for(std::map<const FactoryBase*,int>::iterator it = requests.begin(); it!=requests.end(); it++) {
      std::cout << "ISREQUESTEDBY: " << ename << " gen by " << factory << " requested by " << it->first << " " << it->second << " times " << std::endl;
    }

    if(requests.count(requestedBy) > 0) {
      if(requests.find(requestedBy)->second > 0) return true;
      else return false;
    }
    return false;
#else
    //TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::IsRequestedBy(): " + ename + " not found. Do a request first.");
    if(!variableTable_.IsKey(factory,ename)) return false;
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::IsRequestedBy(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->IsRequested(requestedBy);
#endif
  }

  bool Needs::IsRequestedFactory(const FactoryBase* factory) { //TODO: rename HaveBeenRequested() !!
#if 0
    std::vector<std::string> ekeys = RequestedKeys();
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      std::vector<const FactoryBase*> ehandles = RequestedFactories(*it);
      for (std::vector<const FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
        if (*kt == factory) { // factory is generating factory of requested variable '*it'
          return IsRequested(*it,*kt);
        }
      }
    }
    return false;
#else
    std::vector<const FactoryBase*> factories = RequestedFactories(); // TODO simplify me
    for (std::vector<const FactoryBase*>::iterator kt = factories.begin(); kt != factories.end(); kt++) {
      if(*kt == factory) {
        std::vector<std::string> ekeys = RequestedKeys(factory);
        for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
          if(IsRequested(*it,factory)) return true;
        }
      }
    }
    return false;
#endif
  }

  int Needs::CountRequestedFactory(const FactoryBase* factory) {
#if OLD
    int cnt = 0;
    std::vector<std::string> ekeys = RequestedKeys();
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      std::vector<const FactoryBase*> ehandles = RequestedFactories(*it);
      for (std::vector<const FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
        if (*kt == factory) { // factory is generating factory of requested variable '*it'
          cnt++;
          //cnt += NumRequests(*it,factory);
        }
      }
    }
    return cnt;
#else
    int cnt = 0;
    std::vector<const FactoryBase*> factories = RequestedFactories(); // TODO simplify me
    for (std::vector<const FactoryBase*>::iterator kt = factories.begin(); kt != factories.end(); kt++) {
      if(*kt == factory) {
        std::vector<std::string> ekeys = RequestedKeys(factory);
        for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
          cnt += NumRequests(*it, factory);
        }
      }
    }
    return cnt;
#endif
  }

  bool Needs::IsAvailableFactory(const FactoryBase* factory) {
#if OLD
    std::vector<std::string> ekeys = dataTable_.GetKeyList();
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      std::vector<const FactoryBase*> ehandles = dataTable_.GetKey2List(*it);
      for (std::vector<const FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
        if (*kt == factory) // factory is generating factory of requested variable '*it'
          return true;
      }
    }
    return false;
#else
    std::vector<const FactoryBase*> factories = RequestedFactories(); // TODO simplify me
    for (std::vector<const FactoryBase*>::iterator kt = factories.begin(); kt != factories.end(); kt++) {
      if(*kt == factory) {
        std::vector<std::string> ekeys = RequestedKeys(factory);
        for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
          if(IsAvailable(*it, factory))
            return true;
        }
      }
    }
    return false;
#endif
  }

  int Needs::NumRequests(const std::string & ename, const FactoryBase* factory) const {
#if OLD
    //FIXME should we return 0 instead of throwing an exception?
    TEUCHOS_TEST_FOR_EXCEPTION(!countTable_.IsKey(ename, factory), Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): " + ename + " not found in countTable_");
    return countTable_.Get(ename, factory);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->NumAllRequests();
#endif
  }

  int Needs::NumRequestsBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const {
#if OLD
    if(!IsRequested(ename,factory)) return 0;

    std::map<const FactoryBase*,int> requests = requestTable_.Get(factory,ename);
    if(requests.count(requestedBy)>0) return requests[requestedBy];
    else return 0;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory,ename), Exceptions::RuntimeError, "MueLu::Needs::NumRequestsBy(): " + ename + " not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    TEUCHOS_TEST_FOR_EXCEPTION(var->NumAllRequests() == 0 && var->IsKept() == false, Exceptions::RuntimeError, "MueLu::Needs::NumRequestsBy(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return var->NumRequests(requestedBy);
#endif
  }


#if OLD
  std::vector<std::string> Needs::RequestedKeys() const {
    return countTable_.GetKeyList();
  }
#else
  std::vector<std::string> Needs::RequestedKeys(const FactoryBase* factory) const {
    return variableTable_.GetKey2List(factory);
  }
#endif

#if OLD
  std::vector<const FactoryBase*> Needs::RequestedFactories(const std::string & ename) const {
    return countTable_.GetKey2List(ename);
  }
#else
  std::vector<const FactoryBase*> Needs::RequestedFactories() const {
    return variableTable_.GetKeyList();
  }
#endif

  std::string Needs::GetType(const std::string & ename, const FactoryBase* factory) const {
#if OLD
    return dataTable_.Get(ename, factory).getAny(true).typeName();
#else
    TEUCHOS_TEST_FOR_EXCEPTION(!variableTable_.IsKey(factory, ename), Exceptions::RuntimeError, "MueLu::Needs::Get(): " + ename + " not found in dataTable_");
    const Teuchos::RCP<MueLu::VariableContainer> & var = variableTable_.Get(factory,ename);
    return var->GetData().getAny(true).typeName();
#endif
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
        int reqcount = countTable_.Get(*it, *kt); // request counter
        outputter.outputField(reqcount);
        if (keepTable_.Get(*it, *kt)) outputter.outputField("true");
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
