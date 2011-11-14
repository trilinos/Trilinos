#include "MueLu_Needs.hpp"

namespace MueLu {
  
  Needs::Needs() { }

  Needs::~Needs() { }

  void Needs::Request(const std::string & ename, const FactoryBase* factory) {
    // If it's the first request for 'ename', create a new entry in the map
    if (!countTable_.IsKey(ename, factory))
      countTable_.Set(ename, factory, 0);

    // Increment counter
    int currentCount = countTable_.Get(ename, factory);
    //if (currentCount != -1) { // doesn't matter if counter disabled
    countTable_.Set(ename, factory, ++currentCount);
    //}

  } //Request

  void Needs::Release(const std::string & ename, const FactoryBase* factory) {
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

  } //Release

  void Needs::Keep(const std::string & ename, const FactoryBase* factory, bool keep ) {
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

  }

  bool Needs::IsKept(const std::string & ename, const FactoryBase* factory) const {
    if (keepTable_.IsKey(ename, factory))
      return keepTable_.Get(ename,factory);
    return false;
  }

  bool Needs::IsAvailable(const std::string & ename, const FactoryBase* factory) const {
    return dataTable_.IsKey(ename, factory);
  }

  bool Needs::IsRequested(const std::string & ename, const FactoryBase* factory) const {
    TEUCHOS_TEST_FOR_EXCEPTION(countTable_.IsKey(ename, factory) && (countTable_.Get(ename, factory) == 0) && (IsKept(ename, factory)==false), Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    if (countTable_.IsKey(ename, factory) == false) return false;
    bool bRequested = NumRequests(ename, factory) > 0 ? true : false;
    return bRequested;

    //TEUCHOS_TEST_FOR_EXCEPTION(countTable_.IsKey(ename, factory) && (countTable_.Get(ename, factory) == 0), Exceptions::RuntimeError, "MueLu::Needs::IsRequested(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    //return countTable_.IsKey(ename, factory);
  }

  bool Needs::IsRequestedFactory(const FactoryBase* factory) { //TODO: rename HaveBeenRequested() !!
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
  }

  int Needs::CountRequestedFactory(const FactoryBase* factory) {
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
  }

  bool Needs::IsAvailableFactory(const FactoryBase* factory) {
    std::vector<std::string> ekeys = dataTable_.GetKeyList();
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      std::vector<const FactoryBase*> ehandles = dataTable_.GetKey2List(*it);
      for (std::vector<const FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
        if (*kt == factory) // factory is generating factory of requested variable '*it'
          return true;
      }
    }
    return false;
  }

  int Needs::NumRequests(const std::string & ename, const FactoryBase* factory) const {
    //FIXME should we return 0 instead of throwing an exception?
    TEUCHOS_TEST_FOR_EXCEPTION(!countTable_.IsKey(ename, factory), Exceptions::RuntimeError, "MueLu::Needs::NumRequests(): " + ename + " not found in countTable_");
    return countTable_.Get(ename, factory);
  }

  std::vector<std::string> Needs::RequestedKeys() const {
    return countTable_.GetKeyList();
  }

  std::vector<const FactoryBase*> Needs::RequestedFactories(const std::string & ename) const {
    return countTable_.GetKey2List(ename);
  }

  std::string Needs::GetType(const std::string & ename, const FactoryBase* factory) const {
    return dataTable_.Get(ename, factory).getAny(true).typeName();
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

    std::vector<std::string> ekeys = countTable_.GetKeyList();
    for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
      std::vector<const FactoryBase*> ehandles = countTable_.GetKey2List(*it);
      for (std::vector<const FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
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
