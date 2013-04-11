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

#include "MueLu_Level.hpp"

#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

  Level::Level() : levelID_(-1) { }

  Level::Level(RCP<FactoryManagerBase> & factoryManager) : levelID_(-1), factoryManager_(factoryManager) { }

  RCP<Level> Level::Build() {
    RCP<Level> newLevel = rcp( new Level() );

    // Copy 'keep' status of variables
    // TODO: this only concerns needs_. so a function in Needs class should be provided to do that!
    // TODO: how can i move this to Needs? maybe we need a new constructor for Level which gets a
    // Needs object...

    std::vector<const MueLu::FactoryBase*> ehandles = needs_.RequestedFactories();
    for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
      std::vector<std::string> enames = needs_.RequestedKeys(*kt);
      for (std::vector<std::string>::iterator it = enames.begin(); it != enames.end(); ++it) {
        const std::string & ename = *it;
        const MueLu::FactoryBase* fac = *kt;
        if (IsKept(ename, fac, MueLu::Keep)) { // MueLu::Keep is the only flag propagated
          if (fac == NULL) // TODO: Is this possible?? Throw exception. Not supposed to use the FactoryManager here.
            newLevel->Keep(ename, NoFactory::get());
          else
            newLevel->Keep(ename, fac);
        }
      }
    }

    return newLevel;
  }

  Level::~Level() {}

  int Level::GetLevelID() const { return levelID_; }

  void Level::SetLevelID(int levelID) {
    if (levelID_ != -1 && levelID_ != levelID)
      GetOStream(Warnings1, 0) << "Warning: Level::SetLevelID(): Changing an already defined LevelID (previousID=" << levelID_ << ", newID=" << levelID << ")" << std::endl;

    levelID_ = levelID;
  }

  RCP<Level> & Level::GetPreviousLevel() { return previousLevel_; }

  void Level::SetPreviousLevel(const RCP<Level> & previousLevel) {
    if (previousLevel_ != Teuchos::null && previousLevel_ != previousLevel)
      GetOStream(Warnings1, 0) << "Warning: Level::SetPreviousLevel(): PreviousLevel was already defined" << std::endl;

    previousLevel_ = previousLevel;
  }

  void Level::SetFactoryManager(const RCP<const FactoryManagerBase> & factoryManager) {
    factoryManager_ = factoryManager;
  }

  const RCP<const FactoryManagerBase> Level::GetFactoryManager() {
    return factoryManager_;
  }

  void Level::AddKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keepType) {
    needs_.AddKeepFlag(ename, GetFactory(ename, factory), keepType);
  }

  void Level::RemoveKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keepType) {
    needs_.RemoveKeepFlag(ename, GetFactory(ename, factory), keepType);
  }

  KeepType Level::GetKeepFlag(const std::string& ename, const FactoryBase* factory) const {
    return needs_.GetKeepFlag(ename, GetFactory(ename, factory));
  }

  void Level::Request(const FactoryBase& factory) {
    RequestMode prev = requestMode_;
    requestMode_ = REQUEST;
    factory.CallDeclareInput(*this);
    requestMode_ = prev;
  }

  void Level::Release(const FactoryBase& factory) {
    RequestMode prev = requestMode_;
    requestMode_ = RELEASE;
    factory.CallDeclareInput(*this);
    requestMode_ = prev;
  }

  void Level::DeclareInput(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    if (requestMode_ == REQUEST) {
      Request(ename, factory, requestedBy);
    }
    else if (requestMode_ == RELEASE) {
      Release(ename, factory, requestedBy);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareInput(): requestMode_ undefined.");
  }

  void Level::DeclareDependencies(const FactoryBase* factory, bool bRequestOnly, bool bReleaseOnly) { //TODO: replace bReleaseOnly, bReleaseOnly by one RequestMode enum
    if (bRequestOnly && bReleaseOnly)
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareDependencies(): Both bRequestOnly and bReleaseOnly set to true makes no sense.");

    if (requestMode_ == REQUEST) {

      if (bReleaseOnly == false) Request(*factory);

    } else if (requestMode_ == RELEASE) {

      if (bRequestOnly == false) Release(*factory);

    } else TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareDependencies(): requestMode_ undefined.");
  }

  void Level::Request(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    const FactoryBase* fac = GetFactory(ename, factory);

    // Call request for factory only if the factory has not been requested before and no data has
    // been generated by fact (independent of ename)
    bool test = ((needs_.IsAvailableFactory(fac) == false && needs_.IsRequestedFactory(fac) == false));

    // This request must be done before calling Request(*fac) to avoid circular dependency problems.
    needs_.Request(ename, fac, requestedBy);

    // The value of IsRequestedFactory(fac) is true, due to the above request.
    // That is why a temporary boolean "test" is used!
    TEUCHOS_TEST_FOR_EXCEPTION(needs_.IsRequestedFactory(fac) != true, Exceptions::RuntimeError, "Level::Request(ename, factory): internal logic error.");

    // Call Request for factory dependencies
    if (test) {
      Request(*fac);
    }
  }

  //TODO: finish this
// #define MUELU_LEVEL_ERROR_MESSAGE(function, message)
//   "MueLu::Level[" << levelID_ << "]::" << function << "(" << ename << ", " << factory << " << ): " << message << std::endl
//                                                                                                                  << ((factory == Teuchos::null && factoryManager_ != Teuchos::null) ? (*GetFactory(ename, factory)) : *factory)  << "Generating factory:" << *fac << " NoFactory=" << NoFactory::get()
//   //

  void Level::Release(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    const FactoryBase* fac = GetFactory(ename, factory);

    // Only a factory which has requested (fac,ename) is allowed to release it again.
    // Do not release data if it has not been requested by the factory "requestedBy"
    // Note: when data is released (fac,ename) depends on it often happened that some
    //       of this data has (recursively) been released too often
    if(needs_.IsRequestedBy(fac, ename, requestedBy)) {
      if(needs_.CountRequestedFactory(fac) == 1 &&     // check if factory fac is not requested by another factory
         needs_.IsAvailableFactory(fac) == false ) {   // check if Build function of factory fac has been called
        // In general all data (fac,ename) depends on is released with the Get calls in the Build
        // function of the generating factory fac.
        // Here we have to release the dependencies of some data that has been requested (by factory "requestedBy")
        // but the corresponding Build function of factory "fac" has never been called. Therefore the dependencies
        // have never been released. Do it now.
        Release(*fac);
      }
      needs_.Release(ename,fac,requestedBy);
    }
  }

  bool Level::IsKey(const std::string & ename, const FactoryBase* factory) const {
    return needs_.IsKey(ename, GetFactory(ename, factory));
  }

  bool Level::IsAvailable(const std::string & ename, const FactoryBase* factory) const {
    return needs_.IsAvailable(ename, GetFactory(ename, factory));
  }

  bool Level::IsRequested(const std::string & ename, const FactoryBase* factory) const {
    return needs_.IsRequested(ename, GetFactory(ename, factory));
  }

  std::string Level::description() const {
    std::ostringstream out;
    out << BaseClass::description();
    out << "{ levelID = " << levelID_ << "}";
    return out.str();
  }

  void Level::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    //MUELU_DESCRIBE;
    //out0 << ""; // remove warning

    RCP<Teuchos::FancyOStream> out0 = Teuchos::rcpFromRef(out);
    int previousSetting = out0->getOutputToRootOnly();
    out0->setOutputToRootOnly(0);
    out0->setShowProcRank(true);

    *out0 << "LevelID = " << GetLevelID() << std::endl;

    typedef Teuchos::TabularOutputter TTO;
    TTO outputter(out0);
    outputter.pushFieldSpec("data name",                TTO::STRING, TTO::LEFT, TTO::GENERAL, 20);
    // outputter.pushFieldSpec("generating factory type",  TTO::STRING, TTO::LEFT, TTO::GENERAL, 30);
    outputter.pushFieldSpec("gen. factory addr.",       TTO::STRING, TTO::LEFT, TTO::GENERAL, 18);
    outputter.pushFieldSpec("req",                      TTO::INT,    TTO::LEFT, TTO::GENERAL, 3);
    outputter.pushFieldSpec("keep",                     TTO::STRING, TTO::LEFT, TTO::GENERAL, 5);
    outputter.pushFieldSpec("type",                     TTO::STRING, TTO::LEFT, TTO::GENERAL, 15);
    outputter.pushFieldSpec("data",                     TTO::STRING, TTO::LEFT, TTO::GENERAL, 14);
    outputter.pushFieldSpec("req'd by",                 TTO::STRING, TTO::LEFT, TTO::GENERAL, 20);
    outputter.outputHeader();

    out0->setOutputToRootOnly(-1);
    std::vector<const MueLu::FactoryBase*> ehandles = needs_.RequestedFactories();
    for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
      std::vector<std::string> enames = needs_.RequestedKeys(*kt);
      for (std::vector<std::string>::iterator it = enames.begin(); it != enames.end(); ++it) {
        outputter.outputField(*it);   // variable name

        // NOTE: We cannot dereference the factory pointer and call factory->description() as we do not know
        // if the factory still exist (the factory pointer is a raw pointer by design). Instead, the level
        // should store the factory description internally as a string for debugging purpose (and in debug mode only).
        //         // factory name
        //         std::stringstream ss1;
        //         ss1 << (*kt)->description();
        //         outputter.outputField((ss1.str()).substr(0,30));

        // factory ptr
        outputter.outputField(*kt);

        int reqcount = needs_.NumRequests(*it, *kt); // request counter
        outputter.outputField(reqcount);

        KeepType keepType = needs_.GetKeepFlag(*it, *kt);
        if (keepType != 0) {
          std::stringstream ss;
          if (keepType & MueLu::UserData) { ss << "User";  }
          if (keepType & MueLu::Keep)     { ss << "Keep";  }
          if (keepType & MueLu::Final)    { ss << "Final"; }
          outputter.outputField(ss.str());
        } else {
          outputter.outputField("No");
        }

        if (needs_.IsAvailable(*it, *kt)) {
          std::string strType = needs_.GetType(*it, *kt); // Variable type
          if (strType.find("Xpetra::Matrix") != std::string::npos) {
            outputter.outputField("Matrix" );
            outputter.outputField("available");
          } else if (strType.find("Xpetra::MultiVector") != std::string::npos) {
            outputter.outputField("Vector");
            outputter.outputField("available");
          } else if (strType.find("Xpetra::Map") != std::string::npos) {
            outputter.outputField("Map");
            outputter.outputField("available");
          } else if (strType.find("MueLu::SmootherBase") != std::string::npos) {
            outputter.outputField("SmootherBase");
            outputter.outputField("available");
          } else if (strType.find("MueLu::Aggregates") != std::string::npos) {
            outputter.outputField("Aggregates");
            outputter.outputField("available");
          } else if (strType.find("MueLu::AmalgamationInfo") != std::string::npos) {
            outputter.outputField("AmalgamationInfo");
            outputter.outputField("available");
          } else if (strType.find("MueLu::Graph") != std::string::npos) {
            outputter.outputField("Graph");
            outputter.outputField("available");
          } else if (strType == "int") {
            outputter.outputField(strType);
            int data = needs_.Get<int>(*it, *kt);
            outputter.outputField(data);
          } else if (strType == "double") {
            outputter.outputField(strType);
            double data = needs_.Get<double>(*it, *kt);
            outputter.outputField(data);
          } else if (strType == "string") {
            outputter.outputField(strType);
            std::string data = needs_.Get<std::string>(*it, *kt);
            outputter.outputField(data);
          } else {
            outputter.outputField(strType);
            outputter.outputField("available");
          }
        } else {
          outputter.outputField("unknown");
          outputter.outputField("not available");
        }

        typedef VariableContainer::request_container container_type;
        const container_type& requestedBy = needs_.GetRequests(*it, *kt);
        std::ostringstream ss;
        for (container_type::const_iterator ct = requestedBy.begin(); ct != requestedBy.end(); ct++) {
          if (ct != requestedBy.begin())    ss << ",";
                                            ss << ct->first;
          if (ct->second > 1)               ss << "(" << ct->second << ")";
        }
        outputter.outputField(ss.str());

        outputter.nextRow();
      }
    } //for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++)
    out0->setOutputToRootOnly(previousSetting);
    out0->setShowProcRank(false);

  }

  Level::Level(const Level& source) { }

  // JG Note: should the option IgnoreUserData() moved to the Factory interface or on the specific factories that are using this option? It would simplify the level class.
  const FactoryBase* Level::GetFactory(const std::string& varname, const FactoryBase* factory) const {
    if (factory != NULL)
      return factory;

    // If IgnoreUserData == false and if variable "varname" generated by NoFactory is provided by the user (MueLu::UserData), use user-provided data by default without querying the FactoryManager.
    // When FactoryManager == null, we consider that IgnoreUserData == false.
    if ((factoryManager_ == Teuchos::null || factoryManager_->IgnoreUserData() == false) &&
        (IsAvailable(varname, NoFactory::get()) && IsKept(varname, NoFactory::get(), MueLu::UserData))) {
      return NoFactory::get();
    }

    // Query factory manager
    TEUCHOS_TEST_FOR_EXCEPTION(factoryManager_ == null, Exceptions::RuntimeError, "MueLu::Level("<< levelID_ << ")::GetFactory(" << varname << ", " << factory << "): No FactoryManager");
    const FactoryBase* fac = factoryManager_->GetFactory(varname).get();
    TEUCHOS_TEST_FOR_EXCEPTION(fac == NULL, Exceptions::RuntimeError, "MueLu::Level("<< levelID_ << ")::GetFactory(" << varname << ", " << factory << "): Default factory returned by FactoryManager cannot be NULL");
    return fac;
  }

  Level::RequestMode Level::requestMode_ = UNDEF;

} //namespace MueLu

//TODO: Caps should not matter
