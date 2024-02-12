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
#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <algorithm>  // for swap
#include <map>        // for _Rb_tree_const_iterator, etc
#include <ostream>    // for basic_ostream, etc
#include <string>     // for char_traits, string, etc
#include <utility>    // for pair

#include <Teuchos_Describable.hpp>       // for operator<<
#include <Teuchos_FancyOStream.hpp>      // for FancyOStream
#include <Teuchos_RCPDecl.hpp>           // for RCP
#include <Teuchos_RCP.hpp>               // for RCP::operator->, etc
#include <Teuchos_TestForException.hpp>  // for TEUCHOS_TEST_FOR_EXCEPTION

#include <Xpetra_Map.hpp>  // for UnderlyingLib definition

#include "MueLu_BoostGraphviz.hpp"
#include "MueLu_Exceptions.hpp"  // for RuntimeError
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_KeepType.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_VariableContainer.hpp"
#include "MueLu_VerbosityLevel.hpp"  // for MsgType::Default, VerbLevel

namespace MueLu {

/*!
  @class Level
  @brief Class that holds all level-specific information.

  All data that is stored in the <tt>Level</tt> class need a variable name
  (e.g. "A", "P", ...) and a pointer to the generating factory. Only with
  both the variable name and the generating factory the data can be accessed.

  If no pointer to the generating factory is provided (or it is NULL) then
  the Level class uses the information from a factory manager, which stores
  default factories for different variable names.

  We use a "two key" map for storage of data, with the pointer to the
  generating factory as primary key and the variable name as secondary key.
  The pointer to the generating factory is only used as primary "key".  It
  doesn't matter if the given factory pointer is valid or not. So the NULL
  pointer can also be used.

  The data itself is stored within a VariableContainer object.  A reference
  counter keeps track of the storage and automatically frees the memory if
  the data is not needed any more. In the standard mode, the data first has
  to be requested by calling the Request function. Then the data can be set
  by calling Set.  With Get the user can fetch data when needed. Release
  decrements the reference counter for the given variable.
*/
class Level : public BaseClass {
 public:
  //@{

  //! @name Constructors / Destructors

  Level()
    : lib_(Xpetra::NotSpecified)
    , levelID_(-1) {}

  Level(RCP<FactoryManagerBase>& factoryManager)
    : lib_(Xpetra::UseTpetra)
    , levelID_(-1)
    , factoryManager_(factoryManager) {}

  //! Destructor
  virtual ~Level() {}

  //@}

  //@{
  //! @name Build methods
  //! Builds a new Level object.
  RCP<Level> Build();

  //@}

  //@{
  //! @name Level handling

  //! @brief Return level number.
  int GetLevelID() const;

  //! @brief Set level number.
  void SetLevelID(int levelID);

  //! Previous level
  RCP<Level>& GetPreviousLevel() { return previousLevel_; }

  //! Set previous level object
  //! @\param[in] const RCP<Level>& previousLevel
  void SetPreviousLevel(const RCP<Level>& previousLevel);
  //@}

  //! @name Set/Get factory manager
  //@{
  //! Set default factories (used internally by Hierarchy::SetLevel()).
  // Users should not use this method.
  void SetFactoryManager(const RCP<const FactoryManagerBase>& factoryManager);

  //! returns the current factory manager
  // Users should not use this method
  const RCP<const FactoryManagerBase> GetFactoryManager();
  //@}

  //@{
  //! @name Set functions

  //! Store need label and its associated data. This does not increment the storage counter.
  //! - If factory is not specified or factory == NoFactory::get(), mark data as user-defined by default
  //! - If factory is == NULL, use defaultFactory (FactoryManager have to be available).
  template <class T>
  void Set(const std::string& ename, const T& entry, const FactoryBase* factory = NoFactory::get()) {
    const FactoryBase* fac = GetFactory(ename, factory);

    if (fac == NoFactory::get()) {
      // Any data set with a NoFactory gets UserData keep flag by default
      AddKeepFlag(ename, NoFactory::get(), MueLu::UserData);
    }

    // Store entry only if data have been requested (or any keep flag)
    if (IsRequested(ename, factory) || GetKeepFlag(ename, factory) != 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(!IsKey(factory, ename), Exceptions::RuntimeError, "" + ename + " not found in");
      map_[factory][ename]->SetData(entry);

    } else {
      GetOStream(Runtime1) << "Level::Set: Not storing \"" << ename << "\" generated by factory " << factory->ShortClassName() << "[" << factory->GetID() << "]"
                           << " on level " << toString(GetLevelID()) << ", as it has not been requested and no keep flags were set for it" << std::endl;
    }
  }  // Set

  //@}

  //! @name Get functions
  //! @brief Get functions for accessing stored data

  //@{
  /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).
   *   Usage: Level->Get< RCP<Matrix> >("A", factory)
   *   if factory == NULL => use default factory
   *
   *  @param[in] const std::string& ename
   *  @param[in] const FactoryBase* factory
   *  @return data (templated)
   * */
  template <class T>
  T& Get(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
    const FactoryBase* fac = GetFactory(ename, factory);
    /*      printf("(l=%d)                                               getting    \"%20s\" generated by %10p  [actually, generated by %p (%43s)]\n",
            levelID_, ename.c_str(), factory, fac, fac->description().c_str());*/

    TEUCHOS_TEST_FOR_EXCEPTION(!IsKey(fac, ename), Exceptions::RuntimeError, "\"" + ename + "\" generated by factory \"" + fac->description() + "\" not found on level " + toString(GetLevelID()) + ".");

    if (!IsAvailable(ename, fac)) {
      TEUCHOS_TEST_FOR_EXCEPTION(NumRequests(fac, ename) < 1 && GetKeepFlag(ename, fac) == 0, Exceptions::RuntimeError,
                                 "\"" << ename << "\" has not been requested (counter = " << NumRequests(fac, ename) << ", "
                                                                                                                        "KeepFlag = "
                                      << GetKeepFlag(ename, fac) << "). " << std::endl
                                      << "Generating factory:" << *fac << " NoFactory = " << NoFactory::get());
      fac->CallBuild(*this);
      Release(*fac);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(!IsAvailable(ename, fac), Exceptions::RuntimeError,
                               "MueLu::Level::Get(): factory did not produce expected output on level " << GetLevelID()
                                                                                                        << ". The data \"" << ename << "\" has not been generated by " << *fac);

    return map_[fac][ename]->template GetData<T>();
  }

  /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).*/
  template <class T>
  void Get(const std::string& ename, T& rValue, const FactoryBase* factory = NoFactory::get()) {
    rValue = Get<T>(ename, factory);
  }

  template <class T>
  bool IsType(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
    const FactoryBase* fac = GetFactory(ename, factory);

    TEUCHOS_TEST_FOR_EXCEPTION(!IsKey(fac, ename), Exceptions::RuntimeError, "\"" + ename + "\" generated by factory \"" + fac->description() + "\" not found on level " + toString(GetLevelID()) + ".");

    if (!IsAvailable(ename, fac)) {
      TEUCHOS_TEST_FOR_EXCEPTION(NumRequests(fac, ename) < 1 && GetKeepFlag(ename, fac) == 0, Exceptions::RuntimeError,
                                 "\"" << ename << "\" has not been requested (counter = " << NumRequests(fac, ename) << ", "
                                                                                                                        "KeepFlag = "
                                      << GetKeepFlag(ename, fac) << "). " << std::endl
                                      << "Generating factory:" << *fac << " NoFactory = " << NoFactory::get());
      fac->CallBuild(*this);
      Release(*fac);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(!IsAvailable(ename, fac), Exceptions::RuntimeError,
                               "MueLu::Level::Get(): factory did not produce expected output on level " << GetLevelID()
                                                                                                        << ". The data \"" << ename << "\" has not been generated by " << *fac);

    return map_[fac][ename]->template CheckType<T>();
  }

  /*! @brief GetTypeName returns type string of variable stored using ename and factory
   *
   *  @param[in] const std::string& ename
   *  @param[in] const FactoryBase* factory
   *  @return std::string with type information (e.g., "int")
   * */
  std::string GetTypeName(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
    const FactoryBase* fac = GetFactory(ename, factory);
    TEUCHOS_TEST_FOR_EXCEPTION(!IsKey(fac, ename), Exceptions::RuntimeError, "\"" + ename + "\" not found");

    TEUCHOS_TEST_FOR_EXCEPTION(!IsAvailable(ename, fac), Exceptions::RuntimeError,
                               "MueLu::Level::GetTypeString(): Data "
                               "\"" << ename
                                    << "\" generated by " << *fac << " is not available.");

    return map_[fac][ename]->GetTypeName();
  }

  //@}

  //! @name Permanent storage
  //@{

  //! Request to keep variable 'ename' generated by 'factory' after the setup phase.
  //  This method is intented to be used by user drivers for printing, debugging or to keep some computed data for a next run of the setup phase.
  //
  //  This method is an alias for: AddKeepFlag(ename, factory, MueLu::Keep)
  //  See also the description of KeepEnum for more information.
  //
  //  To undo a keep request, one can use:
  //  - Delete(ename, factory) to delete the data and remove the "Keep" flag
  //  - or RemoveKeepFlag(ename, factory, MueLu::Keep) to go back to previous condition (data are kept only if internal MueLu logic need so).
  //
  // Note: Level variables tagged using this methods are also keep on the levels built using this.Build().
  //       This means that is you request to keep a specific variable on a fine level, all the coarser level that are created automatically during the setup phase will also retain the same variable.
  void Keep(const std::string& ename, const FactoryBase* factory) { AddKeepFlag(ename, factory, MueLu::Keep); }  // Note: do not add default value for input parameter 'factory'

  //! Delete data that have been retained after the setup phase (using Keep(), AddKeepFlag(), or internal MueLu logic).
  // Special cases:
  // - If entry (ename, factory) does not exist, nothing is done.
  // - If entry exists but counter !=, entry cannot be desallocated before counter set to 0 (using Release()) so an exeption is thrown.
  void Delete(const std::string& ename, const FactoryBase* factory) {  // Note: do not add default value for input parameter 'factory'
    if (!IsKey(factory, ename))
      return;

    // Precondition:
    // Delete() should only be called if counter == 0
    // Note: It better to throw an exception rather than deleting the data if counter != 0 because users are not supposed to manipulate data with counter != 0
    TEUCHOS_TEST_FOR_EXCEPTION(IsRequested(ename, factory) == true, Exceptions::RuntimeError, "MueLu::Level::Delete(): IsRequested() == true. Ref counter != 0. You are not allowed to delete data that are still in use.");
    // If counter == 0 and entry exists, this means that a keep flag is set. Or there is an internal logic problem.
    TEUCHOS_TEST_FOR_EXCEPTION(GetKeepFlag(ename, factory) == 0, Exceptions::RuntimeError, "MueLu::Level::Delete(), Keep flag == 0?");

    RemoveKeepFlag(ename, factory, MueLu::All);  // will delete the data if counter == 0

    // Post condition: data must have been deleted
    TEUCHOS_TEST_FOR_EXCEPTION(IsAvailable(ename, factory) == true, Exceptions::RuntimeError, "MueLu::Level::Delete(): Internal error (Post condition). Data have not been deleted.");
  }

  //! Delete all data that have been retained after the setup phase using Final flag
  void Clear();

  //! Delete all data from level that has no Keep flag set.
  //! This is a function for experts only
  void ExpertClear();

  //! Test if a keep flag is set for variable 'ename' generated by 'factory'
  //! The input parameter keep can be a combination of flags. IsKept() will then return true if at least one of the flag is set.
  //! Note: There is no default parameter for IsKept() because it might be confusing (user generally wants to test IsKept with keep=MueLu::Keep but classes Level and Needs generally use keep = All)
  bool IsKept(const std::string& ename, const FactoryBase* factory, KeepType keep) const { return GetKeepFlag(ename, factory) & keep; }

  //! Add a keep flag for variable 'ename' generated by 'factory'
  //! A variable can cumulate several keep flags (UserData+Final for example). This function just add a flag to the current flag combination.
  //! By default, the flag 'Keep' is added. So Keep(ename, factory) is the same as AddKeepFlag(ename, factory).
  //! See also the description of KeepEnum for more information.
  void AddKeepFlag(const std::string& ename, const FactoryBase* factory = NoFactory::get(), KeepType keep = MueLu::Keep);  // TODO: remove default value for input parameter 'factory'?

  //! Remove a keep flag for variable 'ename' generated by 'factory'
  //! A variable can cumulate several keep flags (UserData+Final for example). This function just add a flag to the current flag combination.
  //! By default, all the flags are removed.
  void RemoveKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keep = MueLu::All);

  //! Get the flag combination set for variable 'ename' generated by 'factory'
  KeepType GetKeepFlag(const std::string& ename, const FactoryBase* factory) const;

  //@}

  //! @name Request/Release functions
  //! @brief Request and Release for incrementing/decrementing the reference count pointer for a specific variable.
  //@{

  //! Increment the storage counter for all the inputs of a factory
  void Request(const FactoryBase& factory);

  //! Decrement the storage counter for all the inputs of a factory
  void Release(const FactoryBase& factory);

  //! Callback from FactoryBase::CallDeclareInput() and FactoryBase::DeclareInput()
  void DeclareInput(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy = NoFactory::get());

  //! Callback from FactoryBase::CallDeclareInput() and FactoryBase::DeclareInput() to declare factory dependencies
  void DeclareDependencies(const FactoryBase* factory, bool bRequestOnly = false, bool bReleaseOnly = false);

  //! Indicate that an object is needed. This increments the storage counter.
  void Request(const std::string& ename, const FactoryBase* factory = NoFactory::get(), const FactoryBase* requestedBy = NoFactory::get());

  //! Decrement the storage counter.
  void Release(const std::string& ename, const FactoryBase* factory = NoFactory::get(), const FactoryBase* requestedBy = NoFactory::get());

  //@}

  //! @name Utility functions
  //@{

  //! Test whether a need's value has been saved.
  bool IsAvailable(const std::string& ename, const FactoryBase* factory = NoFactory::get()) const {
    if (!IsKey(factory, ename))
      return false;
    try {
      return Get(factory, ename)->IsAvailable();
    } catch (...) {
      return false;
    }
  }

  //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
  bool IsRequested(const std::string& ename, const FactoryBase* factory = NoFactory::get()) const {
    if (!IsKey(factory, ename))
      return false;
    try {
      return IsRequested(Get(factory, ename));
    } catch (...) {
      return false;
    }
  }

  //@}
  //! @name I/O Functions
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Printing method
  // TODO: print only shows requested variables. check if we also list kept factories with ref counter=0?
  void print(std::ostream& out, const VerbLevel verbLevel = Default) const;

#if defined(HAVE_MUELU_BOOST) && defined(HAVE_MUELU_BOOST_FOR_REAL) && defined(BOOST_VERSION) && (BOOST_VERSION >= 104400)
  void UpdateGraph(std::map<const FactoryBase*, BoostVertex>& vindices,
                   std::map<std::pair<BoostVertex, BoostVertex>, std::string>& edges,
                   BoostProperties& dp,
                   BoostGraph& graph) const;
#endif

  //@}

  enum RequestMode { REQUEST,
                     RELEASE,
                     UNDEF };
  RequestMode GetRequestMode() const { return requestMode_; }

  void setlib(Xpetra::UnderlyingLib lib2) { lib_ = lib2; }
  Xpetra::UnderlyingLib lib() { return lib_; }

  void SetComm(RCP<const Teuchos::Comm<int> > const& comm) { comm_ = comm; }
  RCP<const Teuchos::Comm<int> > GetComm() const { return comm_; }

 private:
  //! Copy constructor.
  Level(const Level& source);

  //! If input factory == NULL, returns the default factory. Else, return input factory.
  //
  //  If factory == NULL, the default factory is defined as follow:
  // - If user data is available, it is considered as the default and the factory manager is ignored.
  //   => The default factory is then NoFactory.
  // - Else, the factory manager is used to get the default factory.
  //
  // This strategy allows to use the same factory manager on the fine and coarse level without any trouble.
  //     Example :
  //
  //     FineLevel:
  //     ----------
  //     A          -> User provided
  //     Nullspace  -> User provided
  //
  //     CoarseLevel:
  //     ------------
  //     A         -> RAPFactory
  //     NullSpace -> NullspaceFactory
  //
  const FactoryBase* GetFactory(const std::string& varname, const FactoryBase* factory) const;

  static RequestMode requestMode_;
  Xpetra::UnderlyingLib lib_;
  RCP<const Teuchos::Comm<int> > comm_;

  typedef const FactoryBase* Key1;
  typedef const std::string Key2;
  typedef RCP<VariableContainer> Value;
  typedef Teuchos::map<Key2, Value> SubMap;      //! Sub-map container (Key2 -> Value)
  typedef Teuchos::map<Key1, SubMap> TwoKeyMap;  //! Map of a map (Key1 -> SubMap)

  int levelID_;  // id number associated with level
  RCP<const FactoryManagerBase> factoryManager_;
  RCP<Level> previousLevel_;  // linked list of Level
  TwoKeyMap map_;

  //! @name Utility functions
  //@{

  //! Test whether some information about (ename, factory) are stored
  bool IsKey(const FactoryBase* factory, const std::string& ename) const {
    TwoKeyMap::const_iterator it = map_.find(factory);
    return (it != map_.end()) ? (it->second).count(ename) : false;
  }

  bool IsAvailableFactory(const FactoryBase* factory) const {
    TwoKeyMap::const_iterator it = map_.find(factory);
    if (it == map_.end())
      return false;
    for (SubMap::const_iterator sit = it->second.begin(); sit != it->second.end(); sit++) {
      if (sit->second->IsAvailable())
        return true;
    }
    return false;
  }

  bool IsRequested(const Value& v) const {
    TEUCHOS_TEST_FOR_EXCEPTION(v->NumAllRequests() == 0 && v->GetKeepFlag() == 0, Exceptions::RuntimeError,
                               "Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return v->IsRequested();
  }

  bool IsRequestedBy(const FactoryBase* factory, const std::string& ename, const FactoryBase* requestedBy) const {
    if (!IsKey(factory, ename))
      return false;

    return IsRequestedBy(Get(factory, ename), requestedBy);
  }

  bool IsRequestedBy(const Value& v, const FactoryBase* requestedBy) const {
    TEUCHOS_TEST_FOR_EXCEPTION(v->NumAllRequests() == 0 && v->GetKeepFlag() == 0, Exceptions::RuntimeError,
                               "Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return v->IsRequested(requestedBy);
  }

  bool IsRequestedFactory(const FactoryBase* factory) const {
    TwoKeyMap::const_iterator it = map_.find(factory);
    if (it == map_.end())
      return false;
    for (SubMap::const_iterator sit = it->second.begin(); sit != it->second.end(); sit++)
      if (IsRequested(sit->second))
        return true;
    return false;
  }

  const Value& Get(const FactoryBase* factory, const std::string& ename) const {
    TwoKeyMap::const_iterator it = map_.find(factory);
    TEUCHOS_TEST_FOR_EXCEPTION(it == map_.end(), Exceptions::RuntimeError, "Key (" << factory << ", *) does not exist.");

    SubMap::const_iterator sit = it->second.find(ename);
    TEUCHOS_TEST_FOR_EXCEPTION(sit == it->second.end(), Exceptions::RuntimeError, "Key (" << factory << ", " << ename << ") does not exist.");

    return sit->second;
  }

  int NumRequests(const FactoryBase* factory, const std::string& ename) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!IsKey(factory, ename), Exceptions::RuntimeError, "\"" + ename + "\" not found. Do a request first.");
    const Teuchos::RCP<MueLu::VariableContainer>& v = Get(factory, ename);
    TEUCHOS_TEST_FOR_EXCEPTION(v->NumAllRequests() == 0 && v->GetKeepFlag() == 0, Exceptions::RuntimeError,
                               "NumRequests(): Internal logic error: if counter == 0, the entry in countTable_ should have been deleted");
    return v->NumAllRequests();
  }

  int CountRequestedFactory(const FactoryBase* factory) const {
    TwoKeyMap::const_iterator it = map_.find(factory);
    if (it == map_.end())
      return 0;

    int cnt = 0;
    for (SubMap::const_iterator sit = it->second.begin(); sit != it->second.end(); sit++)
      cnt += sit->second->NumAllRequests();

    return cnt;
  }

  //@}

};  // class Level

}  // namespace MueLu

// TODO: Caps should not matter

#endif  // MUELU_LEVEL_HPP
