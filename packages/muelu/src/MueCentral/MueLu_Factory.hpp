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
#ifndef MUELU_FACTORY_HPP
#define MUELU_FACTORY_HPP

#include <string>
#include <deque>                         // for _Deque_iterator, operator!=
#include <ostream>                       // for operator<<, etc
#include "Teuchos_ENull.hpp"             // for ENull::null
#include "Teuchos_FilteredIterator.hpp"  // for FilteredIterator, etc
#include "Teuchos_ParameterEntry.hpp"    // for ParameterEntry
#include "Teuchos_ParameterList.hpp"     // for ParameterList, etc
#include "Teuchos_RCPDecl.hpp"           // for RCP
#include "Teuchos_RCPNode.hpp"           // for operator<<
#include "Teuchos_StringIndexedOrderedValueObjectContainer.hpp"
#include "Teuchos_RCP.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryAcceptor.hpp"
#include "MueLu_ParameterListAcceptor.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

class Factory : public FactoryBase, public FactoryAcceptor, public ParameterListAcceptorImpl {
 public:
  //@{ Constructors/Destructors.

  //! Constructor.
  Factory()
#ifdef HAVE_MUELU_DEBUG
    : multipleCallCheck_(FIRSTCALL)
    , lastLevelID_(-1)
#endif
  {
  }

  //! Destructor.
  virtual ~Factory() {}
  //@}

  //@{
  //! Configuration

  //! SetFactory is for expert users only. To change configuration of the preconditioner, use a factory manager.
  virtual void SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
    RCP<const FactoryBase> f = factory;
    SetParameter(varName, ParameterEntry(f));  // parameter validation done in ParameterListAcceptorImpl
  }

  //! Default implementation of FactoryAcceptor::GetFactory()
  const RCP<const FactoryBase> GetFactory(const std::string& varName) const {
    // Special treatment for "NoFactory"
    if (varName == "NoFactory")
      return MueLu::NoFactory::getRCP();

    if (!GetParameterList().isParameter(varName) && GetValidParameterList() == Teuchos::null) {
      // If the parameter is not on the list and there is not validator, the defaults values for 'varName' is not set.
      // Failback by using directly the FactoryManager
      // NOTE: call to GetValidParameterList() can be costly for classes that validate parameters.
      // But it get called only (lazy '&&' operator) if the parameter 'varName' is not on the paramlist and
      // the parameter 'varName' is always on the list when validator is present and 'varName' is valid (at least the default value is set).
      return Teuchos::null;
    }

    return GetParameterList().get<RCP<const FactoryBase> >(varName);
  }

  RCP<ParameterList> RemoveFactoriesFromList(const ParameterList& list) const {
    RCP<ParameterList> paramList = rcp(new ParameterList(list));
    // Remove FactoryBase entries from the list
    // The solution would be much more elegant if ParameterList support std::list like operations
    // In that case, we could simply write:
    //   for (ParameterList::ConstIterator it = paramList.begin(); it != paramList.end(); it++)
    //     if (paramList.isType<RCP<const FactoryBase> >(it->first))
    //       it = paramList.erase(it);
    //     else
    //       it++;
    ParameterList::ConstIterator it = paramList->begin();
    while (it != paramList->end()) {
      it = paramList->begin();

      for (; it != paramList->end(); it++)
        if (paramList->isType<RCP<const FactoryBase> >(it->first))
          paramList->remove(it->first);
    }
    return paramList;
  }

  // SetParameterList(...);

  // GetParameterList(...);

  //@}

  virtual RCP<const ParameterList> GetValidParameterList() const {
    return Teuchos::null;  // Teuchos::null == GetValidParameterList() not implemented == skip validation and no default values (dangerous)
  }

 protected:
  void Input(Level& level, const std::string& varName) const {
    level.DeclareInput(varName, GetFactory(varName).get(), this);
  }
  // Similar to the other Input, but we have an alias (varParamName) to the generated data name (varName)
  void Input(Level& level, const std::string& varName, const std::string& varParamName) const {
    level.DeclareInput(varName, GetFactory(varParamName).get(), this);
  }

  template <class T>
  T Get(Level& level, const std::string& varName) const {
    return level.Get<T>(varName, GetFactory(varName).get());
  }
  // Similar to the other Get, but we have an alias (varParamName) to the generated data name (varName)
  template <class T>
  T Get(Level& level, const std::string& varName, const std::string& varParamName) const {
    return level.Get<T>(varName, GetFactory(varParamName).get());
  }

  template <class T>
  void Set(Level& level, const std::string& varName, const T& data) const {
    return level.Set<T>(varName, data, this);
  }

  template <class T>
  bool IsType(Level& level, const std::string& varName) const {
    return level.IsType<T>(varName, GetFactory(varName).get());
  }

  bool IsAvailable(Level& level, const std::string& varName) const {
    return level.IsAvailable(varName, GetFactory(varName).get());
  }

 public:
  static void EnableTimerSync() { timerSync_ = true; }
  static void DisableTimerSync() { timerSync_ = false; }

 protected:
  static bool timerSync_;

#ifdef HAVE_MUELU_DEBUG
 public:
  enum multipleCallCheckEnum{ENABLED, DISABLED, FIRSTCALL};

  void EnableMultipleCallCheck() const { multipleCallCheck_ = ENABLED; }
  void DisableMultipleCallCheck() const { multipleCallCheck_ = DISABLED; }
  void ResetDebugData() const {
    if (multipleCallCheck_ == FIRSTCALL && lastLevelID_ == -1)
      return;

    multipleCallCheck_ = FIRSTCALL;
    lastLevelID_       = -1;

    const ParameterList& paramList = GetParameterList();

    // We cannot use just FactoryManager to specify which factories call ResetDebugData().
    // The problem is that some factories are not present in the manager, but
    // instead are only accessible through a parameter list of some factory.
    // For instance, FilteredAFactory is only accessible from SaPFactory but
    // nowhere else. So we miss those, and do not reset the data, resulting
    // in problems.
    // Therefore, for each factory we need to go through its dependent
    // factories, and call reset on them.
    for (ParameterList::ConstIterator it = paramList.begin(); it != paramList.end(); it++)
      if (paramList.isType<RCP<const FactoryBase> >(it->first)) {
        RCP<const Factory> fact = rcp_dynamic_cast<const Factory>(paramList.get<RCP<const FactoryBase> >(it->first));
        if (fact != Teuchos::null && fact != NoFactory::getRCP())
          fact->ResetDebugData();
      }
  }

  static void EnableMultipleCheckGlobally() { multipleCallCheckGlobal_ = ENABLED; }
  static void DisableMultipleCheckGlobally() { multipleCallCheckGlobal_ = DISABLED; }

 protected:
  mutable multipleCallCheckEnum multipleCallCheck_;
  static multipleCallCheckEnum multipleCallCheckGlobal_;
  mutable int lastLevelID_;
#else
 public:
  void EnableMultipleCallCheck() const {}
  void DisableMultipleCallCheck() const {}
  void ResetDebugData() const {}
  static void EnableMultipleCheckGlobally() {}
  static void DisableMultipleCheckGlobally() {}
#endif
};  // class Factory

}  // namespace MueLu

#define MUELU_FACTORY_SHORT
#endif  // ifndef MUELU_FACTORY_HPP
