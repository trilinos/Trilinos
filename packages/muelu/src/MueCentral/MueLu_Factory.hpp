// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  Factory();

  //! Destructor.
  virtual ~Factory();
  //@}

  //@{
  //! Configuration

  //! SetFactory is for expert users only. To change configuration of the preconditioner, use a factory manager.
  virtual void SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory);

  //! Default implementation of FactoryAcceptor::GetFactory()
  const RCP<const FactoryBase> GetFactory(const std::string& varName) const;

  RCP<ParameterList> RemoveFactoriesFromList(const ParameterList& list) const;

  // SetParameterList(...);

  // GetParameterList(...);

  //@}

  virtual RCP<const ParameterList> GetValidParameterList() const;

 protected:
  void Input(Level& level, const std::string& varName) const;
  // Similar to the other Input, but we have an alias (varParamName) to the generated data name (varName)
  void Input(Level& level, const std::string& varName, const std::string& varParamName) const;

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

  bool IsAvailable(Level& level, const std::string& varName) const;

 public:
  static void EnableTimerSync();
  static void DisableTimerSync();

 protected:
  static bool timerSync_;

#ifdef HAVE_MUELU_DEBUG
 public:
  enum multipleCallCheckEnum{ENABLED, DISABLED, FIRSTCALL};

  void EnableMultipleCallCheck() const;
  void DisableMultipleCallCheck() const;
  void ResetDebugData() const;

  static void EnableMultipleCheckGlobally();
  static void DisableMultipleCheckGlobally();

 protected:
  mutable multipleCallCheckEnum multipleCallCheck_;
  static multipleCallCheckEnum multipleCallCheckGlobal_;
  mutable int lastLevelID_;
#else
 public:
  void EnableMultipleCallCheck() const;
  void DisableMultipleCallCheck() const;
  void ResetDebugData() const;
  static void EnableMultipleCheckGlobally();
  static void DisableMultipleCheckGlobally();
#endif
};  // class Factory

}  // namespace MueLu

#define MUELU_FACTORY_SHORT
#endif  // ifndef MUELU_FACTORY_HPP
