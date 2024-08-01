// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_Factory.hpp"

namespace MueLu {

Factory::Factory()
#ifdef HAVE_MUELU_DEBUG
  : multipleCallCheck_(FIRSTCALL)
  , lastLevelID_(-1)
#endif
{
}

//! Destructor.
Factory::~Factory() = default;

void Factory::SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
  RCP<const FactoryBase> f = factory;
  SetParameter(varName, ParameterEntry(f));  // parameter validation done in ParameterListAcceptorImpl
}

const RCP<const FactoryBase> Factory::GetFactory(const std::string& varName) const {
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

RCP<ParameterList> Factory::RemoveFactoriesFromList(const ParameterList& list) const {
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

RCP<const ParameterList> Factory::GetValidParameterList() const {
  return Teuchos::null;  // Teuchos::null == GetValidParameterList() not implemented == skip validation and no default values (dangerous)
}

void Factory::Input(Level& level, const std::string& varName) const {
  level.DeclareInput(varName, GetFactory(varName).get(), this);
}
// Similar to the other Input, but we have an alias (varParamName) to the generated data name (varName)
void Factory::Input(Level& level, const std::string& varName, const std::string& varParamName) const {
  level.DeclareInput(varName, GetFactory(varParamName).get(), this);
}

bool Factory::IsAvailable(Level& level, const std::string& varName) const {
  return level.IsAvailable(varName, GetFactory(varName).get());
}

void Factory::EnableTimerSync() { timerSync_ = true; }

void Factory::DisableTimerSync() { timerSync_ = false; }

#ifdef HAVE_MUELU_DEBUG
void Factory::EnableMultipleCallCheck() const { multipleCallCheck_ = ENABLED; }
void Factory::DisableMultipleCallCheck() const { multipleCallCheck_ = DISABLED; }
void Factory::ResetDebugData() const {
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

void Factory::EnableMultipleCheckGlobally() { multipleCallCheckGlobal_ = ENABLED; }
void Factory::DisableMultipleCheckGlobally() { multipleCallCheckGlobal_ = DISABLED; }

#else
void Factory::EnableMultipleCallCheck() const {}
void Factory::DisableMultipleCallCheck() const {}
void Factory::ResetDebugData() const {}
void Factory::EnableMultipleCheckGlobally() {}
void Factory::DisableMultipleCheckGlobally() {}
#endif

bool Factory::timerSync_ = false;
#ifdef HAVE_MUELU_DEBUG
Factory::multipleCallCheckEnum Factory::multipleCallCheckGlobal_ = ENABLED;
#endif

}  // namespace MueLu
