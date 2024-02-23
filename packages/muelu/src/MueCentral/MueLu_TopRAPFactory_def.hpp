/*
 * MueLu_TopRAPFactory_def.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: tawiesn
 */

#ifndef PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPRAPFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPRAPFACTORY_DEF_HPP_

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_FactoryManagerBase.hpp"
//#include "MueLu_HierarchyUtils_fwd.hpp"
#include "MueLu_TopRAPFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager)
  : PFact_(parentFactoryManager->GetFactory("P"))
  , RFact_(parentFactoryManager->GetFactory("R"))
  , AcFact_(parentFactoryManager->GetFactory("A")) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TopRAPFactory(RCP<const FactoryManagerBase> /* parentFactoryManagerFine */, RCP<const FactoryManagerBase> parentFactoryManagerCoarse)
  : PFact_(parentFactoryManagerCoarse->GetFactory("P"))
  , RFact_(parentFactoryManagerCoarse->GetFactory("R"))
  , AcFact_(parentFactoryManagerCoarse->GetFactory("A")) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TopRAPFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* fineLevel */, Level& coarseLevel) const {
  if (PFact_ != Teuchos::null) coarseLevel.DeclareInput("P", PFact_.get());
  if (RFact_ != Teuchos::null) coarseLevel.DeclareInput("R", RFact_.get());
  if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) coarseLevel.DeclareInput("A", AcFact_.get());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& /* fineLevel */, Level& coarseLevel) const {
  if ((PFact_ != Teuchos::null) && (PFact_ != NoFactory::getRCP())) {
    RCP<Operator> oP = coarseLevel.Get<RCP<Operator> >("P", PFact_.get());
    // Don't have a valid operator (e.g., # global aggregates is 0) so we just bail out
    // This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()j
    if (oP == Teuchos::null) return;
    RCP<Matrix> P = rcp_dynamic_cast<Matrix>(oP);
    if (!P.is_null())
      coarseLevel.Set("P", P, NoFactory::get());
    else
      coarseLevel.Set("P", oP, NoFactory::get());
    coarseLevel.AddKeepFlag("P", NoFactory::get(), MueLu::Final);        // FIXME2: Order of Remove/Add matter (data removed otherwise). Should do something about this
    coarseLevel.RemoveKeepFlag("P", NoFactory::get(), MueLu::UserData);  // FIXME: This is a hack, I should change behavior of Level::Set() instead. FIXME3: Should not be removed if flag was there already
  }

  if ((RFact_ != Teuchos::null) && (RFact_ != NoFactory::getRCP())) {
    RCP<Operator> oR = coarseLevel.Get<RCP<Operator> >("R", RFact_.get());
    RCP<Matrix> R    = rcp_dynamic_cast<Matrix>(oR);
    if (!R.is_null())
      coarseLevel.Set("R", R, NoFactory::get());
    else
      coarseLevel.Set("R", oR, NoFactory::get());
    coarseLevel.AddKeepFlag("R", NoFactory::get(), MueLu::Final);
    coarseLevel.RemoveKeepFlag("R", NoFactory::get(), MueLu::UserData);  // FIXME: This is a hack
  }

  if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) {
    RCP<Operator> oA = coarseLevel.Get<RCP<Operator> >("A", AcFact_.get());
    RCP<Matrix> A    = rcp_dynamic_cast<Matrix>(oA);
    if (!A.is_null())
      coarseLevel.Set("A", A, NoFactory::get());
    else
      coarseLevel.Set("A", oA, NoFactory::get());
    coarseLevel.AddKeepFlag("A", NoFactory::get(), MueLu::Final);
    coarseLevel.RemoveKeepFlag("A", NoFactory::get(), MueLu::UserData);  // FIXME: This is a hack
  }
}
}  // namespace MueLu

#endif /* PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPRAPFACTORY_DEF_HPP_ */
