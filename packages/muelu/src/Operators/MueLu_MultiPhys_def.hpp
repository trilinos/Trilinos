// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MULTIPHYS_DEF_HPP
#define MUELU_MULTIPHYS_DEF_HPP

#include <sstream>
#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "MueLu_MultiPhys_decl.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_HierarchyManager.hpp"
#include <MueLu_HierarchyUtils.hpp>
#include "MueLu_VerbosityLevel.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_ML2MueLuParameterTranslator.hpp>

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  return AmatMultiphysics_->getDomainMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  return AmatMultiphysics_->getRangeMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameters(Teuchos::ParameterList& list) {
  // this operator only makes sense for the Combo when using TransP for R

  list.set("multigrid algorithm", "combine");
  list.set("combine: numBlks", nBlks_);

  // Make sure verbosity gets passed to the sublists
  std::string verbosity = list.get("verbosity", "high");
  VerboseObject::SetDefaultVerbLevel(toVerbLevel(verbosity));

  arrayOfParamLists_.resize(nBlks_);
  for (int i = 0; i < nBlks_; i++) {
    std::string listName = "subblockList" + Teuchos::toString(i);
    if (list.isSublist(listName)) {
      arrayOfParamLists_[i] = Teuchos::rcpFromRef(list.sublist(listName));
    } else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Must provide sublist " + listName);

    arrayOfParamLists_[i]->set("verbosity", arrayOfParamLists_[i]->get("verbosity", verbosity));
    arrayOfParamLists_[i]->set("smoother: pre or post", "none");
    arrayOfParamLists_[i]->set("smoother: type", "none");
  }

  // Are we using Kokkos?
  useKokkos_ = !Node::is_serial;
  useKokkos_ = list.get("use kokkos refactor", useKokkos_);

  paramListMultiphysics_ = Teuchos::rcpFromRef(list);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::compute(bool reuse) {
  /*

     Create a set of AMG hierarchies whose interpolation matrices are used to build on combined
     AMG hierarchy for a multiphysics problem

   */

  //#ifdef HAVE_MUELU_CUDA
  //   if (paramListMultiphysics_.get<bool>("multiphysics: cuda profile setup", false)) cudaProfilerStart();
  //#endif

  std::string timerLabel;
  if (reuse)
    timerLabel = "MueLu MultiPhys: compute (reuse)";
  else
    timerLabel = "MueLu MultiPhys: compute";
  RCP<Teuchos::TimeMonitor> tmCompute = getTimer(timerLabel);

  ////////////////////////////////////////////////////////////////////////////////
  // Generate the (iii,iii) Hierarchy

  for (int iii = 0; iii < nBlks_; iii++) {
    if (arrayOfCoords_ != Teuchos::null) {
      arrayOfParamLists_[iii]->sublist("user data").set("Coordinates", arrayOfCoords_[iii]);
    }

    bool wantToRepartition = false;
    if (paramListMultiphysics_->isParameter("repartition: enable"))
      wantToRepartition = paramListMultiphysics_->get<bool>("repartition: enable");

    arrayOfParamLists_[iii]->set("repartition: enable", wantToRepartition);
    arrayOfParamLists_[iii]->set("repartition: rebalance P and R", true);
    arrayOfParamLists_[iii]->set("repartition: explicit via new copy rebalance P and R", true);

    if (paramListMultiphysics_->isParameter("repartition: use subcommunicators"))
      arrayOfParamLists_[iii]->set("repartition: use subcommunicators", paramListMultiphysics_->isParameter("repartition: use subcommunicators"));
    else
      arrayOfParamLists_[iii]->set("repartition: use subcommunicators", true);
  }
  // repartitioning should only happen when createing the individual P's , not
  // when combiing them

  paramListMultiphysics_->set<bool>("repartition: enable", false);

  LO maxLevels = 9999;
  for (int i = 0; i < nBlks_; i++) {
    std::string operatorLabel = "MultiPhys (" + Teuchos::toString(i) + "," + Teuchos::toString(i) + ")";
    arrayOfAuxMatrices_[i]->setObjectLabel(operatorLabel);
    arrayOfHierarchies_[i] = MueLu::CreateXpetraPreconditioner(arrayOfAuxMatrices_[i], *arrayOfParamLists_[i]);
    LO tempNlevels         = arrayOfHierarchies_[i]->GetGlobalNumLevels();
    if (tempNlevels < maxLevels) maxLevels = tempNlevels;
  }

  hierarchyMultiphysics_ = rcp(new Hierarchy("Combo"));
  for (LO i = 0; i < maxLevels; i++) {
    hierarchyMultiphysics_->AddNewLevel();
  }
  for (int i = 0; i < nBlks_; i++) {
    std::string subblkName = "Psubblock" + Teuchos::toString(i);
    MueLu::HierarchyUtils<SC, LO, GO, NO>::CopyBetweenHierarchies(*(arrayOfHierarchies_[i]), *(hierarchyMultiphysics_), "P", subblkName, "RCP<Matrix>");
  }
  paramListMultiphysics_->set("coarse: max size", 1);
  paramListMultiphysics_->set("max levels", maxLevels);

  AmatMultiphysics_->setObjectLabel("A: block " + Teuchos::toString(nBlks_) + " x " + Teuchos::toString(nBlks_) + "multiphysics matrix");

  // Rip off non-serializable data before validation
  Teuchos::ParameterList nonSerialListMultiphysics, processedListMultiphysics;
  MueLu::ExtractNonSerializableData(*paramListMultiphysics_, processedListMultiphysics, nonSerialListMultiphysics);

  // Rip off the subblock List stuff as we  don't need it any more and I think it messes up validator

  Teuchos::ParameterList stripped;
  for (ParameterList::ConstIterator inListEntry = processedListMultiphysics.begin(); inListEntry != processedListMultiphysics.end(); inListEntry++) {
    const std::string& levelName = inListEntry->first;
    if (levelName.find("subblockList") != 0) stripped.setEntry(inListEntry->first, inListEntry->second);
  }

  RCP<HierarchyManager<SC, LO, GO, NO>> mueLuFactory = rcp(new ParameterListInterpreter<SC, LO, GO, NO>(stripped, AmatMultiphysics_->getDomainMap()->getComm()));
  hierarchyMultiphysics_->setlib(AmatMultiphysics_->getDomainMap()->lib());
  hierarchyMultiphysics_->SetProcRankVerbose(AmatMultiphysics_->getDomainMap()->getComm()->getRank());

  // We don't need nullspace or coordinates, since we don't use them when just combining prolongators that have been already created
  hierarchyMultiphysics_->GetLevel(0)->Set("A", AmatMultiphysics_);

  // Stick the non-serializible data on the hierarchy.
  // Not sure that we need this, since we don't use it in building the multiphysics hierarchy
  HierarchyUtils<SC, LO, GO, NO>::AddNonSerializableDataToHierarchy(*mueLuFactory, *hierarchyMultiphysics_, nonSerialListMultiphysics);
  mueLuFactory->SetupHierarchy(*hierarchyMultiphysics_);

  describe(GetOStream(Runtime0));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Teuchos::TimeMonitor> MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTimer(std::string name, RCP<const Teuchos::Comm<int>> comm) const {
  if (IsPrint(Timings)) {
    if (!syncTimers_)
      return Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));
    else {
      if (comm.is_null())
        return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name), AmatMultiphysics_->getRowMap()->getComm().ptr()));
      else
        return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name), comm.ptr()));
    }
  } else
    return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::resetMatrix(RCP<Matrix> AmatMultiphysics_new, bool ComputePrec) {
  bool reuse        = !AmatMultiphysics_.is_null();
  AmatMultiphysics_ = AmatMultiphysics_new;
  if (ComputePrec) compute(reuse);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::applyInverse(const MultiVector& RHS, MultiVector& X) const {
  hierarchyMultiphysics_->Iterate(RHS, X, 1, true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector& RHS, MultiVector& X,
                                                                 Teuchos::ETransp /* mode */,
                                                                 Scalar /* alpha */,
                                                                 Scalar /* beta */) const {
  RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu MultiPhys: solve");
  hierarchyMultiphysics_->Iterate(RHS, X, 1, true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasTransposeApply() const {
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    initialize(const Teuchos::RCP<Matrix>& AmatMultiPhysics,
               const Teuchos::ArrayRCP<RCP<Matrix>> arrayOfAuxMatrices,
               const Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces,
               const Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords,
               const int nBlks,
               Teuchos::ParameterList& List) {
  arrayOfHierarchies_.resize(nBlks_);
  for (int i = 0; i < nBlks_; i++) arrayOfHierarchies_[i] = Teuchos::null;

  // Default settings
  useKokkos_    = false;
  enable_reuse_ = false;
  syncTimers_   = false;

  // set parameters
  setParameters(List);

}  // initialize

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel /* verbLevel */) const {
  std::ostringstream oss;

  RCP<const Teuchos::Comm<int>> comm = AmatMultiphysics_->getDomainMap()->getComm();

  oss << "\n--------------------------------------------------------------------------------\n"
      << "---                            MultiPhysics Summary                            ---\n"
         "--------------------------------------------------------------------------------"
      << std::endl;
  oss << std::endl;

  GlobalOrdinal numRows;
  GlobalOrdinal nnz;

  AmatMultiphysics_->getRowMap()->getComm()->barrier();

  for (int i = 0; i < nBlks_; i++) {
    numRows                  = arrayOfAuxMatrices_[i]->getGlobalNumRows();
    nnz                      = arrayOfAuxMatrices_[i]->getGlobalNumEntries();
    Xpetra::global_size_t tt = numRows;
    int rowspacer            = 3;
    while (tt != 0) {
      tt /= 10;
      rowspacer++;
    }
    tt            = nnz;
    int nnzspacer = 2;
    while (tt != 0) {
      tt /= 10;
      nnzspacer++;
    }

    oss << "block " << std::setw(rowspacer) << " rows " << std::setw(nnzspacer) << " nnz " << std::setw(9) << " nnz/row" << std::endl;
    oss << "(" << Teuchos::toString(i) << ", " << Teuchos::toString(i) << ")" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;
  }
  oss << std::endl;

  out << oss.str();

  for (int i = 0; i < nBlks_; i++) {
    arrayOfHierarchies_[i]->describe(out, GetVerbLevel());
  }

}  // describe

}  // namespace MueLu

#define MUELU_MULTIPHYS_SHORT
#endif  // ifdef MUELU_MULTIPHYS_DEF_HPP
