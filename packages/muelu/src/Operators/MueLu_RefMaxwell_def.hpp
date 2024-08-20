// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REFMAXWELL_DEF_HPP
#define MUELU_REFMAXWELL_DEF_HPP

#include <sstream>

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_TripleMatrixMultiply.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "MueLu_RefMaxwell_decl.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Maxwell_Utils.hpp"

#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"

#include "MueLu_VerbosityLevel.hpp"

#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_ML2MueLuParameterTranslator.hpp>

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

// Stratimikos
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
#include <Xpetra_ThyraLinearOp.hpp>
#endif

namespace MueLu {

template <typename T>
T pop(Teuchos::ParameterList &pl, std::string const &name_in) {
  T result = pl.get<T>(name_in);
  pl.remove(name_in, true);
  return result;
}

template <typename T>
T pop(Teuchos::ParameterList &pl, std::string const &name_in, T def_value) {
  T result = pl.get<T>(name_in, def_value);
  pl.remove(name_in, false);
  return result;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
Matrix2CrsMatrix(Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &matrix) {
  return Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(matrix, true)->getCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  return SM_Matrix_->getDomainMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  return SM_Matrix_->getRangeMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Teuchos::ParameterList>
RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getValidParamterList() {
  bool useKokkosDefault = !Node::is_serial;

  RCP<ParameterList> params = rcp(new ParameterList("RefMaxwell"));

  params->set<RCP<Matrix> >("Dk_1", Teuchos::null);
  params->set<RCP<Matrix> >("Dk_2", Teuchos::null);
  params->set<RCP<Matrix> >("D0", Teuchos::null);

  params->set<RCP<Matrix> >("M1_beta", Teuchos::null);
  params->set<RCP<Matrix> >("M1_alpha", Teuchos::null);
  // for backwards compatibility
  params->set<RCP<Matrix> >("Ms", Teuchos::null);

  params->set<RCP<Matrix> >("Mk_one", Teuchos::null);
  params->set<RCP<Matrix> >("Mk_1_one", Teuchos::null);
  // for backwards compatibility
  params->set<RCP<Matrix> >("M1", Teuchos::null);

  params->set<RCP<Matrix> >("invMk_1_invBeta", Teuchos::null);
  params->set<RCP<Matrix> >("invMk_2_invAlpha", Teuchos::null);
  // for backwards compatibility
  params->set<RCP<Matrix> >("M0inv", Teuchos::null);

  params->set<RCP<MultiVector> >("Nullspace", Teuchos::null);
  params->set<RCP<RealValuedMultiVector> >("Coordinates", Teuchos::null);

  auto spaceValidator = rcp(new Teuchos::EnhancedNumberValidator<int>(1, 2));
  params->set("refmaxwell: space number", 1, "", spaceValidator);
  params->set("verbosity", MasterList::getDefault<std::string>("verbosity"));
  params->set("use kokkos refactor", useKokkosDefault);
  params->set("half precision", false);
  params->set("parameterlist: syntax", MasterList::getDefault<std::string>("parameterlist: syntax"));
  params->set("output filename", MasterList::getDefault<std::string>("output filename"));
  params->set("print initial parameters", MasterList::getDefault<bool>("print initial parameters"));
  params->set("refmaxwell: disable addon", MasterList::getDefault<bool>("refmaxwell: disable addon"));
  params->set("refmaxwell: disable addon 22", true);
  params->set("refmaxwell: mode", MasterList::getDefault<std::string>("refmaxwell: mode"));
  params->set("refmaxwell: use as preconditioner", MasterList::getDefault<bool>("refmaxwell: use as preconditioner"));
  params->set("refmaxwell: dump matrices", MasterList::getDefault<bool>("refmaxwell: dump matrices"));
  params->set("refmaxwell: enable reuse", MasterList::getDefault<bool>("refmaxwell: enable reuse"));
  params->set("refmaxwell: skip first (1,1) level", MasterList::getDefault<bool>("refmaxwell: skip first (1,1) level"));
  params->set("refmaxwell: skip first (2,2) level", false);
  params->set("multigrid algorithm", "Unsmoothed");
  params->set("transpose: use implicit", MasterList::getDefault<bool>("transpose: use implicit"));
  params->set("rap: triple product", MasterList::getDefault<bool>("rap: triple product"));
  params->set("rap: fix zero diagonals", true);
  params->set("rap: fix zero diagonals threshold", MasterList::getDefault<double>("rap: fix zero diagonals threshold"));
  params->set("fuse prolongation and update", MasterList::getDefault<bool>("fuse prolongation and update"));
  params->set("refmaxwell: subsolves on subcommunicators", MasterList::getDefault<bool>("refmaxwell: subsolves on subcommunicators"));
  params->set("refmaxwell: subsolves striding", 1);
  params->set("refmaxwell: row sum drop tol (1,1)", MasterList::getDefault<double>("aggregation: row sum drop tol"));
  params->set("sync timers", false);
  params->set("refmaxwell: num iters coarse 11", 1);
  params->set("refmaxwell: num iters 22", 1);
  params->set("refmaxwell: apply BCs to Anodal", false);
  params->set("refmaxwell: apply BCs to coarse 11", true);
  params->set("refmaxwell: apply BCs to 22", true);
  params->set("refmaxwell: max coarse size", 1);

  ParameterList &precList11 = params->sublist("refmaxwell: 11list");
  precList11.disableRecursiveValidation();
  ParameterList &precList22 = params->sublist("refmaxwell: 22list");
  precList22.disableRecursiveValidation();

  params->set("smoother: type", "CHEBYSHEV");
  ParameterList &smootherList = params->sublist("smoother: params");
  smootherList.disableRecursiveValidation();
  params->set("smoother: pre type", "NONE");
  ParameterList &preSmootherList = params->sublist("smoother: pre params");
  preSmootherList.disableRecursiveValidation();
  params->set("smoother: post type", "NONE");
  ParameterList &postSmootherList = params->sublist("smoother: post params");
  postSmootherList.disableRecursiveValidation();

  ParameterList &matvecParams = params->sublist("matvec params");
  matvecParams.disableRecursiveValidation();

  params->set("multigrid algorithm", "unsmoothed");
  params->set("aggregation: type", MasterList::getDefault<std::string>("aggregation: type"));
  params->set("aggregation: drop tol", MasterList::getDefault<double>("aggregation: drop tol"));
  params->set("aggregation: drop scheme", MasterList::getDefault<std::string>("aggregation: drop scheme"));
  params->set("aggregation: distance laplacian algo", MasterList::getDefault<std::string>("aggregation: distance laplacian algo"));
  params->set("aggregation: min agg size", MasterList::getDefault<int>("aggregation: min agg size"));
  params->set("aggregation: max agg size", MasterList::getDefault<int>("aggregation: max agg size"));
  params->set("aggregation: match ML phase1", MasterList::getDefault<bool>("aggregation: match ML phase1"));
  params->set("aggregation: match ML phase2a", MasterList::getDefault<bool>("aggregation: match ML phase2a"));
  params->set("aggregation: match ML phase2b", MasterList::getDefault<bool>("aggregation: match ML phase2b"));
  params->set("aggregation: export visualization data", MasterList::getDefault<bool>("aggregation: export visualization data"));

  return params;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameters(Teuchos::ParameterList &list) {
  if (list.isType<std::string>("parameterlist: syntax") && list.get<std::string>("parameterlist: syntax") == "ml") {
    Teuchos::ParameterList newList;
    {
      Teuchos::ParameterList newList2                = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list, "refmaxwell"));
      RCP<Teuchos::ParameterList> validateParameters = getValidParamterList();
      for (auto it = newList2.begin(); it != newList2.end(); ++it) {
        const std::string &entry_name = it->first;
        if (validateParameters->isParameter(entry_name)) {
          ParameterEntry theEntry = newList2.getEntry(entry_name);
          newList.setEntry(entry_name, theEntry);
        }
      }
    }

    if (list.isSublist("refmaxwell: 11list") && list.sublist("refmaxwell: 11list").isSublist("edge matrix free: coarse"))
      newList.sublist("refmaxwell: 11list") = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse"), "SA"));
    if (list.isSublist("refmaxwell: 22list"))
      newList.sublist("refmaxwell: 22list") = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list.sublist("refmaxwell: 22list"), "SA"));
    list = newList;
  }

  parameterList_ = list;
  parameterList_.validateParametersAndSetDefaults(*getValidParamterList());
  std::string verbosityLevel = parameterList_.get<std::string>("verbosity");
  VerboseObject::SetDefaultVerbLevel(toVerbLevel(verbosityLevel));
  std::string outputFilename = parameterList_.get<std::string>("output filename");
  if (outputFilename != "")
    VerboseObject::SetMueLuOFileStream(outputFilename);
  if (parameterList_.isType<Teuchos::RCP<Teuchos::FancyOStream> >("output stream"))
    VerboseObject::SetMueLuOStream(parameterList_.get<Teuchos::RCP<Teuchos::FancyOStream> >("output stream"));

  if (parameterList_.get<bool>("print initial parameters"))
    GetOStream(static_cast<MsgType>(Runtime1), 0) << parameterList_ << std::endl;
  disable_addon_             = parameterList_.get<bool>("refmaxwell: disable addon");
  disable_addon_22_          = parameterList_.get<bool>("refmaxwell: disable addon 22");
  mode_                      = parameterList_.get<std::string>("refmaxwell: mode");
  use_as_preconditioner_     = parameterList_.get<bool>("refmaxwell: use as preconditioner");
  dump_matrices_             = parameterList_.get<bool>("refmaxwell: dump matrices");
  enable_reuse_              = parameterList_.get<bool>("refmaxwell: enable reuse");
  implicitTranspose_         = parameterList_.get<bool>("transpose: use implicit");
  fuseProlongationAndUpdate_ = parameterList_.get<bool>("fuse prolongation and update");
  skipFirst11Level_          = parameterList_.get<bool>("refmaxwell: skip first (1,1) level");
  skipFirst22Level_          = parameterList_.get<bool>("refmaxwell: skip first (2,2) level");
  if (spaceNumber_ == 1)
    skipFirst22Level_ = false;
  syncTimers_         = parameterList_.get<bool>("sync timers");
  useKokkos_          = parameterList_.get<bool>("use kokkos refactor");
  numItersCoarse11_   = parameterList_.get<int>("refmaxwell: num iters coarse 11");
  numIters22_         = parameterList_.get<int>("refmaxwell: num iters 22");
  applyBCsToAnodal_   = parameterList_.get<bool>("refmaxwell: apply BCs to Anodal");
  applyBCsToCoarse11_ = parameterList_.get<bool>("refmaxwell: apply BCs to coarse 11");
  applyBCsTo22_       = parameterList_.get<bool>("refmaxwell: apply BCs to 22");

  precList11_ = parameterList_.sublist("refmaxwell: 11list");
  if (!precList11_.isType<std::string>("Preconditioner Type") &&
      !precList11_.isType<std::string>("smoother: type") &&
      !precList11_.isType<std::string>("smoother: pre type") &&
      !precList11_.isType<std::string>("smoother: post type")) {
    precList11_.set("smoother: type", "CHEBYSHEV");
    precList11_.sublist("smoother: params").set("chebyshev: degree", 2);
    precList11_.sublist("smoother: params").set("chebyshev: ratio eigenvalue", 5.4);
    precList11_.sublist("smoother: params").set("chebyshev: eigenvalue max iterations", 30);
  }

  precList22_ = parameterList_.sublist("refmaxwell: 22list");
  if (!precList22_.isType<std::string>("Preconditioner Type") &&
      !precList22_.isType<std::string>("smoother: type") &&
      !precList22_.isType<std::string>("smoother: pre type") &&
      !precList22_.isType<std::string>("smoother: post type")) {
    precList22_.set("smoother: type", "CHEBYSHEV");
    precList22_.sublist("smoother: params").set("chebyshev: degree", 2);
    precList22_.sublist("smoother: params").set("chebyshev: ratio eigenvalue", 7.0);
    precList22_.sublist("smoother: params").set("chebyshev: eigenvalue max iterations", 30);
  }

  if (!parameterList_.isType<std::string>("smoother: type") && !parameterList_.isType<std::string>("smoother: pre type") && !parameterList_.isType<std::string>("smoother: post type")) {
    list.set("smoother: type", "CHEBYSHEV");
    list.sublist("smoother: params").set("chebyshev: degree", 2);
    list.sublist("smoother: params").set("chebyshev: ratio eigenvalue", 20.0);
    list.sublist("smoother: params").set("chebyshev: eigenvalue max iterations", 30);
  }

  if (enable_reuse_ &&
      !precList11_.isType<std::string>("Preconditioner Type") &&
      !precList11_.isParameter("reuse: type"))
    precList11_.set("reuse: type", "full");
  if (enable_reuse_ &&
      !precList22_.isType<std::string>("Preconditioner Type") &&
      !precList22_.isParameter("reuse: type"))
    precList22_.set("reuse: type", "full");

  // This should be taken out again as soon as
  // CoalesceDropFactory_kokkos supports BlockSize > 1 and
  // drop tol != 0.0
  if (useKokkos_ && precList11_.isParameter("aggregation: drop tol") && precList11_.get<double>("aggregation: drop tol") != 0.0) {
    GetOStream(Warnings0) << solverName_ + "::compute(): Setting \"aggregation: drop tol\". to 0.0, since CoalesceDropFactory_kokkos does not "
                          << "support BlockSize > 1 and drop tol != 0.0" << std::endl;
    precList11_.set("aggregation: drop tol", 0.0);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::compute(bool reuse) {
  using memory_space = typename Node::device_type::memory_space;

#ifdef HAVE_MUELU_CUDA
  if (parameterList_.get<bool>("refmaxwell: cuda profile setup", false)) cudaProfilerStart();
#endif

  std::string timerLabel;
  if (reuse)
    timerLabel = "compute (reuse)";
  else
    timerLabel = "compute";
  RCP<Teuchos::TimeMonitor> tmCompute = getTimer(timerLabel);

  ////////////////////////////////////////////////////////////////////////////////
  // COMMENTED OUT SINCE WE SHOULD NOT NEED THIS ANYMORE.
  // Remove explicit zeros from matrices
  // Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(parameterList_,D0_,SM_Matrix_,Mk_one_,M1_beta_);
  // if (!Dk_1_.is_null())
  //   Dk_1_ = Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(Dk_1_, 1e-10, false);

  if (IsPrint(Statistics2)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    params->set("printCommInfo", true);
    GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*SM_Matrix_, "SM_Matrix", params);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Detect Dirichlet boundary conditions
  if (!reuse) {
    magnitudeType rowSumTol = parameterList_.get<double>("refmaxwell: row sum drop tol (1,1)");
    Maxwell_Utils<SC, LO, GO, NO>::detectBoundaryConditionsSM(SM_Matrix_, Dk_1_, rowSumTol,
                                                              BCrows11_, BCcols22_, BCdomain22_,
                                                              globalNumberBoundaryUnknowns11_,
                                                              globalNumberBoundaryUnknowns22_,
                                                              onlyBoundary11_, onlyBoundary22_);
    if (spaceNumber_ == 2) {
      Kokkos::View<bool *, memory_space> BCcolsEdge   = Kokkos::View<bool *, memory_space>(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), Dk_1_->getColMap()->getLocalNumElements());
      Kokkos::View<bool *, memory_space> BCdomainEdge = Kokkos::View<bool *, memory_space>(Kokkos::ViewAllocateWithoutInitializing("dirichletDomains"), Dk_1_->getDomainMap()->getLocalNumElements());
      Utilities::DetectDirichletColsAndDomains(*Dk_1_, BCrows11_, BCcolsEdge, BCdomainEdge);

      Kokkos::View<bool *, memory_space> BCcolsNode   = Kokkos::View<bool *, memory_space>(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), D0_->getColMap()->getLocalNumElements());
      Kokkos::View<bool *, memory_space> BCdomainNode = Kokkos::View<bool *, memory_space>(Kokkos::ViewAllocateWithoutInitializing("dirichletDomains"), D0_->getDomainMap()->getLocalNumElements());
      Utilities::DetectDirichletColsAndDomains(*D0_, BCdomainEdge, BCcolsNode, BCdomainNode);
      BCdomain22_ = BCdomainNode;
    }
    if (IsPrint(Statistics2)) {
      GetOStream(Statistics2) << solverName_ + "::compute(): Detected " << globalNumberBoundaryUnknowns11_ << " BC rows and " << globalNumberBoundaryUnknowns22_ << " BC columns." << std::endl;
    }
    dump(BCrows11_, "BCrows11.m");
    dump(BCcols22_, "BCcols22.m");
    dump(BCdomain22_, "BCdomain22.m");
  }

  if (onlyBoundary11_) {
    // All unknowns of the (1,1) block have been detected as boundary unknowns.
    // Do not attempt to construct sub-hierarchies, but just set up a single level preconditioner.
    GetOStream(Warnings0) << "All unknowns of the (1,1) block have been detected as boundary unknowns!" << std::endl;
    mode_ = "none";
    setFineLevelSmoother11();
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////

  dim_ = NodalCoords_->getNumVectors();

  ////////////////////////////////////////////////////////////////////////////////
  // build special prolongators
  if (!reuse) {
    //////////////////////////////////////////////////////////////////////////////
    // build nullspace for (1,1)-block (if necessary)
    if (Nullspace11_ != null) {  // no need to do anything - nullspace is built
      TEUCHOS_ASSERT(Nullspace11_->getMap()->isCompatible(*(SM_Matrix_->getRowMap())));
    } else if (NodalCoords_ != null) {
      Nullspace11_ = buildNullspace(spaceNumber_, BCrows11_, skipFirst11Level_);
    } else {
      GetOStream(Errors) << solverName_ + "::compute(): either the nullspace or the nodal coordinates must be provided." << std::endl;
    }

    // build special prolongator for (1,1)-block
    {
      RCP<Matrix> A11_nodal;
      if (skipFirst11Level_) {
        // Form A11_nodal = D0^T * M1_beta * D0  (aka TMT_agg)
        std::string label("D0^T*M1_beta*D0");
        A11_nodal = Maxwell_Utils<SC, LO, GO, NO>::PtAPWrapper(M1_beta_, D0_, parameterList_, label);

        if (applyBCsToAnodal_) {
          // Apply boundary conditions to A11_nodal
          Utilities::ApplyOAZToMatrixRows(A11_nodal, BCdomain22_);
        }
        A11_nodal->setObjectLabel(solverName_ + " (1,1) A_nodal");
        dump(A11_nodal, "A11_nodal.m");
      }
      // release it because we won't need it anymore
      M1_beta_ = Teuchos::null;

      // build special prolongator
      buildProlongator(spaceNumber_, A11_nodal, Nullspace11_, P11_, NullspaceCoarse11_, CoordsCoarse11_);

      dump(P11_, "P11.m");
    }

    //////////////////////////////////////////////////////////////////////////////
    // build nullspace for (2,2)-block (if necessary)
    if (Nullspace22_ != null) {
      TEUCHOS_ASSERT(Nullspace22_->getMap()->isCompatible(*(Dk_1_->getDomainMap())));
    } else if (NodalCoords_ != null)
      Nullspace22_ = buildNullspace(spaceNumber_ - 1, BCdomain22_, skipFirst22Level_);
    else {
      GetOStream(Errors) << solverName_ + "::compute(): either the nullspace or the nodal coordinates must be provided." << std::endl;
    }

    // build special prolongator for (2,2)-block
    {
      RCP<Matrix> A22_nodal;
      if (skipFirst22Level_) {
        // Form A22_nodal = D0^T * M1_alpha * D0
        std::string label("D0^T*M1_alpha*D0");
        A22_nodal = Maxwell_Utils<SC, LO, GO, NO>::PtAPWrapper(M1_alpha_, D0_, parameterList_, label);

        if (applyBCsToAnodal_) {
          // Apply boundary conditions to A22_nodal
          Utilities::ApplyOAZToMatrixRows(A22_nodal, BCdomain22_);
        }
        A22_nodal->setObjectLabel(solverName_ + " (2,2) A_nodal");
        dump(A22_nodal, "A22_nodal.m");
      }
      // release it because we won't need it anymore
      M1_alpha_ = Teuchos::null;

      // build special prolongator
      buildProlongator(spaceNumber_ - 1, A22_nodal, Nullspace22_, P22_, CoarseNullspace22_, Coords22_);

      dump(P22_, "P22.m");
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // build coarse grid operator for (1,1)-block
  buildCoarse11Matrix();

  ////////////////////////////////////////////////////////////////////////////////
  // determine the communicator sizes for (1,1)- and (2,2)-blocks
  bool doRebalancing;
  int rebalanceStriding, numProcsCoarseA11, numProcsA22;
  if (!reuse)
    this->determineSubHierarchyCommSizes(doRebalancing, rebalanceStriding, numProcsCoarseA11, numProcsA22);
  else
    doRebalancing = false;

  // rebalance the coarse A11 matrix, as well as P11, CoordsCoarse11 and Addon11
  if (!reuse && doRebalancing)
    rebalanceCoarse11Matrix(rebalanceStriding, numProcsCoarseA11);
  if (!coarseA11_.is_null()) {
    dump(coarseA11_, "coarseA11.m");
    if (!reuse) {
      dumpCoords(CoordsCoarse11_, "CoordsCoarse11.m");
      dump(NullspaceCoarse11_, "NullspaceCoarse11.m");
    }
  }

  if (!reuse) {
    if (!implicitTranspose_) {
      R11_ = Utilities::Transpose(*P11_);
      dump(R11_, "R11.m");
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  // build multigrid for coarse (1,1)-block
  if (!coarseA11_.is_null()) {
    VerbLevel verbosityLevel = VerboseObject::GetDefaultVerbLevel();
    std::string label("coarseA11");
    setupSubSolve(HierarchyCoarse11_, thyraPrecOpH_, coarseA11_, NullspaceCoarse11_, CoordsCoarse11_, precList11_, label, reuse);
    VerboseObject::SetDefaultVerbLevel(verbosityLevel);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Apply BCs to columns of Dk_1
  if (!reuse && applyBCsTo22_) {
    GetOStream(Runtime0) << solverName_ + "::compute(): nuking BC columns of Dk_1" << std::endl;

    Dk_1_->resumeFill();
    Scalar replaceWith = (Dk_1_->getRowMap()->lib() == Xpetra::UseEpetra) ? Teuchos::ScalarTraits<SC>::eps() : Teuchos::ScalarTraits<SC>::zero();
    Utilities::ZeroDirichletCols(Dk_1_, BCcols22_, replaceWith);
    Dk_1_->fillComplete(Dk_1_->getDomainMap(), Dk_1_->getRangeMap());
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Build A22 = Dk_1^T SM Dk_1 and hierarchy for A22
  if (!onlyBoundary22_) {
    GetOStream(Runtime0) << solverName_ + "::compute(): building MG for (2,2)-block" << std::endl;

    // Build A22 = Dk_1^T * SM * Dk_1 and rebalance it, as well as Dk_1_ and P22_ and Coords22_
    build22Matrix(reuse, doRebalancing, rebalanceStriding, numProcsA22);

    if (!P22_.is_null()) {
      std::string label("P22^T*A22*P22");
      coarseA22_ = Maxwell_Utils<SC, LO, GO, NO>::PtAPWrapper(A22_, P22_, parameterList_, label);
      coarseA22_->SetFixedBlockSize(A22_->GetFixedBlockSize());
      coarseA22_->setObjectLabel(solverName_ + " coarse (2, 2)");
      dump(coarseA22_, "coarseA22.m");
    }

    if (!reuse && !implicitTranspose_) {
      Dk_1_T_ = Utilities::Transpose(*Dk_1_);
      if (!P22_.is_null())
        R22_ = Utilities::Transpose(*P22_);
    }

    if (!A22_.is_null()) {
      VerbLevel verbosityLevel = VerboseObject::GetDefaultVerbLevel();
      std::string label("A22");
      if (!P22_.is_null()) {
        precList22_.sublist("level 1 user data").set("A", coarseA22_);
        precList22_.sublist("level 1 user data").set("P", P22_);
        if (!implicitTranspose_)
          precList22_.sublist("level 1 user data").set("R", R22_);
        precList22_.sublist("level 1 user data").set("Nullspace", CoarseNullspace22_);
        precList22_.sublist("level 1 user data").set("Coordinates", Coords22_);
        // A22 is singular, we want to coarsen at least once.
        // So we make sure coarseA22 is not just ignored.
        int maxCoarseSize = precList22_.get("coarse: max size", MasterList::getDefault<int>("coarse: max size"));
        int numRows       = Teuchos::as<int>(coarseA22_->getGlobalNumRows());
        if (maxCoarseSize > numRows)
          precList22_.set("coarse: max size", numRows);
        int maxLevels = precList22_.get("max levels", MasterList::getDefault<int>("max levels"));
        if (maxLevels < 2)
          precList22_.set("max levels", 2);
        setupSubSolve(Hierarchy22_, thyraPrecOp22_, A22_, Teuchos::null, Teuchos::null, precList22_, label, reuse, /*isSingular=*/globalNumberBoundaryUnknowns11_ == 0);
      } else
        setupSubSolve(Hierarchy22_, thyraPrecOp22_, A22_, CoarseNullspace22_, Coords22_, precList22_, label, reuse, /*isSingular=*/globalNumberBoundaryUnknowns11_ == 0);

      VerboseObject::SetDefaultVerbLevel(verbosityLevel);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Apply BCs to rows of Dk_1
  if (!reuse && !onlyBoundary22_ && applyBCsTo22_) {
    GetOStream(Runtime0) << solverName_ + "::compute(): nuking BC rows of Dk_1" << std::endl;

    Dk_1_->resumeFill();
    Scalar replaceWith = (Dk_1_->getRowMap()->lib() == Xpetra::UseEpetra) ? Teuchos::ScalarTraits<SC>::eps() : Teuchos::ScalarTraits<SC>::zero();
    Utilities::ZeroDirichletRows(Dk_1_, BCrows11_, replaceWith);
    Dk_1_->fillComplete(Dk_1_->getDomainMap(), Dk_1_->getRangeMap());
    dump(Dk_1_, "Dk_1_nuked.m");
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Set up the smoother on the finest level
  setFineLevelSmoother11();

  if (!reuse) {
    if (!ImporterCoarse11_.is_null()) {
      RCP<const Import> ImporterP11 = ImportFactory::Build(ImporterCoarse11_->getTargetMap(), P11_->getColMap());
      rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix()->replaceDomainMapAndImporter(ImporterCoarse11_->getTargetMap(), ImporterP11);
    }

    if (!Importer22_.is_null()) {
      if (enable_reuse_) {
        DorigDomainMap_ = Dk_1_->getDomainMap();
        DorigImporter_  = rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_)->getCrsMatrix()->getCrsGraph()->getImporter();
      }
      RCP<const Import> ImporterD = ImportFactory::Build(Importer22_->getTargetMap(), Dk_1_->getColMap());
      rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_)->getCrsMatrix()->replaceDomainMapAndImporter(Importer22_->getTargetMap(), ImporterD);
    }

#ifdef HAVE_MUELU_TPETRA
    if ((!Dk_1_T_.is_null()) &&
        (!R11_.is_null()) &&
        (!rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_T_)->getCrsMatrix()->getCrsGraph()->getImporter().is_null()) &&
        (!rcp_dynamic_cast<CrsMatrixWrap>(R11_)->getCrsMatrix()->getCrsGraph()->getImporter().is_null()) &&
        (Dk_1_T_->getColMap()->lib() == Xpetra::UseTpetra) &&
        (R11_->getColMap()->lib() == Xpetra::UseTpetra))
      Dk_1_T_R11_colMapsMatch_ = Dk_1_T_->getColMap()->isSameAs(*R11_->getColMap());
    else
#endif
      Dk_1_T_R11_colMapsMatch_ = false;
    if (Dk_1_T_R11_colMapsMatch_)
      GetOStream(Runtime0) << solverName_ + "::compute(): Dk_1_T and R11 have matching colMaps" << std::endl;

    // Allocate MultiVectors for solve
    allocateMemory(1);

    // apply matvec params
    if (parameterList_.isSublist("matvec params")) {
      RCP<ParameterList> matvecParams = rcpFromRef(parameterList_.sublist("matvec params"));
      Maxwell_Utils<SC, LO, GO, NO>::setMatvecParams(*SM_Matrix_, matvecParams);
      Maxwell_Utils<SC, LO, GO, NO>::setMatvecParams(*Dk_1_, matvecParams);
      Maxwell_Utils<SC, LO, GO, NO>::setMatvecParams(*P11_, matvecParams);
      if (!Dk_1_T_.is_null()) Maxwell_Utils<SC, LO, GO, NO>::setMatvecParams(*Dk_1_T_, matvecParams);
      if (!R11_.is_null()) Maxwell_Utils<SC, LO, GO, NO>::setMatvecParams(*R11_, matvecParams);
      if (!ImporterCoarse11_.is_null()) ImporterCoarse11_->setDistributorParameters(matvecParams);
      if (!Importer22_.is_null()) Importer22_->setDistributorParameters(matvecParams);
    }
    if (!ImporterCoarse11_.is_null() && parameterList_.isSublist("refmaxwell: ImporterCoarse11 params")) {
      RCP<ParameterList> importerParams = rcpFromRef(parameterList_.sublist("refmaxwell: ImporterCoarse11 params"));
      ImporterCoarse11_->setDistributorParameters(importerParams);
    }
    if (!Importer22_.is_null() && parameterList_.isSublist("refmaxwell: Importer22 params")) {
      RCP<ParameterList> importerParams = rcpFromRef(parameterList_.sublist("refmaxwell: Importer22 params"));
      Importer22_->setDistributorParameters(importerParams);
    }
  }

  describe(GetOStream(Runtime0));

#ifdef HAVE_MUELU_CUDA
  if (parameterList_.get<bool>("refmaxwell: cuda profile setup", false)) cudaProfilerStop();
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    determineSubHierarchyCommSizes(bool &doRebalancing, int &rebalanceStriding, int &numProcsCoarseA11, int &numProcsA22) {
  doRebalancing     = parameterList_.get<bool>("refmaxwell: subsolves on subcommunicators");
  rebalanceStriding = parameterList_.get<int>("refmaxwell: subsolves striding", -1);
  int numProcs      = SM_Matrix_->getDomainMap()->getComm()->getSize();
  if (numProcs == 1) {
    doRebalancing = false;
    return;
  }

#ifdef HAVE_MPI
  if (doRebalancing) {
    {
      // decide on number of ranks for coarse (1, 1) problem

      Level level;
      level.SetFactoryManager(null);
      level.SetLevelID(0);
      level.Set("A", coarseA11_);

      auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
      ParameterList repartheurParams;
      repartheurParams.set("repartition: start level", 0);
      // Setting min == target on purpose.
      int defaultTargetRows = 10000;
      repartheurParams.set("repartition: min rows per proc", precList11_.get<int>("repartition: target rows per proc", defaultTargetRows));
      repartheurParams.set("repartition: target rows per proc", precList11_.get<int>("repartition: target rows per proc", defaultTargetRows));
      repartheurParams.set("repartition: min rows per thread", precList11_.get<int>("repartition: target rows per thread", defaultTargetRows));
      repartheurParams.set("repartition: target rows per thread", precList11_.get<int>("repartition: target rows per thread", defaultTargetRows));
      repartheurParams.set("repartition: max imbalance", precList11_.get<double>("repartition: max imbalance", 1.1));
      repartheurFactory->SetParameterList(repartheurParams);

      level.Request("number of partitions", repartheurFactory.get());
      repartheurFactory->Build(level);
      numProcsCoarseA11 = level.Get<int>("number of partitions", repartheurFactory.get());
      numProcsCoarseA11 = std::min(numProcsCoarseA11, numProcs);
    }

    {
      // decide on number of ranks for (2, 2) problem

      Level level;
      level.SetFactoryManager(null);
      level.SetLevelID(0);

      level.Set("Map", Dk_1_->getDomainMap());

      auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
      ParameterList repartheurParams;
      repartheurParams.set("repartition: start level", 0);
      repartheurParams.set("repartition: use map", true);
      // Setting min == target on purpose.
      int defaultTargetRows = 10000;
      repartheurParams.set("repartition: min rows per proc", precList22_.get<int>("repartition: target rows per proc", defaultTargetRows));
      repartheurParams.set("repartition: target rows per proc", precList22_.get<int>("repartition: target rows per proc", defaultTargetRows));
      repartheurParams.set("repartition: min rows per thread", precList22_.get<int>("repartition: target rows per thread", defaultTargetRows));
      repartheurParams.set("repartition: target rows per thread", precList22_.get<int>("repartition: target rows per thread", defaultTargetRows));
      // repartheurParams.set("repartition: max imbalance",        precList22_.get<double>("repartition: max imbalance", 1.1));
      repartheurFactory->SetParameterList(repartheurParams);

      level.Request("number of partitions", repartheurFactory.get());
      repartheurFactory->Build(level);
      numProcsA22 = level.Get<int>("number of partitions", repartheurFactory.get());
      numProcsA22 = std::min(numProcsA22, numProcs);
    }

    if (rebalanceStriding >= 1) {
      TEUCHOS_ASSERT(rebalanceStriding * numProcsCoarseA11 <= numProcs);
      TEUCHOS_ASSERT(rebalanceStriding * numProcsA22 <= numProcs);
      if (rebalanceStriding * (numProcsCoarseA11 + numProcsA22) > numProcs) {
        GetOStream(Warnings0) << solverName_ + "::compute(): Disabling striding = " << rebalanceStriding << ", since coarseA11 needs " << numProcsCoarseA11
                              << " procs and A22 needs " << numProcsA22 << " procs." << std::endl;
        rebalanceStriding = -1;
      }
      int lclBadMatrixDistribution = (coarseA11_->getLocalNumEntries() == 0) || (Dk_1_->getDomainMap()->getLocalNumElements() == 0);
      int gblBadMatrixDistribution = false;
      MueLu_maxAll(SM_Matrix_->getDomainMap()->getComm(), lclBadMatrixDistribution, gblBadMatrixDistribution);
      if (gblBadMatrixDistribution) {
        GetOStream(Warnings0) << solverName_ + "::compute(): Disabling striding = " << rebalanceStriding << ", since coarseA11 has no entries on at least one rank or Dk_1's domain map has no entries on at least one rank." << std::endl;
        rebalanceStriding = -1;
      }
    }

    if ((numProcsCoarseA11 < 0) || (numProcsA22 < 0) || (numProcsCoarseA11 + numProcsA22 > numProcs)) {
      GetOStream(Warnings0) << solverName_ + "::compute(): Disabling rebalancing of subsolves, since partition heuristic resulted "
                            << "in undesirable number of partitions: " << numProcsCoarseA11 << ", " << numProcsA22 << std::endl;
      doRebalancing = false;
    }
  }
#else
  doRebalancing = false;
#endif  // HAVE_MPI
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    buildAddon(const int spaceNumber) {
  if (spaceNumber == 0)
    return Teuchos::null;

  std::string timerLabel;
  if (spaceNumber == spaceNumber_) {
    if (skipFirst11Level_)
      timerLabel = "Build coarse addon matrix 11";
    else
      timerLabel = "Build addon matrix 11";
  } else
    timerLabel = "Build addon matrix 22";

  RCP<Teuchos::TimeMonitor> tmAddon = getTimer(timerLabel);

  RCP<Matrix> addon;
  RCP<Matrix> Z;
  RCP<Matrix> lumpedInverse;
  if (spaceNumber == spaceNumber_) {
    // catch a failure
    TEUCHOS_TEST_FOR_EXCEPTION(invMk_1_invBeta_ == Teuchos::null, std::invalid_argument,
                               solverName_ +
                                   "::buildCoarse11Matrix(): Inverse of "
                                   "lumped mass matrix required for add-on (i.e. invMk_1_invBeta_ is null)");
    lumpedInverse = invMk_1_invBeta_;

    if (skipFirst11Level_) {
      // construct Zaux = M1 P11
      RCP<Matrix> Zaux;
      Zaux = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Mk_one_, false, *P11_, false, Zaux, GetOStream(Runtime0), true, true);
      // construct Z = D* M1 P11 = D^T Zaux
      Z = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_1_, true, *Zaux, false, Z, GetOStream(Runtime0), true, true);
    } else {
      // construct Z = D* M1 P11 = D^T Zaux
      Z = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_1_, true, *Mk_one_, false, Z, GetOStream(Runtime0), true, true);
    }

  } else if (spaceNumber == spaceNumber_ - 1) {
    // catch a failure
    TEUCHOS_TEST_FOR_EXCEPTION(invMk_2_invAlpha_ == Teuchos::null, std::invalid_argument,
                               solverName_ +
                                   "::buildCoarse11Matrix(): Inverse of "
                                   "lumped mass matrix required for add-on (i.e. invMk_2_invAlpha_ is null)");
    lumpedInverse = invMk_2_invAlpha_;

    // construct Z = Dk_2^T Mk_1_one
    Z = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_2_, true, *Mk_1_one_, false, Z, GetOStream(Runtime0), true, true);
  }

  // construct Z^T lumpedInverse Z
  if (lumpedInverse->getGlobalMaxNumRowEntries() <= 1) {
    // We assume that if lumpedInverse has at most one entry per row then
    // these are all diagonal entries.
    RCP<Vector> diag = VectorFactory::Build(lumpedInverse->getRowMap());
    lumpedInverse->getLocalDiagCopy(*diag);
    {
      ArrayRCP<Scalar> diagVals = diag->getDataNonConst(0);
      for (size_t j = 0; j < diag->getMap()->getLocalNumElements(); j++) {
        diagVals[j] = Teuchos::ScalarTraits<Scalar>::squareroot(diagVals[j]);
      }
    }
    if (Z->getRowMap()->isSameAs(*(diag->getMap())))
      Z->leftScale(*diag);
    else {
      RCP<Import> importer = ImportFactory::Build(diag->getMap(), Z->getRowMap());
      RCP<Vector> diag2    = VectorFactory::Build(Z->getRowMap());
      diag2->doImport(*diag, *importer, Xpetra::INSERT);
      Z->leftScale(*diag2);
    }
    addon = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Z, true, *Z, false, addon, GetOStream(Runtime0), true, true);
  } else if (parameterList_.get<bool>("rap: triple product", false) == false) {
    RCP<Matrix> C2;
    // construct C2 = lumpedInverse Z
    C2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*lumpedInverse, false, *Z, false, C2, GetOStream(Runtime0), true, true);
    // construct Matrix2 = Z* M0inv Z = Z* C2
    addon = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Z, true, *C2, false, addon, GetOStream(Runtime0), true, true);
  } else {
    addon = MatrixFactory::Build(Z->getDomainMap());
    // construct Matrix2 = Z^T lumpedInverse Z
    Xpetra::TripleMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
        MultiplyRAP(*Z, true, *lumpedInverse, false, *Z, false, *addon, true, true);
  }
  return addon;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildCoarse11Matrix() {
  RCP<Teuchos::TimeMonitor> tm = getTimer("Build coarse (1,1) matrix");

  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  // coarse matrix for P11* (M1 + D1* M2 D1) P11
  RCP<Matrix> temp;
  temp = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*SM_Matrix_, false, *P11_, false, temp, GetOStream(Runtime0), true, true);
  if (ImporterCoarse11_.is_null())
    coarseA11_ = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P11_, true, *temp, false, coarseA11_, GetOStream(Runtime0), true, true);
  else {
    RCP<Matrix> temp2;
    temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P11_, true, *temp, false, temp2, GetOStream(Runtime0), true, true);

    RCP<const Map> map = ImporterCoarse11_->getTargetMap()->removeEmptyProcesses();
    temp2->removeEmptyProcessesInPlace(map);
    if (!temp2.is_null() && temp2->getRowMap().is_null())
      temp2 = Teuchos::null;
    coarseA11_ = temp2;
  }

  if (!disable_addon_) {
    RCP<Matrix> addon;

    if (!coarseA11_.is_null() && Addon11_.is_null()) {
      addon = buildAddon(spaceNumber_);
      // Should we keep the addon for next setup?
      if (enable_reuse_)
        Addon11_ = addon;
    } else
      addon = Addon11_;

    if (!coarseA11_.is_null()) {
      // add matrices together
      RCP<Matrix> newCoarseA11;
      Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*coarseA11_, false, one, *addon, false, one, newCoarseA11, GetOStream(Runtime0));
      newCoarseA11->fillComplete();
      coarseA11_ = newCoarseA11;
    }
  }

  if (!coarseA11_.is_null() && !skipFirst11Level_) {
    ArrayRCP<bool> coarseA11BCrows;
    coarseA11BCrows.resize(coarseA11_->getRowMap()->getLocalNumElements());
    for (size_t i = 0; i < BCdomain22_.size(); i++)
      for (size_t k = 0; k < dim_; k++)
        coarseA11BCrows[i * dim_ + k] = BCdomain22_(i);
    magnitudeType rowSumTol = parameterList_.get<double>("refmaxwell: row sum drop tol (1,1)");
    if (rowSumTol > 0.)
      Utilities::ApplyRowSumCriterion(*coarseA11_, rowSumTol, coarseA11BCrows);
    if (applyBCsToCoarse11_)
      Utilities::ApplyOAZToMatrixRows(coarseA11_, coarseA11BCrows);
  }

  if (!coarseA11_.is_null()) {
    // If we already applied BCs to A_nodal, we likely do not need
    // to fix up coarseA11.
    // If we did not apply BCs to A_nodal, we now need to correct
    // the zero diagonals of coarseA11, since we did nuke the nullspace.

    bool fixZeroDiagonal = !applyBCsToAnodal_;
    if (precList11_.isParameter("rap: fix zero diagonals"))
      fixZeroDiagonal = precList11_.get<bool>("rap: fix zero diagonals");

    if (fixZeroDiagonal) {
      magnitudeType threshold = 1e-16;
      Scalar replacement      = 1.0;
      if (precList11_.isType<magnitudeType>("rap: fix zero diagonals threshold"))
        threshold = precList11_.get<magnitudeType>("rap: fix zero diagonals threshold");
      else if (precList11_.isType<double>("rap: fix zero diagonals threshold"))
        threshold = Teuchos::as<magnitudeType>(precList11_.get<double>("rap: fix zero diagonals threshold"));
      if (precList11_.isType<double>("rap: fix zero diagonals replacement"))
        replacement = Teuchos::as<Scalar>(precList11_.get<double>("rap: fix zero diagonals replacement"));
      Xpetra::MatrixUtils<SC, LO, GO, NO>::CheckRepairMainDiagonal(coarseA11_, true, GetOStream(Warnings1), threshold, replacement);
    }

    // Set block size
    coarseA11_->SetFixedBlockSize(dim_);
    if (skipFirst11Level_)
      coarseA11_->setObjectLabel(solverName_ + " coarse (1,1)");
    else
      coarseA11_->setObjectLabel(solverName_ + " (1,1)");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    rebalanceCoarse11Matrix(const int rebalanceStriding, const int numProcsCoarseA11) {
#ifdef HAVE_MPI
  // rebalance coarseA11
  RCP<Teuchos::TimeMonitor> tm = getTimer("Rebalance coarseA11");

  Level fineLevel, coarseLevel;
  fineLevel.SetFactoryManager(null);
  coarseLevel.SetFactoryManager(null);
  coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
  fineLevel.SetLevelID(0);
  coarseLevel.SetLevelID(1);
  coarseLevel.Set("A", coarseA11_);
  coarseLevel.Set("P", P11_);
  coarseLevel.Set("Coordinates", CoordsCoarse11_);
  if (!NullspaceCoarse11_.is_null())
    coarseLevel.Set("Nullspace", NullspaceCoarse11_);
  coarseLevel.Set("number of partitions", numProcsCoarseA11);
  coarseLevel.Set("repartition: heuristic target rows per process", 1000);

  coarseLevel.setlib(coarseA11_->getDomainMap()->lib());
  fineLevel.setlib(coarseA11_->getDomainMap()->lib());
  coarseLevel.setObjectLabel(solverName_ + " coarse (1,1)");
  fineLevel.setObjectLabel(solverName_ + " coarse (1,1)");

  std::string partName = precList11_.get<std::string>("repartition: partitioner", "zoltan2");
  RCP<Factory> partitioner;
  if (partName == "zoltan") {
#ifdef HAVE_MUELU_ZOLTAN
    partitioner = rcp(new ZoltanInterface());
    // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
    // partitioner->SetFactory("number of partitions", repartheurFactory);
#else
    throw Exceptions::RuntimeError("Zoltan interface is not available");
#endif
  } else if (partName == "zoltan2") {
#ifdef HAVE_MUELU_ZOLTAN2
    partitioner = rcp(new Zoltan2Interface());
    ParameterList partParams;
    RCP<const ParameterList> partpartParams = rcp(new ParameterList(precList11_.sublist("repartition: params", false)));
    partParams.set("ParameterList", partpartParams);
    partitioner->SetParameterList(partParams);
    // partitioner->SetFactory("number of partitions", repartheurFactory);
#else
    throw Exceptions::RuntimeError("Zoltan2 interface is not available");
#endif
  }

  auto repartFactory = rcp(new RepartitionFactory());
  ParameterList repartParams;
  repartParams.set("repartition: print partition distribution", precList11_.get<bool>("repartition: print partition distribution", false));
  repartParams.set("repartition: remap parts", precList11_.get<bool>("repartition: remap parts", true));
  if (rebalanceStriding >= 1) {
    bool acceptPart = (SM_Matrix_->getDomainMap()->getComm()->getRank() % rebalanceStriding) == 0;
    if (SM_Matrix_->getDomainMap()->getComm()->getRank() >= numProcsCoarseA11 * rebalanceStriding)
      acceptPart = false;
    repartParams.set("repartition: remap accept partition", acceptPart);
  }
  repartFactory->SetParameterList(repartParams);
  // repartFactory->SetFactory("number of partitions", repartheurFactory);
  repartFactory->SetFactory("Partition", partitioner);

  auto newP = rcp(new RebalanceTransferFactory());
  ParameterList newPparams;
  newPparams.set("type", "Interpolation");
  newPparams.set("repartition: rebalance P and R", precList11_.get<bool>("repartition: rebalance P and R", false));
  newPparams.set("repartition: use subcommunicators", true);
  newPparams.set("repartition: rebalance Nullspace", !NullspaceCoarse11_.is_null());
  newP->SetFactory("Coordinates", NoFactory::getRCP());
  if (!NullspaceCoarse11_.is_null())
    newP->SetFactory("Nullspace", NoFactory::getRCP());
  newP->SetParameterList(newPparams);
  newP->SetFactory("Importer", repartFactory);

  auto newA = rcp(new RebalanceAcFactory());
  ParameterList rebAcParams;
  rebAcParams.set("repartition: use subcommunicators", true);
  newA->SetParameterList(rebAcParams);
  newA->SetFactory("Importer", repartFactory);

  coarseLevel.Request("P", newP.get());
  coarseLevel.Request("Importer", repartFactory.get());
  coarseLevel.Request("A", newA.get());
  coarseLevel.Request("Coordinates", newP.get());
  if (!NullspaceCoarse11_.is_null())
    coarseLevel.Request("Nullspace", newP.get());
  repartFactory->Build(coarseLevel);

  if (!precList11_.get<bool>("repartition: rebalance P and R", false))
    ImporterCoarse11_ = coarseLevel.Get<RCP<const Import> >("Importer", repartFactory.get());
  P11_            = coarseLevel.Get<RCP<Matrix> >("P", newP.get());
  coarseA11_      = coarseLevel.Get<RCP<Matrix> >("A", newA.get());
  CoordsCoarse11_ = coarseLevel.Get<RCP<RealValuedMultiVector> >("Coordinates", newP.get());
  if (!NullspaceCoarse11_.is_null())
    NullspaceCoarse11_ = coarseLevel.Get<RCP<MultiVector> >("Nullspace", newP.get());

  if (!coarseA11_.is_null()) {
    // Set block size
    coarseA11_->SetFixedBlockSize(dim_);
    if (skipFirst11Level_)
      coarseA11_->setObjectLabel(solverName_ + " coarse (1,1)");
    else
      coarseA11_->setObjectLabel(solverName_ + " (1,1)");
  }

  coarseA11_AP_reuse_data_  = Teuchos::null;
  coarseA11_RAP_reuse_data_ = Teuchos::null;

  if (!disable_addon_ && enable_reuse_) {
    // Rebalance the addon for next setup
    RCP<const Import> ImporterCoarse11 = coarseLevel.Get<RCP<const Import> >("Importer", repartFactory.get());
    RCP<const Map> targetMap           = ImporterCoarse11->getTargetMap();
    ParameterList XpetraList;
    XpetraList.set("Restrict Communicator", true);
    Addon11_ = MatrixFactory::Build(Addon11_, *ImporterCoarse11, *ImporterCoarse11, targetMap, targetMap, rcp(&XpetraList, false));
  }
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::build22Matrix(const bool reuse, const bool doRebalancing, const int rebalanceStriding, const int numProcsA22) {
  if (!reuse) {  // build fine grid operator for (2,2)-block, Dk_1^T SM Dk_1  (aka TMT)
    RCP<Teuchos::TimeMonitor> tm = getTimer("Build A22");

    Level fineLevel, coarseLevel;
    fineLevel.SetFactoryManager(null);
    coarseLevel.SetFactoryManager(null);
    coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
    fineLevel.SetLevelID(0);
    coarseLevel.SetLevelID(1);
    fineLevel.Set("A", SM_Matrix_);
    coarseLevel.Set("P", Dk_1_);
    coarseLevel.Set("Coordinates", Coords22_);

    coarseLevel.setlib(SM_Matrix_->getDomainMap()->lib());
    fineLevel.setlib(SM_Matrix_->getDomainMap()->lib());
    coarseLevel.setObjectLabel(solverName_ + " (2,2)");
    fineLevel.setObjectLabel(solverName_ + " (2,2)");

    RCP<RAPFactory> rapFact = rcp(new RAPFactory());
    ParameterList rapList   = *(rapFact->GetValidParameterList());
    rapList.set("transpose: use implicit", true);
    rapList.set("rap: fix zero diagonals", parameterList_.get<bool>("rap: fix zero diagonals", true));
    rapList.set("rap: fix zero diagonals threshold", parameterList_.get<double>("rap: fix zero diagonals threshold", Teuchos::ScalarTraits<double>::eps()));
    rapList.set("rap: triple product", parameterList_.get<bool>("rap: triple product", false));
    rapFact->SetParameterList(rapList);

    if (!A22_AP_reuse_data_.is_null()) {
      coarseLevel.AddKeepFlag("AP reuse data", rapFact.get());
      coarseLevel.Set<Teuchos::RCP<Teuchos::ParameterList> >("AP reuse data", A22_AP_reuse_data_, rapFact.get());
    }
    if (!A22_RAP_reuse_data_.is_null()) {
      coarseLevel.AddKeepFlag("RAP reuse data", rapFact.get());
      coarseLevel.Set<Teuchos::RCP<Teuchos::ParameterList> >("RAP reuse data", A22_RAP_reuse_data_, rapFact.get());
    }

#ifdef HAVE_MPI
    if (doRebalancing) {
      coarseLevel.Set("number of partitions", numProcsA22);
      coarseLevel.Set("repartition: heuristic target rows per process", 1000);

      std::string partName = precList22_.get<std::string>("repartition: partitioner", "zoltan2");
      RCP<Factory> partitioner;
      if (partName == "zoltan") {
#ifdef HAVE_MUELU_ZOLTAN
        partitioner = rcp(new ZoltanInterface());
        partitioner->SetFactory("A", rapFact);
        // partitioner->SetFactory("number of partitions", repartheurFactory);
        // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
#else
        throw Exceptions::RuntimeError("Zoltan interface is not available");
#endif
      } else if (partName == "zoltan2") {
#ifdef HAVE_MUELU_ZOLTAN2
        partitioner = rcp(new Zoltan2Interface());
        ParameterList partParams;
        RCP<const ParameterList> partpartParams = rcp(new ParameterList(precList22_.sublist("repartition: params", false)));
        partParams.set("ParameterList", partpartParams);
        partitioner->SetParameterList(partParams);
        partitioner->SetFactory("A", rapFact);
        // partitioner->SetFactory("number of partitions", repartheurFactory);
#else
        throw Exceptions::RuntimeError("Zoltan2 interface is not available");
#endif
      }

      auto repartFactory = rcp(new RepartitionFactory());
      ParameterList repartParams;
      repartParams.set("repartition: print partition distribution", precList22_.get<bool>("repartition: print partition distribution", false));
      repartParams.set("repartition: remap parts", precList22_.get<bool>("repartition: remap parts", true));
      if (rebalanceStriding >= 1) {
        bool acceptPart = ((SM_Matrix_->getDomainMap()->getComm()->getSize() - 1 - SM_Matrix_->getDomainMap()->getComm()->getRank()) % rebalanceStriding) == 0;
        if (SM_Matrix_->getDomainMap()->getComm()->getSize() - 1 - SM_Matrix_->getDomainMap()->getComm()->getRank() >= numProcsA22 * rebalanceStriding)
          acceptPart = false;
        if (acceptPart)
          TEUCHOS_ASSERT(coarseA11_.is_null());
        repartParams.set("repartition: remap accept partition", acceptPart);
      } else
        repartParams.set("repartition: remap accept partition", coarseA11_.is_null());
      repartFactory->SetParameterList(repartParams);
      repartFactory->SetFactory("A", rapFact);
      // repartFactory->SetFactory("number of partitions", repartheurFactory);
      repartFactory->SetFactory("Partition", partitioner);

      auto newP = rcp(new RebalanceTransferFactory());
      ParameterList newPparams;
      newPparams.set("type", "Interpolation");
      newPparams.set("repartition: rebalance P and R", precList22_.get<bool>("repartition: rebalance P and R", false));
      newPparams.set("repartition: use subcommunicators", true);
      newPparams.set("repartition: rebalance Nullspace", false);
      newP->SetFactory("Coordinates", NoFactory::getRCP());
      newP->SetParameterList(newPparams);
      newP->SetFactory("Importer", repartFactory);

      auto newA = rcp(new RebalanceAcFactory());
      ParameterList rebAcParams;
      rebAcParams.set("repartition: use subcommunicators", true);
      newA->SetParameterList(rebAcParams);
      newA->SetFactory("A", rapFact);
      newA->SetFactory("Importer", repartFactory);

      coarseLevel.Request("P", newP.get());
      coarseLevel.Request("Importer", repartFactory.get());
      coarseLevel.Request("A", newA.get());
      coarseLevel.Request("Coordinates", newP.get());
      rapFact->Build(fineLevel, coarseLevel);
      repartFactory->Build(coarseLevel);

      if (!precList22_.get<bool>("repartition: rebalance P and R", false))
        Importer22_ = coarseLevel.Get<RCP<const Import> >("Importer", repartFactory.get());
      Dk_1_     = coarseLevel.Get<RCP<Matrix> >("P", newP.get());
      A22_      = coarseLevel.Get<RCP<Matrix> >("A", newA.get());
      Coords22_ = coarseLevel.Get<RCP<RealValuedMultiVector> >("Coordinates", newP.get());

      if (!P22_.is_null()) {
        // Todo
      }

    } else
#endif  // HAVE_MPI
    {
      coarseLevel.Request("A", rapFact.get());
      if (enable_reuse_) {
        coarseLevel.Request("AP reuse data", rapFact.get());
        coarseLevel.Request("RAP reuse data", rapFact.get());
      }

      A22_ = coarseLevel.Get<RCP<Matrix> >("A", rapFact.get());

      if (enable_reuse_) {
        if (coarseLevel.IsAvailable("AP reuse data", rapFact.get()))
          A22_AP_reuse_data_ = coarseLevel.Get<RCP<ParameterList> >("AP reuse data", rapFact.get());
        if (coarseLevel.IsAvailable("RAP reuse data", rapFact.get()))
          A22_RAP_reuse_data_ = coarseLevel.Get<RCP<ParameterList> >("RAP reuse data", rapFact.get());
      }
    }
  } else {
    RCP<Teuchos::TimeMonitor> tm = getTimer("Build A22");
    if (Importer22_.is_null()) {
      RCP<Matrix> temp;
      temp = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*SM_Matrix_, false, *Dk_1_, false, temp, GetOStream(Runtime0), true, true);
      if (!implicitTranspose_)
        A22_ = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_1_T_, false, *temp, false, A22_, GetOStream(Runtime0), true, true);
      else
        A22_ = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_1_, true, *temp, false, A22_, GetOStream(Runtime0), true, true);
    } else {
      // we replaced domain map and importer on D, reverse that
      RCP<const Import> Dimporter = rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_)->getCrsMatrix()->getCrsGraph()->getImporter();
      rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_)->getCrsMatrix()->replaceDomainMapAndImporter(DorigDomainMap_, DorigImporter_);

      RCP<Matrix> temp, temp2;
      temp = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*SM_Matrix_, false, *Dk_1_, false, temp, GetOStream(Runtime0), true, true);
      if (!implicitTranspose_)
        temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_1_T_, false, *temp, false, temp2, GetOStream(Runtime0), true, true);
      else
        temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Dk_1_, true, *temp, false, temp2, GetOStream(Runtime0), true, true);

      // and back again
      rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_)->getCrsMatrix()->replaceDomainMapAndImporter(Importer22_->getTargetMap(), Dimporter);

      ParameterList XpetraList;
      XpetraList.set("Restrict Communicator", true);
      XpetraList.set("Timer Label", "MueLu::RebalanceA22");
      RCP<const Map> targetMap = Importer22_->getTargetMap();
      A22_                     = MatrixFactory::Build(temp2, *Importer22_, *Importer22_, targetMap, targetMap, rcp(&XpetraList, false));
    }
  }

  if (not A22_.is_null() and not disable_addon_22_ and spaceNumber_ > 1) {
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    RCP<Matrix> addon22 = buildAddon(spaceNumber_ - 1);

    // add matrices together
    RCP<Matrix> newA22;
    Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*A22_, false, one, *addon22, false, one, newA22, GetOStream(Runtime0));
    newA22->fillComplete();
    A22_ = newA22;
  }

  if (!A22_.is_null()) {
    dump(A22_, "A22.m");
    A22_->setObjectLabel(solverName_ + " (2,2)");
    // Set block size
    if (spaceNumber_ - 1 == 0)
      A22_->SetFixedBlockSize(1);
    else
      A22_->SetFixedBlockSize(dim_);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setFineLevelSmoother11() {
  Level level;
  RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
  level.SetFactoryManager(factoryHandler);
  level.SetLevelID(0);
  level.setObjectLabel(solverName_ + " (1,1)");
  level.Set("A", SM_Matrix_);
  level.setlib(SM_Matrix_->getDomainMap()->lib());
  // For Hiptmair
  level.Set("NodeMatrix", A22_);
  level.Set("D0", Dk_1_);

  if ((parameterList_.get<std::string>("smoother: pre type") != "NONE") && (parameterList_.get<std::string>("smoother: post type") != "NONE")) {
    std::string preSmootherType  = parameterList_.get<std::string>("smoother: pre type");
    std::string postSmootherType = parameterList_.get<std::string>("smoother: post type");

    ParameterList preSmootherList, postSmootherList;
    if (parameterList_.isSublist("smoother: pre params"))
      preSmootherList = parameterList_.sublist("smoother: pre params");
    if (parameterList_.isSublist("smoother: post params"))
      postSmootherList = parameterList_.sublist("smoother: post params");

    RCP<SmootherPrototype> preSmootherPrototype  = rcp(new TrilinosSmoother(preSmootherType, preSmootherList));
    RCP<SmootherPrototype> postSmootherPrototype = rcp(new TrilinosSmoother(postSmootherType, postSmootherList));
    RCP<SmootherFactory> smootherFact            = rcp(new SmootherFactory(preSmootherPrototype, postSmootherPrototype));

    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());
    if (enable_reuse_) {
      ParameterList smootherFactoryParams;
      smootherFactoryParams.set("keep smoother data", true);
      smootherFact->SetParameterList(smootherFactoryParams);
      level.Request("PreSmoother data", smootherFact.get());
      level.Request("PostSmoother data", smootherFact.get());
      if (!PreSmootherData11_.is_null())
        level.Set("PreSmoother data", PreSmootherData11_, smootherFact.get());
      if (!PostSmootherData11_.is_null())
        level.Set("PostSmoother data", PostSmootherData11_, smootherFact.get());
    }
    smootherFact->Build(level);
    PreSmoother11_  = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());
    PostSmoother11_ = level.Get<RCP<SmootherBase> >("PostSmoother", smootherFact.get());
    if (enable_reuse_) {
      PreSmootherData11_  = level.Get<RCP<SmootherPrototype> >("PreSmoother data", smootherFact.get());
      PostSmootherData11_ = level.Get<RCP<SmootherPrototype> >("PostSmoother data", smootherFact.get());
    }
  } else {
    std::string smootherType = parameterList_.get<std::string>("smoother: type");

    ParameterList smootherList;
    if (parameterList_.isSublist("smoother: params"))
      smootherList = parameterList_.sublist("smoother: params");

    RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(smootherType, smootherList));
    RCP<SmootherFactory> smootherFact        = rcp(new SmootherFactory(smootherPrototype));
    level.Request("PreSmoother", smootherFact.get());
    if (enable_reuse_) {
      ParameterList smootherFactoryParams;
      smootherFactoryParams.set("keep smoother data", true);
      smootherFact->SetParameterList(smootherFactoryParams);
      level.Request("PreSmoother data", smootherFact.get());
      if (!PreSmootherData11_.is_null())
        level.Set("PreSmoother data", PreSmootherData11_, smootherFact.get());
    }
    smootherFact->Build(level);
    PreSmoother11_  = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());
    PostSmoother11_ = PreSmoother11_;
    if (enable_reuse_)
      PreSmootherData11_ = level.Get<RCP<SmootherPrototype> >("PreSmoother data", smootherFact.get());
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::allocateMemory(int numVectors) const {
  RCP<Teuchos::TimeMonitor> tmAlloc = getTimer("Allocate MVs");

  // 11 block
  if (!R11_.is_null())
    P11res_ = MultiVectorFactory::Build(R11_->getRangeMap(), numVectors);
  else
    P11res_ = MultiVectorFactory::Build(P11_->getDomainMap(), numVectors);
  P11res_->setObjectLabel("P11res");

  if (Dk_1_T_R11_colMapsMatch_) {
    DTR11Tmp_ = MultiVectorFactory::Build(R11_->getColMap(), numVectors);
    DTR11Tmp_->setObjectLabel("DTR11Tmp");
  }
  if (!ImporterCoarse11_.is_null()) {
    P11resTmp_ = MultiVectorFactory::Build(ImporterCoarse11_->getTargetMap(), numVectors);
    P11resTmp_->setObjectLabel("P11resTmp");
    P11x_ = MultiVectorFactory::Build(ImporterCoarse11_->getTargetMap(), numVectors);
  } else
    P11x_ = MultiVectorFactory::Build(P11_->getDomainMap(), numVectors);
  P11x_->setObjectLabel("P11x");

  // 22 block
  if (!Dk_1_T_.is_null())
    Dres_ = MultiVectorFactory::Build(Dk_1_T_->getRangeMap(), numVectors);
  else
    Dres_ = MultiVectorFactory::Build(Dk_1_->getDomainMap(), numVectors);
  Dres_->setObjectLabel("Dres");

  if (!Importer22_.is_null()) {
    DresTmp_ = MultiVectorFactory::Build(Importer22_->getTargetMap(), numVectors);
    DresTmp_->setObjectLabel("DresTmp");
    Dx_ = MultiVectorFactory::Build(Importer22_->getTargetMap(), numVectors);
  } else if (!onlyBoundary22_)
    Dx_ = MultiVectorFactory::Build(A22_->getDomainMap(), numVectors);
  if (!Dx_.is_null())
    Dx_->setObjectLabel("Dx");

  if (!coarseA11_.is_null()) {
    if (!ImporterCoarse11_.is_null() && !implicitTranspose_)
      P11resSubComm_ = MultiVectorFactory::Build(P11resTmp_, Teuchos::View);
    else
      P11resSubComm_ = MultiVectorFactory::Build(P11res_, Teuchos::View);
    P11resSubComm_->replaceMap(coarseA11_->getRangeMap());
    P11resSubComm_->setObjectLabel("P11resSubComm");

    P11xSubComm_ = MultiVectorFactory::Build(P11x_, Teuchos::View);
    P11xSubComm_->replaceMap(coarseA11_->getDomainMap());
    P11xSubComm_->setObjectLabel("P11xSubComm");
  }

  if (!A22_.is_null()) {
    if (!Importer22_.is_null() && !implicitTranspose_)
      DresSubComm_ = MultiVectorFactory::Build(DresTmp_, Teuchos::View);
    else
      DresSubComm_ = MultiVectorFactory::Build(Dres_, Teuchos::View);
    DresSubComm_->replaceMap(A22_->getRangeMap());
    DresSubComm_->setObjectLabel("DresSubComm");

    DxSubComm_ = MultiVectorFactory::Build(Dx_, Teuchos::View);
    DxSubComm_->replaceMap(A22_->getDomainMap());
    DxSubComm_->setObjectLabel("DxSubComm");
  }

  residual_ = MultiVectorFactory::Build(SM_Matrix_->getDomainMap(), numVectors);
  residual_->setObjectLabel("residual");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dump(const RCP<Matrix> &A, std::string name) const {
  if (dump_matrices_ && !A.is_null()) {
    GetOStream(Runtime0) << "Dumping to " << name << std::endl;
    Xpetra::IO<SC, LO, GO, NO>::Write(name, *A);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dump(const RCP<MultiVector> &X, std::string name) const {
  if (dump_matrices_ && !X.is_null()) {
    GetOStream(Runtime0) << "Dumping to " << name << std::endl;
    Xpetra::IO<SC, LO, GO, NO>::Write(name, *X);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dumpCoords(const RCP<RealValuedMultiVector> &X, std::string name) const {
  if (dump_matrices_ && !X.is_null()) {
    GetOStream(Runtime0) << "Dumping to " << name << std::endl;
    Xpetra::IO<coordinateType, LO, GO, NO>::Write(name, *X);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dump(const Teuchos::ArrayRCP<bool> &v, std::string name) const {
  if (dump_matrices_) {
    GetOStream(Runtime0) << "Dumping to " << name << std::endl;
    std::ofstream out(name);
    for (size_t i = 0; i < Teuchos::as<size_t>(v.size()); i++)
      out << v[i] << "\n";
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dump(const Kokkos::View<bool *, typename Node::device_type> &v, std::string name) const {
  if (dump_matrices_) {
    GetOStream(Runtime0) << "Dumping to " << name << std::endl;
    std::ofstream out(name);
    auto vH = Kokkos::create_mirror_view(v);
    Kokkos::deep_copy(vH, v);
    out << "%%MatrixMarket matrix array real general\n"
        << vH.extent(0) << " 1\n";
    for (size_t i = 0; i < vH.size(); i++)
      out << vH[i] << "\n";
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Teuchos::TimeMonitor> RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTimer(std::string name, RCP<const Teuchos::Comm<int> > comm) const {
  if (IsPrint(Timings)) {
    if (!syncTimers_)
      return Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu " + solverName_ + ": " + name)));
    else {
      if (comm.is_null())
        return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu " + solverName_ + ": " + name), SM_Matrix_->getRowMap()->getComm().ptr()));
      else
        return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu " + solverName_ + ": " + name), comm.ptr()));
    }
  } else
    return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    buildNullspace(const int spaceNumber, const Kokkos::View<bool *, typename Node::device_type> &bcs, const bool applyBCs) {
  std::string spaceLabel;
  if (spaceNumber == 0)
    spaceLabel = "nodal";
  else if (spaceNumber == 1)
    spaceLabel = "edge";
  else if (spaceNumber == 2)
    spaceLabel = "face";
  else
    TEUCHOS_ASSERT(false);

  RCP<Teuchos::TimeMonitor> tm;
  if (spaceNumber > 0) {
    tm = getTimer("nullspace " + spaceLabel);
    GetOStream(Runtime0) << solverName_ + "::compute(): building " + spaceLabel + " nullspace" << std::endl;
  }

  if (spaceNumber == 0) {
    return Teuchos::null;

  } else if (spaceNumber == 1) {
    RCP<MultiVector> CoordsSC;
    CoordsSC                   = Utilities::RealValuedToScalarMultiVector(NodalCoords_);
    RCP<MultiVector> Nullspace = MultiVectorFactory::Build(D0_->getRowMap(), NodalCoords_->getNumVectors());
    D0_->apply(*CoordsSC, *Nullspace);

    bool normalize = parameterList_.get<bool>("refmaxwell: normalize nullspace", MasterList::getDefault<bool>("refmaxwell: normalize nullspace"));

    coordinateType minLen, maxLen, meanLen;
    if (IsPrint(Statistics2) || normalize) {
      // compute edge lengths
      ArrayRCP<ArrayRCP<const Scalar> > localNullspace(dim_);
      for (size_t i = 0; i < dim_; i++)
        localNullspace[i] = Nullspace->getData(i);
      coordinateType localMinLen  = Teuchos::ScalarTraits<coordinateType>::rmax();
      coordinateType localMeanLen = Teuchos::ScalarTraits<coordinateType>::zero();
      coordinateType localMaxLen  = Teuchos::ScalarTraits<coordinateType>::zero();
      for (size_t j = 0; j < Nullspace->getMap()->getLocalNumElements(); j++) {
        Scalar lenSC = Teuchos::ScalarTraits<Scalar>::zero();
        for (size_t i = 0; i < dim_; i++)
          lenSC += localNullspace[i][j] * localNullspace[i][j];
        coordinateType len = Teuchos::as<coordinateType>(Teuchos::ScalarTraits<Scalar>::real(Teuchos::ScalarTraits<Scalar>::squareroot(lenSC)));
        localMinLen        = std::min(localMinLen, len);
        localMaxLen        = std::max(localMaxLen, len);
        localMeanLen += len;
      }

      RCP<const Teuchos::Comm<int> > comm = Nullspace->getMap()->getComm();
      MueLu_minAll(comm, localMinLen, minLen);
      MueLu_sumAll(comm, localMeanLen, meanLen);
      MueLu_maxAll(comm, localMaxLen, maxLen);
      meanLen /= Nullspace->getMap()->getGlobalNumElements();
    }

    if (IsPrint(Statistics2)) {
      GetOStream(Statistics2) << "Edge length (min/mean/max): " << minLen << " / " << meanLen << " / " << maxLen << std::endl;
    }

    if (normalize) {
      // normalize the nullspace
      GetOStream(Runtime0) << solverName_ + "::compute(): normalizing nullspace" << std::endl;

      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

      Array<Scalar> normsSC(NodalCoords_->getNumVectors(), one / Teuchos::as<Scalar>(meanLen));
      Nullspace->scale(normsSC());
    }

    if (applyBCs) {
      // Nuke the BC edges in nullspace
      Utilities::ZeroDirichletRows(Nullspace, bcs);
    }
    dump(Nullspace, "nullspaceEdge.m");

    return Nullspace;

  } else if (spaceNumber == 2) {
    using ATS         = Kokkos::ArithTraits<Scalar>;
    using impl_Scalar = typename ATS::val_type;
    using impl_ATS    = Kokkos::ArithTraits<impl_Scalar>;
    using range_type  = Kokkos::RangePolicy<LO, typename NO::execution_space>;

    RCP<Matrix> facesToNodes;
    {
      RCP<Matrix> edgesToNodes = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(D0_);
      Maxwell_Utils<SC, LO, GO, NO>::thresholdedAbs(edgesToNodes, 1e-3);

      // dump(edgesToNodes, "edgesToNodes.m");

      RCP<Matrix> facesToEdges = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(Dk_1_);
      Maxwell_Utils<SC, LO, GO, NO>::thresholdedAbs(facesToEdges, 1e-3);
      facesToEdges = Maxwell_Utils<SC, LO, GO, NO>::removeExplicitZeros(facesToEdges, 1e-3, false);

      // dump(facesToEdges, "facesToEdges.m");

      facesToNodes = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*facesToEdges, false, *edgesToNodes, false, facesToNodes, GetOStream(Runtime0), true, true);
      Maxwell_Utils<SC, LO, GO, NO>::thresholdedAbs(facesToNodes, 1e-3);
      facesToNodes = Maxwell_Utils<SC, LO, GO, NO>::removeExplicitZeros(facesToNodes, 1e-3, false);
    }

    // dump(facesToNodes, "facesToNodes.m");

    RCP<RealValuedMultiVector> ghostedNodalCoordinates;
    auto importer = facesToNodes->getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      ghostedNodalCoordinates = Xpetra::MultiVectorFactory<coordinateType, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), dim_);
      ghostedNodalCoordinates->doImport(*NodalCoords_, *importer, Xpetra::INSERT);
    } else
      ghostedNodalCoordinates = NodalCoords_;

    RCP<MultiVector> Nullspace = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(facesToNodes->getRangeMap(), dim_);
    {
      auto facesToNodesLocal     = facesToNodes->getLocalMatrixDevice();
      auto localNodalCoordinates = ghostedNodalCoordinates->getDeviceLocalView(Xpetra::Access::ReadOnly);
      auto localFaceNullspace    = Nullspace->getDeviceLocalView(Xpetra::Access::ReadWrite);

      // enter values
      Kokkos::parallel_for(
          solverName_ + "::buildFaceProjection_nullspace",
          range_type(0, Nullspace->getMap()->getLocalNumElements()),
          KOKKOS_LAMBDA(const size_t f) {
            size_t n0                      = facesToNodesLocal.graph.entries(facesToNodesLocal.graph.row_map(f));
            size_t n1                      = facesToNodesLocal.graph.entries(facesToNodesLocal.graph.row_map(f) + 1);
            size_t n2                      = facesToNodesLocal.graph.entries(facesToNodesLocal.graph.row_map(f) + 2);
            impl_Scalar elementNullspace00 = localNodalCoordinates(n1, 0) - localNodalCoordinates(n0, 0);
            impl_Scalar elementNullspace10 = localNodalCoordinates(n2, 0) - localNodalCoordinates(n0, 0);
            impl_Scalar elementNullspace01 = localNodalCoordinates(n1, 1) - localNodalCoordinates(n0, 1);
            impl_Scalar elementNullspace11 = localNodalCoordinates(n2, 1) - localNodalCoordinates(n0, 1);
            impl_Scalar elementNullspace02 = localNodalCoordinates(n1, 2) - localNodalCoordinates(n0, 2);
            impl_Scalar elementNullspace12 = localNodalCoordinates(n2, 2) - localNodalCoordinates(n0, 2);

            localFaceNullspace(f, 0) = impl_ATS::magnitude(elementNullspace01 * elementNullspace12 - elementNullspace02 * elementNullspace11) / 6.0;
            localFaceNullspace(f, 1) = impl_ATS::magnitude(elementNullspace02 * elementNullspace10 - elementNullspace00 * elementNullspace12) / 6.0;
            localFaceNullspace(f, 2) = impl_ATS::magnitude(elementNullspace00 * elementNullspace11 - elementNullspace01 * elementNullspace10) / 6.0;
          });
    }

    if (applyBCs) {
      // Nuke the BC faces in nullspace
      Utilities::ZeroDirichletRows(Nullspace, bcs);
    }

    dump(Nullspace, "nullspaceFace.m");

    return Nullspace;

  } else
    TEUCHOS_ASSERT(false);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildProjection(const int spaceNumber, const RCP<MultiVector> &Nullspace) const {
  using ATS         = Kokkos::ArithTraits<Scalar>;
  using impl_Scalar = typename ATS::val_type;
  using impl_ATS    = Kokkos::ArithTraits<impl_Scalar>;
  using range_type  = Kokkos::RangePolicy<LO, typename NO::execution_space>;

  typedef typename Matrix::local_matrix_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  const impl_Scalar impl_SC_ONE  = impl_ATS::one();
  const impl_Scalar impl_SC_ZERO = impl_ATS::zero();
  const impl_Scalar impl_half    = impl_SC_ONE / (impl_SC_ONE + impl_SC_ONE);

  std::string spaceLabel;
  if (spaceNumber == 0)
    spaceLabel = "nodal";
  else if (spaceNumber == 1)
    spaceLabel = "edge";
  else if (spaceNumber == 2)
    spaceLabel = "face";
  else
    TEUCHOS_ASSERT(false);

  RCP<Teuchos::TimeMonitor> tm;
  if (spaceNumber > 0) {
    tm = getTimer("projection " + spaceLabel);
    GetOStream(Runtime0) << solverName_ + "::compute(): building " + spaceLabel + " projection" << std::endl;
  }

  RCP<Matrix> incidence;
  if (spaceNumber == 0) {
    // identity projection
    return Teuchos::null;

  } else if (spaceNumber == 1) {
    // D0 is incidence from nodes to edges
    incidence = D0_;

  } else if (spaceNumber == 2) {
    // get incidence from nodes to faces by multiplying D0 and D1

    TEUCHOS_ASSERT(spaceNumber_ == 2);

    RCP<Matrix> facesToNodes;
    {
      RCP<Matrix> edgesToNodes = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(D0_);
      Maxwell_Utils<SC, LO, GO, NO>::thresholdedAbs(edgesToNodes, 1e-10);

      dump(edgesToNodes, "edgesToNodes.m");

      RCP<Matrix> facesToEdges = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(Dk_1_);
      Maxwell_Utils<SC, LO, GO, NO>::thresholdedAbs(facesToEdges, 1e-10);
      // facesToEdges = Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(facesToEdges, 1e-2, false);

      dump(facesToEdges, "facesToEdges.m");

      facesToNodes = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*facesToEdges, false, *edgesToNodes, false, facesToNodes, GetOStream(Runtime0), true, true);
      Maxwell_Utils<SC, LO, GO, NO>::thresholdedAbs(facesToNodes, 1e-10);
      facesToNodes = Maxwell_Utils<SC, LO, GO, NO>::removeExplicitZeros(facesToNodes, 1e-2, false);
    }

    dump(facesToNodes, "facesToNodes.m");

    incidence = facesToNodes;

  } else
    TEUCHOS_ASSERT(false);

  size_t dim = dim_;

  // Create maps
  RCP<const Map> rowMap         = incidence->getRowMap();
  RCP<const Map> blockColMap    = MapFactory::Build(incidence->getColMap(), dim);
  RCP<const Map> blockDomainMap = MapFactory::Build(incidence->getDomainMap(), dim);

  auto localIncidence    = incidence->getLocalMatrixDevice();
  size_t numLocalRows    = rowMap->getLocalNumElements();
  size_t numLocalColumns = dim * incidence->getColMap()->getLocalNumElements();
  size_t nnzEstimate     = dim * localIncidence.graph.entries.size();
  lno_view_t rowptr(Kokkos::ViewAllocateWithoutInitializing("projection_rowptr_" + spaceLabel), numLocalRows + 1);
  lno_nnz_view_t colind(Kokkos::ViewAllocateWithoutInitializing("projection_colind_" + spaceLabel), nnzEstimate);
  scalar_view_t vals("projection_vals_" + spaceLabel, nnzEstimate);

  // set rowpointer
  Kokkos::parallel_for(
      solverName_ + "::buildProjection_adjustRowptr_" + spaceLabel,
      range_type(0, numLocalRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        rowptr(i) = dim * localIncidence.graph.row_map(i);
      });

  auto localNullspace = Nullspace->getDeviceLocalView(Xpetra::Access::ReadOnly);

  // set column indices and values
  magnitudeType tol = 1e-5;
  Kokkos::parallel_for(
      solverName_ + "::buildProjection_enterValues_" + spaceLabel,
      range_type(0, numLocalRows),
      KOKKOS_LAMBDA(const size_t f) {
        for (size_t jj = localIncidence.graph.row_map(f); jj < localIncidence.graph.row_map(f + 1); jj++) {
          for (size_t k = 0; k < dim; k++) {
            colind(dim * jj + k) = dim * localIncidence.graph.entries(jj) + k;
            if (impl_ATS::magnitude(localIncidence.values(jj)) > tol)
              vals(dim * jj + k) = impl_half * localNullspace(f, k);
            else
              vals(dim * jj + k) = impl_SC_ZERO;
          }
        }
      });

  // Create matrix
  typename CrsMatrix::local_matrix_type lclProjection("local projection " + spaceLabel,
                                                      numLocalRows, numLocalColumns, nnzEstimate,
                                                      vals, rowptr, colind);
  RCP<Matrix> projection = MatrixFactory::Build(lclProjection,
                                                rowMap, blockColMap,
                                                blockDomainMap, rowMap);

  return projection;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildNodalProlongator(const Teuchos::RCP<Matrix> &A_nodal,
                                                                                  Teuchos::RCP<Matrix> &P_nodal,
                                                                                  Teuchos::RCP<MultiVector> &Nullspace_nodal,
                                                                                  Teuchos::RCP<RealValuedMultiVector> &CoarseCoords_nodal) const {
  RCP<Teuchos::TimeMonitor> tm = getTimer("nodal prolongator");
  GetOStream(Runtime0) << solverName_ + "::compute(): building nodal prolongator" << std::endl;

  // build prolongator: algorithm 1 in the reference paper
  // First, build nodal unsmoothed prolongator using the matrix A_nodal

  const SC SC_ONE = Teuchos::ScalarTraits<SC>::one();

  {
    Level fineLevel, coarseLevel;
    fineLevel.SetFactoryManager(null);
    coarseLevel.SetFactoryManager(null);
    coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
    fineLevel.SetLevelID(0);
    coarseLevel.SetLevelID(1);
    fineLevel.Set("A", A_nodal);
    fineLevel.Set("Coordinates", NodalCoords_);
    fineLevel.Set("DofsPerNode", 1);
    coarseLevel.setlib(A_nodal->getDomainMap()->lib());
    fineLevel.setlib(A_nodal->getDomainMap()->lib());
    coarseLevel.setObjectLabel(A_nodal->getObjectLabel());
    fineLevel.setObjectLabel(A_nodal->getObjectLabel());

    LocalOrdinal NSdim         = 1;
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A_nodal->getRowMap(), NSdim);
    nullSpace->putScalar(SC_ONE);
    fineLevel.Set("Nullspace", nullSpace);

    std::string algo = parameterList_.get<std::string>("multigrid algorithm");

    RCP<Factory> amalgFact, dropFact, UncoupledAggFact, coarseMapFact, TentativePFact, Tfact, SaPFact;
    amalgFact        = rcp(new AmalgamationFactory());
    coarseMapFact    = rcp(new CoarseMapFactory());
    Tfact            = rcp(new CoordinatesTransferFactory());
    UncoupledAggFact = rcp(new UncoupledAggregationFactory());
    if (useKokkos_) {
      dropFact       = rcp(new CoalesceDropFactory_kokkos());
      TentativePFact = rcp(new TentativePFactory_kokkos());
    } else {
      dropFact       = rcp(new CoalesceDropFactory());
      TentativePFact = rcp(new TentativePFactory());
    }
    if (algo == "sa")
      SaPFact = rcp(new SaPFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    double dropTol           = parameterList_.get<double>("aggregation: drop tol");
    std::string dropScheme   = parameterList_.get<std::string>("aggregation: drop scheme");
    std::string distLaplAlgo = parameterList_.get<std::string>("aggregation: distance laplacian algo");
    dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
    dropFact->SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(dropScheme));
    if (!useKokkos_)
      dropFact->SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(distLaplAlgo));

    UncoupledAggFact->SetFactory("Graph", dropFact);
    int minAggSize = parameterList_.get<int>("aggregation: min agg size");
    UncoupledAggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(minAggSize));
    int maxAggSize = parameterList_.get<int>("aggregation: max agg size");
    UncoupledAggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(maxAggSize));
    bool matchMLbehavior1 = parameterList_.get<bool>("aggregation: match ML phase1");
    UncoupledAggFact->SetParameter("aggregation: match ML phase1", Teuchos::ParameterEntry(matchMLbehavior1));
    bool matchMLbehavior2a = parameterList_.get<bool>("aggregation: match ML phase2a");
    UncoupledAggFact->SetParameter("aggregation: match ML phase2a", Teuchos::ParameterEntry(matchMLbehavior2a));
    bool matchMLbehavior2b = parameterList_.get<bool>("aggregation: match ML phase2b");
    UncoupledAggFact->SetParameter("aggregation: match ML phase2b", Teuchos::ParameterEntry(matchMLbehavior2b));

    coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);

    TentativePFact->SetFactory("Aggregates", UncoupledAggFact);
    TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
    TentativePFact->SetFactory("CoarseMap", coarseMapFact);

    Tfact->SetFactory("Aggregates", UncoupledAggFact);
    Tfact->SetFactory("CoarseMap", coarseMapFact);

    if (algo == "sa") {
      SaPFact->SetFactory("P", TentativePFact);
      coarseLevel.Request("P", SaPFact.get());
    } else
      coarseLevel.Request("P", TentativePFact.get());
    coarseLevel.Request("Nullspace", TentativePFact.get());
    coarseLevel.Request("Coordinates", Tfact.get());

    RCP<AggregationExportFactory> aggExport;
    bool exportVizData = parameterList_.get<bool>("aggregation: export visualization data");
    if (exportVizData) {
      aggExport = rcp(new AggregationExportFactory());
      ParameterList aggExportParams;
      aggExportParams.set("aggregation: output filename", "aggs.vtk");
      aggExportParams.set("aggregation: output file: agg style", "Jacks");
      aggExport->SetParameterList(aggExportParams);

      aggExport->SetFactory("Aggregates", UncoupledAggFact);
      aggExport->SetFactory("UnAmalgamationInfo", amalgFact);
      fineLevel.Request("Aggregates", UncoupledAggFact.get());
      fineLevel.Request("UnAmalgamationInfo", amalgFact.get());
    }

    if (algo == "sa")
      coarseLevel.Get("P", P_nodal, SaPFact.get());
    else
      coarseLevel.Get("P", P_nodal, TentativePFact.get());
    coarseLevel.Get("Nullspace", Nullspace_nodal, TentativePFact.get());
    coarseLevel.Get("Coordinates", CoarseCoords_nodal, Tfact.get());

    if (exportVizData)
      aggExport->Build(fineLevel, coarseLevel);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildVectorNodalProlongator(const Teuchos::RCP<Matrix> &P_nodal) const {
  RCP<Teuchos::TimeMonitor> tm = getTimer("vectorial nodal prolongator");
  GetOStream(Runtime0) << solverName_ + "::compute(): building vectorial nodal prolongator" << std::endl;

  using range_type = Kokkos::RangePolicy<LO, typename NO::execution_space>;

  typedef typename Matrix::local_matrix_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  size_t dim = dim_;

  // Create the matrix object
  RCP<Map> blockRowMap    = MapFactory::Build(P_nodal->getRowMap(), dim);
  RCP<Map> blockColMap    = MapFactory::Build(P_nodal->getColMap(), dim);
  RCP<Map> blockDomainMap = MapFactory::Build(P_nodal->getDomainMap(), dim);

  // Get data out of P_nodal.
  auto localP_nodal = P_nodal->getLocalMatrixDevice();

  size_t numLocalRows    = blockRowMap->getLocalNumElements();
  size_t numLocalColumns = blockColMap->getLocalNumElements();
  size_t nnzEstimate     = dim * localP_nodal.graph.entries.size();
  lno_view_t rowptr(Kokkos::ViewAllocateWithoutInitializing("vectorPNodal_rowptr"), numLocalRows + 1);
  lno_nnz_view_t colind(Kokkos::ViewAllocateWithoutInitializing("vectorPNodal_colind"), nnzEstimate);
  scalar_view_t vals(Kokkos::ViewAllocateWithoutInitializing("vectorPNodal_vals"), nnzEstimate);

  // fill rowpointer
  Kokkos::parallel_for(
      solverName_ + "::buildVectorNodalProlongator_adjustRowptr",
      range_type(0, localP_nodal.numRows() + 1),
      KOKKOS_LAMBDA(const LocalOrdinal i) {
        if (i < localP_nodal.numRows()) {
          for (size_t k = 0; k < dim; k++) {
            rowptr(dim * i + k) = dim * localP_nodal.graph.row_map(i) + k;
          }
        } else
          rowptr(dim * localP_nodal.numRows()) = dim * localP_nodal.graph.row_map(i);
      });

  // fill column indices and values
  Kokkos::parallel_for(
      solverName_ + "::buildVectorNodalProlongator_adjustColind",
      range_type(0, localP_nodal.graph.entries.size()),
      KOKKOS_LAMBDA(const size_t jj) {
        for (size_t k = 0; k < dim; k++) {
          colind(dim * jj + k) = dim * localP_nodal.graph.entries(jj) + k;
          // vals(dim*jj+k) = localP_nodal.values(jj);
          vals(dim * jj + k) = 1.;
        }
      });

  typename CrsMatrix::local_matrix_type lclVectorNodalP("local vector nodal prolongator",
                                                        numLocalRows, numLocalColumns, nnzEstimate,
                                                        vals, rowptr, colind);
  RCP<Matrix> vectorNodalP = MatrixFactory::Build(lclVectorNodalP,
                                                  blockRowMap, blockColMap,
                                                  blockDomainMap, blockRowMap);

  return vectorNodalP;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    buildProlongator(const int spaceNumber,
                     const Teuchos::RCP<Matrix> &A_nodal,
                     const Teuchos::RCP<MultiVector> &Nullspace,
                     Teuchos::RCP<Matrix> &Prolongator,
                     Teuchos::RCP<MultiVector> &coarseNullspace,
                     Teuchos::RCP<RealValuedMultiVector> &coarseNodalCoords) const {
  using ATS         = Kokkos::ArithTraits<Scalar>;
  using impl_Scalar = typename ATS::val_type;
  using range_type  = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  std::string typeStr;
  switch (spaceNumber) {
    case 0:
      typeStr = "node";
      TEUCHOS_ASSERT(A_nodal.is_null());
      break;
    case 1:
      typeStr = "edge";
      break;
    case 2:
      typeStr = "face";
      break;
    default:
      TEUCHOS_ASSERT(false);
  }

  const bool skipFirstLevel = !A_nodal.is_null();

  RCP<Teuchos::TimeMonitor> tm;
  if (spaceNumber > 0) {
    tm = getTimer("special prolongator " + typeStr);
    GetOStream(Runtime0) << solverName_ + "::compute(): building special " + typeStr + " prolongator" << std::endl;
  }

  RCP<Matrix> projection = buildProjection(spaceNumber, Nullspace);
  dump(projection, typeStr + "Projection.m");

  if (skipFirstLevel) {
    RCP<Matrix> P_nodal;
    RCP<MultiVector> coarseNodalNullspace;

    buildNodalProlongator(A_nodal, P_nodal, coarseNodalNullspace, coarseNodalCoords);

    dump(P_nodal, "P_nodal_" + typeStr + ".m");
    dump(coarseNodalNullspace, "coarseNullspace_nodal_" + typeStr + ".m");

    RCP<Matrix> vectorP_nodal = buildVectorNodalProlongator(P_nodal);

    dump(vectorP_nodal, "vectorP_nodal_" + typeStr + ".m");

    Prolongator = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*projection, false, *vectorP_nodal, false, Prolongator, GetOStream(Runtime0), true, true);

    // This is how ML computes P22 for Darcy.
    // The difference is the scaling by nonzeros. I don't think that that is actually needed.
    //
    // if (spaceNumber==2) {

    //   RCP<Matrix> facesToNodes, aggsToFaces;
    //   {
    //     RCP<Matrix> edgesToNodes = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::BuildCopy(D0_);
    //     Maxwell_Utils<SC,LO,GO,NO>::thresholdedAbs(edgesToNodes, 1e-10);

    //     dump(edgesToNodes, "edgesToNodes.m");

    //     RCP<Matrix> facesToEdges = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::BuildCopy(Dk_1_);
    //     Maxwell_Utils<SC,LO,GO,NO>::thresholdedAbs(facesToEdges, 1e-10);
    //     // facesToEdges = Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(facesToEdges, 1e-2, false);

    //     dump(facesToEdges, "facesToEdges.m");

    //     facesToNodes = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*facesToEdges,false,*edgesToNodes,false,facesToNodes,GetOStream(Runtime0),true,true);
    //     Maxwell_Utils<SC,LO,GO,NO>::thresholdedAbs(facesToNodes, 1e-10);
    //     facesToNodes = Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(facesToNodes, 1e-2, false);
    //   }
    //   aggsToFaces = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*facesToNodes,false,*P_nodal,false,aggsToFaces,GetOStream(Runtime0),true,true);

    //   auto localP = Prolongator->getLocalMatrixDevice();
    //   auto localAggsToFaces = aggsToFaces->getLocalMatrixDevice();
    //   auto localNullspace = Nullspace->getDeviceLocalView(Xpetra::Access::ReadOnly);

    //   size_t dim = dim_;
    //   Kokkos::parallel_for(solverName_+"::buildVectorNodalProlongator_adjustRowptr",
    //                        range_type(0,localP.numRows()),
    //                        KOKKOS_LAMBDA(const LocalOrdinal i) {
    //                          LocalOrdinal nonzeros = localAggsToFaces.graph.row_map(i+1)-localAggsToFaces.graph.row_map(i);
    //                          for (LocalOrdinal jj = localAggsToFaces.graph.row_map(i); jj < localAggsToFaces.graph.row_map(i+1); jj++ ) {
    //                            LocalOrdinal j = localAggsToFaces.graph.entries(jj);
    //                            for (LocalOrdinal k = 0; k<dim; k++)
    //                              for (LocalOrdinal kk = localP.graph.row_map(i); kk < localP.graph.row_map(i+1); kk++)
    //                                if (localP.graph.entries(kk) == (dim * j+k)) {
    //                                  localP.values(kk) = localNullspace(i, k) / nonzeros;
    //                                  break;
    //                                }
    //                          }
    //                        });
    // }
    //

    size_t dim      = dim_;
    coarseNullspace = MultiVectorFactory::Build(vectorP_nodal->getDomainMap(), dim);

    auto localNullspace_nodal  = coarseNodalNullspace->getDeviceLocalView(Xpetra::Access::ReadOnly);
    auto localNullspace_coarse = coarseNullspace->getDeviceLocalView(Xpetra::Access::ReadWrite);
    Kokkos::parallel_for(
        solverName_ + "::buildProlongator_nullspace_" + typeStr,
        range_type(0, coarseNodalNullspace->getLocalLength()),
        KOKKOS_LAMBDA(const size_t i) {
          impl_Scalar val = localNullspace_nodal(i, 0);
          for (size_t j = 0; j < dim; j++)
            localNullspace_coarse(dim * i + j, j) = val;
        });

  } else {
    Prolongator       = projection;
    coarseNodalCoords = NodalCoords_;

    if (spaceNumber == 0) {
      // nothing, just use the default constant vector
    } else if (spaceNumber >= 1) {
      size_t dim                 = dim_;
      coarseNullspace            = MultiVectorFactory::Build(projection->getDomainMap(), dim);
      auto localNullspace_coarse = coarseNullspace->getDeviceLocalView(Xpetra::Access::ReadWrite);
      Kokkos::parallel_for(
          solverName_ + "::buildProlongator_nullspace_" + typeStr,
          range_type(0, coarseNullspace->getLocalLength() / dim),
          KOKKOS_LAMBDA(const size_t i) {
            for (size_t j = 0; j < dim; j++)
              localNullspace_coarse(dim * i + j, j) = 1.0;
          });
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setupSubSolve(Teuchos::RCP<Hierarchy> &hierarchy,
                                                                          Teuchos::RCP<Operator> &thyraPrecOp,
                                                                          const Teuchos::RCP<Matrix> &A,
                                                                          const Teuchos::RCP<MultiVector> &Nullspace,
                                                                          const Teuchos::RCP<RealValuedMultiVector> &Coords,
                                                                          Teuchos::ParameterList &params,
                                                                          std::string &label,
                                                                          const bool reuse,
                                                                          const bool isSingular) {
  int oldRank = SetProcRankVerbose(A->getDomainMap()->getComm()->getRank());
  if (IsPrint(Statistics2)) {
    RCP<ParameterList> pl = rcp(new ParameterList());
    pl->set("printLoadBalancingInfo", true);
    pl->set("printCommInfo", true);
    GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*A, label, pl);
  }
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  if (params.isType<std::string>("Preconditioner Type")) {
    TEUCHOS_ASSERT(!reuse);
    // build a Stratimikos preconditioner
    if (params.get<std::string>("Preconditioner Type") == "MueLu") {
      ParameterList &userParamList = params.sublist("Preconditioner Types").sublist("MueLu").sublist("user data");
      if (!Nullspace.is_null())
        userParamList.set<RCP<MultiVector> >("Nullspace", Nullspace);
      userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", Coords);
    }
    thyraPrecOp = rcp(new XpetraThyraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(coarseA11_, rcp(&params, false)));
  } else
#endif
  {
    // build a MueLu hierarchy

    if (!reuse) {
      ParameterList &userParamList = params.sublist("user data");
      if (!Coords.is_null())
        userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", Coords);
      if (!Nullspace.is_null())
        userParamList.set<RCP<MultiVector> >("Nullspace", Nullspace);

      if (isSingular) {
        std::string coarseType = "";
        if (params.isParameter("coarse: type")) {
          coarseType = params.get<std::string>("coarse: type");
          // Transform string to "Abcde" notation
          std::transform(coarseType.begin(), coarseType.end(), coarseType.begin(), ::tolower);
          std::transform(coarseType.begin(), ++coarseType.begin(), coarseType.begin(), ::toupper);
        }
        if ((coarseType == "" ||
             coarseType == "Klu" ||
             coarseType == "Klu2" ||
             coarseType == "Superlu" ||
             coarseType == "Superlu_dist" ||
             coarseType == "Superludist" ||
             coarseType == "Basker" ||
             coarseType == "Cusolver" ||
             coarseType == "Tacho") &&
            (!params.isSublist("coarse: params") ||
             !params.sublist("coarse: params").isParameter("fix nullspace")))
          params.sublist("coarse: params").set("fix nullspace", true);
      }

      hierarchy = MueLu::CreateXpetraPreconditioner(A, params);
    } else {
      RCP<MueLu::Level> level0 = hierarchy->GetLevel(0);
      level0->Set("A", A);
      hierarchy->SetupRe();
    }
  }
  SetProcRankVerbose(oldRank);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::resetMatrix(RCP<Matrix> SM_Matrix_new, bool ComputePrec) {
  bool reuse = !SM_Matrix_.is_null();
  SM_Matrix_ = SM_Matrix_new;
  dump(SM_Matrix_, "SM.m");
  if (ComputePrec)
    compute(reuse);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::applyInverseAdditive(const MultiVector &RHS, MultiVector &X) const {
  // residual(SM_Matrix_, X, RHS, residual_)
  //
  // P11res_ = P11_^T*residual_ or P11res_ = R11_*residual_
  //
  // Dres_ = Dk_1_^T*residual or Dres_ = Dk_1_T_*residual
  //
  // if ImporterCoarse11_ is not null
  //     ImporterCoarse11: P11res_ -> P11resTmp_
  // if Importer22_ is not null
  //     Importer22: Dres_ -> DresTmp_
  //
  // if coarseA11 is not null
  //
  //     Hierarchy11(P11resSubComm, P11xSubComm)   P11resSubComm aliases P11res or P11resTmp
  //                                               P11xSubComm aliases P11x
  //
  // if A22 is not null
  //
  //     Hierarchy22(DresSubComm, DxSubComm)       DresSubComm aliases Dres or DresTmp
  //                                               DxSubComm aliases Dx
  //
  // if ImporterCoarse11_ is not null
  //     ImporterCoarse11: P11xTmp_ -> P11x
  // if Importer22_ is not null
  //     Importer22: DxTmp_ -> Dx_
  //
  // if fuse
  //     X += P11*P11x
  //     X += P11*Dx
  // else
  //     residual = P11*P11x
  //     residual += Dk_1*Dx
  //     X += residual

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  {  // compute residual

    RCP<Teuchos::TimeMonitor> tmRes = getTimer("residual calculation");
    Utilities::Residual(*SM_Matrix_, X, RHS, *residual_);
  }

  {  // restrict residual to sub-hierarchies

    if (implicitTranspose_) {
      {
        RCP<Teuchos::TimeMonitor> tmRes = getTimer("restriction coarse (1,1) (implicit)");
        P11_->apply(*residual_, *P11res_, Teuchos::TRANS);
      }
      if (!onlyBoundary22_) {
        RCP<Teuchos::TimeMonitor> tmD = getTimer("restriction (2,2) (implicit)");
        Dk_1_->apply(*residual_, *Dres_, Teuchos::TRANS);
      }
    } else {
      if (Dk_1_T_R11_colMapsMatch_) {
        // Column maps of D_T and R11 match, and we're running Tpetra
        {
          RCP<Teuchos::TimeMonitor> tmD = getTimer("restrictions import");
          DTR11Tmp_->doImport(*residual_, *rcp_dynamic_cast<CrsMatrixWrap>(R11_)->getCrsMatrix()->getCrsGraph()->getImporter(), Xpetra::INSERT);
        }
        if (!onlyBoundary22_) {
          RCP<Teuchos::TimeMonitor> tmD = getTimer("restriction (2,2) (explicit)");
          rcp_dynamic_cast<TpetraCrsMatrix>(rcp_dynamic_cast<CrsMatrixWrap>(Dk_1_T_)->getCrsMatrix())->getTpetra_CrsMatrix()->localApply(toTpetra(*DTR11Tmp_), toTpetra(*Dres_), Teuchos::NO_TRANS);
        }
        {
          RCP<Teuchos::TimeMonitor> tmP11 = getTimer("restriction coarse (1,1) (explicit)");
          rcp_dynamic_cast<TpetraCrsMatrix>(rcp_dynamic_cast<CrsMatrixWrap>(R11_)->getCrsMatrix())->getTpetra_CrsMatrix()->localApply(toTpetra(*DTR11Tmp_), toTpetra(*P11res_), Teuchos::NO_TRANS);
        }
      } else {
        {
          RCP<Teuchos::TimeMonitor> tmP11 = getTimer("restriction coarse (1,1) (explicit)");
          R11_->apply(*residual_, *P11res_, Teuchos::NO_TRANS);
        }
        if (!onlyBoundary22_) {
          RCP<Teuchos::TimeMonitor> tmD = getTimer("restriction (2,2) (explicit)");
          Dk_1_T_->apply(*residual_, *Dres_, Teuchos::NO_TRANS);
        }
      }
    }
  }

  {
    RCP<Teuchos::TimeMonitor> tmSubSolves = getTimer("subsolves");

    // block diagonal preconditioner on 2x2 (V-cycle for diagonal blocks)

    if (!ImporterCoarse11_.is_null() && !implicitTranspose_) {
      RCP<Teuchos::TimeMonitor> tmH = getTimer("import coarse (1,1)");
      P11resTmp_->beginImport(*P11res_, *ImporterCoarse11_, Xpetra::INSERT);
    }
    if (!onlyBoundary22_ && !Importer22_.is_null() && !implicitTranspose_) {
      RCP<Teuchos::TimeMonitor> tm22 = getTimer("import (2,2)");
      DresTmp_->beginImport(*Dres_, *Importer22_, Xpetra::INSERT);
    }

    // iterate on coarse (1, 1) block
    if (!coarseA11_.is_null()) {
      if (!ImporterCoarse11_.is_null() && !implicitTranspose_)
        P11resTmp_->endImport(*P11res_, *ImporterCoarse11_, Xpetra::INSERT);

      RCP<Teuchos::TimeMonitor> tmH = getTimer("solve coarse (1,1)", coarseA11_->getRowMap()->getComm());

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
      if (!thyraPrecOpH_.is_null()) {
        Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
        thyraPrecOpH_->apply(*P11resSubComm_, *P11xSubComm_, Teuchos::NO_TRANS, one, zero);
      } else
#endif
        HierarchyCoarse11_->Iterate(*P11resSubComm_, *P11xSubComm_, numItersCoarse11_, true);
    }

    // iterate on (2, 2) block
    if (!A22_.is_null()) {
      if (!onlyBoundary22_ && !Importer22_.is_null() && !implicitTranspose_)
        DresTmp_->endImport(*Dres_, *Importer22_, Xpetra::INSERT);

      RCP<Teuchos::TimeMonitor> tm22 = getTimer("solve (2,2)", A22_->getRowMap()->getComm());
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
      if (!thyraPrecOp22_.is_null()) {
        Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
        thyraPrecOp22_->apply(*DresSubComm_, *DxSubComm_, Teuchos::NO_TRANS, one, zero);
      } else
#endif
        Hierarchy22_->Iterate(*DresSubComm_, *DxSubComm_, numIters22_, true);
    }

    if (coarseA11_.is_null() && !ImporterCoarse11_.is_null() && !implicitTranspose_)
      P11resTmp_->endImport(*P11res_, *ImporterCoarse11_, Xpetra::INSERT);
    if (A22_.is_null() && !onlyBoundary22_ && !Importer22_.is_null() && !implicitTranspose_)
      DresTmp_->endImport(*Dres_, *Importer22_, Xpetra::INSERT);
  }

  if (fuseProlongationAndUpdate_) {
    {  // prolongate (1,1) block
      RCP<Teuchos::TimeMonitor> tmP11 = getTimer("prolongation coarse (1,1) (fused)");
      P11_->apply(*P11x_, X, Teuchos::NO_TRANS, one, one);
    }

    if (!onlyBoundary22_) {  // prolongate (2,2) block
      RCP<Teuchos::TimeMonitor> tmD = getTimer("prolongation (2,2) (fused)");
      Dk_1_->apply(*Dx_, X, Teuchos::NO_TRANS, one, one);
    }
  } else {
    {  // prolongate (1,1) block
      RCP<Teuchos::TimeMonitor> tmP11 = getTimer("prolongation coarse (1,1) (unfused)");
      P11_->apply(*P11x_, *residual_, Teuchos::NO_TRANS);
    }

    if (!onlyBoundary22_) {  // prolongate (2,2) block
      RCP<Teuchos::TimeMonitor> tmD = getTimer("prolongation (2,2) (unfused)");
      Dk_1_->apply(*Dx_, *residual_, Teuchos::NO_TRANS, one, one);
    }

    {  // update current solution
      RCP<Teuchos::TimeMonitor> tmUpdate = getTimer("update");
      X.update(one, *residual_, one);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::solveH(const MultiVector &RHS, MultiVector &X) const {
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  {  // compute residual
    RCP<Teuchos::TimeMonitor> tmRes = getTimer("residual calculation");
    Utilities::Residual(*SM_Matrix_, X, RHS, *residual_);
    if (implicitTranspose_)
      P11_->apply(*residual_, *P11res_, Teuchos::TRANS);
    else
      R11_->apply(*residual_, *P11res_, Teuchos::NO_TRANS);
  }

  {  // solve coarse (1,1) block
    if (!ImporterCoarse11_.is_null() && !implicitTranspose_) {
      RCP<Teuchos::TimeMonitor> tmH = getTimer("import coarse (1,1)");
      P11resTmp_->doImport(*P11res_, *ImporterCoarse11_, Xpetra::INSERT);
    }
    if (!coarseA11_.is_null()) {
      RCP<Teuchos::TimeMonitor> tmH = getTimer("solve coarse (1,1)", coarseA11_->getRowMap()->getComm());
      HierarchyCoarse11_->Iterate(*P11resSubComm_, *P11xSubComm_, numItersCoarse11_, true);
    }
  }

  {  // update current solution
    RCP<Teuchos::TimeMonitor> tmUp = getTimer("update");
    P11_->apply(*P11x_, *residual_, Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::solve22(const MultiVector &RHS, MultiVector &X) const {
  if (onlyBoundary22_)
    return;

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  {  // compute residual
    RCP<Teuchos::TimeMonitor> tmRes = getTimer("residual calculation");
    Utilities::Residual(*SM_Matrix_, X, RHS, *residual_);
    if (implicitTranspose_)
      Dk_1_->apply(*residual_, *Dres_, Teuchos::TRANS);
    else
      Dk_1_T_->apply(*residual_, *Dres_, Teuchos::NO_TRANS);
  }

  {  // solve (2,2) block
    if (!Importer22_.is_null() && !implicitTranspose_) {
      RCP<Teuchos::TimeMonitor> tm22 = getTimer("import (2,2)");
      DresTmp_->doImport(*Dres_, *Importer22_, Xpetra::INSERT);
    }
    if (!A22_.is_null()) {
      RCP<Teuchos::TimeMonitor> tm22 = getTimer("solve (2,2)", A22_->getRowMap()->getComm());
      Hierarchy22_->Iterate(*DresSubComm_, *DxSubComm_, numIters22_, true);
    }
  }

  {  // update current solution
    RCP<Teuchos::TimeMonitor> tmUp = getTimer("update");
    Dk_1_->apply(*Dx_, *residual_, Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector &RHS, MultiVector &X,
                                                                  Teuchos::ETransp /* mode */,
                                                                  Scalar /* alpha */,
                                                                  Scalar /* beta */) const {
  RCP<Teuchos::TimeMonitor> tm = getTimer("solve");

  // make sure that we have enough temporary memory
  if (!onlyBoundary11_ && X.getNumVectors() != P11res_->getNumVectors())
    allocateMemory(X.getNumVectors());

  {  // apply pre-smoothing

    RCP<Teuchos::TimeMonitor> tmSm = getTimer("smoothing");

    PreSmoother11_->Apply(X, RHS, use_as_preconditioner_);
  }

  // do solve for the 2x2 block system
  if (mode_ == "additive")
    applyInverseAdditive(RHS, X);
  else if (mode_ == "121") {
    solveH(RHS, X);
    solve22(RHS, X);
    solveH(RHS, X);
  } else if (mode_ == "212") {
    solve22(RHS, X);
    solveH(RHS, X);
    solve22(RHS, X);
  } else if (mode_ == "1")
    solveH(RHS, X);
  else if (mode_ == "2")
    solve22(RHS, X);
  else if (mode_ == "7") {
    solveH(RHS, X);
    {  // apply pre-smoothing

      RCP<Teuchos::TimeMonitor> tmSm = getTimer("smoothing");

      PreSmoother11_->Apply(X, RHS, false);
    }
    solve22(RHS, X);
    {  // apply post-smoothing

      RCP<Teuchos::TimeMonitor> tmSm = getTimer("smoothing");

      PostSmoother11_->Apply(X, RHS, false);
    }
    solveH(RHS, X);
  } else if (mode_ == "none") {
    // do nothing
  } else
    applyInverseAdditive(RHS, X);

  {  // apply post-smoothing

    RCP<Teuchos::TimeMonitor> tmSm = getTimer("smoothing");

    PostSmoother11_->Apply(X, RHS, false);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasTransposeApply() const {
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    RefMaxwell(const Teuchos::RCP<Matrix> &SM_Matrix,
               Teuchos::ParameterList &List,
               bool ComputePrec) {
  int spaceNumber = List.get<int>("refmaxwell: space number", 1);

  RCP<Matrix> Dk_1, Dk_2, D0;
  RCP<Matrix> M1_beta, M1_alpha;
  RCP<Matrix> Mk_one, Mk_1_one;
  RCP<Matrix> invMk_1_invBeta, invMk_2_invAlpha;
  RCP<MultiVector> Nullspace11, Nullspace22;
  RCP<RealValuedMultiVector> NodalCoords;

  Dk_1 = pop(List, "Dk_1", Dk_1);
  Dk_2 = pop<RCP<Matrix> >(List, "Dk_2", Dk_2);
  D0   = pop<RCP<Matrix> >(List, "D0", D0);

  M1_beta  = pop<RCP<Matrix> >(List, "M1_beta", M1_beta);
  M1_alpha = pop<RCP<Matrix> >(List, "M1_alpha", M1_alpha);

  Mk_one   = pop<RCP<Matrix> >(List, "Mk_one", Mk_one);
  Mk_1_one = pop<RCP<Matrix> >(List, "Mk_1_one", Mk_1_one);

  invMk_1_invBeta  = pop<RCP<Matrix> >(List, "invMk_1_invBeta", invMk_1_invBeta);
  invMk_2_invAlpha = pop<RCP<Matrix> >(List, "invMk_2_invAlpha", invMk_2_invAlpha);

  Nullspace11 = pop<RCP<MultiVector> >(List, "Nullspace11", Nullspace11);
  Nullspace22 = pop<RCP<MultiVector> >(List, "Nullspace22", Nullspace22);
  NodalCoords = pop<RCP<RealValuedMultiVector> >(List, "Coordinates", NodalCoords);

  // old parameter names
  if (List.isType<RCP<Matrix> >("Ms")) {
    if (M1_beta.is_null())
      M1_beta = pop<RCP<Matrix> >(List, "Ms");
    else
      TEUCHOS_ASSERT(false);
  }
  if (List.isType<RCP<Matrix> >("M1")) {
    if (Mk_one.is_null())
      Mk_one = pop<RCP<Matrix> >(List, "M1");
    else
      TEUCHOS_ASSERT(false);
  }
  if (List.isType<RCP<Matrix> >("M0inv")) {
    if (invMk_1_invBeta.is_null())
      invMk_1_invBeta = pop<RCP<Matrix> >(List, "M0inv");
    else
      TEUCHOS_ASSERT(false);
  }
  if (List.isType<RCP<MultiVector> >("Nullspace")) {
    if (Nullspace11.is_null())
      Nullspace11 = pop<RCP<MultiVector> >(List, "Nullspace");
    else
      TEUCHOS_ASSERT(false);
  }

  if (spaceNumber == 1) {
    if (Dk_1.is_null())
      Dk_1 = D0;
    else if (D0.is_null())
      D0 = Dk_1;
    if (M1_beta.is_null())
      M1_beta = Mk_one;
  } else if (spaceNumber == 2) {
    if (Dk_2.is_null())
      Dk_2 = D0;
    else if (D0.is_null())
      D0 = Dk_2;
  }

  initialize(spaceNumber,
             Dk_1, Dk_2, D0,
             M1_beta, M1_alpha,
             Mk_one, Mk_1_one,
             invMk_1_invBeta, invMk_2_invAlpha,
             Nullspace11, Nullspace22,
             NodalCoords,
             List);

  if (SM_Matrix != Teuchos::null)
    resetMatrix(SM_Matrix, ComputePrec);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    initialize(const Teuchos::RCP<Matrix> &D0_Matrix,
               const Teuchos::RCP<Matrix> &Ms_Matrix,
               const Teuchos::RCP<Matrix> &M0inv_Matrix,
               const Teuchos::RCP<Matrix> &M1_Matrix,
               const Teuchos::RCP<MultiVector> &Nullspace11,
               const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
               Teuchos::ParameterList &List) {
  initialize(1,
             D0_Matrix, Teuchos::null, D0_Matrix,
             Ms_Matrix, Teuchos::null,
             M1_Matrix, Teuchos::null,
             M0inv_Matrix, Teuchos::null,
             Nullspace11, Teuchos::null,
             NodalCoords,
             List);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    initialize(const int k,
               const Teuchos::RCP<Matrix> &Dk_1,
               const Teuchos::RCP<Matrix> &Dk_2,
               const Teuchos::RCP<Matrix> &D0,
               const Teuchos::RCP<Matrix> &M1_beta,
               const Teuchos::RCP<Matrix> &M1_alpha,
               const Teuchos::RCP<Matrix> &Mk_one,
               const Teuchos::RCP<Matrix> &Mk_1_one,
               const Teuchos::RCP<Matrix> &invMk_1_invBeta,
               const Teuchos::RCP<Matrix> &invMk_2_invAlpha,
               const Teuchos::RCP<MultiVector> &Nullspace11,
               const Teuchos::RCP<MultiVector> &Nullspace22,
               const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
               Teuchos::ParameterList &List) {
  spaceNumber_ = k;
  if (spaceNumber_ == 1)
    solverName_ = "RefMaxwell";
  else if (spaceNumber_ == 2)
    solverName_ = "RefDarcy";
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                               "spaceNumber needs to be 1 (HCurl) or 2 (HDiv)");
  HierarchyCoarse11_ = Teuchos::null;
  Hierarchy22_       = Teuchos::null;
  PreSmoother11_     = Teuchos::null;
  PostSmoother11_    = Teuchos::null;
  disable_addon_     = false;
  disable_addon_22_  = true;
  mode_              = "additive";

  // set parameters
  setParameters(List);

  // some pre-conditions
  TEUCHOS_ASSERT((k == 1) || (k == 2));
  // Need Dk_1
  TEUCHOS_ASSERT(Dk_1 != Teuchos::null);
  // Need D0 for aggregation
  TEUCHOS_ASSERT(D0 != Teuchos::null);

  // Need M1_beta for aggregation
  TEUCHOS_ASSERT(M1_beta != Teuchos::null);
  // Need M1_alpha for aggregation if k>=1
  if (k >= 2)
    TEUCHOS_ASSERT(M1_alpha != Teuchos::null);

  if (!disable_addon_) {
    // Need Mk_one and invMk_1_invBeta for addon11
    TEUCHOS_ASSERT(Mk_one != Teuchos::null);
    TEUCHOS_ASSERT(invMk_1_invBeta != Teuchos::null);
  }

  if ((k >= 2) && !disable_addon_22_) {
    // Need Dk_2, Mk_1_one and invMk_2_invAlpha for addon22
    TEUCHOS_ASSERT(Dk_2 != Teuchos::null);
    TEUCHOS_ASSERT(Mk_1_one != Teuchos::null);
    TEUCHOS_ASSERT(invMk_2_invAlpha != Teuchos::null);
  }

#ifdef HAVE_MUELU_DEBUG

  TEUCHOS_ASSERT(D0->getRangeMap()->isSameAs(*D0->getRowMap()));

  // M1_beta is square
  TEUCHOS_ASSERT(M1_beta->getDomainMap()->isSameAs(*M1_beta->getRangeMap()));
  TEUCHOS_ASSERT(M1_beta->getDomainMap()->isSameAs(*M1_beta->getRowMap()));

  // M1_beta is consistent with D0
  TEUCHOS_ASSERT(M1_beta->getDomainMap()->isSameAs(*D0->getRangeMap()));

  if (k >= 2) {
    // M1_alpha is square
    TEUCHOS_ASSERT(M1_alpha->getDomainMap()->isSameAs(*M1_alpha->getRangeMap()));
    TEUCHOS_ASSERT(M1_alpha->getDomainMap()->isSameAs(*M1_alpha->getRowMap()));

    // M1_alpha is consistent with D0
    TEUCHOS_ASSERT(M1_alpha->getDomainMap()->isSameAs(*D0->getRangeMap()))
  }

  if (!disable_addon_) {
    // Mk_one is square
    TEUCHOS_ASSERT(Mk_one->getDomainMap()->isSameAs(*Mk_one->getRangeMap()));
    TEUCHOS_ASSERT(Mk_one->getDomainMap()->isSameAs(*Mk_one->getRowMap()));

    // Mk_one is consistent with Dk_1
    TEUCHOS_ASSERT(Mk_one->getDomainMap()->isSameAs(*Dk_1->getRangeMap()));

    // invMk_1_invBeta is square
    TEUCHOS_ASSERT(invMk_1_invBeta->getDomainMap()->isSameAs(*invMk_1_invBeta->getRangeMap()));
    TEUCHOS_ASSERT(invMk_1_invBeta->getDomainMap()->isSameAs(*invMk_1_invBeta->getRowMap()));

    // invMk_1_invBeta is consistent with Dk_1
    TEUCHOS_ASSERT(Mk_one->getDomainMap()->isSameAs(*Dk_1->getRangeMap()));
  }

  if ((k >= 2) && !disable_addon_22_) {
    // Mk_1_one is square
    TEUCHOS_ASSERT(Mk_1_one->getDomainMap()->isSameAs(*Mk_1_one->getRangeMap()));
    TEUCHOS_ASSERT(Mk_1_one->getDomainMap()->isSameAs(*Mk_1_one->getRowMap()));

    // Mk_1_one is consistent with Dk_1
    TEUCHOS_ASSERT(Mk_1_one->getDomainMap()->isSameAs(*Dk_1->getDomainMap()));

    // Mk_1_one is consistent with Dk_2
    TEUCHOS_ASSERT(Mk_1_one->getDomainMap()->isSameAs(*Dk_2->getRangeMap()));

    // invMk_2_invAlpha is square
    TEUCHOS_ASSERT(invMk_2_invAlpha->getDomainMap()->isSameAs(*invMk_2_invAlpha->getRangeMap()));
    TEUCHOS_ASSERT(invMk_2_invAlpha->getDomainMap()->isSameAs(*invMk_2_invAlpha->getRowMap()));

    // invMk_2_invAlpha is consistent with Dk_2
    TEUCHOS_ASSERT(invMk_2_invAlpha->getDomainMap()->isSameAs(*Dk_2->getDomainMap()));
  }
#endif

  D0_ = D0;
  if (Dk_1->getRowMap()->lib() == Xpetra::UseTpetra) {
    // We will remove boundary conditions from Dk_1, and potentially change maps, so we copy the input.
    // Fortunately, Dk_1 is quite sparse.
    // We cannot use the Tpetra copy constructor, since it does not copy the graph.

    RCP<Matrix> Dk_1copy       = MatrixFactory::Build(Dk_1->getRowMap(), Dk_1->getColMap(), 0);
    RCP<CrsMatrix> Dk_1copyCrs = rcp_dynamic_cast<CrsMatrixWrap>(Dk_1copy, true)->getCrsMatrix();
    ArrayRCP<const size_t> Dk_1rowptr_RCP;
    ArrayRCP<const LO> Dk_1colind_RCP;
    ArrayRCP<const SC> Dk_1vals_RCP;
    rcp_dynamic_cast<CrsMatrixWrap>(Dk_1, true)->getCrsMatrix()->getAllValues(Dk_1rowptr_RCP, Dk_1colind_RCP, Dk_1vals_RCP);

    ArrayRCP<size_t> Dk_1copyrowptr_RCP;
    ArrayRCP<LO> Dk_1copycolind_RCP;
    ArrayRCP<SC> Dk_1copyvals_RCP;
    Dk_1copyCrs->allocateAllValues(Dk_1vals_RCP.size(), Dk_1copyrowptr_RCP, Dk_1copycolind_RCP, Dk_1copyvals_RCP);
    Dk_1copyrowptr_RCP.deepCopy(Dk_1rowptr_RCP());
    Dk_1copycolind_RCP.deepCopy(Dk_1colind_RCP());
    Dk_1copyvals_RCP.deepCopy(Dk_1vals_RCP());
    Dk_1copyCrs->setAllValues(Dk_1copyrowptr_RCP,
                              Dk_1copycolind_RCP,
                              Dk_1copyvals_RCP);
    Dk_1copyCrs->expertStaticFillComplete(Dk_1->getDomainMap(), Dk_1->getRangeMap(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(Dk_1, true)->getCrsMatrix()->getCrsGraph()->getImporter(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(Dk_1, true)->getCrsMatrix()->getCrsGraph()->getExporter());
    Dk_1_ = Dk_1copy;
  } else
    Dk_1_ = MatrixFactory::BuildCopy(Dk_1);

  if ((!Dk_2.is_null()) && (Dk_2->getRowMap()->lib() == Xpetra::UseTpetra)) {
    // We will remove boundary conditions from Dk_2, and potentially change maps, so we copy the input.
    // Fortunately, Dk_2 is quite sparse.
    // We cannot use the Tpetra copy constructor, since it does not copy the graph.

    RCP<Matrix> Dk_2copy       = MatrixFactory::Build(Dk_2->getRowMap(), Dk_2->getColMap(), 0);
    RCP<CrsMatrix> Dk_2copyCrs = rcp_dynamic_cast<CrsMatrixWrap>(Dk_2copy, true)->getCrsMatrix();
    ArrayRCP<const size_t> Dk_2rowptr_RCP;
    ArrayRCP<const LO> Dk_2colind_RCP;
    ArrayRCP<const SC> Dk_2vals_RCP;
    rcp_dynamic_cast<CrsMatrixWrap>(Dk_2, true)->getCrsMatrix()->getAllValues(Dk_2rowptr_RCP, Dk_2colind_RCP, Dk_2vals_RCP);

    ArrayRCP<size_t> Dk_2copyrowptr_RCP;
    ArrayRCP<LO> Dk_2copycolind_RCP;
    ArrayRCP<SC> Dk_2copyvals_RCP;
    Dk_2copyCrs->allocateAllValues(Dk_2vals_RCP.size(), Dk_2copyrowptr_RCP, Dk_2copycolind_RCP, Dk_2copyvals_RCP);
    Dk_2copyrowptr_RCP.deepCopy(Dk_2rowptr_RCP());
    Dk_2copycolind_RCP.deepCopy(Dk_2colind_RCP());
    Dk_2copyvals_RCP.deepCopy(Dk_2vals_RCP());
    Dk_2copyCrs->setAllValues(Dk_2copyrowptr_RCP,
                              Dk_2copycolind_RCP,
                              Dk_2copyvals_RCP);
    Dk_2copyCrs->expertStaticFillComplete(Dk_2->getDomainMap(), Dk_2->getRangeMap(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(Dk_2, true)->getCrsMatrix()->getCrsGraph()->getImporter(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(Dk_2, true)->getCrsMatrix()->getCrsGraph()->getExporter());
    Dk_2_ = Dk_2copy;
  } else if (!Dk_2.is_null())
    Dk_2_ = MatrixFactory::BuildCopy(Dk_2);

  M1_beta_  = M1_beta;
  M1_alpha_ = M1_alpha;

  Mk_one_   = Mk_one;
  Mk_1_one_ = Mk_1_one;

  invMk_1_invBeta_  = invMk_1_invBeta;
  invMk_2_invAlpha_ = invMk_2_invAlpha;

  NodalCoords_ = NodalCoords;
  Nullspace11_ = Nullspace11;
  Nullspace22_ = Nullspace22;

  dump(D0_, "D0.m");
  dump(Dk_1_, "Dk_1_clean.m");
  dump(Dk_2_, "Dk_2_clean.m");

  dump(M1_beta_, "M1_beta.m");
  dump(M1_alpha_, "M1_alpha.m");

  dump(Mk_one_, "Mk_one.m");
  dump(Mk_1_one_, "Mk_1_one.m");

  dump(invMk_1_invBeta_, "invMk_1_invBeta.m");
  dump(invMk_2_invAlpha_, "invMk_2_invAlpha.m");

  dumpCoords(NodalCoords_, "coords.m");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel /* verbLevel */) const {
  std::ostringstream oss;

  RCP<const Teuchos::Comm<int> > comm = SM_Matrix_->getDomainMap()->getComm();

#ifdef HAVE_MPI
  int root;
  if (!coarseA11_.is_null())
    root = comm->getRank();
  else
    root = -1;

  int actualRoot;
  reduceAll(*comm, Teuchos::REDUCE_MAX, root, Teuchos::ptr(&actualRoot));
  root = actualRoot;
#endif

  oss << "\n--------------------------------------------------------------------------------\n"
      << "---                            " + solverName_ +
             " Summary                            ---\n"
             "--------------------------------------------------------------------------------"
      << std::endl;
  oss << std::endl;

  GlobalOrdinal numRows;
  GlobalOrdinal nnz;

  SM_Matrix_->getRowMap()->getComm()->barrier();

  numRows = SM_Matrix_->getGlobalNumRows();
  nnz     = SM_Matrix_->getGlobalNumEntries();

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
  oss << "(1, 1)" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;

  if (!A22_.is_null()) {
    numRows = A22_->getGlobalNumRows();
    nnz     = A22_->getGlobalNumEntries();

    oss << "(2, 2)" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;
  }

  oss << std::endl;

  {
    if (PreSmoother11_ != null && PreSmoother11_ == PostSmoother11_)
      oss << "Smoother 11 both : " << PreSmoother11_->description() << std::endl;
    else {
      oss << "Smoother 11 pre  : "
          << (PreSmoother11_ != null ? PreSmoother11_->description() : "no smoother") << std::endl;
      oss << "Smoother 11 post : "
          << (PostSmoother11_ != null ? PostSmoother11_->description() : "no smoother") << std::endl;
    }
  }
  oss << std::endl;

  std::string outstr = oss.str();

#ifdef HAVE_MPI
  RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
  MPI_Comm rawComm                          = (*mpiComm->getRawMpiComm())();

  int strLength = outstr.size();
  MPI_Bcast(&strLength, 1, MPI_INT, root, rawComm);
  if (comm->getRank() != root)
    outstr.resize(strLength);
  MPI_Bcast(&outstr[0], strLength, MPI_CHAR, root, rawComm);
#endif

  out << outstr;

  if (!HierarchyCoarse11_.is_null())
    HierarchyCoarse11_->describe(out, GetVerbLevel());

  if (!Hierarchy22_.is_null())
    Hierarchy22_->describe(out, GetVerbLevel());

  if (IsPrint(Statistics2)) {
    // Print the grid of processors
    std::ostringstream oss2;

    oss2 << "Sub-solver distribution over ranks" << std::endl;
    oss2 << "( (1,1) block only is indicated by '1', (2,2) block only by '2', and both blocks by 'B' and none by '.')" << std::endl;

    int numProcs = comm->getSize();
#ifdef HAVE_MPI
    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();
#endif

    char status = 0;
    if (!coarseA11_.is_null())
      status += 1;
    if (!A22_.is_null())
      status += 2;
    std::vector<char> states(numProcs, 0);
#ifdef HAVE_MPI
    MPI_Gather(&status, 1, MPI_CHAR, &states[0], 1, MPI_CHAR, 0, *rawMpiComm);
#else
    states.push_back(status);
#endif

    int rowWidth = std::min(Teuchos::as<int>(ceil(sqrt(numProcs))), 100);
    for (int proc = 0; proc < numProcs; proc += rowWidth) {
      for (int j = 0; j < rowWidth; j++)
        if (proc + j < numProcs)
          if (states[proc + j] == 0)
            oss2 << ".";
          else if (states[proc + j] == 1)
            oss2 << "1";
          else if (states[proc + j] == 2)
            oss2 << "2";
          else
            oss2 << "B";
        else
          oss2 << " ";

      oss2 << "      " << proc << ":" << std::min(proc + rowWidth, numProcs) - 1 << std::endl;
    }
    oss2 << std::endl;
    GetOStream(Statistics2) << oss2.str();
  }
}

}  // namespace MueLu

#define MUELU_REFMAXWELL_SHORT
#endif  // ifdef MUELU_REFMAXWELL_DEF_HPP
