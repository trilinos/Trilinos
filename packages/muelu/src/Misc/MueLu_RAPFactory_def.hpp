// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RAPFACTORY_DEF_HPP
#define MUELU_RAPFACTORY_DEF_HPP

#include <sstream>

#include <Xpetra_Matrix.hpp>

#include "MueLu_RAPFactory_decl.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RAPFactory()
  : hasDeclaredInput_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RAPFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("transpose: use implicit");
  SET_VALID_ENTRY("rap: triple product");
  SET_VALID_ENTRY("rap: fix zero diagonals");
  SET_VALID_ENTRY("rap: fix zero diagonals threshold");
  SET_VALID_ENTRY("rap: fix zero diagonals replacement");
  SET_VALID_ENTRY("rap: relative diagonal floor");
#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("P", null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase> >("R", null, "Restrictor factory");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  if (pL.get<bool>("transpose: use implicit") == false)
    Input(coarseLevel, "R");

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);

  hasDeclaredInput_ = true;
}

// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
// void TripleMatrixProduct(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& R,
//                          const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
//                          const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& P,
//                          Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ac,
//                          const Teuchos::ParameterList& pL,
//                          const MueLu::BaseClass& verbObj,
//                          Teuchos::RCP<Teuchos::ParameterList>& APparams = Teuchos::null,
//                          Teuchos::RCP<Teuchos::ParameterList>& RAPparams = Teuchos::null,
//                          Level* coarseLevel = nullptr) {
//   using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//   using MatrixMatrix = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//   using MatrixFactory = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//   using PerfUtils = MueLu::PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

//   const bool doTranspose       = true;
//   const bool doFillComplete    = true;
//   const bool doOptimizeStorage = true;

//   std::string levelstr, labelstr;
//   if (coarseLevel != nullptr) {
//     std::ostringstream levelss;
//     levelss << coarseLevel->GetLevelID();
//     levelstr = levelss.str();
//     labelstr = FormattingHelper::getColonLabel(coarseLevel->getObjectLabel());
//   }

//   bool isGPU = Node::is_gpu;

//   // Reuse coarse matrix memory if available (multiple solve)
//   if (RAPparams.is_null())
//     RAPparams = rcp(new ParameterList);
//   if (pL.isSublist("matrixmatrix: kernel params"))
//     RAPparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));

//   if (RAPparams->isParameter("graph")) {
//     Ac = RAPparams->get<RCP<Matrix> >("graph");

//     // Some eigenvalue may have been cached with the matrix in the previous run.
//     // As the matrix values will be updated, we need to reset the eigenvalue.
//     Ac->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<Scalar>::one());
//   }

//   // We *always* need global constants for the RAP, but not for the temps
//   RAPparams->set("compute global constants: temporaries", RAPparams->get("compute global constants: temporaries", false));
//   RAPparams->set("compute global constants", true);

//   if (pL.get<bool>("rap: triple product") == false || isGPU) {
//     if (pL.get<bool>("rap: triple product") && isGPU)
//       verbObj.GetOStream(Warnings1) << "Switching from triple product to R x (A x P) since triple product has not been implemented for "
//       << Node::execution_space::name() << std::endl;

//     RCP<Matrix> AP;

//     // Reuse pattern if available (multiple solve)
//     if (APparams.is_null())
//       APparams = rcp(new ParameterList);
//     if (pL.isSublist("matrixmatrix: kernel params"))
//       APparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));

//     // By default, we don't need global constants for A*P
//     APparams->set("compute global constants: temporaries", APparams->get("compute global constants: temporaries", false));
//     APparams->set("compute global constants", APparams->get("compute global constants", false));

//     if (APparams->isParameter("graph"))
//       AP = APparams->get<RCP<Matrix> >("graph");

//     std::string timerstrAP = "MueLu::A*P";
//     if (!labelstr.empty())
//       timerstrAP = labelstr+timerstrAP;
//     if (!levelstr.empty())
//       timerstrAP = timerstrAP+"-"+levelstr;

//     {
//       SubFactoryMonitor subM(verbObj, "MxM: A x P", *coarseLevel);

//       AP = MatrixMatrix::Multiply(*A, !doTranspose, *P, !doTranspose, AP, verbObj.GetOStream(Statistics2),
//                                   doFillComplete, doOptimizeStorage, timerstrAP, APparams);
//     }

//     // Allow optimization of storage.
//     // This is necessary for new faster Epetra MM kernels.
//     // Seems to work with matrix modifications to repair diagonal entries.

//     std::string timerstrRAP;
//     if (pL.get<bool>("transpose: use implicit") == true)
//       timerstrRAP = "MueLu::R*(AP)-implicit";
//     else
//       timerstrRAP = "MueLu::R*(AP)-explicit";
//     if (!labelstr.empty())
//       timerstrRAP = labelstr+timerstrRAP;
//     if (!levelstr.empty())
//       timerstrRAP = timerstrRAP+"-"+levelstr;

//     if (pL.get<bool>("transpose: use implicit") == true) {
//       SubFactoryMonitor m2(verbObj, "MxM: P' x (AP) (implicit)", *coarseLevel);

//       Ac = MatrixMatrix::Multiply(*P, doTranspose, *AP, !doTranspose, Ac, verbObj.GetOStream(Statistics2),
//                                   doFillComplete, doOptimizeStorage, timerstrRAP, RAPparams);

//     } else {
//       SubFactoryMonitor m2(verbObj, "MxM: R x (AP) (explicit)", *coarseLevel);

//       Ac = MatrixMatrix::Multiply(*R, !doTranspose, *AP, !doTranspose, Ac, verbObj.GetOStream(Statistics2),
//                                   doFillComplete, doOptimizeStorage, timerstrRAP, RAPparams);
//     }

//     if (!isGPU) {
//       APparams->set("graph", AP);
//     }

//   } else {

//     std::string timerstrRAP;
//     if (pL.get<bool>("transpose: use implicit") == true)
//       timerstrRAP = "MueLu::R*A*P-implicit";
//     else
//       timerstrRAP = "MueLu::R*A*P-explicit";
//     if (!labelstr.empty())
//       timerstrRAP = labelstr+timerstrRAP;
//     if (!levelstr.empty())
//       timerstrRAP = timerstrRAP+"-"+levelstr;

//     if (pL.get<bool>("transpose: use implicit") == true) {
//       Ac = MatrixFactory::Build(P->getDomainMap(), Teuchos::as<LocalOrdinal>(0));

//       SubFactoryMonitor m2(verbObj, "MxMxM: R x A x P (implicit)", *coarseLevel);

//       Xpetra::TripleMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
//             MultiplyRAP(*P, doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
//                         doOptimizeStorage, timerstrRAP,
//                         RAPparams);
//     } else {
//       Ac = MatrixFactory::Build(R->getRowMap(), Teuchos::as<LocalOrdinal>(0));

//       SubFactoryMonitor m2(verbObj, "MxMxM: R x A x P (explicit)", *coarseLevel);

//       Xpetra::TripleMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
//             MultiplyRAP(*R, !doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
//                         doOptimizeStorage, timerstrRAP,
//                         RAPparams);
//     }
//   }

//   Teuchos::ArrayView<const double> relativeFloor = pL.get<Teuchos::Array<double> >("rap: relative diagonal floor")();
//   if (relativeFloor.size() > 0) {
//     Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RelativeDiagonalBoost(Ac, relativeFloor, verbObj.GetOStream(Statistics2));
//   }

//   bool repairZeroDiagonals = pL.get<bool>("RepairMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
//   bool checkAc             = pL.get<bool>("CheckMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");

//   if (checkAc || repairZeroDiagonals) {
//     using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
//     magnitudeType threshold;
//     if (pL.isType<magnitudeType>("rap: fix zero diagonals threshold"))
//       threshold = pL.get<magnitudeType>("rap: fix zero diagonals threshold");
//     else
//       threshold = Teuchos::as<magnitudeType>(pL.get<double>("rap: fix zero diagonals threshold"));
//     Scalar replacement = Teuchos::as<Scalar>(pL.get<double>("rap: fix zero diagonals replacement"));
//     Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckRepairMainDiagonal(Ac, repairZeroDiagonals, verbObj.GetOStream(Warnings1), threshold, replacement);
//   }

//   if (verbObj.IsPrint(Statistics2)) {
//     RCP<ParameterList> params = rcp(new ParameterList());
//     params->set("printLoadBalancingInfo", true);
//     params->set("printCommInfo", true);
//     verbObj.GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);
//   }

//   if (!isGPU) {
//     RAPparams->set("graph", Ac);
//   }
// }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  RCP<Matrix> Ac;

  {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                               "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

    RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");
    RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P"), AP, R;
    // We don't have a valid P (e.g., # global aggregates = 0) so we bail.
    // This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
    if (P.is_null()) {
      Ac = Teuchos::null;
      Set(coarseLevel, "A", Ac);
      return;
    }

    const Teuchos::ParameterList& pL = GetParameterList();
    const bool useImplicit           = pL.get<bool>("transpose: use implicit");
    bool isGPU                       = Node::is_gpu;

    Teuchos::RCP<Teuchos::ParameterList> APparams;
    Teuchos::RCP<Teuchos::ParameterList> RAPparams;
    if (coarseLevel.IsAvailable("AP reuse data", this))
      APparams = coarseLevel.Get<RCP<ParameterList> >("AP reuse data", this);
    if (coarseLevel.IsAvailable("RAP reuse data", this))
      RAPparams = coarseLevel.Get<RCP<ParameterList> >("RAP reuse data", this);
    if (!useImplicit)
      R = Get<RCP<Matrix> >(coarseLevel, "R");
    Utilities::TripleMatrixProduct(R, A, P, Ac, pL, *this, APparams, RAPparams, &coarseLevel);

    if (!Ac.is_null()) {
      std::ostringstream oss;
      oss << "A_" << coarseLevel.GetLevelID();
      Ac->setObjectLabel(oss.str());
    }
    Set(coarseLevel, "A", Ac);
    if (!isGPU) {
      if (!useImplicit)
        Set(coarseLevel, "AP reuse data", APparams);
      Set(coarseLevel, "RAP reuse data", RAPparams);
    }
  }

#ifdef HAVE_MUELU_DEBUG
  MatrixUtils::checkLocalRowMapMatchesColMap(*Ac);
#endif  // HAVE_MUELU_DEBUG

  if (transferFacts_.begin() != transferFacts_.end()) {
    SubFactoryMonitor m(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
      RCP<const FactoryBase> fac = *it;
      GetOStream(Runtime0) << "RAPFactory: call transfer factory: " << fac->description() << std::endl;
      fac->CallBuild(coarseLevel);
      // Coordinates transfer is marginally different from all other operations
      // because it is *optional*, and not required. For instance, we may need
      // coordinates only on level 4 if we start repartitioning from that level,
      // but we don't need them on level 1,2,3. As our current Hierarchy setup
      // assumes propagation of dependencies only through three levels, this
      // means that we need to rely on other methods to propagate optional data.
      //
      // The method currently used is through RAP transfer factories, which are
      // simply factories which are called at the end of RAP with a single goal:
      // transfer some fine data to coarser level. Because these factories are
      // kind of outside of the mainline factories, they behave different. In
      // particular, we call their Build method explicitly, rather than through
      // Get calls. This difference is significant, as the Get call is smart
      // enough to know when to release all factory dependencies, and Build is
      // dumb. This led to the following CoordinatesTransferFactory sequence:
      // 1. Request level 0
      // 2. Request level 1
      // 3. Request level 0
      // 4. Release level 0
      // 5. Release level 1
      //
      // The problem is missing "6. Release level 0". Because it was missing,
      // we had outstanding request on "Coordinates", "Aggregates" and
      // "CoarseMap" on level 0.
      //
      // This was fixed by explicitly calling Release on transfer factories in
      // RAPFactory. I am still unsure how exactly it works, but now we have
      // clear data requests for all levels.
      coarseLevel.Release(*fac);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");
  transferFacts_.push_back(factory);
}

}  // namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif  // MUELU_RAPFACTORY_DEF_HPP
