// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RAPSHIFTFACTORY_DEF_HPP
#define MUELU_RAPSHIFTFACTORY_DEF_HPP

#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixFactory2.hpp>

#include "MueLu_RAPShiftFactory_decl.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RAPShiftFactory()
  : implicitTranspose_(false) {}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("transpose: use implicit");
  SET_VALID_ENTRY("rap: fix zero diagonals");
  SET_VALID_ENTRY("rap: shift");
  SET_VALID_ENTRY("rap: shift array");
  SET_VALID_ENTRY("rap: cfl array");
  SET_VALID_ENTRY("rap: shift diagonal M");
  SET_VALID_ENTRY("rap: shift low storage");
  SET_VALID_ENTRY("rap: relative diagonal floor");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("M", Teuchos::null, "Generating factory of the matrix M used during the non-Galerkin RAP");
  validParamList->set<RCP<const FactoryBase> >("Mdiag", Teuchos::null, "Generating factory of the matrix Mdiag used during the non-Galerkin RAP");
  validParamList->set<RCP<const FactoryBase> >("K", Teuchos::null, "Generating factory of the matrix K used during the non-Galerkin RAP");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase> >("R", Teuchos::null, "Restrictor factory");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  validParamList->set<RCP<const FactoryBase> >("deltaT", Teuchos::null, "user deltaT");
  validParamList->set<RCP<const FactoryBase> >("cfl", Teuchos::null, "user cfl");
  validParamList->set<RCP<const FactoryBase> >("cfl-based shift array", Teuchos::null, "MueLu-generated shift array for CFL-based shifting");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  const Teuchos::ParameterList &pL = GetParameterList();

  bool use_mdiag = false;
  if (pL.isParameter("rap: shift diagonal M"))
    use_mdiag = pL.get<bool>("rap: shift diagonal M");

  // The low storage version requires mdiag
  bool use_low_storage = false;
  if (pL.isParameter("rap: shift low storage")) {
    use_low_storage = pL.get<bool>("rap: shift low storage");
    use_mdiag       = use_low_storage ? true : use_mdiag;
  }

  if (implicitTranspose_ == false) {
    Input(coarseLevel, "R");
  }

  if (!use_low_storage)
    Input(fineLevel, "K");
  else
    Input(fineLevel, "A");
  Input(coarseLevel, "P");

  if (!use_mdiag)
    Input(fineLevel, "M");
  else
    Input(fineLevel, "Mdiag");

  // CFL array stuff
  if (pL.isParameter("rap: cfl array") && pL.get<Teuchos::Array<double> >("rap: cfl array").size() > 0) {
    if (fineLevel.GetLevelID() == 0) {
      if (fineLevel.IsAvailable("deltaT", NoFactory::get())) {
        fineLevel.DeclareInput("deltaT", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("fine deltaT", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "deltaT was not provided by the user on level0!");
      }

      if (fineLevel.IsAvailable("cfl", NoFactory::get())) {
        fineLevel.DeclareInput("cfl", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("fine cfl", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "cfl was not provided by the user on level0!");
      }
    } else {
      Input(fineLevel, "cfl-based shift array");
    }
  }

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
    (*it)->CallDeclareInput(coarseLevel);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {  // FIXME make fineLevel const
  {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);
    const Teuchos::ParameterList &pL = GetParameterList();

    bool M_is_diagonal = false;
    if (pL.isParameter("rap: shift diagonal M"))
      M_is_diagonal = pL.get<bool>("rap: shift diagonal M");

    // The low storage version requires mdiag
    bool use_low_storage = false;
    if (pL.isParameter("rap: shift low storage")) {
      use_low_storage = pL.get<bool>("rap: shift low storage");
      M_is_diagonal   = use_low_storage ? true : M_is_diagonal;
    }

    Teuchos::ArrayView<const double> doubleShifts;
    Teuchos::ArrayRCP<double> myshifts;
    if (pL.isParameter("rap: shift array") && pL.get<Teuchos::Array<double> >("rap: shift array").size() > 0) {
      // Do we have an array of shifts?  If so, we set doubleShifts_
      doubleShifts = pL.get<Teuchos::Array<double> >("rap: shift array")();
    }
    if (pL.isParameter("rap: cfl array") && pL.get<Teuchos::Array<double> >("rap: cfl array").size() > 0) {
      // Do we have an array of CFLs?  If so, we calculated the shifts from them.
      Teuchos::ArrayView<const double> CFLs = pL.get<Teuchos::Array<double> >("rap: cfl array")();
      if (fineLevel.GetLevelID() == 0) {
        double dt         = Get<double>(fineLevel, "deltaT");
        double cfl        = Get<double>(fineLevel, "cfl");
        double ts_at_cfl1 = dt / cfl;
        myshifts.resize(CFLs.size());
        Teuchos::Array<double> myCFLs(CFLs.size());
        myCFLs[0] = cfl;

        // Never make the CFL bigger
        for (int i = 1; i < (int)CFLs.size(); i++)
          myCFLs[i] = (CFLs[i] > cfl) ? cfl : CFLs[i];

        {
          std::ostringstream ofs;
          ofs << "RAPShiftFactory: CFL schedule = ";
          for (int i = 0; i < (int)CFLs.size(); i++)
            ofs << " " << myCFLs[i];
          GetOStream(Statistics0) << ofs.str() << std::endl;
        }
        GetOStream(Statistics0) << "RAPShiftFactory: Timestep at CFL=1 is " << ts_at_cfl1 << "    " << std::endl;

        // The shift array needs to be 1/dt
        for (int i = 0; i < (int)myshifts.size(); i++)
          myshifts[i] = 1.0 / (ts_at_cfl1 * myCFLs[i]);
        doubleShifts = myshifts();

        {
          std::ostringstream ofs;
          ofs << "RAPShiftFactory: shift schedule = ";
          for (int i = 0; i < (int)doubleShifts.size(); i++)
            ofs << " " << doubleShifts[i];
          GetOStream(Statistics0) << ofs.str() << std::endl;
        }
        Set(coarseLevel, "cfl-based shift array", myshifts);
      } else {
        myshifts     = Get<Teuchos::ArrayRCP<double> >(fineLevel, "cfl-based shift array");
        doubleShifts = myshifts();
        Set(coarseLevel, "cfl-based shift array", myshifts);
        // NOTE: If we're not on level zero, then we should have a shift array
      }
    }

    // Inputs: K, M, P
    // Note: In the low-storage case we do not keep a separate "K", we just use A
    RCP<Matrix> K;
    RCP<Matrix> M;
    RCP<Vector> Mdiag;

    if (use_low_storage)
      K = Get<RCP<Matrix> >(fineLevel, "A");
    else
      K = Get<RCP<Matrix> >(fineLevel, "K");
    if (!M_is_diagonal)
      M = Get<RCP<Matrix> >(fineLevel, "M");
    else
      Mdiag = Get<RCP<Vector> >(fineLevel, "Mdiag");

    RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P");

    // Build Kc = RKP, Mc = RMP
    RCP<Matrix> KP, MP;

    // Reuse pattern if available (multiple solve)
    // FIXME: Old style reuse doesn't work any more
    //      if (IsAvailable(coarseLevel, "AP Pattern")) {
    //        KP = Get< RCP<Matrix> >(coarseLevel, "AP Pattern");
    //        MP = Get< RCP<Matrix> >(coarseLevel, "AP Pattern");
    //      }

    {
      SubFactoryMonitor subM(*this, "MxM: K x P", coarseLevel);
      KP = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*K, false, *P, false, KP, GetOStream(Statistics2));
      if (!M_is_diagonal) {
        MP = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*M, false, *P, false, MP, GetOStream(Statistics2));
      } else {
        MP = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(P);
        MP->leftScale(*Mdiag);
      }

      Set(coarseLevel, "AP Pattern", KP);
    }

    bool doOptimizedStorage = true;

    RCP<Matrix> Ac, Kc, Mc;

    // Reuse pattern if available (multiple solve)
    //     if (IsAvailable(coarseLevel, "RAP Pattern"))
    // Ac = Get< RCP<Matrix> >(coarseLevel, "RAP Pattern");

    bool doFillComplete = true;
    if (implicitTranspose_) {
      SubFactoryMonitor m2(*this, "MxM: P' x (KP) (implicit)", coarseLevel);
      Kc = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P, true, *KP, false, Kc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
      Mc = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P, true, *MP, false, Mc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
    } else {
      RCP<Matrix> R = Get<RCP<Matrix> >(coarseLevel, "R");
      SubFactoryMonitor m2(*this, "MxM: R x (KP) (explicit)", coarseLevel);
      Kc = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*R, false, *KP, false, Kc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
      Mc = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*R, false, *MP, false, Mc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
    }

    // Get the shift
    // FIXME - We should really get rid of the shifts array and drive this the same way everything else works
    // If we're using the recursive "low storage" version, we need to shift by ( \prod_{i=1}^k shift[i] - \prod_{i=1}^{k-1} shift[i]) to
    // get the recursive relationships correct
    int level    = coarseLevel.GetLevelID();
    Scalar shift = Teuchos::ScalarTraits<Scalar>::zero();
    if (!use_low_storage) {
      // High Storage version
      if (level < (int)shifts_.size())
        shift = shifts_[level];
      else
        shift = Teuchos::as<Scalar>(pL.get<double>("rap: shift"));
    } else {
      // Low Storage Version
      if (level < (int)shifts_.size()) {
        if (level == 1)
          shift = shifts_[level];
        else {
          Scalar prod1 = Teuchos::ScalarTraits<Scalar>::one();
          for (int i = 1; i < level - 1; i++) {
            prod1 *= shifts_[i];
          }
          shift = (prod1 * shifts_[level] - prod1);
        }
      } else if (doubleShifts.size() != 0) {
        double d_shift = 0.0;
        if (level < doubleShifts.size())
          d_shift = doubleShifts[level] - doubleShifts[level - 1];

        if (d_shift < 0.0)
          GetOStream(Warnings1) << "WARNING: RAPShiftFactory has detected a negative shift... This implies a less stable coarse grid." << std::endl;
        shift = Teuchos::as<Scalar>(d_shift);
      } else {
        double base_shift = pL.get<double>("rap: shift");
        if (level == 1)
          shift = Teuchos::as<Scalar>(base_shift);
        else
          shift = Teuchos::as<Scalar>(pow(base_shift, level) - pow(base_shift, level - 1));
      }
    }
    GetOStream(Runtime0) << "RAPShiftFactory: Using shift " << shift << std::endl;

    // recombine to get K+shift*M
    {
      SubFactoryMonitor m2(*this, "Add: RKP + s*RMP", coarseLevel);
      Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*Kc, false, Teuchos::ScalarTraits<Scalar>::one(), *Mc, false, shift, Ac, GetOStream(Statistics2));
      Ac->fillComplete();
    }

    Teuchos::ArrayView<const double> relativeFloor = pL.get<Teuchos::Array<double> >("rap: relative diagonal floor")();
    if (relativeFloor.size() > 0)
      Xpetra::MatrixUtils<SC, LO, GO, NO>::RelativeDiagonalBoost(Ac, relativeFloor, GetOStream(Statistics2));

    bool repairZeroDiagonals = pL.get<bool>("RepairMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
    bool checkAc             = pL.get<bool>("CheckMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
    ;
    if (checkAc || repairZeroDiagonals)
      Xpetra::MatrixUtils<SC, LO, GO, NO>::CheckRepairMainDiagonal(Ac, repairZeroDiagonals, GetOStream(Warnings1));

    RCP<ParameterList> params = rcp(new ParameterList());
    ;
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);

    Set(coarseLevel, "A", Ac);
    // We only need K in the 'high storage' mode
    if (!use_low_storage)
      Set(coarseLevel, "K", Kc);

    if (!M_is_diagonal) {
      Set(coarseLevel, "M", Mc);
    } else {
      // If M is diagonal, then we only pass that part down the hierarchy
      // NOTE: Should we be doing some kind of rowsum instead?
      RCP<Vector> Mcv = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Mc->getRowMap(), false);
      Mc->getLocalDiagCopy(*Mcv);
      Set(coarseLevel, "Mdiag", Mcv);
    }

    //      Set(coarseLevel, "RAP Pattern", Ac);
  }

  if (transferFacts_.begin() != transferFacts_.end()) {
    SubFactoryMonitor m(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
      RCP<const FactoryBase> fac = *it;
      GetOStream(Runtime0) << "RAPShiftFactory: call transfer factory: " << fac->description() << std::endl;
      fac->CallBuild(coarseLevel);
      // AP (11/11/13): I am not sure exactly why we need to call Release, but we do need it to get rid
      // of dangling data for CoordinatesTransferFactory
      coarseLevel.Release(*fac);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase> &factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "MueLu::RAPShiftFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
  transferFacts_.push_back(factory);
}

}  // namespace MueLu

#define MUELU_RAPSHIFTFACTORY_SHORT
#endif  // MUELU_RAPSHIFTFACTORY_DEF_HPP
