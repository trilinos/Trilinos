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
#ifndef MUELU_SAPFACTORY_KOKKOS_DEF_HPP
#define MUELU_SAPFACTORY_KOKKOS_DEF_HPP

#ifdef out
#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_spmv.hpp"
#endif
#include "MueLu_SaPFactory_kokkos_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_IteratorOps.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities.hpp"

#include <sstream>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
RCP<const ParameterList> SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("sa: damping factor");
  SET_VALID_ENTRY("sa: calculate eigenvalue estimate");
  SET_VALID_ENTRY("sa: eigenvalue estimate num iterations");
  SET_VALID_ENTRY("sa: use rowsumabs diagonal scaling");
  SET_VALID_ENTRY("sa: enforce constraints");
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("sa: max eigenvalue");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase>>("P", Teuchos::null, "Tentative prolongator factory");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
void SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "A");

  // Get default tentative prolongator factory
  // Getting it that way ensure that the same factory instance will be used for both SaPFactory_kokkos and NullspaceFactory.
  RCP<const FactoryBase> initialPFact = GetFactory("P");
  if (initialPFact == Teuchos::null) {
    initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent");
  }
  coarseLevel.DeclareInput("P", initialPFact.get(), this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
void SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
void SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

  // Get default tentative prolongator factory
  // Getting it that way ensure that the same factory instance will be used for both SaPFactory_kokkos and NullspaceFactory.
  // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
  RCP<const FactoryBase> initialPFact = GetFactory("P");
  if (initialPFact == Teuchos::null) {
    initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent");
  }
  const ParameterList& pL = GetParameterList();

  // Level Get
  RCP<Matrix> A     = Get<RCP<Matrix>>(fineLevel, "A");
  RCP<Matrix> Ptent = coarseLevel.Get<RCP<Matrix>>("P", initialPFact.get());
  RCP<Matrix> finalP;
  // If Tentative facctory bailed out (e.g., number of global aggregates is 0), then SaPFactory bails
  //  This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
  if (Ptent == Teuchos::null) {
    finalP = Teuchos::null;
    Set(coarseLevel, "P", finalP);
    return;
  }

  if (restrictionMode_) {
    SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);

    A = Utilities::Transpose(*A, true);  // build transpose of A explicitly
  }

  // Build final prolongator

  // Reuse pattern if available
  RCP<ParameterList> APparams;
  if (pL.isSublist("matrixmatrix: kernel params"))
    APparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
  else
    APparams = rcp(new ParameterList);
  if (coarseLevel.IsAvailable("AP reuse data", this)) {
    GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous AP data" << std::endl;

    APparams = coarseLevel.Get<RCP<ParameterList>>("AP reuse data", this);

    if (APparams->isParameter("graph"))
      finalP = APparams->get<RCP<Matrix>>("graph");
  }
  // By default, we don't need global constants for SaP
  APparams->set("compute global constants: temporaries", APparams->get("compute global constants: temporaries", false));
  APparams->set("compute global constants", APparams->get("compute global constants", false));

  const SC dampingFactor        = as<SC>(pL.get<double>("sa: damping factor"));
  const LO maxEigenIterations   = as<LO>(pL.get<int>("sa: eigenvalue estimate num iterations"));
  const bool estimateMaxEigen   = pL.get<bool>("sa: calculate eigenvalue estimate");
  const bool useAbsValueRowSum  = pL.get<bool>("sa: use rowsumabs diagonal scaling");
  const bool doQRStep           = pL.get<bool>("tentative: calculate qr");
  const bool enforceConstraints = pL.get<bool>("sa: enforce constraints");
  const SC userDefinedMaxEigen  = as<SC>(pL.get<double>("sa: max eigenvalue"));

  // Sanity checking
  TEUCHOS_TEST_FOR_EXCEPTION(doQRStep && enforceConstraints, Exceptions::RuntimeError,
                             "MueLu::TentativePFactory::MakeTentative: cannot use 'enforce constraints' and 'calculate qr' at the same time");

  if (dampingFactor != Teuchos::ScalarTraits<SC>::zero()) {
    SC lambdaMax;
    RCP<Vector> invDiag;
    if (Teuchos::ScalarTraits<SC>::real(userDefinedMaxEigen) == Teuchos::ScalarTraits<SC>::real(-1.0)) {
      SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
      lambdaMax = A->GetMaxEigenvalueEstimate();
      if (lambdaMax == -Teuchos::ScalarTraits<SC>::one() || estimateMaxEigen) {
        GetOStream(Statistics1) << "Calculating max eigenvalue estimate now (max iters = " << maxEigenIterations << ((useAbsValueRowSum) ? ", use rowSumAbs diagonal)" : ", use point diagonal)") << std::endl;
        Magnitude stopTol = 1e-4;
        invDiag           = Utilities::GetMatrixDiagonalInverse(*A, Teuchos::ScalarTraits<SC>::eps() * 100, Teuchos::ScalarTraits<SC>::zero(), useAbsValueRowSum);
        if (useAbsValueRowSum)
          lambdaMax = Utilities::PowerMethod(*A, invDiag, maxEigenIterations, stopTol);
        else
          lambdaMax = Utilities::PowerMethod(*A, true, maxEigenIterations, stopTol);
        A->SetMaxEigenvalueEstimate(lambdaMax);
      } else {
        GetOStream(Statistics1) << "Using cached max eigenvalue estimate" << std::endl;
      }
    } else {
      lambdaMax = userDefinedMaxEigen;
      A->SetMaxEigenvalueEstimate(lambdaMax);
    }
    GetOStream(Statistics0) << "Prolongator damping factor = " << dampingFactor / lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;

    {
      SubFactoryMonitor m2(*this, "Fused (I-omega*D^{-1} A)*Ptent", coarseLevel);
      {
        SubFactoryMonitor m3(*this, "Diagonal Extraction", coarseLevel);
        if (useAbsValueRowSum)
          GetOStream(Runtime0) << "Using rowSumAbs diagonal" << std::endl;
        if (invDiag == Teuchos::null)
          invDiag = Utilities::GetMatrixDiagonalInverse(*A, Teuchos::ScalarTraits<SC>::eps() * 100, Teuchos::ScalarTraits<SC>::zero(), useAbsValueRowSum);
      }
      SC omega = dampingFactor / lambdaMax;
      TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(Teuchos::ScalarTraits<SC>::magnitude(omega)), Exceptions::RuntimeError, "Prolongator damping factor needs to be finite.");

      {
        SubFactoryMonitor m3(*this, "Xpetra::IteratorOps::Jacobi", coarseLevel);
        // finalP = Ptent + (I - \omega D^{-1}A) Ptent
        finalP = Xpetra::IteratorOps<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Jacobi(omega, *invDiag, *A, *Ptent, finalP, GetOStream(Statistics2), std::string("MueLu::SaP-") + toString(coarseLevel.GetLevelID()), APparams);
        if (enforceConstraints) SatisfyPConstraints(A, finalP);
      }
    }

  } else {
    finalP = Ptent;
  }

  // Level Set
  if (!restrictionMode_) {
    // prolongation factory is in prolongation mode
    if (!finalP.is_null()) {
      std::ostringstream oss;
      oss << "P_" << coarseLevel.GetLevelID();
      finalP->setObjectLabel(oss.str());
    }
    Set(coarseLevel, "P", finalP);

    // NOTE: EXPERIMENTAL
    if (Ptent->IsView("stridedMaps"))
      finalP->CreateView("stridedMaps", Ptent);

  } else {
    // prolongation factory is in restriction mode
    RCP<Matrix> R = Utilities::Transpose(*finalP, true);
    Set(coarseLevel, "R", R);
    if (!R.is_null()) {
      std::ostringstream oss;
      oss << "R_" << coarseLevel.GetLevelID();
      R->setObjectLabel(oss.str());
    }

    // NOTE: EXPERIMENTAL
    if (Ptent->IsView("stridedMaps"))
      R->CreateView("stridedMaps", Ptent, true);
  }

  if (IsPrint(Statistics2)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    params->set("printCommInfo", true);
    GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*finalP, (!restrictionMode_ ? "P" : "R"), params);
  }

}  // Build()

// Analyze the grid transfer produced by smoothed aggregation and make
// modifications if it does not look right. In particular, if there are
// negative entries or entries larger than 1, modify P's rows.
//
// Note: this kind of evaluation probably only makes sense if not doing QR
// when constructing tentative P.
//
// For entries that do not satisfy the two constraints (>= 0  or <=1) we set
// these entries to the constraint value and modify the rest of the row
// so that the row sum remains the same as before by adding an equal
// amount to each remaining entry. However, if the original row sum value
// violates the constraints, we set the row sum back to 1 (the row sum of
// tentative P). After doing the modification to a row, we need to check
// again the entire row to make sure that the modified row does not violate
// the constraints.

template <typename local_matrix_type>
struct constraintKernel {
  using Scalar      = typename local_matrix_type::non_const_value_type;
  using SC          = Scalar;
  using LO          = typename local_matrix_type::non_const_ordinal_type;
  using Device      = typename local_matrix_type::device_type;
  const Scalar zero = Kokkos::ArithTraits<SC>::zero();
  const Scalar one  = Kokkos::ArithTraits<SC>::one();
  LO nPDEs;
  local_matrix_type localP;
  Kokkos::View<SC**, Device> ConstraintViolationSum;
  Kokkos::View<SC**, Device> Rsum;
  Kokkos::View<size_t**, Device> nPositive;

  constraintKernel(LO nPDEs_, local_matrix_type localP_)
    : nPDEs(nPDEs_)
    , localP(localP_) {
    ConstraintViolationSum = Kokkos::View<SC**, Device>("ConstraintViolationSum", localP_.numRows(), nPDEs);
    Rsum                   = Kokkos::View<SC**, Device>("Rsum", localP_.numRows(), nPDEs);
    nPositive              = Kokkos::View<size_t**, Device>("nPositive", localP_.numRows(), nPDEs);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t rowIdx) const {
    auto rowPtr = localP.graph.row_map;
    auto values = localP.values;

    bool checkRow = true;

    if (rowPtr(rowIdx + 1) == rowPtr(rowIdx)) checkRow = false;

    while (checkRow) {
      // check constraints and compute the row sum

      for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++) {
        Rsum(rowIdx, entryIdx % nPDEs) += values(entryIdx);
        if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) < Kokkos::ArithTraits<SC>::real(zero)) {
          ConstraintViolationSum(rowIdx, entryIdx % nPDEs) += values(entryIdx);
          values(entryIdx) = zero;
        } else {
          if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) != Kokkos::ArithTraits<SC>::real(zero))
            nPositive(rowIdx, entryIdx % nPDEs) = nPositive(rowIdx, entryIdx % nPDEs) + 1;

          if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) > Kokkos::ArithTraits<SC>::real(1.00001)) {
            ConstraintViolationSum(rowIdx, entryIdx % nPDEs) += (values(entryIdx) - one);
            values(entryIdx) = one;
          }
        }
      }

      checkRow = false;

      // take into account any row sum that violates the contraints

      for (size_t k = 0; k < (size_t)nPDEs; k++) {
        if (Kokkos::ArithTraits<SC>::real(Rsum(rowIdx, k)) < Kokkos::ArithTraits<SC>::magnitude(zero)) {
          ConstraintViolationSum(rowIdx, k) = ConstraintViolationSum(rowIdx, k) - Rsum(rowIdx, k);  // rstumin
        } else if (Kokkos::ArithTraits<SC>::real(Rsum(rowIdx, k)) > Kokkos::ArithTraits<SC>::magnitude(1.00001)) {
          ConstraintViolationSum(rowIdx, k) = ConstraintViolationSum(rowIdx, k) + (one - Rsum(rowIdx, k));  // rstumin
        }
      }

      // check if row need modification
      for (size_t k = 0; k < (size_t)nPDEs; k++) {
        if (Kokkos::ArithTraits<SC>::magnitude(ConstraintViolationSum(rowIdx, k)) != Kokkos::ArithTraits<SC>::magnitude(zero))
          checkRow = true;
      }
      // modify row
      if (checkRow) {
        for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++) {
          if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) > Kokkos::ArithTraits<SC>::real(zero)) {
            values(entryIdx) = values(entryIdx) +
                               (ConstraintViolationSum(rowIdx, entryIdx % nPDEs) / (Scalar(nPositive(rowIdx, entryIdx % nPDEs)) != zero ? Scalar(nPositive(rowIdx, entryIdx % nPDEs)) : one));
          }
        }
        for (size_t k = 0; k < (size_t)nPDEs; k++) ConstraintViolationSum(rowIdx, k) = zero;
      }
      for (size_t k = 0; k < (size_t)nPDEs; k++) Rsum(rowIdx, k) = zero;
      for (size_t k = 0; k < (size_t)nPDEs; k++) nPositive(rowIdx, k) = 0;
    }  // while (checkRow) ...
  }
};

template <typename local_matrix_type>
struct optimalSatisfyConstraintsForScalarPDEsKernel {
  using Scalar      = typename local_matrix_type::non_const_value_type;
  using SC          = Scalar;
  using LO          = typename local_matrix_type::non_const_ordinal_type;
  using Device      = typename local_matrix_type::device_type;
  using KAT         = Kokkos::ArithTraits<SC>;
  const Scalar zero = Kokkos::ArithTraits<SC>::zero();
  const Scalar one  = Kokkos::ArithTraits<SC>::one();
  LO nPDEs;
  local_matrix_type localP;
  Kokkos::View<SC**, Device> origSorted;
  Kokkos::View<SC**, Device> fixedSorted;
  Kokkos::View<LO**, Device> inds;

  optimalSatisfyConstraintsForScalarPDEsKernel(LO nPDEs_, LO maxRowEntries_, local_matrix_type localP_)
    : nPDEs(nPDEs_)
    , localP(localP_) {
    origSorted  = Kokkos::View<SC**, Device>("origSorted", localP_.numRows(), maxRowEntries_);
    fixedSorted = Kokkos::View<SC**, Device>("fixedSorted", localP_.numRows(), maxRowEntries_);
    inds        = Kokkos::View<LO**, Device>("inds", localP_.numRows(), maxRowEntries_);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t rowIdx) const {
    auto rowPtr = localP.graph.row_map;
    auto values = localP.values;

    auto nnz = rowPtr(rowIdx + 1) - rowPtr(rowIdx);

    if (nnz != 0) {
      Scalar rsumTarget = zero;
      for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++) rsumTarget += values(entryIdx);

      {
        SC aBigNumber;
        SC rowSumDeviation, temp, delta;
        SC closestToLeftBoundDist, closestToRghtBoundDist;
        LO closestToLeftBound, closestToRghtBound;
        bool flipped;

        SC leftBound = zero;
        SC rghtBound = one;
        if ((KAT::real(rsumTarget) >= KAT::real(leftBound * (static_cast<SC>(nnz)))) &&
            (KAT::real(rsumTarget) <= KAT::real(rghtBound * (static_cast<SC>(nnz))))) {  // has Feasible solution

          flipped = false;
          // compute aBigNumber to handle some corner cases where we need
          // something large so that an if statement will be false
          aBigNumber = KAT::zero();
          for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++) {
            if (KAT::magnitude(values(entryIdx)) > KAT::magnitude(aBigNumber))
              aBigNumber = KAT::magnitude(values(entryIdx));
          }
          aBigNumber = aBigNumber + (KAT::magnitude(leftBound) + KAT::magnitude(rghtBound)) * (100 * one);

          LO ind = 0;
          for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++) {
            origSorted(rowIdx, ind) = values(entryIdx);
            inds(rowIdx, ind)       = ind;
            ind++;
          }

          auto sortVals = Kokkos::subview(origSorted, rowIdx, Kokkos::ALL());
          auto sortInds = Kokkos::subview(inds, rowIdx, Kokkos::make_pair(0, ind));
          // Need permutation indices to sort row values from smallest to largest,
          // and unsort once row constraints are applied.
          // serial insertion sort workaround from https://github.com/kokkos/kokkos/issues/645
          for (LO i = 1; i < static_cast<LO>(nnz); ++i) {
            ind  = sortInds(i);
            LO j = i;

            if (KAT::real(sortVals(sortInds(i))) < KAT::real(sortVals(sortInds(0)))) {
              for (; j > 0; --j) sortInds(j) = sortInds(j - 1);

              sortInds(0) = ind;
            } else {
              for (; KAT::real(sortVals(ind)) < KAT::real(sortVals(sortInds(j - 1))); --j) sortInds(j) = sortInds(j - 1);

              sortInds(j) = ind;
            }
          }

          for (LO i = 0; i < static_cast<LO>(nnz); i++) origSorted(rowIdx, i) = values(rowPtr(rowIdx) + inds(rowIdx, i));  // values is no longer used
          // find entry in origSorted just to the right of the leftBound
          closestToLeftBound = 0;
          while ((closestToLeftBound < static_cast<LO>(nnz)) && (KAT::real(origSorted(rowIdx, closestToLeftBound)) <= KAT::real(leftBound))) closestToLeftBound++;

          // find entry in origSorted just to the right of the rghtBound
          closestToRghtBound = closestToLeftBound;
          while ((closestToRghtBound < static_cast<LO>(nnz)) && (KAT::real(origSorted(rowIdx, closestToRghtBound)) <= KAT::real(rghtBound))) closestToRghtBound++;

          // compute distance between closestToLeftBound and the left bound and the
          // distance between closestToRghtBound and the right bound.

          closestToLeftBoundDist = origSorted(rowIdx, closestToLeftBound) - leftBound;
          if (closestToRghtBound == static_cast<LO>(nnz))
            closestToRghtBoundDist = aBigNumber;
          else
            closestToRghtBoundDist = origSorted(rowIdx, closestToRghtBound) - rghtBound;

          // compute how far the rowSum is off from the target row sum taking into account
          // numbers that have been shifted to satisfy bound constraint

          rowSumDeviation = leftBound * (static_cast<SC>(closestToLeftBound)) + (static_cast<SC>(nnz - closestToRghtBound)) * rghtBound - rsumTarget;
          for (LO i = closestToLeftBound; i < closestToRghtBound; i++) rowSumDeviation += origSorted(rowIdx, i);

          // the code that follow after this if statement assumes that rowSumDeviation is positive. If this
          // is not the case, flip the signs of everything so that rowSumDeviation is now positive.
          // Later we will flip the data back to its original form.
          if (KAT::real(rowSumDeviation) < KAT::real(KAT::zero())) {
            flipped   = true;
            temp      = leftBound;
            leftBound = -rghtBound;
            rghtBound = temp;

            /* flip sign of origSorted and reverse ordering so that the negative version is sorted */

            // TODO: the following is bad for GPU performance. Switch to bit shifting if brave.
            if ((nnz % 2) == 1) origSorted(rowIdx, (nnz / 2)) = -origSorted(rowIdx, (nnz / 2));
            for (LO i = 0; i < static_cast<LO>(nnz / 2); i++) {
              temp                            = origSorted(rowIdx, i);
              origSorted(rowIdx, i)           = -origSorted(rowIdx, nnz - 1 - i);
              origSorted(rowIdx, nnz - i - 1) = -temp;
            }

            /* reverse bounds */

            LO itemp               = closestToLeftBound;
            closestToLeftBound     = nnz - closestToRghtBound;
            closestToRghtBound     = nnz - itemp;
            closestToLeftBoundDist = origSorted(rowIdx, closestToLeftBound) - leftBound;
            if (closestToRghtBound == static_cast<LO>(nnz))
              closestToRghtBoundDist = aBigNumber;
            else
              closestToRghtBoundDist = origSorted(rowIdx, closestToRghtBound) - rghtBound;

            rowSumDeviation = -rowSumDeviation;
          }

          // initial fixedSorted so bounds are satisfied and interiors correspond to origSorted

          for (LO i = 0; i < closestToLeftBound; i++) fixedSorted(rowIdx, i) = leftBound;
          for (LO i = closestToLeftBound; i < closestToRghtBound; i++) fixedSorted(rowIdx, i) = origSorted(rowIdx, i);
          for (LO i = closestToRghtBound; i < static_cast<LO>(nnz); i++) fixedSorted(rowIdx, i) = rghtBound;

          while ((KAT::magnitude(rowSumDeviation) > KAT::magnitude((one * 1.e-10) * rsumTarget))) {  // && ( (closestToLeftBound < nEntries ) || (closestToRghtBound < nEntries))) {
            if (closestToRghtBound != closestToLeftBound)
              delta = rowSumDeviation / static_cast<SC>(closestToRghtBound - closestToLeftBound);
            else
              delta = aBigNumber;

            if (KAT::magnitude(closestToLeftBoundDist) <= KAT::magnitude(closestToRghtBoundDist)) {
              if (KAT::magnitude(delta) <= KAT::magnitude(closestToLeftBoundDist)) {
                rowSumDeviation = zero;
                for (LO i = closestToLeftBound; i < closestToRghtBound; i++) fixedSorted(rowIdx, i) = origSorted(rowIdx, i) - delta;
              } else {
                rowSumDeviation                         = rowSumDeviation - closestToLeftBoundDist;
                fixedSorted(rowIdx, closestToLeftBound) = leftBound;
                closestToLeftBound++;
                if (closestToLeftBound < static_cast<LO>(nnz))
                  closestToLeftBoundDist = origSorted(rowIdx, closestToLeftBound) - leftBound;
                else
                  closestToLeftBoundDist = aBigNumber;
              }
            } else {
              if (KAT::magnitude(delta) <= KAT::magnitude(closestToRghtBoundDist)) {
                rowSumDeviation = 0;
                for (LO i = closestToLeftBound; i < closestToRghtBound; i++) fixedSorted(rowIdx, i) = origSorted(rowIdx, i) - delta;
              } else {
                rowSumDeviation = rowSumDeviation + closestToRghtBoundDist;
                //        if (closestToRghtBound < nEntries) {
                fixedSorted(rowIdx, closestToRghtBound) = origSorted(rowIdx, closestToRghtBound);
                closestToRghtBound++;
                //       }
                if (closestToRghtBound >= static_cast<LO>(nnz))
                  closestToRghtBoundDist = aBigNumber;
                else
                  closestToRghtBoundDist = origSorted(rowIdx, closestToRghtBound) - rghtBound;
              }
            }
          }

          auto rowStart = rowPtr(rowIdx);
          if (flipped) {
            /* flip sign of fixedSorted and reverse ordering so that the positve version is sorted */

            if ((nnz % 2) == 1) fixedSorted(rowIdx, (nnz / 2)) = -fixedSorted(rowIdx, (nnz / 2));
            for (LO i = 0; i < static_cast<LO>(nnz / 2); i++) {
              temp                             = fixedSorted(rowIdx, i);
              fixedSorted(rowIdx, i)           = -fixedSorted(rowIdx, nnz - 1 - i);
              fixedSorted(rowIdx, nnz - i - 1) = -temp;
            }
          }
          // unsort and update row values with new values
          for (LO i = 0; i < static_cast<LO>(nnz); i++) values(rowStart + inds(rowIdx, i)) = fixedSorted(rowIdx, i);

        } else {  // row does not have feasible solution to match constraint
          // just set all entries to the same value giving a row sum of 1
          for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++) values(entryIdx) = one / (static_cast<SC>(nnz));
        }
      }
    }
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
void SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::SatisfyPConstraints(const RCP<Matrix> A, RCP<Matrix>& P) const {
  using Device = typename Matrix::local_matrix_type::device_type;
  LO nPDEs     = A->GetFixedBlockSize();

  using local_mat_type = typename Matrix::local_matrix_type;
  constraintKernel<local_mat_type> myKernel(nPDEs, P->getLocalMatrixDevice());
  Kokkos::parallel_for("enforce constraint", Kokkos::RangePolicy<typename Device::execution_space>(0, P->getRowMap()->getLocalNumElements()),
                       myKernel);

}  // SatsifyPConstraints()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
void SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::optimalSatisfyPConstraintsForScalarPDEs(RCP<Matrix>& P) const {
  using Device = typename Matrix::local_matrix_type::device_type;
  LO nPDEs     = 1;  // A->GetFixedBlockSize();

  using local_mat_type = typename Matrix::local_matrix_type;
  optimalSatisfyConstraintsForScalarPDEsKernel<local_mat_type> myKernel(nPDEs, P->getLocalMaxNumRowEntries(), P->getLocalMatrixDevice());
  Kokkos::parallel_for("enforce constraint", Kokkos::RangePolicy<typename Device::execution_space>(0, P->getLocalNumRows()),
                       myKernel);

}  // SatsifyPConstraints()

}  // namespace MueLu

#endif  // MUELU_SAPFACTORY_KOKKOS_DEF_HPP

// TODO: restrictionMode_ should use the parameter list.
