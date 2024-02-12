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
#ifndef MUELU_SAPFACTORY_DEF_HPP
#define MUELU_SAPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_IteratorOps.hpp>  // containing routines to generate Jacobi iterator
#include <Xpetra_IO.hpp>
#include <sstream>

#include "MueLu_SaPFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("sa: damping factor");
  SET_VALID_ENTRY("sa: calculate eigenvalue estimate");
  SET_VALID_ENTRY("sa: eigenvalue estimate num iterations");
  SET_VALID_ENTRY("sa: use rowsumabs diagonal scaling");
  SET_VALID_ENTRY("sa: enforce constraints");
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("sa: max eigenvalue");
  SET_VALID_ENTRY("sa: rowsumabs diagonal replacement tolerance");
  SET_VALID_ENTRY("sa: rowsumabs diagonal replacement value");
  SET_VALID_ENTRY("sa: rowsumabs replace single entry row with zero");
  SET_VALID_ENTRY("sa: rowsumabs use automatic diagonal tolerance");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Tentative prolongator factory");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(fineLevel, "A");

  // Get default tentative prolongator factory
  // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
  RCP<const FactoryBase> initialPFact = GetFactory("P");
  if (initialPFact == Teuchos::null) {
    initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent");
  }
  coarseLevel.DeclareInput("P", initialPFact.get(), this);  // --
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

  std::string levelIDs = toString(coarseLevel.GetLevelID());

  const std::string prefix = "MueLu::SaPFactory(" + levelIDs + "): ";

  typedef typename Teuchos::ScalarTraits<SC>::coordinateType Coordinate;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;

  // Get default tentative prolongator factory
  // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
  // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
  RCP<const FactoryBase> initialPFact = GetFactory("P");
  if (initialPFact == Teuchos::null) {
    initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent");
  }
  const ParameterList &pL = GetParameterList();

  // Level Get
  RCP<Matrix> A     = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<Matrix> Ptent = coarseLevel.Get<RCP<Matrix> >("P", initialPFact.get());
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
    A = Utilities::Transpose(*A, true);  // build transpose of A explicitely
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

    APparams = coarseLevel.Get<RCP<ParameterList> >("AP reuse data", this);

    if (APparams->isParameter("graph"))
      finalP = APparams->get<RCP<Matrix> >("graph");
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
  const MT userDefinedMaxEigen  = as<MT>(pL.get<double>("sa: max eigenvalue"));
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  double dTol                                  = pL.get<double>("sa: rowsumabs diagonal replacement tolerance");
  const Magnitude diagonalReplacementTolerance = (dTol == as<double>(-1) ? Teuchos::ScalarTraits<Scalar>::eps() * 100 : as<Magnitude>(pL.get<double>("sa: rowsumabs diagonal replacement tolerance")));
  const SC diagonalReplacementValue            = as<SC>(pL.get<double>("sa: rowsumabs diagonal replacement value"));
  const bool replaceSingleEntryRowWithZero     = pL.get<bool>("sa: rowsumabs replace single entry row with zero");
  const bool useAutomaticDiagTol               = pL.get<bool>("sa: rowsumabs use automatic diagonal tolerance");

  // Sanity checking
  TEUCHOS_TEST_FOR_EXCEPTION(doQRStep && enforceConstraints, Exceptions::RuntimeError,
                             "MueLu::TentativePFactory::MakeTentative: cannot use 'enforce constraints' and 'calculate qr' at the same time");

  if (dampingFactor != Teuchos::ScalarTraits<SC>::zero()) {
    Scalar lambdaMax;
    Teuchos::RCP<Vector> invDiag;
    if (userDefinedMaxEigen == -1.) {
      SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
      lambdaMax = A->GetMaxEigenvalueEstimate();
      if (lambdaMax == -Teuchos::ScalarTraits<SC>::one() || estimateMaxEigen) {
        GetOStream(Statistics1) << "Calculating max eigenvalue estimate now (max iters = " << maxEigenIterations << ((useAbsValueRowSum) ? ", use rowSumAbs diagonal)" : ", use point diagonal)") << std::endl;
        Coordinate stopTol = 1e-4;
        if (useAbsValueRowSum) {
          const bool returnReciprocal = true;
          invDiag                     = Utilities::GetLumpedMatrixDiagonal(*A, returnReciprocal,
                                                                           diagonalReplacementTolerance,
                                                                           diagonalReplacementValue,
                                                                           replaceSingleEntryRowWithZero,
                                                                           useAutomaticDiagTol);
          TEUCHOS_TEST_FOR_EXCEPTION(invDiag.is_null(), Exceptions::RuntimeError,
                                     "SaPFactory: eigenvalue estimate: diagonal reciprocal is null.");
          lambdaMax = Utilities::PowerMethod(*A, invDiag, maxEigenIterations, stopTol);
        } else
          lambdaMax = Utilities::PowerMethod(*A, true, maxEigenIterations, stopTol);
        A->SetMaxEigenvalueEstimate(lambdaMax);
      } else {
        GetOStream(Statistics1) << "Using cached max eigenvalue estimate" << std::endl;
      }
    } else {
      lambdaMax = userDefinedMaxEigen;
      A->SetMaxEigenvalueEstimate(lambdaMax);
    }
    GetOStream(Statistics1) << "Prolongator damping factor = " << dampingFactor / lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;

    {
      SubFactoryMonitor m2(*this, "Fused (I-omega*D^{-1} A)*Ptent", coarseLevel);
      if (!useAbsValueRowSum)
        invDiag = Utilities::GetMatrixDiagonalInverse(*A);  // default
      else if (invDiag == Teuchos::null) {
        GetOStream(Runtime0) << "Using rowsumabs diagonal" << std::endl;
        const bool returnReciprocal = true;
        invDiag                     = Utilities::GetLumpedMatrixDiagonal(*A, returnReciprocal,
                                                                         diagonalReplacementTolerance,
                                                                         diagonalReplacementValue,
                                                                         replaceSingleEntryRowWithZero,
                                                                         useAutomaticDiagTol);
        TEUCHOS_TEST_FOR_EXCEPTION(invDiag.is_null(), Exceptions::RuntimeError, "SaPFactory: diagonal reciprocal is null.");
      }

      SC omega = dampingFactor / lambdaMax;
      TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(Teuchos::ScalarTraits<SC>::magnitude(omega)), Exceptions::RuntimeError, "Prolongator damping factor needs to be finite.");

      // finalP = Ptent + (I - \omega D^{-1}A) Ptent
      finalP = Xpetra::IteratorOps<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Jacobi(omega, *invDiag, *A, *Ptent, finalP,
                                                                                      GetOStream(Statistics2), std::string("MueLu::SaP-") + levelIDs, APparams);
      if (enforceConstraints) {
        if (A->GetFixedBlockSize() == 1)
          optimalSatisfyPConstraintsForScalarPDEs(finalP);
        else
          SatisfyPConstraints(A, finalP);
      }
    }

  } else {
    finalP = Ptent;
  }

  // Level Set
  if (!restrictionMode_) {
    // The factory is in prolongation mode
    if (!finalP.is_null()) {
      std::ostringstream oss;
      oss << "P_" << coarseLevel.GetLevelID();
      finalP->setObjectLabel(oss.str());
    }
    Set(coarseLevel, "P", finalP);

    APparams->set("graph", finalP);
    Set(coarseLevel, "AP reuse data", APparams);

    // NOTE: EXPERIMENTAL
    if (Ptent->IsView("stridedMaps"))
      finalP->CreateView("stridedMaps", Ptent);

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*finalP, "P", params);
    }

  } else {
    // The factory is in restriction mode
    RCP<Matrix> R;
    {
      SubFactoryMonitor m2(*this, "Transpose P", coarseLevel);
      R = Utilities::Transpose(*finalP, true);
      if (!R.is_null()) {
        std::ostringstream oss;
        oss << "R_" << coarseLevel.GetLevelID();
        R->setObjectLabel(oss.str());
      }
    }

    Set(coarseLevel, "R", R);

    // NOTE: EXPERIMENTAL
    if (Ptent->IsView("stridedMaps"))
      R->CreateView("stridedMaps", Ptent, true /*transposeA*/);

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*R, "R", params);
    }
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SatisfyPConstraints(const RCP<Matrix> A, RCP<Matrix> &P) const {
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
  LO nPDEs          = A->GetFixedBlockSize();
  Teuchos::ArrayRCP<Scalar> ConstraintViolationSum(nPDEs);
  Teuchos::ArrayRCP<Scalar> Rsum(nPDEs);
  Teuchos::ArrayRCP<size_t> nPositive(nPDEs);
  for (size_t k = 0; k < (size_t)nPDEs; k++) ConstraintViolationSum[k] = zero;
  for (size_t k = 0; k < (size_t)nPDEs; k++) Rsum[k] = zero;
  for (size_t k = 0; k < (size_t)nPDEs; k++) nPositive[k] = 0;

  for (size_t i = 0; i < as<size_t>(P->getRowMap()->getLocalNumElements()); i++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals1;
    Teuchos::ArrayView<Scalar> vals;
    P->getLocalRowView((LO)i, indices, vals1);
    size_t nnz = indices.size();
    if (nnz == 0) continue;

    vals = ArrayView<Scalar>(const_cast<SC *>(vals1.getRawPtr()), nnz);

    bool checkRow = true;

    while (checkRow) {
      // check constraints and compute the row sum

      for (LO j = 0; j < indices.size(); j++) {
        Rsum[j % nPDEs] += vals[j];
        if (Teuchos::ScalarTraits<SC>::real(vals[j]) < Teuchos::ScalarTraits<SC>::real(zero)) {
          ConstraintViolationSum[j % nPDEs] += vals[j];
          vals[j] = zero;
        } else {
          if (Teuchos::ScalarTraits<SC>::real(vals[j]) != Teuchos::ScalarTraits<SC>::real(zero))
            (nPositive[j % nPDEs])++;

          if (Teuchos::ScalarTraits<SC>::real(vals[j]) > Teuchos::ScalarTraits<SC>::real(1.00001)) {
            ConstraintViolationSum[j % nPDEs] += (vals[j] - one);
            vals[j] = one;
          }
        }
      }

      checkRow = false;

      // take into account any row sum that violates the contraints

      for (size_t k = 0; k < (size_t)nPDEs; k++) {
        if (Teuchos::ScalarTraits<SC>::real(Rsum[k]) < Teuchos::ScalarTraits<SC>::magnitude(zero)) {
          ConstraintViolationSum[k] += (-Rsum[k]);  // rstumin
        } else if (Teuchos::ScalarTraits<SC>::real(Rsum[k]) > Teuchos::ScalarTraits<SC>::magnitude(1.00001)) {
          ConstraintViolationSum[k] += (one - Rsum[k]);  // rstumin
        }
      }

      // check if row need modification
      for (size_t k = 0; k < (size_t)nPDEs; k++) {
        if (Teuchos::ScalarTraits<SC>::magnitude(ConstraintViolationSum[k]) != Teuchos::ScalarTraits<SC>::magnitude(zero))
          checkRow = true;
      }
      // modify row
      if (checkRow) {
        for (LO j = 0; j < indices.size(); j++) {
          if (Teuchos::ScalarTraits<SC>::real(vals[j]) > Teuchos::ScalarTraits<SC>::real(zero)) {
            vals[j] += (ConstraintViolationSum[j % nPDEs] / (as<Scalar>(nPositive[j % nPDEs])));
          }
        }
        for (size_t k = 0; k < (size_t)nPDEs; k++) ConstraintViolationSum[k] = zero;
      }
      for (size_t k = 0; k < (size_t)nPDEs; k++) Rsum[k] = zero;
      for (size_t k = 0; k < (size_t)nPDEs; k++) nPositive[k] = 0;
    }  // while (checkRow) ...
  }    // for (size_t i = 0; i < as<size_t>(P->getRowMap()->getNumNodeElements()); i++) ...
}  // SatsifyPConstraints()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::optimalSatisfyPConstraintsForScalarPDEs(RCP<Matrix> &P) const {
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();

  LocalOrdinal maxEntriesPerRow = 100;  // increased later if needed
  Teuchos::ArrayRCP<Scalar> scalarData(3 * maxEntriesPerRow);
  bool hasFeasible;

  for (size_t i = 0; i < as<size_t>(P->getRowMap()->getLocalNumElements()); i++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals1;
    Teuchos::ArrayView<Scalar> vals;
    P->getLocalRowView((LocalOrdinal)i, indices, vals1);
    size_t nnz = indices.size();
    if (nnz != 0) {
      vals              = ArrayView<Scalar>(const_cast<SC *>(vals1.getRawPtr()), nnz);
      Scalar rsumTarget = zero;
      for (size_t j = 0; j < nnz; j++) rsumTarget += vals[j];

      if (nnz > as<size_t>(maxEntriesPerRow)) {
        maxEntriesPerRow = nnz * 3;
        scalarData.resize(3 * maxEntriesPerRow);
      }
      hasFeasible = constrainRow(vals.getRawPtr(), as<LocalOrdinal>(nnz), zero, one, rsumTarget, vals.getRawPtr(), scalarData.getRawPtr());

      if (!hasFeasible) {  // just set all entries to the same value giving a row sum of 1
        for (size_t j = 0; j < nnz; j++) vals[j] = one / as<Scalar>(nnz);
      }
    }

  }  // for (size_t i = 0; i < as<size_t>(P->getRowMap()->getNumNodeElements()); i++) ...
}  // SatsifyPConstraints()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::constrainRow(Scalar *orig, LocalOrdinal nEntries, Scalar leftBound, Scalar rghtBound, Scalar rsumTarget, Scalar *fixedUnsorted, Scalar *scalarData) const {
  /*
     Input
       orig          data that should be adjusted to satisfy bound constraints and
                     row sum constraint. orig is not modified by this function.

       nEntries      length or 'orig'

       leftBound,    define bound constraints for the results.
       rghtBound

       rsumTarget    defines an equality constraint for the row sum

       fixedUnsorted  on output, if a feasible solutuion exists then
                      || orig - fixedUnsorted || = min  when also
                      leftBound <= fixedUnsorted[i] <= rghtBound for all i
                      and  sum(fixedUnsorted) = rsumTarget.

                      Note: it is possible to use the same pointer for
                      fixedUnsorted and orig. In this case, orig gets
                      overwritten with the new constraint satisfying values.

        scalarData    a work array that should be 3x nEntries.

     On return constrain() indicates whether or not a feasible solution exists.
  */

  /*
    Given a sequence of numbers    o1  ... on, fix these so that they are as
    close as possible to the original but satisfy bound constraints and also
    have the same row sum as the oi's. If we know who is going to lie on a
    bound, then the "best" answer (i.e., || o - f ||_2 = min)  perturbs
    each element that doesn't lie on a bound by the same amount.

    We can represent the oi's by considering scattered points on a number line

                   |                                                   |
                   |                                                   |
     o      o    o |  o        o     o           o       o     o       |o       o
                   |                                                   |

                   \____/                                              \____/
                   <----                                               <----
                   delta                                               delta

    Bounds are shown by vertical lines. The fi's must lie within the bounds, so
    the 3 leftmost points must be shifted to the right and the 2 rightmost must
    be shifted to the left. If these fi points are all shifted to the bounds
    while the others remain in the same location, the row sum constraint is
    likely not satisfied and so more shifting is necessary.  In the figure, the f
    rowsum is too large and so there must be more shifting to the left.

    To minimize || o - f ||_2, we basically shift all "interiors" by the same
    amount, denoted delta. The only trick is that some points near bounds are
    still affected by the bounds and so these points might be shifted more or less
    than delta. In the example,t he 3 rightmost points are shifted in the opposite
    direction as delta to the bound. The 4th point is shifted by something less
    than delta so it does not violate the lower bound. The rightmost point is
    shifted to the bound by some amount larger than delta. However, the 2nd point
    is shifted by delta (i.e., it lies inside the two bounds).

    If we know delta, we can figure out everything. If we know which points
    are special (not shifted by delta), we can also figure out everything.
    The problem is these two things (delta and the special points) are
    inter-connected. An algorithm for computing follows.

    1) move exterior points to the bounds and compute how much the row sum is off
       (rowSumDeviation). We assume now that the new row sum is high, so interior
       points must be shifted left.

    2) Mark closest point just left of the leftmost bound, closestToLeftBound,
       and compute its distance to the leftmost bound.  Mark closest point to the
       left of the rightmost bound, closestToRghtBound, and compute its distance to
       right bound. There are two cases to consider.

    3) Case 1: closestToLeftBound is closer than closestToRghtBound.
       Assume that shifting by delta does not move closestToLeftBound past the
       left bound. This means that it will be shifted by delta. However,
       closestToRghtBound will be shifted by more than delta.  So the total
       number of points shifted by delta (|interiorToBounds|) includes
       closestToLeftBound up to and including the point just to the left of
       closestToRghtBound. So

                    delta = rowSumDeviation/ |interiorToBounds| .

       Recall that rowSumDeviation already accounts for the non-delta shift of
       of closestToRightBound.  Now check whether our assumption is valid.

       If delta <= closestToLeftBoundDist, assumption is true so delta can be
       applied to interiorToBounds ... and we are done.
       Else assumption is false. Shift closestToLeftBound to the left bound.
       Update rowSumDeviation, interiorToBounds, and identify new
       closestToLeftBound. Repeat step 3).

       Case 2: closestToRghtBound is closer than closestToLeftBound.
       Assume that shifting by delta does not move closestToRghtBound past right
       bound. This means that it must be shifted by more than delta to right
       bound.  So the total number of points shifted by delta again includes
       closestToLeftBound up to and including the point just to the left of
       closestToRghtBound. So again compute

                    delta = rowSumDeviation/ |interiorToBounds| .

       If delta <= closestToRghtBoundDist, assumption is true so delta is
       can be applied to interiorToBounds ... and we are done
       Else  assumption is false. Put closestToRghtBound in the
       interiorToBounds set. Remove it's contribution to rowSumDeviation,
       identify new closestToRghtBound. Repeat step 3)


    To implement, sort the oi's so things like closestToLeftBound is just index
    into sorted array.  Updaing closestToLeftBound or closestToRghtBound requires
    increment by 1.  |interiorToBounds|= closestToRghtBound -  closestToLeftBound
    To handle the case when the rowsum is low (requiring right interior shifts),
    just flip the signs on data and use the left-shift code (and then flip back
    before exiting function.
  */
  bool hasFeasibleSol;
  Scalar notFlippedLeftBound, notFlippedRghtBound, aBigNumber, *origSorted;
  Scalar rowSumDeviation, temp, *fixedSorted, delta;
  Scalar closestToLeftBoundDist, closestToRghtBoundDist;
  LocalOrdinal closestToLeftBound, closestToRghtBound;
  LocalOrdinal *inds;
  bool flipped;

  notFlippedLeftBound = leftBound;
  notFlippedRghtBound = rghtBound;

  if ((Teuchos::ScalarTraits<SC>::real(rsumTarget) >= Teuchos::ScalarTraits<SC>::real(leftBound * as<Scalar>(nEntries))) &&
      (Teuchos::ScalarTraits<SC>::real(rsumTarget) <= Teuchos::ScalarTraits<SC>::real(rghtBound * as<Scalar>(nEntries))))
    hasFeasibleSol = true;
  else {
    hasFeasibleSol = false;
    return hasFeasibleSol;
  }
  flipped = false;
  // compute aBigNumber to handle some corner cases where we need
  // something large so that an if statement will be false
  aBigNumber = Teuchos::ScalarTraits<SC>::zero();
  for (LocalOrdinal i = 0; i < nEntries; i++) {
    if (Teuchos::ScalarTraits<SC>::magnitude(orig[i]) > Teuchos::ScalarTraits<SC>::magnitude(aBigNumber))
      aBigNumber = Teuchos::ScalarTraits<SC>::magnitude(orig[i]);
  }
  aBigNumber = aBigNumber + (Teuchos::ScalarTraits<SC>::magnitude(leftBound) + Teuchos::ScalarTraits<SC>::magnitude(rghtBound)) * as<Scalar>(100.0);

  origSorted  = &scalarData[0];
  fixedSorted = &(scalarData[nEntries]);
  inds        = (LocalOrdinal *)&(scalarData[2 * nEntries]);

  for (LocalOrdinal i = 0; i < nEntries; i++) inds[i] = i;
  for (LocalOrdinal i = 0; i < nEntries; i++) origSorted[i] = orig[i]; /* orig no longer used */

  // sort so that  orig[inds] is sorted.
  std::sort(inds, inds + nEntries,
            [origSorted](LocalOrdinal leftIndex, LocalOrdinal rightIndex) { return Teuchos::ScalarTraits<SC>::real(origSorted[leftIndex]) < Teuchos::ScalarTraits<SC>::real(origSorted[rightIndex]); });

  for (LocalOrdinal i = 0; i < nEntries; i++) origSorted[i] = orig[inds[i]];
  // find entry in origSorted just to the right of the leftBound
  closestToLeftBound = 0;
  while ((closestToLeftBound < nEntries) && (Teuchos::ScalarTraits<SC>::real(origSorted[closestToLeftBound]) <= Teuchos::ScalarTraits<SC>::real(leftBound))) closestToLeftBound++;

  // find entry in origSorted just to the right of the rghtBound
  closestToRghtBound = closestToLeftBound;
  while ((closestToRghtBound < nEntries) && (Teuchos::ScalarTraits<SC>::real(origSorted[closestToRghtBound]) <= Teuchos::ScalarTraits<SC>::real(rghtBound))) closestToRghtBound++;

  // compute distance between closestToLeftBound and the left bound and the
  // distance between closestToRghtBound and the right bound.

  closestToLeftBoundDist = origSorted[closestToLeftBound] - leftBound;
  if (closestToRghtBound == nEntries)
    closestToRghtBoundDist = aBigNumber;
  else
    closestToRghtBoundDist = origSorted[closestToRghtBound] - rghtBound;

  // compute how far the rowSum is off from the target row sum taking into account
  // numbers that have been shifted to satisfy bound constraint

  rowSumDeviation = leftBound * as<Scalar>(closestToLeftBound) + as<Scalar>((nEntries - closestToRghtBound)) * rghtBound - rsumTarget;
  for (LocalOrdinal i = closestToLeftBound; i < closestToRghtBound; i++) rowSumDeviation += origSorted[i];

  // the code that follow after this if statement assumes that rowSumDeviation is positive. If this
  // is not the case, flip the signs of everything so that rowSumDeviation is now positive.
  // Later we will flip the data back to its original form.
  if (Teuchos::ScalarTraits<SC>::real(rowSumDeviation) < Teuchos::ScalarTraits<SC>::real(Teuchos::ScalarTraits<SC>::zero())) {
    flipped   = true;
    temp      = leftBound;
    leftBound = -rghtBound;
    rghtBound = temp;

    /* flip sign of origSorted and reverse ordering so that the negative version is sorted */

    if ((nEntries % 2) == 1) origSorted[(nEntries / 2)] = -origSorted[(nEntries / 2)];
    for (LocalOrdinal i = 0; i < nEntries / 2; i++) {
      temp                         = origSorted[i];
      origSorted[i]                = -origSorted[nEntries - 1 - i];
      origSorted[nEntries - i - 1] = -temp;
    }

    /* reverse bounds */

    LocalOrdinal itemp     = closestToLeftBound;
    closestToLeftBound     = nEntries - closestToRghtBound;
    closestToRghtBound     = nEntries - itemp;
    closestToLeftBoundDist = origSorted[closestToLeftBound] - leftBound;
    if (closestToRghtBound == nEntries)
      closestToRghtBoundDist = aBigNumber;
    else
      closestToRghtBoundDist = origSorted[closestToRghtBound] - rghtBound;

    rowSumDeviation = -rowSumDeviation;
  }

  // initial fixedSorted so bounds are satisfied and interiors correspond to origSorted

  for (LocalOrdinal i = 0; i < closestToLeftBound; i++) fixedSorted[i] = leftBound;
  for (LocalOrdinal i = closestToLeftBound; i < closestToRghtBound; i++) fixedSorted[i] = origSorted[i];
  for (LocalOrdinal i = closestToRghtBound; i < nEntries; i++) fixedSorted[i] = rghtBound;

  while ((Teuchos::ScalarTraits<SC>::magnitude(rowSumDeviation) > Teuchos::ScalarTraits<SC>::magnitude(as<Scalar>(1.e-10) * rsumTarget))) {  // && ( (closestToLeftBound < nEntries ) || (closestToRghtBound < nEntries))) {
    if (closestToRghtBound != closestToLeftBound)
      delta = rowSumDeviation / as<Scalar>(closestToRghtBound - closestToLeftBound);
    else
      delta = aBigNumber;

    if (Teuchos::ScalarTraits<SC>::magnitude(closestToLeftBoundDist) <= Teuchos::ScalarTraits<SC>::magnitude(closestToRghtBoundDist)) {
      if (Teuchos::ScalarTraits<SC>::magnitude(delta) <= Teuchos::ScalarTraits<SC>::magnitude(closestToLeftBoundDist)) {
        rowSumDeviation = Teuchos::ScalarTraits<SC>::zero();
        for (LocalOrdinal i = closestToLeftBound; i < closestToRghtBound; i++) fixedSorted[i] = origSorted[i] - delta;
      } else {
        rowSumDeviation                 = rowSumDeviation - closestToLeftBoundDist;
        fixedSorted[closestToLeftBound] = leftBound;
        closestToLeftBound++;
        if (closestToLeftBound < nEntries)
          closestToLeftBoundDist = origSorted[closestToLeftBound] - leftBound;
        else
          closestToLeftBoundDist = aBigNumber;
      }
    } else {
      if (Teuchos::ScalarTraits<SC>::magnitude(delta) <= Teuchos::ScalarTraits<SC>::magnitude(closestToRghtBoundDist)) {
        rowSumDeviation = 0;
        for (LocalOrdinal i = closestToLeftBound; i < closestToRghtBound; i++) fixedSorted[i] = origSorted[i] - delta;
      } else {
        rowSumDeviation = rowSumDeviation + closestToRghtBoundDist;
        //        if (closestToRghtBound < nEntries) {
        fixedSorted[closestToRghtBound] = origSorted[closestToRghtBound];
        closestToRghtBound++;
        //       }
        if (closestToRghtBound >= nEntries)
          closestToRghtBoundDist = aBigNumber;
        else
          closestToRghtBoundDist = origSorted[closestToRghtBound] - rghtBound;
      }
    }
  }

  if (flipped) {
    /* flip sign of fixedSorted and reverse ordering so that the positve version is sorted */

    if ((nEntries % 2) == 1) fixedSorted[(nEntries / 2)] = -fixedSorted[(nEntries / 2)];
    for (LocalOrdinal i = 0; i < nEntries / 2; i++) {
      temp                          = fixedSorted[i];
      fixedSorted[i]                = -fixedSorted[nEntries - 1 - i];
      fixedSorted[nEntries - i - 1] = -temp;
    }
  }
  for (LocalOrdinal i = 0; i < nEntries; i++) fixedUnsorted[inds[i]] = fixedSorted[i];

  /* check that no constraints are violated */

  bool lowerViolation = false;
  bool upperViolation = false;
  bool sumViolation   = false;
  using TST           = Teuchos::ScalarTraits<SC>;
  temp                = TST::zero();
  for (LocalOrdinal i = 0; i < nEntries; i++) {
    if (TST::real(fixedUnsorted[i]) < TST::real(notFlippedLeftBound)) lowerViolation = true;
    if (TST::real(fixedUnsorted[i]) > TST::real(notFlippedRghtBound)) upperViolation = true;
    temp += fixedUnsorted[i];
  }
  SC tol = as<Scalar>(std::max(1.0e-8, as<double>(100 * TST::eps())));
  if (TST::magnitude(temp - rsumTarget) > TST::magnitude(tol * rsumTarget)) sumViolation = true;

  TEUCHOS_TEST_FOR_EXCEPTION(lowerViolation, Exceptions::RuntimeError, "MueLu::SaPFactory::constrainRow: feasible solution but computation resulted in a lower bound violation??? ");
  TEUCHOS_TEST_FOR_EXCEPTION(upperViolation, Exceptions::RuntimeError, "MueLu::SaPFactory::constrainRow: feasible solution but computation resulted in an upper bound violation??? ");
  TEUCHOS_TEST_FOR_EXCEPTION(sumViolation, Exceptions::RuntimeError, "MueLu::SaPFactory::constrainRow: feasible solution but computation resulted in a row sum violation??? ");

  return hasFeasibleSol;
}

}  // namespace MueLu

#endif  // MUELU_SAPFACTORY_DEF_HPP

// TODO: restrictionMode_ should use the parameter list.
