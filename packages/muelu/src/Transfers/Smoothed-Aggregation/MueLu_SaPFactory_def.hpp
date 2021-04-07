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
#include <Xpetra_IteratorOps.hpp> // containing routines to generate Jacobi iterator
#include <Xpetra_IO.hpp>
#include <sstream>

#include "MueLu_SaPFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
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
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Tentative prolongator factory");

    // Make sure we don't recursively validate options for the matrixmatrix kernels
    ParameterList norecurse;
    norecurse.disableRecursiveValidation();
    validParamList->set<ParameterList> ("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");


    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P", initialPFact.get(), this); // --
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

    std::string levelIDs = toString(coarseLevel.GetLevelID());

    const std::string prefix = "MueLu::SaPFactory(" + levelIDs + "): ";

    typedef typename Teuchos::ScalarTraits<SC>::coordinateType Coordinate;

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    const ParameterList& pL = GetParameterList();

    // Level Get
    RCP<Matrix> A     = Get< RCP<Matrix> >(fineLevel, "A");
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    if (restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
      A = Utilities::Transpose(*A, true); // build transpose of A explicitely
    }

    // Build final prolongator
    RCP<Matrix> finalP;

    // Reuse pattern if available
    RCP<ParameterList> APparams = rcp(new ParameterList);;
    if(pL.isSublist("matrixmatrix: kernel params"))
      APparams->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

    if (coarseLevel.IsAvailable("AP reuse data", this)) {
      GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous AP data" << std::endl;

      APparams = coarseLevel.Get< RCP<ParameterList> >("AP reuse data", this);

      if (APparams->isParameter("graph"))
        finalP = APparams->get< RCP<Matrix> >("graph");
    }
    // By default, we don't need global constants for SaP
    APparams->set("compute global constants: temporaries",APparams->get("compute global constants: temporaries",false));
    APparams->set("compute global constants",APparams->get("compute global constants",false));

    const SC dampingFactor      = as<SC>(pL.get<double>("sa: damping factor"));
    const LO maxEigenIterations = as<LO>(pL.get<int>   ("sa: eigenvalue estimate num iterations"));
    const bool estimateMaxEigen =        pL.get<bool>  ("sa: calculate eigenvalue estimate");
    const bool useAbsValueRowSum=        pL.get<bool>  ("sa: use rowsumabs diagonal scaling");
    const bool doQRStep         =        pL.get<bool>("tentative: calculate qr");
    const bool enforceConstraints=       pL.get<bool>("sa: enforce constraints");
    const SC   userDefinedMaxEigen =  as<SC>(pL.get<double>("sa: max eigenvalue"));
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    double dTol = pL.get<double>("sa: rowsumabs diagonal replacement tolerance");
    const Magnitude diagonalReplacementTolerance = (dTol == as<double>(-1) ? Teuchos::ScalarTraits<Scalar>::eps()*100 : as<Magnitude>(pL.get<double>("sa: rowsumabs diagonal replacement tolerance")));
    const SC diagonalReplacementValue =  as<SC>(pL.get<double>("sa: rowsumabs diagonal replacement value"));

    // Sanity checking
    TEUCHOS_TEST_FOR_EXCEPTION(doQRStep && enforceConstraints,Exceptions::RuntimeError,
                               "MueLu::TentativePFactory::MakeTentative: cannot use 'enforce constraints' and 'calculate qr' at the same time");

    if (dampingFactor != Teuchos::ScalarTraits<SC>::zero()) {

      Scalar lambdaMax;
      Teuchos::RCP<Vector> invDiag;
      if (userDefinedMaxEigen == -1.)
      {
        SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
        lambdaMax = A->GetMaxEigenvalueEstimate();
        if (lambdaMax == -Teuchos::ScalarTraits<SC>::one() || estimateMaxEigen) {
          GetOStream(Statistics1) << "Calculating max eigenvalue estimate now (max iters = "<< maxEigenIterations <<
          ( (useAbsValueRowSum) ?  ", use rowSumAbs diagonal)" :  ", use point diagonal)") << std::endl;
          Coordinate stopTol = 1e-4;
          if (useAbsValueRowSum) {
            const bool returnReciprocal=true;
            const bool replaceSingleEntryRowWithZero=true;
            invDiag = Utilities::GetLumpedMatrixDiagonal(*A,returnReciprocal,
                                                        diagonalReplacementTolerance,
                                                        diagonalReplacementValue,
                                                        replaceSingleEntryRowWithZero);
            TEUCHOS_TEST_FOR_EXCEPTION(invDiag.is_null(), Exceptions::RuntimeError,
                                       "SaPFactory: eigenvalue estimate: diagonal reciprocal is null.");
            lambdaMax = Utilities::PowerMethod(*A, invDiag, maxEigenIterations, stopTol);
          } else
            lambdaMax = Utilities::PowerMethod(*A, true, maxEigenIterations, stopTol);
          A->SetMaxEigenvalueEstimate(lambdaMax);
        } else {
          GetOStream(Statistics1) << "Using cached max eigenvalue estimate" << std::endl;
        }
      }
      else {
          lambdaMax = userDefinedMaxEigen;
          A->SetMaxEigenvalueEstimate(lambdaMax);
      }
      GetOStream(Statistics1) << "Prolongator damping factor = " << dampingFactor/lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;

      {
        SubFactoryMonitor m2(*this, "Fused (I-omega*D^{-1} A)*Ptent", coarseLevel);
        if (!useAbsValueRowSum)
          invDiag = Utilities::GetMatrixDiagonalInverse(*A); //default
        else if (invDiag == Teuchos::null) {
          GetOStream(Runtime0) << "Using rowsumabs diagonal" << std::endl;
          const bool returnReciprocal=true;
          const bool replaceSingleEntryRowWithZero=true;
          invDiag = Utilities::GetLumpedMatrixDiagonal(*A,returnReciprocal,
                                                        diagonalReplacementTolerance,
                                                        diagonalReplacementValue,
                                                        replaceSingleEntryRowWithZero);
          TEUCHOS_TEST_FOR_EXCEPTION(invDiag.is_null(), Exceptions::RuntimeError, "SaPFactory: diagonal reciprocal is null.");
        }	
	
	      SC omega = dampingFactor / lambdaMax;
	      TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(Teuchos::ScalarTraits<SC>::magnitude(omega)), Exceptions::RuntimeError, "Prolongator damping factor needs to be finite.");

        // finalP = Ptent + (I - \omega D^{-1}A) Ptent
        finalP = Xpetra::IteratorOps<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Jacobi(omega, *invDiag, *A, *Ptent, finalP,
                    GetOStream(Statistics2), std::string("MueLu::SaP-")+levelIDs, APparams);
        if (enforceConstraints) SatisfyPConstraints( A, finalP);
      }

    } else {
      finalP = Ptent;
    }

    // Level Set
    if (!restrictionMode_) {
      // The factory is in prolongation mode
      if(!finalP.is_null()) {std::ostringstream oss; oss << "P_" << coarseLevel.GetLevelID(); finalP->setObjectLabel(oss.str());}
      Set(coarseLevel, "P",             finalP);
    
      APparams->set("graph", finalP);
      Set(coarseLevel, "AP reuse data", APparams);

      // NOTE: EXPERIMENTAL
      if (Ptent->IsView("stridedMaps"))
        finalP->CreateView("stridedMaps", Ptent);

      if (IsPrint(Statistics2)) {
        RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        params->set("printCommInfo",          true);
        GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*finalP, "P", params);
      }

    } else {
      // The factory is in restriction mode
      RCP<Matrix> R;
      {
        SubFactoryMonitor m2(*this, "Transpose P", coarseLevel);
        R = Utilities::Transpose(*finalP, true);
        if(!R.is_null()) {std::ostringstream oss; oss << "R_" << coarseLevel.GetLevelID(); R->setObjectLabel(oss.str());}
      }

      Set(coarseLevel, "R", R);

      // NOTE: EXPERIMENTAL
      if (Ptent->IsView("stridedMaps"))
        R->CreateView("stridedMaps", Ptent, true/*transposeA*/);

      if (IsPrint(Statistics2)) {
        RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        params->set("printCommInfo",          true);
        GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*R, "R", params);
      }
    }

  } //Build()

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
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SatisfyPConstraints(const RCP<Matrix> A, RCP<Matrix>& P) const {

    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
    LO nPDEs = A->GetFixedBlockSize();
    Teuchos::ArrayRCP<Scalar> ConstraintViolationSum(nPDEs);
    Teuchos::ArrayRCP<Scalar> Rsum(nPDEs);
    Teuchos::ArrayRCP<size_t> nPositive(nPDEs);
    for (size_t k=0; k < (size_t) nPDEs; k++) ConstraintViolationSum[k] = zero;
    for (size_t k=0; k < (size_t) nPDEs; k++) Rsum[k] = zero;
    for (size_t k=0; k < (size_t) nPDEs; k++) nPositive[k] = 0;


    for (size_t i = 0; i < as<size_t>(P->getRowMap()->getNodeNumElements()); i++) {

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals1;
      Teuchos::ArrayView<      Scalar> vals;
      P->getLocalRowView((LO) i, indices, vals1);
      size_t nnz = indices.size();
      if (nnz == 0) continue;

      vals = ArrayView<Scalar>(const_cast<SC*>(vals1.getRawPtr()), nnz);


      bool checkRow = true;

      while (checkRow) { 

        // check constraints and compute the row sum
    
        for (LO j = 0; j < indices.size(); j++)  {
          Rsum[ j%nPDEs ] += vals[j]; 
          if (Teuchos::ScalarTraits<SC>::real(vals[j]) < Teuchos::ScalarTraits<SC>::real(zero)) { 
            ConstraintViolationSum[ j%nPDEs ] += vals[j]; 
            vals[j] = zero;
          }
          else {
            if (Teuchos::ScalarTraits<SC>::real(vals[j]) != Teuchos::ScalarTraits<SC>::real(zero))
              (nPositive[ j%nPDEs])++;
    
            if (Teuchos::ScalarTraits<SC>::real(vals[j]) > Teuchos::ScalarTraits<SC>::real(1.00001  )) { 
              ConstraintViolationSum[ j%nPDEs ] += (vals[j] - one); 
              vals[j] = one;
            }
          }
        }
    
        checkRow = false;

        // take into account any row sum that violates the contraints
    
        for (size_t k=0; k < (size_t) nPDEs; k++) {

          if (Teuchos::ScalarTraits<SC>::real(Rsum[ k ]) < Teuchos::ScalarTraits<SC>::magnitude(zero)) {
              ConstraintViolationSum[k] +=  (-Rsum[k]);  // rstumin 
          }
          else if (Teuchos::ScalarTraits<SC>::real(Rsum[ k ]) > Teuchos::ScalarTraits<SC>::magnitude(1.00001)) {
              ConstraintViolationSum[k] += (one - Rsum[k]);  // rstumin
          }
        }

        // check if row need modification 
        for (size_t k=0; k < (size_t) nPDEs; k++) {
          if (Teuchos::ScalarTraits<SC>::magnitude(ConstraintViolationSum[ k ]) != Teuchos::ScalarTraits<SC>::magnitude(zero))
             checkRow = true;
        }
        // modify row
        if (checkRow) {
           for (LO j = 0; j < indices.size(); j++)  {
             if (Teuchos::ScalarTraits<SC>::real(vals[j]) > Teuchos::ScalarTraits<SC>::real(zero)) { 
                vals[j] += (ConstraintViolationSum[j%nPDEs]/ (as<Scalar>(nPositive[j%nPDEs])));
             }
           }
           for (size_t k=0; k < (size_t) nPDEs; k++) ConstraintViolationSum[k] = zero; 
        }
        for (size_t k=0; k < (size_t) nPDEs; k++) Rsum[k] = zero; 
        for (size_t k=0; k < (size_t) nPDEs; k++) nPositive[k] = 0;
      } // while (checkRow) ...
    } // for (size_t i = 0; i < as<size_t>(P->getRowMap()->getNumNodeElements()); i++) ...
  } //SatsifyPConstraints()

} //namespace MueLu

#endif // MUELU_SAPFACTORY_DEF_HPP

//TODO: restrictionMode_ should use the parameter list.
