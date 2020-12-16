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

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

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
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities_kokkos.hpp"

#include <sstream>

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> SaPFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("sa: damping factor");
    SET_VALID_ENTRY("sa: calculate eigenvalue estimate");
    SET_VALID_ENTRY("sa: eigenvalue estimate num iterations");
    SET_VALID_ENTRY("sa: use rowsumabs diagonal scaling");
    SET_VALID_ENTRY("sa: enforce constraints");
    SET_VALID_ENTRY("tentative: calculate qr");
    SET_VALID_ENTRY("sa: max eigenvalue");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    validParamList->set< RCP<const FactoryBase> >("P", Teuchos::null, "Tentative prolongator factory");

    // Make sure we don't recursively validate options for the matrixmatrix kernels
    ParameterList norecurse;
    norecurse.disableRecursiveValidation();
    validParamList->set<ParameterList> ("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");


    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void SaPFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory_kokkos and NullspaceFactory.
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P", initialPFact.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void SaPFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void SaPFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

    // Add debugging information
    DeviceType::execution_space::print_configuration(GetOStream(Runtime1));

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory_kokkos and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    const ParameterList& pL = GetParameterList();

    // Level Get
    RCP<Matrix> A     = Get< RCP<Matrix> >(fineLevel, "A");
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    if(restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);

      A = Utilities_kokkos::Transpose(*A, true); // build transpose of A explicitly
    }

    //Build final prolongator
    RCP<Matrix> finalP; // output

    // Reuse pattern if available
    RCP<ParameterList> APparams;
    if(pL.isSublist("matrixmatrix: kernel params")) 
      APparams=rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
    else 
      APparams= rcp(new ParameterList);
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
    const LO maxEigenIterations = as<LO>(pL.get<int>("sa: eigenvalue estimate num iterations"));
    const bool estimateMaxEigen = pL.get<bool>("sa: calculate eigenvalue estimate");
    const bool useAbsValueRowSum = pL.get<bool>  ("sa: use rowsumabs diagonal scaling");
    const bool doQRStep         =        pL.get<bool>("tentative: calculate qr");
    const bool enforceConstraints=       pL.get<bool>("sa: enforce constraints");
    const SC   userDefinedMaxEigen =  as<SC>(pL.get<double>("sa: max eigenvalue"));

    // Sanity checking
    TEUCHOS_TEST_FOR_EXCEPTION(doQRStep && enforceConstraints,Exceptions::RuntimeError,
                               "MueLu::TentativePFactory::MakeTentative: cannot use 'enforce constraints' and 'calculate qr' at the same time");

    if (dampingFactor != Teuchos::ScalarTraits<SC>::zero()) {

      SC lambdaMax;
      RCP<Vector> invDiag;
      if (Teuchos::ScalarTraits<SC>::real(userDefinedMaxEigen) == Teuchos::ScalarTraits<SC>::real(-1.0))
      {
        SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
        lambdaMax = A->GetMaxEigenvalueEstimate();
        if (lambdaMax == -Teuchos::ScalarTraits<SC>::one() || estimateMaxEigen) {
          GetOStream(Statistics1) << "Calculating max eigenvalue estimate now (max iters = "<< maxEigenIterations <<
          ( (useAbsValueRowSum) ?  ", use rowSumAbs diagonal)" :  ", use point diagonal)") << std::endl;
          Magnitude stopTol = 1e-4;
          invDiag = Utilities_kokkos::GetMatrixDiagonalInverse(*A, Teuchos::ScalarTraits<SC>::eps()*100, useAbsValueRowSum);
          if (useAbsValueRowSum)
            lambdaMax = Utilities_kokkos::PowerMethod(*A, invDiag, maxEigenIterations, stopTol);
          else
            lambdaMax = Utilities_kokkos::PowerMethod(*A, true, maxEigenIterations, stopTol);
          A->SetMaxEigenvalueEstimate(lambdaMax);
        } else {
          GetOStream(Statistics1) << "Using cached max eigenvalue estimate" << std::endl;
        }
      }
      else {
          lambdaMax = userDefinedMaxEigen;
          A->SetMaxEigenvalueEstimate(lambdaMax);
      }
      GetOStream(Statistics0) << "Prolongator damping factor = " << dampingFactor/lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;

      {
        SubFactoryMonitor m2(*this, "Fused (I-omega*D^{-1} A)*Ptent", coarseLevel);
        {
          SubFactoryMonitor m3(*this, "Diagonal Extraction", coarseLevel);
          if (useAbsValueRowSum)
            GetOStream(Runtime0) << "Using rowSumAbs diagonal" << std::endl;
          if (invDiag == Teuchos::null)
            invDiag = Utilities_kokkos::GetMatrixDiagonalInverse(*A, Teuchos::ScalarTraits<SC>::eps()*100, useAbsValueRowSum);
        }
        SC omega = dampingFactor / lambdaMax;
        TEUCHOS_TEST_FOR_EXCEPTION(!std::isfinite(Teuchos::ScalarTraits<SC>::magnitude(omega)), Exceptions::RuntimeError, "Prolongator damping factor needs to be finite.");

        {
          SubFactoryMonitor m3(*this, "Xpetra::IteratorOps::Jacobi", coarseLevel);
          // finalP = Ptent + (I - \omega D^{-1}A) Ptent
          finalP = Xpetra::IteratorOps<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Jacobi(omega, *invDiag, *A, *Ptent, finalP, GetOStream(Statistics2), std::string("MueLu::SaP-") + toString(coarseLevel.GetLevelID()), APparams);
        if (enforceConstraints) SatisfyPConstraints( A, finalP);
        }
      }

    } else {
      finalP = Ptent;
    }

    // Level Set
    if (!restrictionMode_) {
      // prolongation factory is in prolongation mode
      if(!finalP.is_null()) {std::ostringstream oss; oss << "P_" << coarseLevel.GetLevelID(); finalP->setObjectLabel(oss.str());}
      Set(coarseLevel, "P", finalP);

      // NOTE: EXPERIMENTAL
      if (Ptent->IsView("stridedMaps"))
        finalP->CreateView("stridedMaps", Ptent);

    } else {
      // prolongation factory is in restriction mode
      RCP<Matrix> R = Utilities_kokkos::Transpose(*finalP, true);
      Set(coarseLevel, "R", R);
      if(!R.is_null()) {std::ostringstream oss; oss << "R_" << coarseLevel.GetLevelID(); R->setObjectLabel(oss.str());}

      // NOTE: EXPERIMENTAL
      if (Ptent->IsView("stridedMaps"))
        R->CreateView("stridedMaps", Ptent, true);
    }

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo",          true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*finalP, (!restrictionMode_ ? "P" : "R"), params);
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

template<typename local_matrix_type>
struct constraintKernel {

   using Scalar= typename local_matrix_type::non_const_value_type;
   using SC=Scalar;
   using LO=typename local_matrix_type::non_const_ordinal_type;
   using Device= typename local_matrix_type::device_type;
   const Scalar zero = Kokkos::ArithTraits<SC>::zero();
   const Scalar one  = Kokkos::ArithTraits<SC>::one();
   LO     nPDEs;
   local_matrix_type localP;
   Kokkos::View<SC**, Device> ConstraintViolationSum;
   Kokkos::View<SC**, Device> Rsum;
   Kokkos::View<size_t**, Device> nPositive;

  constraintKernel(LO nPDEs_,local_matrix_type localP_) : nPDEs(nPDEs_), localP(localP_)
   {
     ConstraintViolationSum = Kokkos::View<SC**, Device>("ConstraintViolationSum", localP_.numRows(), nPDEs);
     Rsum = Kokkos::View<SC**, Device>("Rsum", localP_.numRows(), nPDEs);
     nPositive = Kokkos::View<size_t**, Device>("nPositive", localP_.numRows(), nPDEs);
   }

   KOKKOS_INLINE_FUNCTION
   void operator() (const size_t rowIdx) const {

      auto rowPtr = localP.graph.row_map;
      auto values = localP.values;

      bool checkRow = true;

      if (rowPtr(rowIdx + 1) == rowPtr(rowIdx)) checkRow = false;


      while (checkRow) { 

        // check constraints and compute the row sum
    
        for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++)  {
          Rsum(rowIdx, entryIdx%nPDEs) += values(entryIdx); 
          if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) < Kokkos::ArithTraits<SC>::real(zero)) { 

            ConstraintViolationSum(rowIdx, entryIdx%nPDEs) += values(entryIdx); 
            values(entryIdx) = zero;
          }
          else {
            if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) != Kokkos::ArithTraits<SC>::real(zero))
              nPositive(rowIdx, entryIdx%nPDEs) = nPositive(rowIdx, entryIdx%nPDEs) + 1;
    
            if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) > Kokkos::ArithTraits<SC>::real(1.00001  )) { 
              ConstraintViolationSum(rowIdx, entryIdx%nPDEs) += (values(entryIdx) - one); 
              values(entryIdx) =  one;
            }
          }
        }
    
        checkRow = false;

        // take into account any row sum that violates the contraints
    
        for (size_t k=0; k < (size_t) nPDEs; k++) {

          if (Kokkos::ArithTraits<SC>::real(Rsum(rowIdx, k)) < Kokkos::ArithTraits<SC>::magnitude(zero)) {
              ConstraintViolationSum(rowIdx, k) = ConstraintViolationSum(rowIdx, k) - Rsum(rowIdx, k);  // rstumin 
          }
          else if (Kokkos::ArithTraits<SC>::real(Rsum(rowIdx, k)) > Kokkos::ArithTraits<SC>::magnitude(1.00001)) {
              ConstraintViolationSum(rowIdx, k) = ConstraintViolationSum(rowIdx, k)+ (one - Rsum(rowIdx, k));  // rstumin
          }
        }

        // check if row need modification 
        for (size_t k=0; k < (size_t) nPDEs; k++) {
          if (Kokkos::ArithTraits<SC>::magnitude(ConstraintViolationSum(rowIdx, k)) != Kokkos::ArithTraits<SC>::magnitude(zero))
             checkRow = true;
        }
        // modify row
        if (checkRow) {
	   for (auto entryIdx = rowPtr(rowIdx); entryIdx < rowPtr(rowIdx + 1); entryIdx++)  {
             if (Kokkos::ArithTraits<SC>::real(values(entryIdx)) > Kokkos::ArithTraits<SC>::real(zero)) { 
	       values(entryIdx) = values(entryIdx) +
		 (ConstraintViolationSum(rowIdx, entryIdx%nPDEs)/ (Scalar (nPositive(rowIdx, entryIdx%nPDEs)) != zero ? Scalar (nPositive(rowIdx, entryIdx%nPDEs)) : one));
             }
           }
           for (size_t k=0; k < (size_t) nPDEs; k++) ConstraintViolationSum(rowIdx, k) = zero; 
        }
        for (size_t k=0; k < (size_t) nPDEs; k++) Rsum(rowIdx, k) = zero; 
        for (size_t k=0; k < (size_t) nPDEs; k++) nPositive(rowIdx, k) = 0;
      } // while (checkRow) ...

   }
};
   

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void SaPFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::SatisfyPConstraints(const RCP<Matrix> A, RCP<Matrix>& P) const {

    using Device = typename Matrix::local_matrix_type::device_type;
    LO nPDEs = A->GetFixedBlockSize();

    using local_mat_type = typename Matrix::local_matrix_type;
    constraintKernel<local_mat_type> myKernel(nPDEs,P->getLocalMatrix() );
    Kokkos::parallel_for("enforce constraint",Kokkos::RangePolicy<typename Device::execution_space>(0, P->getRowMap()->getNodeNumElements() ),
                        myKernel );

  } //SatsifyPConstraints()

} //namespace MueLu

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_SAPFACTORY_KOKKOS_DEF_HPP

//TODO: restrictionMode_ should use the parameter list.
