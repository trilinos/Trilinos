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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PGPFACTORY_DEF_HPP
#define MUELU_PGPFACTORY_DEF_HPP

#include <vector>

#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_PgPFactory_decl.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PgPFactory()
    : diagonalView_("current"), min_norm_(DINVANORM), bReUseRowBasedOmegas_(false) {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDiagonalView(std::string const& diagView) {
    diagonalView_ = diagView;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMinimizationMode(MinimizationNorm minnorm) { min_norm_ = minnorm; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MueLu::MinimizationNorm PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMinimizationMode() { return min_norm_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDiagonalView() {
    return diagonalView_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P", initialPFact.get(), this);

    /* If PgPFactory is reusing the row based damping parameters omega for
     * restriction, it has to request the data here.
     * we have the following scenarios:
     * 1) Reuse omegas:
     * PgPFactory.DeclareInput for prolongation mode requests A and P0
     * PgPFactory.DeclareInput for restriction mode requests A, P0 and RowBasedOmega (call triggered by GenericRFactory)
     * PgPFactory.Build for prolongation mode calculates RowBasedOmega and stores it (as requested)
     * PgPFactory.Build for restriction mode reuses RowBasedOmega (and Releases the data with the Get call)
     * 2) do not reuse omegas
     * PgPFactory.DeclareInput for prolongation mode requests A and P0
     * PgPFactory.DeclareInput for restriction mode requests A and P0
     * PgPFactory.Build for prolongation mode calculates RowBasedOmega for prolongation operator
     * PgPFactory.Build for restriction mode calculates RowBasedOmega for restriction operator
     */
    if( bReUseRowBasedOmegas_ == true && restrictionMode_ == true ) {
      coarseLevel.DeclareInput("RowBasedOmega", this, this); // RowBasedOmega is calculated by this PgPFactory and requested by this PgPFactory
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator smoothing (PG-AMG)", coarseLevel);

    // Level Get
    RCP<Matrix> A = Get< RCP<Matrix> >(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    /////////////////// switch from A to A^T in restriction mode (necessary as long as implicit transpose not working for Epetra)
    if(restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
      A = Utils2::Transpose(A, true); // build transpose of A explicitely
    }

    /////////////////// calculate D^{-1} A Ptent (needed for smoothing)
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Matrix> DinvAP0 = Utils::TwoMatrixMultiply(A, false, Ptent, false, doFillComplete, optimizeStorage);

    doFillComplete=true;
    optimizeStorage=false;
    Teuchos::ArrayRCP<Scalar> diag = Utils::GetMatrixDiagonal(*A);
    Utils::MyOldScaleMatrix(DinvAP0, diag, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag

    /////////////////// calculate local damping factors omega

    Teuchos::RCP<Vector > RowBasedOmega = Teuchos::null;

    if(restrictionMode_ == false || bReUseRowBasedOmegas_ == false) {
      // if in prolongation mode: calculate row based omegas
      // if in restriction mode: calculate omegas only if row based omegas are not used from prolongation mode
      ComputeRowBasedOmega(fineLevel, coarseLevel, A, Ptent, DinvAP0, RowBasedOmega);
    } // if(bReUseRowBasedOmegas == false)
    else  {
      // reuse row based omegas, calculated by this factory in the run before (with restrictionMode_ == false)
      RowBasedOmega = coarseLevel.Get<Teuchos::RCP<Vector > >("RowBasedOmega", this);

      // RowBasedOmega is now based on row map of A (not transposed)
      // for restriction we use A^T instead of A
      // -> recommunicate row based omega

      // exporter: overlapping row map to nonoverlapping domain map (target map is unique)
      // since A is already transposed we use the RangeMap of A
      Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(RowBasedOmega->getMap(), A->getRangeMap());

      Teuchos::RCP<Vector > noRowBasedOmega =
        VectorFactory::Build(A->getRangeMap());

      noRowBasedOmega->doExport(*RowBasedOmega, *exporter, Xpetra::INSERT);

      // importer: nonoverlapping map to overlapping map

      // importer: source -> target maps
      Teuchos::RCP<const Import > importer =
        ImportFactory::Build(A->getRangeMap(), A->getRowMap());

      // doImport target->doImport(*source, importer, action)
      RowBasedOmega->doImport(*noRowBasedOmega, *importer, Xpetra::INSERT);
    }

    Teuchos::ArrayRCP< Scalar > RowBasedOmega_local = RowBasedOmega->getDataNonConst(0);

    /////////////////// prolongator smoothing using local damping parameters omega
    RCP<Matrix> P_smoothed = Teuchos::null;
    Utils::MyOldScaleMatrix(DinvAP0, RowBasedOmega_local, false, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag

    Utils2::TwoMatrixAdd(Ptent, false, Teuchos::ScalarTraits<Scalar>::one(),
                         DinvAP0, false, -Teuchos::ScalarTraits<Scalar>::one(),
                         P_smoothed);
    P_smoothed->fillComplete(Ptent->getDomainMap(), Ptent->getRangeMap());

    //////////////////// store results in Level

    // Level Set
    if(!restrictionMode_)
      {
        // prolongation factory is in prolongation mode
        Set(coarseLevel, "P", P_smoothed);

        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) P_smoothed->CreateView("stridedMaps", Ptent);
        ///////////////////////// EXPERIMENTAL
      }
    else
      {
        // prolongation factory is in restriction mode
        RCP<Matrix> R = Utils2::Transpose(P_smoothed, true); // use Utils2 -> specialization for double
        Set(coarseLevel, "R", R);

        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) R->CreateView("stridedMaps", Ptent, true);
        ///////////////////////// EXPERIMENTAL
      }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeRowBasedOmega(Level& fineLevel, Level &coarseLevel, const RCP<Matrix>& A, const RCP<Matrix>& P0, const RCP<Matrix>& DinvAP0, RCP<Vector > & RowBasedOmega) const {
    FactoryMonitor m(*this, "PgPFactory::ComputeRowBasedOmega", coarseLevel);

    Teuchos::RCP<Vector > Numerator = Teuchos::null;
    Teuchos::RCP<Vector > Denominator = Teuchos::null;

    switch (min_norm_)
      {
      case ANORM: {
        // MUEMAT mode (=paper)
        // Minimize with respect to the (A)' A norm.
        // Need to be smart here to avoid the construction of A' A
        //
        //                   diag( P0' (A' A) D^{-1} A P0)
        //   omega =   ------------------------------------------
        //             diag( P0' A' D^{-1}' ( A'  A) D^{-1} A P0)
        //
        // expensive, since we have to recalculate AP0 due to the lack of an explicit scaling routine for DinvAP0

        // calculate A * P0
        bool doFillComplete=true;
        bool optimizeStorage=false;
        RCP<Matrix> AP0 = Utils::TwoMatrixMultiply(A, false, P0, false, doFillComplete, optimizeStorage);

        // compute A * D^{-1} * A * P0
        RCP<Matrix> ADinvAP0 = Utils::TwoMatrixMultiply(A, false, DinvAP0, false, doFillComplete, optimizeStorage);

        Numerator =   VectorFactory::Build(ADinvAP0->getColMap(), true);
        Denominator = VectorFactory::Build(ADinvAP0->getColMap(), true);
        MultiplyAll(AP0, ADinvAP0, Numerator);
        MultiplySelfAll(ADinvAP0, Denominator);
      }
        break;
      case L2NORM: {

        // ML mode 1 (cheapest)
        // Minimize with respect to L2 norm
        //                  diag( P0' D^{-1} A P0)
        //   omega =   -----------------------------
        //             diag( P0' A' D^{-1}' D^{-1} A P0)
        //
        Numerator =   VectorFactory::Build(DinvAP0->getColMap(), true);
        Denominator = VectorFactory::Build(DinvAP0->getColMap(), true);
        MultiplyAll(P0, DinvAP0, Numerator);
        MultiplySelfAll(DinvAP0, Denominator);
      }
        break;
      case DINVANORM: {
        // ML mode 2
        // Minimize with respect to the (D^{-1} A)' D^{-1} A norm.
        // Need to be smart here to avoid the construction of A' A
        //
        //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
        //   omega =   --------------------------------------------------------
        //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
        //

        // compute D^{-1} * A * D^{-1} * A * P0
        bool doFillComplete=true;
        bool optimizeStorage=false;
        Teuchos::ArrayRCP<Scalar> diagA = Utils::GetMatrixDiagonal(*A);
        RCP<Matrix> DinvADinvAP0 = Utils::TwoMatrixMultiply(A, false, DinvAP0, false, doFillComplete, optimizeStorage);
        Utils::MyOldScaleMatrix(DinvADinvAP0, diagA, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag

        Numerator =   VectorFactory::Build(DinvADinvAP0->getColMap(), true);
        Denominator = VectorFactory::Build(DinvADinvAP0->getColMap(), true);
        MultiplyAll(DinvAP0, DinvADinvAP0, Numerator);
        MultiplySelfAll(DinvADinvAP0, Denominator);
      }
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::Build: minimization mode not supported. error");
      }


    //////////// build Column based omegas /////////////
    Teuchos::RCP<Vector > ColBasedOmega =
      VectorFactory::Build(Numerator->getMap()/*DinvAP0->getColMap()*/, true);

    ColBasedOmega->putScalar(-666/**Teuchos::ScalarTraits<Scalar>::one()*/);

    Teuchos::ArrayRCP< const Scalar > Numerator_local = Numerator->getData(0);
    Teuchos::ArrayRCP< const Scalar > Denominator_local = Denominator->getData(0);
    Teuchos::ArrayRCP< Scalar >       ColBasedOmega_local = ColBasedOmega->getDataNonConst(0);
    GlobalOrdinal zero_local = 0;  // count negative colbased omegas
    GlobalOrdinal nan_local = 0;   // count NaNs -> set them to zero
    Scalar min_local = 1000000; //Teuchos::ScalarTraits<Scalar>::one() * (Scalar) 1000000;
    Scalar max_local = Teuchos::ScalarTraits<Scalar>::zero();
    for(LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(Numerator->getLocalLength()); i++) {
      if(std::abs(Denominator_local[i]) == 0.0)
        {
          ColBasedOmega_local[i] = 0.0; // fallback: nonsmoothed basis function since denominator == 0.0
          nan_local++;
        }
      else
        {
          ColBasedOmega_local[i] = Numerator_local[i] / Denominator_local[i];  // default case
        }

      if(std::abs(ColBasedOmega_local[i]) < std::abs(Teuchos::ScalarTraits<Scalar>::zero())) { // negative omegas are not valid. set them to zero
        ColBasedOmega_local[i] = Teuchos::ScalarTraits<Scalar>::zero();
        zero_local++; // count zero omegas
      }

      // handle case that Nominator == Denominator -> Dirichlet bcs in A?
      // fallback if ColBasedOmega == 1 -> very strong smoothing may lead to zero rows in P
      // TAW: this is somewhat nonstandard and a rough fallback strategy to avoid problems
      // also avoid "overshooting" with omega > 0.8
      if(std::abs(ColBasedOmega_local[i]) >= 0.8) {
        ColBasedOmega_local[i] = 0.0;
      }

      if(std::abs(ColBasedOmega_local[i]) < std::abs(min_local)) { min_local = std::abs(ColBasedOmega_local[i]); }
      if(std::abs(ColBasedOmega_local[i]) > std::abs(max_local)) { max_local = std::abs(ColBasedOmega_local[i]); }
    }

    { // be verbose
      GlobalOrdinal zero_all;
      GlobalOrdinal nan_all;
      Scalar min_all;
      Scalar max_all;
      sumAll(A->getRowMap()->getComm(), zero_local, zero_all);
      sumAll(A->getRowMap()->getComm(), nan_local, nan_all);
      minAll(A->getRowMap()->getComm(), min_local, min_all);
      maxAll(A->getRowMap()->getComm(), max_local, max_all);

      GetOStream(MueLu::Statistics1, 0) << "PgPFactory: smoothed aggregation (scheme: ";
      switch (min_norm_)
        {
        case ANORM:     { GetOStream(MueLu::Statistics1, 0) << "Anorm)"     << std::endl;   }   break;
        case L2NORM:    { GetOStream(MueLu::Statistics1, 0) << "L2norm)"    << std::endl;   }   break;
        case DINVANORM: { GetOStream(MueLu::Statistics1, 0) << "DinvAnorm)" << std::endl;   }    break;
        default:          GetOStream(MueLu::Statistics1, 0) << "unknown)" << std::endl;
        }
      GetOStream(MueLu::Statistics1, 0) << "Damping parameter: min = " << min_all << ", max = " << max_all << std::endl;
      GetOStream(MueLu::Statistics, 0) << "# negative omegas: " << zero_all << " out of " << ColBasedOmega->getGlobalLength() << " column-based omegas" << std::endl;
      GetOStream(MueLu::Statistics, 0) << "# NaNs: " << nan_all << " out of " << ColBasedOmega->getGlobalLength() << " column-based omegas" << std::endl;
    }

    if(coarseLevel.IsRequested("ColBasedOmega", this)) {
      coarseLevel.Set("ColBasedOmega", ColBasedOmega, this);
    }

    //////////// build Row based omegas /////////////
    // transform column based omegas to row based omegas
    RowBasedOmega =
      VectorFactory::Build(DinvAP0->getRowMap(), true);

    RowBasedOmega->putScalar(-666); // TODO bad programming style

    bool bAtLeastOneDefined = false;
    Teuchos::ArrayRCP< Scalar > RowBasedOmega_local = RowBasedOmega->getDataNonConst(0);
    for(LocalOrdinal row = 0; row<Teuchos::as<LocalOrdinal>(A->getNodeNumRows()); row++) {
      Teuchos::ArrayView<const LocalOrdinal> lindices;
      Teuchos::ArrayView<const Scalar> lvals;
      DinvAP0->getLocalRowView(row, lindices, lvals);
      bAtLeastOneDefined = false;
      for(size_t j=0; j<Teuchos::as<size_t>(lindices.size()); j++) {
        Scalar omega = ColBasedOmega_local[lindices[j]];
        if (std::abs(omega) != -666) { // TODO bad programming style
          bAtLeastOneDefined = true;
          if(std::abs(RowBasedOmega_local[row]) == -666)    RowBasedOmega_local[row] = omega;
          else if(std::abs(omega) < std::abs(RowBasedOmega_local[row])) RowBasedOmega_local[row] = omega;
        }
      }
      if(bAtLeastOneDefined == true) {
        if(std::abs(RowBasedOmega_local[row]) < 0 /*Teuchos::ScalarTraits<Scalar>::zero()*/) RowBasedOmega_local[row] = 0; /* Teuchos::ScalarTraits<Scalar>::zero();*/
      }
    }

    if(coarseLevel.IsRequested("RowBasedOmega", this)) {
      Set(coarseLevel, "RowBasedOmega", RowBasedOmega);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiplySelfAll(const RCP<Matrix>& Op, Teuchos::RCP<Vector >& InnerProdVec) const {

    // note: InnerProdVec is based on column map of Op
    TEUCHOS_TEST_FOR_EXCEPTION(!InnerProdVec->getMap()->isSameAs(*Op->getColMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplySelfAll: map of InnerProdVec must be same as column map of operator. error");

    Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);

    Teuchos::ArrayView<const LocalOrdinal> lindices;
    Teuchos::ArrayView<const Scalar> lvals;

    for(size_t n=0; n<Op->getNodeNumRows(); n++) {
      Op->getLocalRowView(n, lindices, lvals);
      for(size_t i=0; i<Teuchos::as<size_t>(lindices.size()); i++) {
        InnerProd_local[lindices[i]] += lvals[i]*lvals[i];
      }
    }

    // exporter: overlapping map to nonoverlapping map (target map is unique)
    Teuchos::RCP<const Export> exporter =
      ExportFactory::Build(Op->getColMap(), Op->getDomainMap());

    Teuchos::RCP<Vector > nonoverlap =
      VectorFactory::Build(Op->getDomainMap());

    nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Import > importer =
      ImportFactory::Build(Op->getDomainMap(), Op->getColMap());

    // doImport target->doImport(*source, importer, action)
    InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiplyAll(const RCP<Matrix>& left, const RCP<Matrix>& right, Teuchos::RCP<Vector >& InnerProdVec) const {

    TEUCHOS_TEST_FOR_EXCEPTION(!left->getDomainMap()->isSameAs(*right->getDomainMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: domain maps of left and right do not match. Error.");
    TEUCHOS_TEST_FOR_EXCEPTION(!left->getRowMap()->isSameAs(*right->getRowMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: row maps of left and right do not match. Error.");
#if 1 // 1=new "fast code, 0=old "slow", but safe code
#if 0 // not necessary - remove me
    if(InnerProdVec->getMap()->isSameAs(*left->getColMap())) {
      // initialize NewRightLocal vector and assign all entries to
      // left->getColMap()->getNodeNumElements() + 1
      std::vector<LocalOrdinal> NewRightLocal(right->getColMap()->getNodeNumElements(), Teuchos::as<LocalOrdinal>(left->getColMap()->getNodeNumElements()+1));

      LocalOrdinal i = 0;
      for (size_t j=0; j < right->getColMap()->getNodeNumElements(); j++) {
        while ( (i < Teuchos::as<LocalOrdinal>(left->getColMap()->getNodeNumElements())) &&
                (left->getColMap()->getGlobalElement(i) < right->getColMap()->getGlobalElement(j)) ) i++;
        if (left->getColMap()->getGlobalElement(i) == right->getColMap()->getGlobalElement(j)) {
          NewRightLocal[j] = i;
        }
      }

      Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);
      std::vector<Scalar> temp_array(left->getColMap()->getNodeNumElements()+1, 0.0);

      for(size_t n=0; n<right->getNodeNumRows(); n++) {
        Teuchos::ArrayView<const LocalOrdinal> lindices_left;
        Teuchos::ArrayView<const Scalar> lvals_left;
        Teuchos::ArrayView<const LocalOrdinal> lindices_right;
        Teuchos::ArrayView<const Scalar> lvals_right;

        left->getLocalRowView (n, lindices_left,  lvals_left);
        right->getLocalRowView(n, lindices_right, lvals_right);

        for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++) {
          temp_array[NewRightLocal[lindices_right[j] ] ] = lvals_right[j];
        }
        for (size_t j=0; j < Teuchos::as<size_t>(lindices_left.size()); j++) {
          InnerProd_local[lindices_left[j]] += temp_array[lindices_left[j] ]*lvals_left[j];
        }
        for (size_t j=0; j < Teuchos::as<size_t>(lindices_right.size()); j++) {
          temp_array[NewRightLocal[lindices_right[j] ] ] = 0.0;
        }
      }
      // exporter: overlapping map to nonoverlapping map (target map is unique)
      Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(left->getColMap(), left->getDomainMap()); // TODO: change left to right?

      Teuchos::RCP<Vector > nonoverlap =
        VectorFactory::Build(left->getDomainMap()); // TODO: change left to right?

      nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

      // importer: nonoverlapping map to overlapping map

      // importer: source -> target maps
      Teuchos::RCP<const Import > importer =
        ImportFactory::Build(left->getDomainMap(), left->getColMap()); // TODO: change left to right?

      // doImport target->doImport(*source, importer, action)
      InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);


    } else
#endif // end remove me
      if(InnerProdVec->getMap()->isSameAs(*right->getColMap())) {
        size_t szNewLeftLocal = std::max(left->getColMap()->getNodeNumElements(), right->getColMap()->getNodeNumElements());
        Teuchos::RCP<std::vector<LocalOrdinal> > NewLeftLocal = Teuchos::rcp(new std::vector<LocalOrdinal>(szNewLeftLocal, Teuchos::as<LocalOrdinal>(right->getColMap()->getMaxLocalIndex()+1)));

        LocalOrdinal j = 0;
        for (size_t i=0; i < left->getColMap()->getNodeNumElements(); i++) {
          while ( (j < Teuchos::as<LocalOrdinal>(right->getColMap()->getNodeNumElements())) &&
                  (right->getColMap()->getGlobalElement(j) < left->getColMap()->getGlobalElement(i)) ) j++;
          if (right->getColMap()->getGlobalElement(j) == left->getColMap()->getGlobalElement(i)) {
            (*NewLeftLocal)[i] = j;
          }
        }

        /*for (size_t i=0; i < right->getColMap()->getNodeNumElements(); i++) {
          std::cout << "left col map: " << (*NewLeftLocal)[i] << " GID: " << left->getColMap()->getGlobalElement((*NewLeftLocal)[i]) << " GID: " << right->getColMap()->getGlobalElement(i) << " right col map: " << i << std::endl;
          }*/

        Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);
        Teuchos::RCP<std::vector<Scalar> > temp_array = Teuchos::rcp(new std::vector<Scalar>(right->getColMap()->getMaxLocalIndex()+2, 0.0));

        for(size_t n=0; n<left->getNodeNumRows(); n++) {
          Teuchos::ArrayView<const LocalOrdinal> lindices_left;
          Teuchos::ArrayView<const Scalar> lvals_left;
          Teuchos::ArrayView<const LocalOrdinal> lindices_right;
          Teuchos::ArrayView<const Scalar> lvals_right;

          left->getLocalRowView (n, lindices_left,  lvals_left);
          right->getLocalRowView(n, lindices_right, lvals_right);

          for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++) {
            (*temp_array)[(*NewLeftLocal)[lindices_left[i] ] ] = lvals_left[i];
          }
          for (size_t i=0; i < Teuchos::as<size_t>(lindices_right.size()); i++) {
            InnerProd_local[lindices_right[i]] += (*temp_array)[lindices_right[i] ] * lvals_right[i];
          }
          for (size_t i=0; i < Teuchos::as<size_t>(lindices_left.size()); i++) {
            (*temp_array)[(*NewLeftLocal)[lindices_left[i] ] ] = 0.0;
          }
        }

        // exporter: overlapping map to nonoverlapping map (target map is unique)
        Teuchos::RCP<const Export> exporter =
          ExportFactory::Build(right->getColMap(), right->getDomainMap()); // TODO: change left to right?

        Teuchos::RCP<Vector> nonoverlap =
          VectorFactory::Build(right->getDomainMap()); // TODO: change left to right?

        nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

        // importer: nonoverlapping map to overlapping map

        // importer: source -> target maps
        Teuchos::RCP<const Import > importer =
          ImportFactory::Build(right->getDomainMap(), right->getColMap()); // TODO: change left to right?
        // doImport target->doImport(*source, importer, action)
        InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: map of InnerProdVec must be same as column map of left operator? Error.");
      }

#else // old "safe" code
    if(InnerProdVec->getMap()->isSameAs(*left->getColMap())) {

      Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);

      // declare variables
      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;

      for(size_t n=0; n<left->getNodeNumRows(); n++)
        {

          left->getLocalRowView (n, lindices_left,  lvals_left);
          right->getLocalRowView(n, lindices_right, lvals_right);

          for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
            {
              GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
              for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
                {
                  GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
                  if(left_gid == right_gid)
                    {
                      InnerProd_local[lindices_left[i]] += lvals_left[i]*lvals_right[j];
                      break; // skip remaining gids of right operator
                    }
                }
            }
        }

      // exporter: overlapping map to nonoverlapping map (target map is unique)
      Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(left->getColMap(), left->getDomainMap()); // TODO: change left to right?

      Teuchos::RCP<Vector > nonoverlap =
        VectorFactory::Build(left->getDomainMap()); // TODO: change left to right?

      nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

      // importer: nonoverlapping map to overlapping map

      // importer: source -> target maps
      Teuchos::RCP<const Import > importer =
        ImportFactory::Build(left->getDomainMap(), left->getColMap()); // TODO: change left to right?

      // doImport target->doImport(*source, importer, action)
      InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
    }
    else if(InnerProdVec->getMap()->isSameAs(*right->getColMap())) {
      Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);

      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;

      for(size_t n=0; n<left->getNodeNumRows(); n++)
        {
          left->getLocalRowView(n, lindices_left, lvals_left);
          right->getLocalRowView(n, lindices_right, lvals_right);

          for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
            {
              GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
              for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
                {
                  GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
                  if(left_gid == right_gid)
                    {
                      InnerProd_local[lindices_right[j]] += lvals_left[i]*lvals_right[j];
                      break; // skip remaining gids of right operator
                    }
                }
            }
        }

      // exporter: overlapping map to nonoverlapping map (target map is unique)
      Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(right->getColMap(), right->getDomainMap()); // TODO: change left to right?

      Teuchos::RCP<Vector> nonoverlap =
        VectorFactory::Build(right->getDomainMap()); // TODO: change left to right?

      nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

      // importer: nonoverlapping map to overlapping map

      // importer: source -> target maps
      Teuchos::RCP<const Import > importer =
        ImportFactory::Build(right->getDomainMap(), right->getColMap()); // TODO: change left to right?

      // doImport target->doImport(*source, importer, action)
      InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: map of InnerProdVec must be same as column map of left or right operator? Error.");
    }
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    std::cout << "TODO: remove me" << std::endl;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ReUseDampingParameters(bool bReuse) {
    bReUseRowBasedOmegas_ = bReuse;
  }

#if 0
  case ANORM: {
    // MUEMAT mode (=paper)
    // Minimize with respect to the (A)' A norm.
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' (A' A) D^{-1} A P0)
    //   omega =   ------------------------------------------
    //             diag( P0' A' D^{-1}' ( A'  A) D^{-1} A P0)
    //
    // expensive, since we have to recalculate AP0 due to the lack of an explicit scaling routine for DinvAP0

    // calculate A * Ptent
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Matrix> AP0 = Utils::TwoMatrixMultiply(A, false, Ptent, false, doFillComplete, optimizeStorage);

    // compute A * D^{-1} * A * P0
    RCP<Matrix> ADinvAP0 = Utils::TwoMatrixMultiply(A, false, DinvAPtent, false, doFillComplete, optimizeStorage);

    Numerator = MultiplyAll(AP0, ADinvAP0, GID2localgid);
    Denominator = MultiplySelfAll(ADinvAP0, GID2localgid);
  }
  break;
  case L2NORM: {
    // ML mode 1 (cheapest)
    // Minimize with respect to L2 norm
    //                  diag( P0' D^{-1} A P0)
    //   omega =   -----------------------------
    //             diag( P0' A' D^{-1}' D^{-1} A P0)
    //
    Numerator = MultiplyAll(Ptent, DinvAPtent, GID2localgid);
    Denominator = MultiplySelfAll(DinvAPtent, GID2localgid);
  }
  break;
  case DINVANORM: {
    // ML mode 2
    // Minimize with respect to the (D^{-1} A)' D^{-1} A norm.
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //   omega =   --------------------------------------------------------
    //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //

    // compute D^{-1} * A * D^{-1} * A * P0
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Matrix> DinvADinvAP0 = Utils::TwoMatrixMultiply(A, false, DinvAPtent, false, doFillComplete, optimizeStorage);
    Utils::MyOldScaleMatrix(DinvADinvAP0, diagA, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag

    Numerator = MultiplyAll(DinvAPtent, DinvADinvAP0, GID2localgid);
    Denominator = MultiplySelfAll(DinvADinvAP0, GID2localgid);
  }
  break;
  case ATDINVTPLUSDINVANORM: {
    // ML mode 3 (most expensive)
    //             diag( P0' ( A'D' + DA) D A P0)
    //   omega =   -----------------------------
    //             diag( P0'A'D' ( A'D' + DA) D A P0)
    //
    //             diag( DinvAP0'DinvAP0 + P0'DinvADinvAP0)
    //         =   -----------------------------
    //                2*diag( DinvADinvAP0'DinvAP0)
    //
    //

    // compute D^{-1} * A * D^{-1} * A * P0
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Matrix> DinvADinvAP0 = Utils::TwoMatrixMultiply(A, false, DinvAPtent, false, doFillComplete, optimizeStorage);
    Utils::MyOldScaleMatrix(DinvADinvAP0, diagA, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag

    Numerator = MultiplyAll(Ptent, DinvADinvAP0, GID2localgid);
    RCP<Teuchos::Array<Scalar> > Numerator2= MultiplySelfAll(DinvAPtent, GID2localgid);
    TEUCHOS_TEST_FOR_EXCEPTION(Numerator->size() != Numerator2->size(), Exceptions::RuntimeError, "PgPFactory::ComputeRowBasedOmegas: size of Numerator and Numerator2 different. Error");
    for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
      (*Numerator)[i] += (*Numerator2)[i];
    Denominator = MultiplyAll(DinvAPtent, DinvADinvAP0, GID2localgid);
    for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
      (*Denominator)[i] *= 2.;

  }
  break;

  /////////////////// DEBUG: check for zeros in denominator
  size_t zeros_in_denominator = 0;
  for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
    {
      if((*Denominator)[i] == Teuchos::ScalarTraits<Scalar>::zero()) zeros_in_denominator ++;
    }
  if(zeros_in_denominator>Teuchos::ScalarTraits<Scalar>::zero())
    GetOStream(Warnings0, 0) << "Found " << zeros_in_denominator<< " zeros in Denominator. very suspicious!" << std::endl;

  /////////////////// build ColBasedOmegas
  RCP<Teuchos::ArrayRCP<Scalar> > ColBasedOmegas = Teuchos::rcp(new Teuchos::ArrayRCP<Scalar>(Numerator->size(), Teuchos::ScalarTraits<Scalar>::zero()));
  for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
    {
      (*ColBasedOmegas)[i] = (*Numerator)[i]/(*Denominator)[i];
      if((*ColBasedOmegas)[i] < Teuchos::ScalarTraits<Scalar>::zero())
        (*ColBasedOmegas)[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }
#endif // if 0

} //namespace MueLu

#endif /* MUELU_PGPFACTORY_DEF_HPP */
