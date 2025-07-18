// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PGPFACTORY_DEF_HPP
#define MUELU_PGPFACTORY_DEF_HPP

#include <vector>

#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_PgPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Tentative prolongator factory");
  validParamList->set<MinimizationNorm>("Minimization norm", DINVANORM, "Norm to be minimized");
  validParamList->set<bool>("ReUseRowBasedOmegas", false, "Reuse omegas for prolongator for restrictor");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetMinimizationMode(MinimizationNorm minnorm) {
  SetParameter("Minimization norm", ParameterEntry(minnorm));  // revalidate
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MueLu::MinimizationNorm PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMinimizationMode() {
  const ParameterList& pL = GetParameterList();
  return pL.get<MueLu::MinimizationNorm>("Minimization norm");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "A");

  // Get default tentative prolongator factory
  // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
  // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
  RCP<const FactoryBase> initialPFact = GetFactory("P");
  if (initialPFact == Teuchos::null) {
    initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent");
  }
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
  const ParameterList& pL   = GetParameterList();
  bool bReUseRowBasedOmegas = pL.get<bool>("ReUseRowBasedOmegas");
  if (bReUseRowBasedOmegas == true && restrictionMode_ == true) {
    coarseLevel.DeclareInput("RowBasedOmega", this, this);  // RowBasedOmega is calculated by this PgPFactory and requested by this PgPFactory
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Prolongator smoothing (PG-AMG)", coarseLevel);

  // Level Get
  RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");

  // Get default tentative prolongator factory
  // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
  // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
  RCP<const FactoryBase> initialPFact = GetFactory("P");
  if (initialPFact == Teuchos::null) {
    initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent");
  }
  RCP<Matrix> Ptent = coarseLevel.Get<RCP<Matrix> >("P", initialPFact.get());

  /////////////////// switch from A to A^T in restriction mode (necessary as long as implicit transpose not working for Epetra)
  if (restrictionMode_) {
    SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
    A = Utilities::Transpose(*A, true);  // build transpose of A explicitely
  }

  /////////////////// calculate D^{-1} A Ptent (needed for smoothing)
  bool doFillComplete  = true;
  bool optimizeStorage = true;
  RCP<Matrix> DinvAP0  = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *Ptent, false, GetOStream(Statistics2), doFillComplete, optimizeStorage);

  doFillComplete                 = true;
  optimizeStorage                = false;
  Teuchos::ArrayRCP<Scalar> diag = Utilities::GetMatrixDiagonal_arcp(*A);
  Utilities::MyOldScaleMatrix(*DinvAP0, diag, true, doFillComplete, optimizeStorage);  // scale matrix with reciprocal of diag

  /////////////////// calculate local damping factors omega

  Teuchos::RCP<Vector> RowBasedOmega = Teuchos::null;

  const ParameterList& pL   = GetParameterList();
  bool bReUseRowBasedOmegas = pL.get<bool>("ReUseRowBasedOmegas");
  if (restrictionMode_ == false || bReUseRowBasedOmegas == false) {
    // if in prolongation mode: calculate row based omegas
    // if in restriction mode: calculate omegas only if row based omegas are not used from prolongation mode
    ComputeRowBasedOmega(fineLevel, coarseLevel, A, Ptent, DinvAP0, RowBasedOmega);
  }  // if(bReUseRowBasedOmegas == false)
  else {
    // reuse row based omegas, calculated by this factory in the run before (with restrictionMode_ == false)
    RowBasedOmega = coarseLevel.Get<Teuchos::RCP<Vector> >("RowBasedOmega", this);

    // RowBasedOmega is now based on row map of A (not transposed)
    // for restriction we use A^T instead of A
    // -> recommunicate row based omega

    // exporter: overlapping row map to nonoverlapping domain map (target map is unique)
    // since A is already transposed we use the RangeMap of A
    Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(RowBasedOmega->getMap(), A->getRangeMap());

    Teuchos::RCP<Vector> noRowBasedOmega =
        VectorFactory::Build(A->getRangeMap());

    noRowBasedOmega->doExport(*RowBasedOmega, *exporter, Xpetra::INSERT);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Import> importer =
        ImportFactory::Build(A->getRangeMap(), A->getRowMap());

    // doImport target->doImport(*source, importer, action)
    RowBasedOmega->doImport(*noRowBasedOmega, *importer, Xpetra::INSERT);
  }

  Teuchos::ArrayRCP<Scalar> RowBasedOmega_local = RowBasedOmega->getDataNonConst(0);

  /////////////////// prolongator smoothing using local damping parameters omega
  RCP<Matrix> P_smoothed = Teuchos::null;
  Utilities::MyOldScaleMatrix(*DinvAP0, RowBasedOmega_local, false, doFillComplete, optimizeStorage);  // scale matrix with reciprocal of diag

  Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*Ptent, false, Teuchos::ScalarTraits<Scalar>::one(),
                                                                                *DinvAP0, false, -Teuchos::ScalarTraits<Scalar>::one(),
                                                                                P_smoothed, GetOStream(Statistics2));
  P_smoothed->fillComplete(Ptent->getDomainMap(), Ptent->getRangeMap());

  //////////////////// store results in Level

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("printLoadBalancingInfo", true);

  // Level Set
  if (!restrictionMode_) {
    // prolongation factory is in prolongation mode
    Set(coarseLevel, "P", P_smoothed);

    // RfromPFactory used to indicate to TogglePFactory that a factory
    // capable  or producing R can be invoked later. TogglePFactory
    // replaces dummy value with an index into it's array of prolongators
    // pointing to the correct prolongator factory. This is later used by
    // RfromP_Or_TransP to invoke the prolongatorfactory in RestrictionMode
    int dummy = 7;
    Set(coarseLevel, "RfromPfactory", dummy);

    if (IsPrint(Statistics1))
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*P_smoothed, "P", params);

    // NOTE: EXPERIMENTAL
    if (Ptent->IsView("stridedMaps"))
      P_smoothed->CreateView("stridedMaps", Ptent);

  } else {
    // prolongation factory is in restriction mode
    RCP<Matrix> R = Utilities::Transpose(*P_smoothed, true);
    Set(coarseLevel, "R", R);

    if (IsPrint(Statistics1))
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*R, "P", params);

    // NOTE: EXPERIMENTAL
    if (Ptent->IsView("stridedMaps"))
      R->CreateView("stridedMaps", Ptent, true);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeRowBasedOmega(Level& /* fineLevel */, Level& coarseLevel, const RCP<Matrix>& A, const RCP<Matrix>& P0, const RCP<Matrix>& DinvAP0, RCP<Vector>& RowBasedOmega) const {
  FactoryMonitor m(*this, "PgPFactory::ComputeRowBasedOmega", coarseLevel);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  Scalar sZero    = Teuchos::ScalarTraits<Scalar>::zero();
  Magnitude mZero = Teuchos::ScalarTraits<Scalar>::magnitude(sZero);

  Teuchos::RCP<Vector> Numerator   = Teuchos::null;
  Teuchos::RCP<Vector> Denominator = Teuchos::null;

  const ParameterList& pL   = GetParameterList();
  MinimizationNorm min_norm = pL.get<MinimizationNorm>("Minimization norm");

  switch (min_norm) {
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
      bool doFillComplete  = true;
      bool optimizeStorage = false;
      RCP<Matrix> AP0      = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *P0, false, GetOStream(Statistics2), doFillComplete, optimizeStorage);

      // compute A * D^{-1} * A * P0
      RCP<Matrix> ADinvAP0 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *DinvAP0, false, GetOStream(Statistics2), doFillComplete, optimizeStorage);

      Numerator   = VectorFactory::Build(ADinvAP0->getColMap(), true);
      Denominator = VectorFactory::Build(ADinvAP0->getColMap(), true);
      MultiplyAll(AP0, ADinvAP0, Numerator);
      MultiplySelfAll(ADinvAP0, Denominator);
    } break;
    case L2NORM: {
      // ML mode 1 (cheapest)
      // Minimize with respect to L2 norm
      //                  diag( P0' D^{-1} A P0)
      //   omega =   -----------------------------
      //             diag( P0' A' D^{-1}' D^{-1} A P0)
      //
      Numerator   = VectorFactory::Build(DinvAP0->getColMap(), true);
      Denominator = VectorFactory::Build(DinvAP0->getColMap(), true);
      MultiplyAll(P0, DinvAP0, Numerator);
      MultiplySelfAll(DinvAP0, Denominator);
    } break;
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
      bool doFillComplete             = true;
      bool optimizeStorage            = true;
      Teuchos::ArrayRCP<Scalar> diagA = Utilities::GetMatrixDiagonal_arcp(*A);
      RCP<Matrix> DinvADinvAP0        = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *DinvAP0, false, GetOStream(Statistics2), doFillComplete, optimizeStorage);
      Utilities::MyOldScaleMatrix(*DinvADinvAP0, diagA, true, doFillComplete, optimizeStorage);  // scale matrix with reciprocal of diag
      diagA = Teuchos::ArrayRCP<Scalar>();

      Numerator   = VectorFactory::Build(DinvADinvAP0->getColMap(), true);
      Denominator = VectorFactory::Build(DinvADinvAP0->getColMap(), true);
      MultiplyAll(DinvAP0, DinvADinvAP0, Numerator);
      MultiplySelfAll(DinvADinvAP0, Denominator);
    } break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::Build: minimization mode not supported. error");
  }

  //////////// build Column based omegas /////////////
  Teuchos::RCP<Vector> ColBasedOmega =
      VectorFactory::Build(Numerator->getMap() /*DinvAP0->getColMap()*/, true);

  ColBasedOmega->putScalar(-666 /**Teuchos::ScalarTraits<Scalar>::one()*/);

  Teuchos::ArrayRCP<const Scalar> Numerator_local   = Numerator->getData(0);
  Teuchos::ArrayRCP<const Scalar> Denominator_local = Denominator->getData(0);
  Teuchos::ArrayRCP<Scalar> ColBasedOmega_local     = ColBasedOmega->getDataNonConst(0);
  GlobalOrdinal zero_local                          = 0;          // count negative colbased omegas
  GlobalOrdinal nan_local                           = 0;          // count NaNs -> set them to zero
  Magnitude min_local                               = 1000000.0;  // Teuchos::ScalarTraits<Scalar>::one() * (Scalar) 1000000;
  Magnitude max_local                               = 0.0;
  for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(Numerator->getLocalLength()); i++) {
    if (Teuchos::ScalarTraits<Scalar>::magnitude(Denominator_local[i]) == mZero) {
      ColBasedOmega_local[i] = 0.0;  // fallback: nonsmoothed basis function since denominator == 0.0
      nan_local++;
    } else {
      ColBasedOmega_local[i] = Numerator_local[i] / Denominator_local[i];  // default case
    }

    if (Teuchos::ScalarTraits<Scalar>::magnitude(ColBasedOmega_local[i]) < mZero) {  // negative omegas are not valid. set them to zero
      ColBasedOmega_local[i] = Teuchos::ScalarTraits<Scalar>::zero();
      zero_local++;  // count zero omegas
    }

    // handle case that Nominator == Denominator -> Dirichlet bcs in A?
    // fallback if ColBasedOmega == 1 -> very strong smoothing may lead to zero rows in P
    // TAW: this is somewhat nonstandard and a rough fallback strategy to avoid problems
    // also avoid "overshooting" with omega > 0.8
    if (Teuchos::ScalarTraits<Scalar>::magnitude(ColBasedOmega_local[i]) >= 0.8) {
      ColBasedOmega_local[i] = 0.0;
    }

    if (Teuchos::ScalarTraits<Scalar>::magnitude(ColBasedOmega_local[i]) < min_local) {
      min_local = Teuchos::ScalarTraits<Scalar>::magnitude(ColBasedOmega_local[i]);
    }
    if (Teuchos::ScalarTraits<Scalar>::magnitude(ColBasedOmega_local[i]) > max_local) {
      max_local = Teuchos::ScalarTraits<Scalar>::magnitude(ColBasedOmega_local[i]);
    }
  }

  {  // be verbose
    GlobalOrdinal zero_all;
    GlobalOrdinal nan_all;
    Magnitude min_all;
    Magnitude max_all;
    MueLu_sumAll(A->getRowMap()->getComm(), zero_local, zero_all);
    MueLu_sumAll(A->getRowMap()->getComm(), nan_local, nan_all);
    MueLu_minAll(A->getRowMap()->getComm(), min_local, min_all);
    MueLu_maxAll(A->getRowMap()->getComm(), max_local, max_all);

    GetOStream(MueLu::Statistics1, 0) << "PgPFactory: smoothed aggregation (scheme: ";
    switch (min_norm) {
      case ANORM: GetOStream(Statistics1) << "Anorm)" << std::endl; break;
      case L2NORM: GetOStream(Statistics1) << "L2norm)" << std::endl; break;
      case DINVANORM: GetOStream(Statistics1) << "DinvAnorm)" << std::endl; break;
      default: GetOStream(Statistics1) << "unknown)" << std::endl; break;
    }
    GetOStream(Statistics1) << "Damping parameter: min = " << min_all << ", max = " << max_all << std::endl;
    GetOStream(Statistics) << "# negative omegas: " << zero_all << " out of " << ColBasedOmega->getGlobalLength() << " column-based omegas" << std::endl;
    GetOStream(Statistics) << "# NaNs: " << nan_all << " out of " << ColBasedOmega->getGlobalLength() << " column-based omegas" << std::endl;
  }

  if (coarseLevel.IsRequested("ColBasedOmega", this)) {
    coarseLevel.Set("ColBasedOmega", ColBasedOmega, this);
  }

  //////////// build Row based omegas /////////////
  // transform column based omegas to row based omegas
  RowBasedOmega =
      VectorFactory::Build(DinvAP0->getRowMap(), true);

  RowBasedOmega->putScalar(-666);  // TODO bad programming style

  bool bAtLeastOneDefined                       = false;
  Teuchos::ArrayRCP<Scalar> RowBasedOmega_local = RowBasedOmega->getDataNonConst(0);
  for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(A->getLocalNumRows()); row++) {
    Teuchos::ArrayView<const LocalOrdinal> lindices;
    Teuchos::ArrayView<const Scalar> lvals;
    DinvAP0->getLocalRowView(row, lindices, lvals);
    bAtLeastOneDefined = false;
    for (size_t j = 0; j < Teuchos::as<size_t>(lindices.size()); j++) {
      Scalar omega = ColBasedOmega_local[lindices[j]];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(omega) != -666) {  // TODO bad programming style
        bAtLeastOneDefined = true;
        if (Teuchos::ScalarTraits<Scalar>::magnitude(RowBasedOmega_local[row]) == -666)
          RowBasedOmega_local[row] = omega;
        else if (Teuchos::ScalarTraits<Scalar>::magnitude(omega) < Teuchos::ScalarTraits<Scalar>::magnitude(RowBasedOmega_local[row]))
          RowBasedOmega_local[row] = omega;
      }
    }
    if (bAtLeastOneDefined == true) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(RowBasedOmega_local[row]) < mZero)
        RowBasedOmega_local[row] = sZero;
    }
  }

  if (coarseLevel.IsRequested("RowBasedOmega", this)) {
    Set(coarseLevel, "RowBasedOmega", RowBasedOmega);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiplySelfAll(const RCP<Matrix>& Op, Teuchos::RCP<Vector>& InnerProdVec) const {
  // note: InnerProdVec is based on column map of Op
  TEUCHOS_TEST_FOR_EXCEPTION(!InnerProdVec->getMap()->isSameAs(*Op->getColMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplySelfAll: map of InnerProdVec must be same as column map of operator. error");

  Teuchos::ArrayRCP<Scalar> InnerProd_local = InnerProdVec->getDataNonConst(0);

  Teuchos::ArrayView<const LocalOrdinal> lindices;
  Teuchos::ArrayView<const Scalar> lvals;

  for (size_t n = 0; n < Op->getLocalNumRows(); n++) {
    Op->getLocalRowView(n, lindices, lvals);
    for (size_t i = 0; i < Teuchos::as<size_t>(lindices.size()); i++) {
      InnerProd_local[lindices[i]] += lvals[i] * lvals[i];
    }
  }
  InnerProd_local = Teuchos::ArrayRCP<Scalar>();

  // exporter: overlapping map to nonoverlapping map (target map is unique)
  Teuchos::RCP<const Export> exporter =
      ExportFactory::Build(Op->getColMap(), Op->getDomainMap());

  Teuchos::RCP<Vector> nonoverlap =
      VectorFactory::Build(Op->getDomainMap());

  nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

  // importer: nonoverlapping map to overlapping map

  // importer: source -> target maps
  Teuchos::RCP<const Import> importer =
      ImportFactory::Build(Op->getDomainMap(), Op->getColMap());

  // doImport target->doImport(*source, importer, action)
  InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiplyAll(const RCP<Matrix>& left, const RCP<Matrix>& right, Teuchos::RCP<Vector>& InnerProdVec) const {
  TEUCHOS_TEST_FOR_EXCEPTION(!left->getDomainMap()->isSameAs(*right->getDomainMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: domain maps of left and right do not match. Error.");
  TEUCHOS_TEST_FOR_EXCEPTION(!left->getRowMap()->isSameAs(*right->getRowMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: row maps of left and right do not match. Error.");
#if 1   // 1=new "fast code, 0=old "slow", but safe code
#if 0   // not necessary - remove me
    if(InnerProdVec->getMap()->isSameAs(*left->getColMap())) {
      // initialize NewRightLocal vector and assign all entries to
      // left->getColMap()->getLocalNumElements() + 1
      std::vector<LocalOrdinal> NewRightLocal(right->getColMap()->getLocalNumElements(), Teuchos::as<LocalOrdinal>(left->getColMap()->getLocalNumElements()+1));

      LocalOrdinal i = 0;
      for (size_t j=0; j < right->getColMap()->getLocalNumElements(); j++) {
        while ( (i < Teuchos::as<LocalOrdinal>(left->getColMap()->getLocalNumElements())) &&
                (left->getColMap()->getGlobalElement(i) < right->getColMap()->getGlobalElement(j)) ) i++;
        if (left->getColMap()->getGlobalElement(i) == right->getColMap()->getGlobalElement(j)) {
          NewRightLocal[j] = i;
        }
      }

      Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);
      std::vector<Scalar> temp_array(left->getColMap()->getLocalNumElements()+1, 0.0);

      for(size_t n=0; n<right->getLocalNumRows(); n++) {
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
#endif  // end remove me
  if (InnerProdVec->getMap()->isSameAs(*right->getColMap())) {
    size_t szNewLeftLocal                                 = TEUCHOS_MAX(left->getColMap()->getLocalNumElements(), right->getColMap()->getLocalNumElements());
    Teuchos::RCP<std::vector<LocalOrdinal> > NewLeftLocal = Teuchos::rcp(new std::vector<LocalOrdinal>(szNewLeftLocal, Teuchos::as<LocalOrdinal>(right->getColMap()->getMaxLocalIndex() + 1)));

    LocalOrdinal j = 0;
    for (size_t i = 0; i < left->getColMap()->getLocalNumElements(); i++) {
      while ((j < Teuchos::as<LocalOrdinal>(right->getColMap()->getLocalNumElements())) &&
             (right->getColMap()->getGlobalElement(j) < left->getColMap()->getGlobalElement(i))) j++;
      if (right->getColMap()->getGlobalElement(j) == left->getColMap()->getGlobalElement(i)) {
        (*NewLeftLocal)[i] = j;
      }
    }

    /*for (size_t i=0; i < right->getColMap()->getLocalNumElements(); i++) {
      std::cout << "left col map: " << (*NewLeftLocal)[i] << " GID: " << left->getColMap()->getGlobalElement((*NewLeftLocal)[i]) << " GID: " << right->getColMap()->getGlobalElement(i) << " right col map: " << i << std::endl;
      }*/

    Teuchos::ArrayRCP<Scalar> InnerProd_local     = InnerProdVec->getDataNonConst(0);
    Teuchos::RCP<std::vector<Scalar> > temp_array = Teuchos::rcp(new std::vector<Scalar>(right->getColMap()->getMaxLocalIndex() + 2, 0.0));

    for (size_t n = 0; n < left->getLocalNumRows(); n++) {
      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;

      left->getLocalRowView(n, lindices_left, lvals_left);
      right->getLocalRowView(n, lindices_right, lvals_right);

      for (size_t i = 0; i < Teuchos::as<size_t>(lindices_left.size()); i++) {
        (*temp_array)[(*NewLeftLocal)[lindices_left[i]]] = lvals_left[i];
      }
      for (size_t i = 0; i < Teuchos::as<size_t>(lindices_right.size()); i++) {
        InnerProd_local[lindices_right[i]] += (*temp_array)[lindices_right[i]] * lvals_right[i];
      }
      for (size_t i = 0; i < Teuchos::as<size_t>(lindices_left.size()); i++) {
        (*temp_array)[(*NewLeftLocal)[lindices_left[i]]] = 0.0;
      }
    }
    InnerProd_local = Teuchos::ArrayRCP<Scalar>();
    // exporter: overlapping map to nonoverlapping map (target map is unique)
    Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(right->getColMap(), right->getDomainMap());  // TODO: change left to right?

    Teuchos::RCP<Vector> nonoverlap =
        VectorFactory::Build(right->getDomainMap());  // TODO: change left to right?

    nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Import> importer =
        ImportFactory::Build(right->getDomainMap(), right->getColMap());  // TODO: change left to right?
    // doImport target->doImport(*source, importer, action)
    InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: map of InnerProdVec must be same as column map of left operator? Error.");
  }

#else  // old "safe" code
  if (InnerProdVec->getMap()->isSameAs(*left->getColMap())) {
    Teuchos::ArrayRCP<Scalar> InnerProd_local = InnerProdVec->getDataNonConst(0);

    // declare variables
    Teuchos::ArrayView<const LocalOrdinal> lindices_left;
    Teuchos::ArrayView<const Scalar> lvals_left;
    Teuchos::ArrayView<const LocalOrdinal> lindices_right;
    Teuchos::ArrayView<const Scalar> lvals_right;

    for (size_t n = 0; n < left->getLocalNumRows(); n++) {
      left->getLocalRowView(n, lindices_left, lvals_left);
      right->getLocalRowView(n, lindices_right, lvals_right);

      for (size_t i = 0; i < Teuchos::as<size_t>(lindices_left.size()); i++) {
        GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
        for (size_t j = 0; j < Teuchos::as<size_t>(lindices_right.size()); j++) {
          GlobalOrdinal right_gid = right->getColMap()->getGlobalElement(lindices_right[j]);
          if (left_gid == right_gid) {
            InnerProd_local[lindices_left[i]] += lvals_left[i] * lvals_right[j];
            break;  // skip remaining gids of right operator
          }
        }
      }
    }

    // exporter: overlapping map to nonoverlapping map (target map is unique)
    Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(left->getColMap(), left->getDomainMap());  // TODO: change left to right?

    Teuchos::RCP<Vector> nonoverlap =
        VectorFactory::Build(left->getDomainMap());  // TODO: change left to right?

    nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Import> importer =
        ImportFactory::Build(left->getDomainMap(), left->getColMap());  // TODO: change left to right?

    // doImport target->doImport(*source, importer, action)
    InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
  } else if (InnerProdVec->getMap()->isSameAs(*right->getColMap())) {
    Teuchos::ArrayRCP<Scalar> InnerProd_local = InnerProdVec->getDataNonConst(0);

    Teuchos::ArrayView<const LocalOrdinal> lindices_left;
    Teuchos::ArrayView<const Scalar> lvals_left;
    Teuchos::ArrayView<const LocalOrdinal> lindices_right;
    Teuchos::ArrayView<const Scalar> lvals_right;

    for (size_t n = 0; n < left->getLocalNumRows(); n++) {
      left->getLocalRowView(n, lindices_left, lvals_left);
      right->getLocalRowView(n, lindices_right, lvals_right);

      for (size_t i = 0; i < Teuchos::as<size_t>(lindices_left.size()); i++) {
        GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
        for (size_t j = 0; j < Teuchos::as<size_t>(lindices_right.size()); j++) {
          GlobalOrdinal right_gid = right->getColMap()->getGlobalElement(lindices_right[j]);
          if (left_gid == right_gid) {
            InnerProd_local[lindices_right[j]] += lvals_left[i] * lvals_right[j];
            break;  // skip remaining gids of right operator
          }
        }
      }
    }

    // exporter: overlapping map to nonoverlapping map (target map is unique)
    Teuchos::RCP<const Export> exporter =
        ExportFactory::Build(right->getColMap(), right->getDomainMap());  // TODO: change left to right?

    Teuchos::RCP<Vector> nonoverlap =
        VectorFactory::Build(right->getDomainMap());  // TODO: change left to right?

    nonoverlap->doExport(*InnerProdVec, *exporter, Xpetra::ADD);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Import> importer =
        ImportFactory::Build(right->getDomainMap(), right->getColMap());  // TODO: change left to right?

    // doImport target->doImport(*source, importer, action)
    InnerProdVec->doImport(*nonoverlap, *importer, Xpetra::INSERT);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: map of InnerProdVec must be same as column map of left or right operator? Error.");
  }
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& /* fineLevel */, Level& /* coarseLevel */) const {
  std::cout << "TODO: remove me" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReUseDampingParameters(bool bReuse) {
  SetParameter("ReUseRowBasedOmegas", ParameterEntry(bReuse));
}

}  // namespace MueLu

#endif /* MUELU_PGPFACTORY_DEF_HPP */
