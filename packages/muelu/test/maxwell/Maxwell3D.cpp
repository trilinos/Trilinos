// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <map>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_IO.hpp>

// MueLu
#include <MueLu_RefMaxwell.hpp>
#include <MueLu_Maxwell1.hpp>
#include <MueLu_Maxwell_Utils.hpp>
#include <MueLu_TestHelpers_Common.hpp>
#include <MueLu_Exceptions.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

#include <MueLu_TpetraOperator.hpp>

// Belos
#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#endif
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp

// Stratimikos
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
// Stratimikos includes
#include <Stratimikos_LinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>
#endif

// Support for ML interface
#if defined(HAVE_MUELU_ML) and defined(HAVE_MUELU_EPETRA)
#include <Xpetra_EpetraOperator.hpp>
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_RefMaxwell.h"
#endif

// Helper functions for compilation purposes
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct EpetraSolvers_Wrapper {
  static void Generate_ML_MaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& SM,
                                                Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& D0,
                                                Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Kn,
                                                Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& nullspace,
                                                Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& coords,
                                                Teuchos::ParameterList& mueluList,
                                                Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& mlopX) {
    throw std::runtime_error("Template parameter mismatch");
  }

  static void Generate_ML_RefMaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& SM,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& D0,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Ms,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& M0inv,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& M1,
                                                   Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& nullspace,
                                                   Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& node_material,
                                                   Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& coords,
                                                   Teuchos::ParameterList& mueluList,
                                                   Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& mlopX) {
    throw std::runtime_error("Template parameter mismatch");
  }
};

#if defined(HAVE_MUELU_EPETRA)
template <class GlobalOrdinal>
struct EpetraSolvers_Wrapper<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> {
  static void Generate_ML_MaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& SM,
                                                Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& D0,
                                                Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& Kn,
                                                Teuchos::RCP<Xpetra::MultiVector<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& nullspace,
                                                Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<double>::coordinateType, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& coords,
                                                Teuchos::ParameterList& mueluList,
                                                Teuchos::RCP<Xpetra::Operator<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& mlopX) {
#if defined(HAVE_MUELU_ML)
    typedef double SC;
    typedef int LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode NO;
    typedef typename Teuchos::ScalarTraits<SC>::coordinateType coordinate_type;
    typedef typename Xpetra::Matrix<SC, LO, GO, NO> Matrix;

    RCP<const Epetra_CrsMatrix> epetraSM = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(SM);
    RCP<const Epetra_CrsMatrix> epetraD0 = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(D0);
    if (!coords.is_null()) {
      RCP<const Epetra_MultiVector> epetraCoord = MueLu::Utilities<coordinate_type, LO, GO, NO>::MV2EpetraMV(coords);
      if (epetraCoord->NumVectors() > 0) mueluList.set("x-coordinates", (*epetraCoord)[0]);
      if (epetraCoord->NumVectors() > 1) mueluList.set("y-coordinates", (*epetraCoord)[1]);
      if (epetraCoord->NumVectors() > 2) mueluList.set("z-coordinates", (*epetraCoord)[2]);
    }
    if (!nullspace.is_null()) {
      RCP<const Epetra_MultiVector> epetraNullspace = MueLu::Utilities<SC, LO, GO, NO>::MV2EpetraMV(nullspace);
      mueluList.set("null space: dimension", epetraNullspace->NumVectors());
      mueluList.set("null space: vectors", (*epetraNullspace)[0]);
      mueluList.set("null space: type", "pre-computed");
    }
    RCP<const Epetra_CrsMatrix> epetraKn;
    if (Kn.is_null()) {
      RCP<Matrix> temp = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(SM->getRangeMap());
      Xpetra::MatrixMatrix<SC, LO, GO, NO>::Multiply(*SM, false, *D0, false, *temp, true, true);
      RCP<Matrix> Kn2 = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(D0->getDomainMap());
      Xpetra::MatrixMatrix<SC, LO, GO, NO>::Multiply(*D0, true, *temp, false, *Kn2, true, true);
      epetraKn = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(Kn2);
    } else
      epetraKn = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(Kn);

    RCP<Epetra_Operator> mlop = rcp<Epetra_Operator>(new ML_Epetra::MultiLevelPreconditioner(*epetraSM, *epetraD0, *epetraKn, mueluList, true));
#if defined(HAVE_MUELU_BELOS)
    // NOTE: Belos needs the Apply() and AppleInverse() routines of ML swapped.  So...
    mlop = rcp<Belos::EpetraPrecOp>(new Belos::EpetraPrecOp(mlop));
#endif

    mlopX = rcp(new Xpetra::EpetraOperator<GO, NO>(mlop));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "Need ML & Epetra support");
#endif
  }

  static void Generate_ML_RefMaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& SM,
                                                   Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& D0,
                                                   Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& Ms,
                                                   Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& M0inv,
                                                   Teuchos::RCP<Xpetra::Matrix<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& M1,
                                                   Teuchos::RCP<Xpetra::MultiVector<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& nullspace,
                                                   Teuchos::RCP<Xpetra::MultiVector<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& node_material,
                                                   Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<double>::coordinateType, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& coords,
                                                   Teuchos::ParameterList& mueluList,
                                                   Teuchos::RCP<Xpetra::Operator<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& mlopX) {
#if defined(HAVE_MUELU_ML)
    typedef double SC;
    typedef int LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode NO;
    typedef typename Teuchos::ScalarTraits<SC>::coordinateType coordinate_type;

    RCP<const Epetra_CrsMatrix> epetraSM    = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(SM);
    RCP<const Epetra_CrsMatrix> epetraD0    = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(D0);
    RCP<const Epetra_CrsMatrix> epetraM0inv = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(M0inv);
    RCP<const Epetra_CrsMatrix> epetraMs;
    RCP<const Epetra_CrsMatrix> epetraM1 = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(M1);
    if (!Ms.is_null())
      epetraMs = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(Ms);
    else
      epetraMs = epetraM1;
    mueluList.set("D0", epetraD0);
    mueluList.set("Ms", epetraMs);
    mueluList.set("M0inv", epetraM0inv);
    mueluList.set("M1", epetraM1);
    if (!coords.is_null()) {
      RCP<const Epetra_MultiVector> epetraCoord = MueLu::Utilities<coordinate_type, LO, GO, NO>::MV2EpetraMV(coords);
      if (epetraCoord->NumVectors() > 0) mueluList.sublist("refmaxwell: 11list").set("x-coordinates", (*epetraCoord)[0]);
      if (epetraCoord->NumVectors() > 1) mueluList.sublist("refmaxwell: 11list").set("y-coordinates", (*epetraCoord)[1]);
      if (epetraCoord->NumVectors() > 2) mueluList.sublist("refmaxwell: 11list").set("z-coordinates", (*epetraCoord)[2]);
    }
    if (!node_material.is_null()) {
      RCP<const Epetra_MultiVector> epetraMaterial = MueLu::Utilities<coordinate_type, LO, GO, NO>::MV2EpetraMV(node_material);
      mueluList.sublist("refmaxwell: 11list").set("material coordinates", (*epetraMaterial)[0]);
    }
    if (!nullspace.is_null()) {
      RCP<const Epetra_MultiVector> epetraNullspace = MueLu::Utilities<SC, LO, GO, NO>::MV2EpetraMV(nullspace);
      mueluList.sublist("refmaxwell: 11list").set("null space: dimension", epetraNullspace->NumVectors());
      mueluList.sublist("refmaxwell: 11list").set("null space: vectors", (*epetraNullspace)[0]);
      mueluList.sublist("refmaxwell: 11list").set("null space: type", "pre-computed");
    }

    RCP<Epetra_Operator> mlop = rcp<Epetra_Operator>(new ML_Epetra::RefMaxwellPreconditioner(*epetraSM, mueluList, true));
#if defined(HAVE_MUELU_BELOS)
    // NOTE: Belos needs the Apply() and AppleInverse() routines of ML swapped.  So...
    mlop = rcp<Belos::EpetraPrecOp>(new Belos::EpetraPrecOp(mlop));
#endif

    mlopX = rcp(new Xpetra::EpetraOperator<GO, NO>(mlop));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "Need ML & Epetra support");
#endif
  }
};
#endif  // HAVE_MUELU_EPETRA

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool SetupSolve(std::map<std::string, void*> inputs) {
#include <MueLu_UseShortNames.hpp>

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> coordMV;

  RCP<Matrix> SM_Matrix    = *static_cast<RCP<Matrix>*>(inputs["SM"]);
  RCP<Matrix> GmhdA_Matrix = *static_cast<RCP<Matrix>*>(inputs["GmhdA"]);
  RCP<Matrix> D0_Matrix    = *static_cast<RCP<Matrix>*>(inputs["D0"]);
  RCP<Matrix> M1_Matrix    = *static_cast<RCP<Matrix>*>(inputs["M1"]);
  RCP<Matrix> Ms_Matrix    = *static_cast<RCP<Matrix>*>(inputs["Ms"]);
  RCP<Matrix> M0inv_Matrix = *static_cast<RCP<Matrix>*>(inputs["M0inv"]);
  RCP<Matrix> Kn_Matrix    = *static_cast<RCP<Matrix>*>(inputs["Kn"]);

  RCP<coordMV> coords        = *static_cast<RCP<coordMV>*>(inputs["coordinates"]);
  RCP<MultiVector> nullspace = *static_cast<RCP<MultiVector>*>(inputs["nullspace"]);
  RCP<MultiVector> material  = *static_cast<RCP<MultiVector>*>(inputs["material"]);

  RCP<MultiVector> B  = *static_cast<RCP<MultiVector>*>(inputs["B"]);
  RCP<MultiVector> X  = *static_cast<RCP<MultiVector>*>(inputs["X"]);
  RCP<MultiVector> X0 = *static_cast<RCP<MultiVector>*>(inputs["X0"]);

  Teuchos::ParameterList params            = *static_cast<Teuchos::ParameterList*>(inputs["params"]);
  Teuchos::ParameterList belosParams       = *static_cast<Teuchos::ParameterList*>(inputs["belosParams"]);
  Teuchos::ParameterList stratimikosParams = *static_cast<Teuchos::ParameterList*>(inputs["stratimikosParams"]);
  std::string solverName                   = *static_cast<std::string*>(inputs["solverName"]);
  std::string belosSolverType              = *static_cast<std::string*>(inputs["belosSolverType"]);
  std::string precType                     = *static_cast<std::string*>(inputs["precType"]);
  int numResolves                          = *static_cast<int*>(inputs["numResolves"]);
  bool reuse                               = *static_cast<bool*>(inputs["reuse"]);

  RCP<const Teuchos::Comm<int> > comm = *static_cast<RCP<const Teuchos::Comm<int> >*>(inputs["comm"]);
  RCP<Teuchos::FancyOStream> out      = *static_cast<RCP<Teuchos::FancyOStream>*>(inputs["out"]);

  bool success = false;

  auto tm2 = TimeMonitor::getNewTimer("Maxwell: 2 - Build solver and preconditioner");
#ifdef HAVE_MUELU_BELOS
  if (solverName == "Belos") {
    // construct preconditioner
    RCP<Operator> preconditioner;
    if (precType == "MueLu-RefMaxwell") {
      preconditioner = rcp(new MueLu::RefMaxwell<SC, LO, GO, NO>(SM_Matrix, D0_Matrix, Ms_Matrix, M0inv_Matrix,
                                                                 M1_Matrix, nullspace, coords, params));
    } else if (precType == "MueLu-Maxwell1" || precType == "MueLu-Reitzinger") {
      if (GmhdA_Matrix.is_null())  // are we doing MHD as opposed to GMHD?
        preconditioner = rcp(new MueLu::Maxwell1<SC, LO, GO, NO>(SM_Matrix, D0_Matrix, Kn_Matrix, nullspace, coords, params));
      else
        preconditioner = rcp(new MueLu::Maxwell1<SC, LO, GO, NO>(SM_Matrix, D0_Matrix, Kn_Matrix, nullspace, coords, params, GmhdA_Matrix));

    }
#ifdef HAVE_MUELU_EPETRA
    else if (precType == "ML-RefMaxwell") {
      Xpetra::UnderlyingLib lib = *static_cast<Xpetra::UnderlyingLib*>(inputs["lib"]);
      TEUCHOS_ASSERT(lib == Xpetra::UseEpetra);
      EpetraSolvers_Wrapper<SC, LO, GO, NO>::Generate_ML_RefMaxwellPreconditioner(SM_Matrix, D0_Matrix, Ms_Matrix, M0inv_Matrix,
                                                                                  M1_Matrix, nullspace, material, coords, params, preconditioner);
    } else if (precType == "ML-Maxwell") {
      Xpetra::UnderlyingLib lib = *static_cast<Xpetra::UnderlyingLib*>(inputs["lib"]);
      TEUCHOS_ASSERT(lib == Xpetra::UseEpetra);
      EpetraSolvers_Wrapper<SC, LO, GO, NO>::Generate_ML_MaxwellPreconditioner(SM_Matrix, D0_Matrix, Kn_Matrix,
                                                                               nullspace, coords, params, preconditioner);
    }
#endif

    {
      // A test to make sure we can wrap this guy as a MueLu::TpetraOperator
      RCP<Operator> precOp = Teuchos::rcp_dynamic_cast<Operator>(preconditioner);
      MueLu::TpetraOperator<SC, LO, GO, NO> OpT(precOp);
    }

    // Belos linear problem
    typedef MultiVector MV;
    typedef Belos::OperatorT<MV> OP;
    RCP<OP> belosOp;
    if (GmhdA_Matrix.is_null())                                                // are we doing MHD as opposed to GMHD?
      belosOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(SM_Matrix));  // Turns a Xpetra::Matrix object into a Belos operator
    else
      belosOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(GmhdA_Matrix));  // Turns a Xpetra::Matrix object into a Belos operator

    RCP<Belos::LinearProblem<SC, MV, OP> > problem = rcp(new Belos::LinearProblem<SC, MV, OP>());
    problem->setOperator(belosOp);
    RCP<OP> belosPrecOp;
    if (precType != "none") {
      belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner));  // Turns a Xpetra::Matrix object into a Belos operator
      problem->setRightPrec(belosPrecOp);
    }
    problem->setProblem(X, B);

    bool set = problem->setProblem();
    if (set == false) {
      *out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return false;
    }

    // Belos solver
    RCP<Belos::SolverManager<SC, MV, OP> > solver;
    RCP<Belos::SolverFactory<SC, MV, OP> > factory = rcp(new Belos::SolverFactory<SC, MV, OP>());
    solver                                         = factory->create(belosSolverType, Teuchos::rcpFromRef(belosParams.sublist(belosSolverType)));

    comm->barrier();
    tm2 = Teuchos::null;

    auto tm3 = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

    // set problem and solve
    solver->setProblem(problem);
    for (int solveno = 0; solveno <= numResolves; solveno++) {
      if (X0.is_null())
        X->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      else
        X = X0;
      Belos::ReturnType status = solver->solve();
      int iters                = solver->getNumIters();
      success                  = (iters < 50 && status == Belos::Converged);
      if (success)
        *out << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        *out << "FAILURE! Belos did not converge fast enough." << std::endl;
    }

    if (reuse) {
      TEUCHOS_ASSERT(precType == "MueLu-RefMaxwell");
      for (int solveno = 0; solveno < 2; solveno++) {
        // SM_Matrix->resumeFill();
        // SM_Matrix->fillComplete();
        if (X0.is_null())
          X->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
        else
          X = X0;
        problem->setProblem(X, B);
        if (GmhdA_Matrix.is_null())  // are we doing MHD as opposed to GMHD?
          Teuchos::rcp_dynamic_cast<MueLu::RefMaxwell<SC, LO, GO, NO> >(preconditioner)->resetMatrix(SM_Matrix);
        else
          Teuchos::rcp_dynamic_cast<MueLu::RefMaxwell<SC, LO, GO, NO> >(preconditioner)->resetMatrix(GmhdA_Matrix);
        Belos::ReturnType status = solver->solve();
        int iters                = solver->getNumIters();
        success                  = (iters < 50 && status == Belos::Converged);
        if (success)
          *out << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
        else
          *out << "FAILURE! Belos did not converge fast enough." << std::endl;
      }
    }
  }
#endif  // HAVE_MUELU_BELOS
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  if (solverName == "Stratimikos") {
    using Teuchos::ParameterList;
    using Teuchos::rcp_dynamic_cast;

    // Add matrices to parameterlist
    Teuchos::ParameterList& precTypesPL = stratimikosParams.sublist("Preconditioner Types");
    std::list<Teuchos::ParameterList*> sublists;
    sublists.push_back(&precTypesPL);
    while (sublists.size() > 0) {
      Teuchos::ParameterList* sublist = sublists.front();
      sublists.pop_front();
      std::list<std::string> keys;
      for (Teuchos::ParameterList::ConstIterator it = sublist->begin(); it != sublist->end(); ++it) {
        const std::string& entryName            = sublist->name(it);
        const Teuchos::ParameterEntry& theEntry = sublist->entry(it);
        if (theEntry.isList()) {
          Teuchos::ParameterList& sl = sublist->sublist(entryName);
          sublists.push_back(&sl);
        } else if (theEntry.isType<std::string>() && Teuchos::getValue<std::string>(theEntry).find("substitute ") != std::string::npos) {
          keys.push_back(entryName);
        }
      }
      for (auto key_it = keys.begin(); key_it != keys.end(); ++key_it) {
        std::string value = sublist->get<std::string>(*key_it).substr(11, std::string::npos);
        if (value == "D0")
          sublist->set(*key_it, D0_Matrix);
        else if (value == "M0inv")
          sublist->set(*key_it, M0inv_Matrix);
        else if (value == "M1")
          sublist->set(*key_it, M1_Matrix);
        else if (value == "Ms")
          sublist->set(*key_it, Ms_Matrix);
        else if (value == "Kn")
          sublist->set(*key_it, Kn_Matrix);
        else if (value == "Coordinates")
          sublist->set(*key_it, coords);
          // else if (*key_it == "Nullspace")
          //     sublist->set(*key_it, Nullspace);
#ifdef HAVE_MUELU_EPETRA
        else if (value == "eD0") {
          auto eD0 = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix, true)->getCrsMatrix(), true)->getEpetra_CrsMatrix();
          sublist->set(*key_it, eD0);
        } else if (value == "eCoordinates")
          sublist->set(*key_it, Teuchos::rcp_dynamic_cast<Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node> >(coords, true)->getEpetra_MultiVector());
#endif
        else if (value == "tD0") {
          auto tD0 = Teuchos::rcp_dynamic_cast<TpetraCrsMatrix>(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix, true)->getCrsMatrix(), true)->getTpetra_CrsMatrix();
          sublist->set(*key_it, tD0);
        } else if (value == "tCoordinates") {
          sublist->set(*key_it, Teuchos::rcp_dynamic_cast<TpetraMultiVector>(coords, true)->getTpetra_MultiVector());
        }
      }
    }

    // Build Thyra linear algebra objects
    RCP<const Thyra::LinearOpBase<Scalar> > thyraA;
    if (GmhdA_Matrix.is_null())  // are we doing MHD as opposed to GMHD?
      thyraA = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(SM_Matrix)->getCrsMatrix());
    else
      thyraA = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(GmhdA_Matrix)->getCrsMatrix());
    RCP<Thyra::VectorBase<Scalar> > thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraVector(X->getVectorNonConst(0)));
    // TODO: Why do we loose a reference when running this with Epetra?
    RCP<const Thyra::VectorBase<Scalar> > thyraB = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraVector(B->getVector(0));

    // Build Stratimikos solver
    Stratimikos::LinearSolverBuilder<Scalar> linearSolverBuilder;                                        // This is the Stratimikos main class (= factory of solver factory).
    Stratimikos::enableMueLuRefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);  // Register MueLu as a Stratimikos preconditioner strategy.
    Stratimikos::enableMueLuMaxwell1<Scalar, LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);
    linearSolverBuilder.setParameterList(rcp(&stratimikosParams, false));  // Setup solver parameters using a Stratimikos parameter list.

    // Build a new "solver factory" according to the previously specified parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    auto precFactory                                                = solverFactory->getPreconditionerFactory();
    auto prec                                                       = precFactory->createPrec();

    // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver.
    Thyra::initializePrec<Scalar>(*precFactory, thyraA, prec.ptr());
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > thyraInverseA = solverFactory->createOp();
    Thyra::initializePreconditionedOp<Scalar>(*solverFactory, thyraA, prec, thyraInverseA.ptr());

    comm->barrier();
    tm2 = Teuchos::null;

    auto tm5 = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

    // Solve Ax = b.
    Thyra::SolveStatus<Scalar> status = Thyra::solve<Scalar>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());
    comm->barrier();
    tm5 = Teuchos::null;

    std::cout << status << std::endl;

    success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);

    if (reuse) {
      for (int solveno = 0; solveno < 2; solveno++) {
        auto tm6 = TimeMonitor::getNewTimer("Maxwell: 4 - Setup Re");
        Thyra::initializePrec<Scalar>(*precFactory, thyraA, prec.ptr());
        comm->barrier();
        tm6 = Teuchos::null;

        auto tm7 = TimeMonitor::getNewTimer("Maxwell: 5 - Solve Re");
        status   = Thyra::solve<Scalar>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());
        comm->barrier();
        tm7 = Teuchos::null;
      }
    }

  }  // HAVE_MUELU_STRATIMIKOS && HAVE_MUELU_THYRA
#endif
  comm->barrier();

  return success;
}  // SetupSolve

template <typename Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  std::string timingsFormat = "table-fixed";
  clp.setOption("time-format", &timingsFormat, "timings format (table-fixed | table-scientific | yaml)");
  double scaling = 1.0;
  clp.setOption("scaling", &scaling, "scale mass term");
  std::string solverName = "Belos";
  clp.setOption("solverName", &solverName,
                "Name of iterative linear solver "
                "to use for solving the linear system. "
                "(\"Belos\" or \"Stratimikos\")");
  std::string belosSolverType = "Block CG";
  clp.setOption("belosSolverType", &belosSolverType, "Name of the Belos linear solver");
  std::string precType = "MueLu-RefMaxwell";
  clp.setOption("precType", &precType, "preconditioner to use (MueLu-RefMaxwell|ML-RefMaxwell|none)");
  std::string xml = "";
  clp.setOption("xml", &xml, "xml file with solver parameters (default: \"Maxwell.xml\")");
  std::string belos_xml = "Belos.xml";
  clp.setOption("belos_xml", &belos_xml, "xml file with Belos solver parameters (default: \"Belos.xml\")");
  std::string stratimikos_xml = "Stratimikos.xml";
  clp.setOption("stratimikos_xml", &stratimikos_xml, "xml file with Stratimikos solver parameters (default: \"Stratimikos.xml\")");
  int numResolves = 0;
  clp.setOption("resolve", &numResolves, "#times to redo solve");
  bool reuse = false;
  clp.setOption("reuse", "no-reuse", &reuse, "test reuse");
  double tol = 1e-10;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  int maxIts = 200;
  clp.setOption("its", &maxIts, "maximum number of solver iterations");
  bool use_stacked_timer = false;
  clp.setOption("stacked-timer",
                "no-stacked-timer", &use_stacked_timer, "use stacked timer");
  bool ensure_kn = false;
  clp.setOption("ensure-kn",
                "no-ensure-kn", &ensure_kn, "generate a kn matrix if the user doesn't provide one");

  std::string S_file, D0_file, M1_file, M0_file;
  if (!Teuchos::ScalarTraits<Scalar>::isComplex) {
    S_file  = "S.mat";
    D0_file = "D0.mat";
    M1_file = "M1.mat";
    M0_file = "M0.mat";
  } else {
    S_file  = "S_complex.mat";
    D0_file = "D0_complex.mat";
    M1_file = "M1_complex.mat";
    M0_file = "M0_complex.mat";
  }
  clp.setOption("S", &S_file);
  std::string SM_file = "";
  clp.setOption("SM", &SM_file);
  std::string Kn_file = "";
  clp.setOption("Kn", &Kn_file);
  clp.setOption("D0", &D0_file);
  clp.setOption("M1", &M1_file);
  std::string Ms_file = "";
  clp.setOption("Ms", &Ms_file);
  clp.setOption("M0", &M0_file);
  std::string M0inv_file = "";
  clp.setOption("M0inv", &M0inv_file);

  std::string coords_file = "coords.mat";
  clp.setOption("coords", &coords_file);
  std::string nullspace_file = "";
  clp.setOption("nullspace", &nullspace_file);
  std::string material_file = "";
  clp.setOption("material", &material_file);

  std::string rhs_file = "";
  clp.setOption("rhs", &rhs_file);
  std::string x0_file = "";
  clp.setOption("x0", &x0_file);
  std::string GmhdA_file = "";
  clp.setOption("GmhdA", &GmhdA_file);

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  if (precType == "MueLu-Reitzinger" || precType == "MueLu-Maxwell1" || precType == "ML-Maxwell") {
    if (SM_file != "")
      M1_file = "";
    M0_file = "";
  }

  if (xml == "") {
    if (precType == "MueLu-RefMaxwell")
      xml = "Maxwell.xml";
    else if (precType == "MueLu-Reitzinger" || precType == "MueLu-Maxwell1")
      xml = "Maxwell_Reitzinger.xml";
    else if (precType == "ML-RefMaxwell")
      xml = "Maxwell_ML.xml";
    else if (precType == "ML-Maxwell")
      xml = "Maxwell_ML1.xml";
    else if (precType == "Hypre")

      xml = "Hypre.xml";
  }

  RCP<Teuchos::StackedTimer> stacked_timer;
  if (use_stacked_timer)
    stacked_timer = rcp(new Teuchos::StackedTimer("Maxwell Driver"));
  TimeMonitor::setStackedTimer(stacked_timer);

  auto globalTimeMonitor = TimeMonitor::getNewTimer("Maxwell: S - Global Time");
  auto tm                = TimeMonitor::getNewTimer("Maxwell: 1 - Read and Build Matrices");

  // Read matrices in from files
  RCP<Matrix> D0_Matrix, SM_Matrix, M1_Matrix, Ms_Matrix, M0inv_Matrix, Kn_Matrix, GmhdA_Matrix;

  // maps for nodal and edge matrices
  RCP<const Map> node_map;
  RCP<const Map> edge_map;
  RCP<const Map> Gmhd_map;
  RCP<const Map> edgeOrGmhdMap;

  // gradient matrix
  try {
    std::string base         = D0_file.substr(0, D0_file.find_last_of('/') + 1);
    std::string D0_filename  = D0_file.substr(D0_file.find_last_of('/') + 1, std::string::npos);
    std::string edgeMap_file = base + "rowmap_" + D0_filename;
    std::string nodeMap_file = base + "domainmap_" + D0_filename;
    std::string colMap_file  = base + "colmap_" + D0_filename;
    node_map                 = Xpetra::IO<SC, LO, GO, NO>::ReadMap(nodeMap_file, lib, comm);
    edge_map                 = Xpetra::IO<SC, LO, GO, NO>::ReadMap(edgeMap_file, lib, comm);
    RCP<const Map> colMap;
    if (comm->getSize() > 1)
      colMap = Xpetra::IO<SC, LO, GO, NO>::ReadMap(colMap_file, lib, comm);
    D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, edge_map, colMap, node_map, edge_map);
  } catch (const std::exception&) {
    // *out << "Skipping D0 maps, because: " << e.what() << std::endl;
    D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, lib, comm);
    node_map  = D0_Matrix->getDomainMap();
    edge_map  = D0_Matrix->getRangeMap();
  }

  // build stiffness plus mass matrix (SM_Matrix)
  if (SM_file == "") {
    // edge stiffness matrix
    RCP<Matrix> S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(S_file, edge_map);
    M1_Matrix            = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
    Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*S_Matrix, false, (SC)1.0, *M1_Matrix, false, scaling, SM_Matrix, *out);
    SM_Matrix->fillComplete();
  } else {
    SM_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(SM_file, edge_map);
    if (M1_file != "")
      M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
  }
  if (GmhdA_file != "") {
    Teuchos::ArrayView<const GO> edgeMapEleList = edge_map->getLocalElementList();
    Teuchos::Array<GO> elementList(edgeMapEleList.size() * 2);
    for (int i = 0; i < edgeMapEleList.size(); i++) {
      elementList[2 * i]     = 2 * edgeMapEleList[i];
      elementList[2 * i + 1] = 2 * edgeMapEleList[i] + 1;
    }
    edgeOrGmhdMap = MapFactory::Build(edge_map->lib(), edge_map->getGlobalNumElements() * 2, elementList, edge_map->getIndexBase(), edge_map->getComm());
    GmhdA_Matrix  = Xpetra::IO<SC, LO, GO, NO>::Read(GmhdA_file, edgeOrGmhdMap);
  } else
    edgeOrGmhdMap = edge_map;

  if (Ms_file != "")
    Ms_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(Ms_file, edge_map);
  else
    Ms_Matrix = M1_Matrix;
  if (Kn_file != "") {
    *out << "User provided Kn matrix" << std::endl;
    Kn_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(Kn_file, node_map);
  } else if (ensure_kn) {
    // The user requested we generate Kn_Matrix with an SpGEMM
    Teuchos::ParameterList params;
    *out << "Making Kn Matrix as requested by yser" << std::endl;
    Kn_Matrix = MueLu::Maxwell_Utils<SC, LO, GO, NO>::PtAPWrapper(SM_Matrix, D0_Matrix, params, "User Kn");
  } else {
    *out << "NOT using a Kn matrix" << std::endl;
  }

  if ((M0inv_file == "") && (M0_file != "")) {
    // nodal mass matrix
    RCP<Matrix> M0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0_file, node_map);
    // build lumped mass matrix inverse (M0inv_Matrix)
    RCP<Vector> diag                    = Utilities::GetLumpedMatrixDiagonal(*M0_Matrix);
    RCP<CrsMatrixWrap> M0inv_MatrixWrap = rcp(new CrsMatrixWrap(node_map, node_map, 0));
    RCP<CrsMatrix> M0inv_CrsMatrix      = M0inv_MatrixWrap->getCrsMatrix();
    Teuchos::ArrayRCP<size_t> rowPtr;
    Teuchos::ArrayRCP<LO> colInd;
    Teuchos::ArrayRCP<SC> values;
    Teuchos::ArrayRCP<const SC> diags = diag->getData(0);
    size_t nodeNumElements            = node_map->getLocalNumElements();
    M0inv_CrsMatrix->allocateAllValues(nodeNumElements, rowPtr, colInd, values);
    SC ONE = Teuchos::ScalarTraits<Scalar>::one();
    for (size_t i = 0; i < nodeNumElements; i++) {
      rowPtr[i] = i;
      colInd[i] = i;
      values[i] = ONE / diags[i];
    }
    rowPtr[nodeNumElements] = nodeNumElements;
    M0inv_CrsMatrix->setAllValues(rowPtr, colInd, values);
    M0inv_CrsMatrix->expertStaticFillComplete(node_map, node_map);
    M0inv_Matrix = Teuchos::rcp_dynamic_cast<Matrix>(M0inv_MatrixWrap);
  } else if (M0inv_file != "") {
    M0inv_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0inv_file, node_map);
  }

  // coordinates
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::ReadMultiVector(coords_file, node_map);

  RCP<MultiVector> nullspace = Teuchos::null;
  if (nullspace_file != "")
    nullspace = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace_file, edge_map);

  RCP<MultiVector> material = Teuchos::null;
  if (material_file != "") {
    material = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(material_file, node_map);
  }

  // set parameters
  Teuchos::ParameterList params, belosParams, stratimikosParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xml, Teuchos::Ptr<Teuchos::ParameterList>(&params), *comm);
  Teuchos::updateParametersFromXmlFileAndBroadcast(belos_xml, Teuchos::Ptr<Teuchos::ParameterList>(&belosParams), *comm);
  Teuchos::updateParametersFromXmlFileAndBroadcast(stratimikos_xml, Teuchos::Ptr<Teuchos::ParameterList>(&stratimikosParams), *comm);
  belosParams.sublist(belosSolverType).set("Maximum Iterations", maxIts);
  belosParams.sublist(belosSolverType).set("Convergence Tolerance", tol);

  // setup LHS, RHS
  RCP<MultiVector> X, X0, B;
  if (rhs_file == "") {
    B                    = MultiVectorFactory::Build(edgeOrGmhdMap, 1);
    RCP<MultiVector> vec = MultiVectorFactory::Build(edgeOrGmhdMap, 1);
    vec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    if (GmhdA_Matrix.is_null())  // are we doing MHD as opposed to GMHD?
      SM_Matrix->apply(*vec, *B);
    else
      GmhdA_Matrix->apply(*vec, *B);
  } else
    B = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(rhs_file, edgeOrGmhdMap);

  X = MultiVectorFactory::Build(edgeOrGmhdMap, 1);
  X->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
  if (x0_file != "")
    X0 = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(x0_file, edgeOrGmhdMap);

  comm->barrier();
  tm = Teuchos::null;

  std::map<std::string, void*> inputs;
  inputs["SM"]    = &SM_Matrix;
  inputs["GmhdA"] = &GmhdA_Matrix;
  inputs["D0"]    = &D0_Matrix;
  inputs["M1"]    = &M1_Matrix;
  inputs["Ms"]    = &Ms_Matrix;
  inputs["M0inv"] = &M0inv_Matrix;
  inputs["Kn"]    = &Kn_Matrix;

  inputs["coordinates"] = &coords;
  inputs["nullspace"]   = &nullspace;
  inputs["material"]    = &material;

  inputs["B"]  = &B;
  inputs["X"]  = &X;
  inputs["X0"] = &X0;

  inputs["params"]            = &params;
  inputs["solverName"]        = &solverName;
  inputs["belosSolverType"]   = &belosSolverType;
  inputs["precType"]          = &precType;
  inputs["belosParams"]       = &belosParams;
  inputs["stratimikosParams"] = &stratimikosParams;
  inputs["numResolves"]       = &numResolves;
  inputs["reuse"]             = &reuse;

  inputs["lib"]  = &lib;
  inputs["comm"] = &comm;
  inputs["out"]  = &out;

  bool success = SetupSolve<Scalar, LocalOrdinal, GlobalOrdinal, Node>(inputs);

  globalTimeMonitor = Teuchos::null;

  if (printTimings) {
    if (use_stacked_timer) {
      stacked_timer->stop("Maxwell Driver");
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = options.output_histogram = options.output_minmax = true;
      stacked_timer->report(*out, comm, options);
    } else {
      RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
      if (timingsFormat == "yaml") {
        reportParams->set("Report format", "YAML");  // "Table" or "YAML"
        reportParams->set("YAML style", "compact");  // "spacious" or "compact"
      }
      reportParams->set("How to merge timer sets", "Union");
      reportParams->set("alwaysWriteLocal", false);
      reportParams->set("writeGlobalStats", true);
      reportParams->set("writeZeroTimers", false);
      // FIXME: no "ignoreZeroTimers"

      const std::string filter = "";

      std::ios_base::fmtflags ff(out->flags());
      if (timingsFormat == "table-fixed")
        *out << std::fixed;
      else
        *out << std::scientific;
      TimeMonitor::report(comm.ptr(), *out, filter, reportParams);
      *out << std::setiosflags(ff);
    }
  }

  TimeMonitor::clearCounters();

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
