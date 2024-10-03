// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_AdaptiveSolutionManager.hpp"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ParameterList.hpp"


using Teuchos::rcp;
using Teuchos::RCP;

Piro::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(
           const Teuchos::RCP<Teuchos::ParameterList>& appParams,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_) :
  LOCA::Epetra::AdaptiveSolutionManager(appParams->sublist("Discretization").get<int>("Number Of Time Derivatives"),
      map_, overlapMap_, overlapJacGraph_),
  createPrec(false),
  adaptiveMesh(false)
{

  // Create problem PL
  RCP<Teuchos::ParameterList> problemParams =
    Teuchos::sublist(appParams, "Problem", true);

  piroParams =
    sublist(appParams, "Piro", true);

  if(problemParams->isSublist("Adaptation")){ // If the user has specified adaptation on input, grab the sublist

      adaptParams =  sublist(problemParams, "Adaptation", true);
      adaptiveMesh = true;

  }

}

Piro::Epetra::AdaptiveSolutionManager::~AdaptiveSolutionManager(){

}

void
Piro::Epetra::AdaptiveSolutionManager::initialize(
       const Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
       const Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface>& interface_,
       const Teuchos::RCP<LOCA::ParameterVector>& pVector_,
       const Teuchos::RCP<LOCA::GlobalData>& globalData_,
       bool createPrec_){

  model = model_;
  interface = interface_;
  pVector = pVector_;
  globalData = globalData_;
  createPrec = createPrec_;

}

Teuchos::RCP<LOCA::Epetra::Group>
Piro::Epetra::AdaptiveSolutionManager::buildSolutionGroup(){

  // Create the Jacobian matrix
  Teuchos::RCP<Epetra_Operator> A;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
  Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq
     = Teuchos::rcp_dynamic_cast<NOX::Epetra::ModelEvaluatorInterface>(interface);
  Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");
  Teuchos::ParameterList& printParams =
       noxParams.sublist("Printing");
  Teuchos::ParameterList& noxstratlsParams =
    noxParams.sublist("Direction").
    sublist("Newton").sublist("Stratimikos Linear Solver");


  // Inexact Newton must be set in a second sublist when using
  // Stratimikos: This code snippet sets it automatically
  bool inexact = (noxParams.sublist("Direction").sublist("Newton").
                  get("Forcing Term Method", "Constant") != "Constant");
  noxstratlsParams.sublist("NOX Stratimikos Options").
                   set("Use Linear Solve Tolerance From NOX", inexact);

  std::string jacobianSource = piroParams->get("Jacobian Operator", "Have Jacobian");
  bool leanMatrixFree = piroParams->get("Lean Matrix Free",false);

  if (jacobianSource == "Have Jacobian" || jacobianSource == "Matrix-Free") {
    A = model->create_W();
    iJac = interface;
  }
  else if (jacobianSource == "Finite Difference") {
    A = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams,
                                            iReq, *currentSolution));
    iJac = Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(A);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::Epetra::NOXSolver " <<
                 "Invalid value for parameter \" Jacobian Operator\"= " <<
                  jacobianSource << std::endl);


  /* Check flag if separate memory is requested for shifted matrix.
   * This allows 2 Matrix versions of eigensolvers, which save lots
   * of matrix recomputations.  */

  Teuchos::RCP<Epetra_Operator> Ashift=A;
  bool separateMatrixMem = piroParams->get("LOCASolver: Create Second Matrix",false);
  if (separateMatrixMem) {
    Teuchos::RCP<Epetra_CrsMatrix> Acrs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A);
    if (Acrs != Teuchos::null)
     // Ashift = Teuchos::rcp(new Epetra_CrsMatrix(*Acrs));
     Ashift = model->create_W();
    else TEUCHOS_TEST_FOR_EXCEPTION(Acrs == Teuchos::null, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::Epetra::LOCASolver " <<
                 "LOCASolver: Create Second Matrix was requested, but only implemented for CrsMatrix\n");
  }

  // Create separate preconditioner if the model supports it
  /* NOTE: How do we want to decide between using an
   * available preconditioner: (1) If the model supports
   * it, then we use it, or (2) if a parameter list says
   * User_Defined ?  [Below, logic is option (1).]
   */

  Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner> WPrec;
  if (createPrec)
    WPrec = model->create_WPrec();

  // Create the linear system
  // also Build shifted linear system for eigensolver
  Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos> shiftedLinSys;
  if (WPrec != Teuchos::null) {
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

    linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
                  noxstratlsParams, iJac, A, iPrec, WPrec->PrecOp,
                  *currentSolution, WPrec->isAlreadyInverted));
    shiftedLinSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
  		  noxstratlsParams, iJac, Ashift, iPrec, WPrec->PrecOp,
                  *currentSolution, WPrec->isAlreadyInverted));
  }
  else {
     linsys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
                          noxstratlsParams, iJac, A, *currentSolution));
     shiftedLinSys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
                          noxstratlsParams, iJac, Ashift, *currentSolution));
  }

  // Create the LOCA Group
  grp = Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams, iTime,
                                       *currentSolution, linsys,
                                        shiftedLinSys, *pVector));

  grp->setDerivUtils(interface);
  if (separateMatrixMem) grp->declareSeparateMatrixMemory();

  // Saves one resid calculation per solve, but not as safe
  if (leanMatrixFree) grp->disableLinearResidualComputation(true);

  return grp;

}

void
Piro::Epetra::AdaptiveSolutionManager::destroySolutionGroup(){

  // Release the RCPs to prevent a circular reference situation

  model = Teuchos::null;
  interface = Teuchos::null;
  pVector = Teuchos::null;
  globalData = Teuchos::null;
  linsys = Teuchos::null;
  grp = Teuchos::null;

}



