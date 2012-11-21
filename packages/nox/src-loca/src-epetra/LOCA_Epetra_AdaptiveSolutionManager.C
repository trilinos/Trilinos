//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2012) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Glen Hansen (gahanse@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Epetra_AdaptiveSolutionManager.H"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"


using Teuchos::rcp;

LOCA::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(
           const Teuchos::RCP<Teuchos::ParameterList>& params,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_) :
   NOX::Epetra::AdaptiveSolutionManager(params, map_, overlapMap_, overlapJacGraph_)
{
}

void
LOCA::Epetra::AdaptiveSolutionManager::
projectCurrentSolution()
{

  // Epetra_Vector element copy for now

  // Note that this must be called after the currentSolution vector is formed for the new mesh, but prior to
  // building the new solution group (grp).

    *currentSolution = grp->getX();
}

void
LOCA::Epetra::AdaptiveSolutionManager::initialize(
       Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
       Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface_,
       Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
       Teuchos::RCP<LOCA::ParameterVector>& pVector_,
       Teuchos::RCP<LOCA::GlobalData>& globalData_,
       bool createPrec_){

  model = model_;
  interface = interface_;
  piroParams = piroParams_;
  pVector = pVector_;
  globalData = globalData_;
  createPrec = createPrec_;

}

Teuchos::RCP<LOCA::Epetra::Group> 
LOCA::Epetra::AdaptiveSolutionManager::buildSolutionGroup(){

  // Create the Jacobian matrix
  Teuchos::RCP<Epetra_Operator> A;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
  Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq
     = Teuchos::rcp_dynamic_cast<NOX::Epetra::ModelEvaluatorInterface>(interface);
  Teuchos::ParameterList& printParams =
       piroParams->sublist("NOX").sublist("Printing");
  Teuchos::ParameterList& noxstratlsParams =
    piroParams->sublist("NOX").sublist("Direction").
    sublist("Newton").sublist("Stratimikos Linear Solver");
  Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");


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

