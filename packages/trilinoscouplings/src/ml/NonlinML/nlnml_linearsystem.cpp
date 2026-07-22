// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>


// ----------   User Defined Includes   ----------
#include "nlnml_linearsystem.H"
#include "nlnml_preconditioner.H"
#include "nlnml_ConstrainedMultiLevelOperator.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_LinearSystem::NLNML_LinearSystem(
               RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac,
               RefCountPtr<Epetra_Operator> J,
               RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec,
               RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface,
               RefCountPtr<Epetra_Operator> Prec,
               //const Epetra_Vector& soln,
               bool ismatrixfree,
               int level,
               int printlevel) :
iJac_(iJac),
J_(J),
iPrec_(iPrec),
coarseinterface_(coarseinterface),
Prec_(null),
//soln_(soln),
ismatrixfree_(ismatrixfree),
level_(level),
printlevel_(printlevel)
{
  // check corrct type of preconditioner
  Epetra_Operator* tmp = dynamic_cast<Epetra_Operator*>(Prec.get());
  if (!tmp)
  {
     cout << "**ERR**: NLNML::NLNML_LinearSystem::NLNML_LinearSystem:\n"
          << "**ERR**: supplied preconditioner is not an Epetra_Operator!\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  // incoming type are allowed to be
  // NLNML::NLNML_Preconditioner
  // ML_Epetra::MultiLevelOperator
  // NLNML::NLNML_ConstrainedMultiLevelOperator
  NLNML::NLNML_Preconditioner* tmp1 = dynamic_cast<NLNML::NLNML_Preconditioner*>(Prec.get());
  if (tmp1 && printlevel_>9 && tmp->Comm().MyPID()==0)
     cout << "ML (level " << level << "): Preconditioner in linear system is a NLNML::NLNML_Preconditioner\n";

  ML_Epetra::MultiLevelOperator* tmp2 = dynamic_cast<ML_Epetra::MultiLevelOperator*>(Prec.get());
   if (tmp2 && printlevel>9 && tmp->Comm().MyPID()==0)
     cout << "ML (level " << level << "): Preconditioner in linear system is a ML_Epetra::MultiLevelOperator\n";

  NLNML::NLNML_ConstrainedMultiLevelOperator* tmp3 = dynamic_cast<NLNML::NLNML_ConstrainedMultiLevelOperator*>(Prec.get());
  if (tmp3 && printlevel>9 && tmp->Comm().MyPID()==0)
    cout << "ML (level " << level << "): Preconditioner in linear system is a NLNML::NLNML_ConstrainedMultiLevelOperator\n";

  if (!tmp1 && !tmp2 && !tmp3)
  {
    cout << "**ERR**: NLNML::NLNML_LinearSystem::NLNML_LinearSystem:\n"
         << "**ERR**: supplied preconditioner is not of any recognized type\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  Prec_ = Prec;
  return;
}

/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyRightPreconditioning(
                                 bool useTranspose,
         Teuchos::ParameterList& params,
         const NOX::Epetra::Vector& input,
         NOX::Epetra::Vector& result) const
{
  const Epetra_Vector& in   = input.getEpetraVector();
  Epetra_Vector&       out  = result.getEpetraVector();

  int ierr = Prec_->ApplyInverse(in,out);
  if (!ierr) return true;
  else return false;
}


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::computeJacobian(const NOX::Epetra::Vector& x)
{
  const Epetra_Vector& ex = x.getEpetraVector();
  if (!ismatrixfree_)
  {
    bool err = iJac_->computeJacobian(ex,*J_);
    return err;
  }
  else
  {
    NOX::Epetra::MatrixFree* Jac =
                          dynamic_cast<NOX::Epetra::MatrixFree*>(J_.get());
    if (!Jac)
    {
      cout << "**ERR**: NLNML::NLNML_LinearSystem::computeJacobian:\n"
           << "**ERR**: Jacobian is not a NOX::Epetra::MatrixFree operator\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    Jac->computeJacobian(ex,*Jac);
  }
  return true;
}



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::createPreconditioner(
                            const NOX::Epetra::Vector& x,
          Teuchos::ParameterList& p,
          bool recomputeGraph) const
{
  if (level_==0)
  {
      // on the fine level, the preconditioner can either be
      // NLNML::NLNML_Preconditioner or
      // ML_Epetra::MultiLevelOperator or
      // NLNML::NLNML_ConstrainedMultiLevelOperator
      NLNML::NLNML_Preconditioner* Prec1 =
        dynamic_cast<NLNML::NLNML_Preconditioner*>(Prec_.get());
      if (Prec1)
      {
        const Epetra_Vector& ex = x.getEpetraVector();
        bool err = Prec1->computePreconditioner(ex,*Prec1,NULL);
        return err;
      }

      ML_Epetra::MultiLevelOperator* Prec2 =
        dynamic_cast<ML_Epetra::MultiLevelOperator*>(Prec_.get());
      if (Prec2) return true;

      NLNML::NLNML_ConstrainedMultiLevelOperator* Prec3 =
        dynamic_cast<NLNML::NLNML_ConstrainedMultiLevelOperator*>(Prec_.get());
      if (Prec3) return true;
      cout << "**ERR**: NLNML::NLNML_LinearSystem::createPreconditioner:\n"
           << "**ERR**: Preconditioning operator is of unknown type\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  else
  {
      ML_Epetra::MultiLevelOperator* Prec2 =
        dynamic_cast<ML_Epetra::MultiLevelOperator*>(Prec_.get());
      if (Prec2) return true;

      NLNML::NLNML_ConstrainedMultiLevelOperator* Prec3 =
        dynamic_cast<NLNML::NLNML_ConstrainedMultiLevelOperator*>(Prec_.get());
      if (Prec3) return true;

      cout << "**ERR**: NLNML::NLNML_LinearSystem::createPreconditioner:\n"
           << "**ERR**: Preconditioning operator is of unknown type\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  return true;
}


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_LinearSystem::setJacobianOperatorForSolve(
                const Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp)
{
  cout << "**ERR**: NLNML::NLNML_LinearSystem::setJacobianOperatorForSolve:\n"
       << "**ERR**: not supposed to be called???????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return;
}


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_LinearSystem::setPrecOperatorForSolve(
                 const Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp)
{
  cout << "**ERR**: NLNML::NLNML_LinearSystem::setPrecOperatorForSolve:\n"
       << "**ERR**: not supposed to be called???????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return;
}

/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyJacobian(
                                         const NOX::Epetra::Vector& input,
                                   NOX::Epetra::Vector& result) const
{
  cout << "**ERR**: NLNML::NLNML_LinearSystem::applyJacobian:\n"
       << "**ERR**: not supposed to be called???????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return false;
}


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyJacobianTranspose(
                                         const NOX::Epetra::Vector& input,
                             NOX::Epetra::Vector& result) const
{
  cout << "**ERR**: NLNML::NLNML_LinearSystem::applyJacobianTranspose:\n"
       << "**ERR**: not supposed to be called???????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return false;
}

/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyJacobianInverse(
                                         Teuchos::ParameterList &params,
                             const NOX::Epetra::Vector &input,
                             NOX::Epetra::Vector &result)
{
  cout << "**ERR**: NLNML::NLNML_LinearSystem::applyJacobianInverse:\n"
       << "**ERR**: not supposed to be called???????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return false;
}
#endif









