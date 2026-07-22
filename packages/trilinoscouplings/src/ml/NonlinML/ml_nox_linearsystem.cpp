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
#include "ml_nox_linearsystem.H"
#include "ml_nox_preconditioner.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_LinearSystem::Ml_Nox_LinearSystem(
                             NOX::EpetraNew::Interface::Jacobian& iJac,
                             Epetra_Operator& J,
                             NOX::EpetraNew::Interface::Preconditioner& iPrec,
                             ML_NOX::Nox_CoarseProblem_Interface* coarseinterface,
                             Epetra_Operator& Prec,
                             const Epetra_Vector& soln,
                             bool ismatrixfree,
                             int level,
                             int printlevel) :
iJac_(iJac),
J_(J),
iPrec_(iPrec),
soln_(soln),
coarseinterface_(coarseinterface)
{
   level_        = level;
   ismatrixfree_ = ismatrixfree;
   printlevel_   = printlevel;
   Precptr_      = 0;
   
   // check whether the correct type of preconditioner was supplied
   Epetra_Operator* tmp = dynamic_cast<Epetra_Operator*>(&Prec);
   if (!tmp)
   {
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::Ml_Nox_LinearSystem:\n"
           << "**ERR**: supplied preconditioner is not an Epetra_Operator!\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // find out a little bit more detailed, just to know....
   ML_NOX::ML_Nox_Preconditioner* tmp1 = dynamic_cast<ML_NOX::ML_Nox_Preconditioner*>(&Prec);
   if (tmp1 && printlevel>9 && tmp->Comm().MyPID()==0)
      cout << "ML (level " << level << "): Preconditioner is a ML_NOX::ML_Nox_Preconditioner\n";
   
   ML_Epetra::MultiLevelOperator* tmp2 = dynamic_cast<ML_Epetra::MultiLevelOperator*>(&Prec);
   if (tmp2 && printlevel>9 && tmp->Comm().MyPID()==0)
      cout << "ML (level " << level << "): Preconditioner is a ML_Epetra::MultiLevelOperator\n";
   
   ML_NOX::ML_Nox_ConstrainedMultiLevelOperator* tmp3 = dynamic_cast<ML_NOX::ML_Nox_ConstrainedMultiLevelOperator*>(&Prec);
   if (tmp3 && printlevel>9 && tmp->Comm().MyPID()==0)
      cout << "ML (level " << level << "): Preconditioner is a ML_NOX::ML_Nox_ConstrainedMultiLevelOperator\n";
   
   Precptr_ = tmp;

   return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_LinearSystem::~Ml_Nox_LinearSystem() 
{
  coarseinterface_ = NULL;
}


/*----------------------------------------------------------------------*
 |  applyRightPreconditioning (public)                       m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::applyRightPreconditioning(
                                           bool useTranspose,
			                   NOX::Parameter::List& params, 
			                   const NOX::Epetra::Vector& input, 
			                   NOX::Epetra::Vector& result) const 
{
   int ierr=1;
   if (!Precptr_)
   {
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyRightPreconditioning:\n"
           << "**ERR**: ptr to Preconditioner is NULL\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   if (useTranspose)
   {
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyRightPreconditioning:\n"
           << "**ERR**: apply transposed of preconditioner is not a good idea\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   const Epetra_Vector& in   = input.getEpetraVector();
   Epetra_Vector&       out  = result.getEpetraVector();

#if 0
   cout << "ML_NOX::Ml_Nox_LinearSystem::applyRightPreconditioning:\n";
   cout << "in on input:\n";
   cout << in;
   cout << "out on input:\n";
   cout << out;
#endif
   
   // Note that this preconditioner is either a ML_NOX::ML_Nox_Preconditioner
   //                                      or a ML_Epetra::MultiLevelOperator
   // both implement an Epetra_Operator so we are happy here anyway

#if 0
   Epetra_Operator* tmp = dynamic_cast<Epetra_Operator*>(Precptr_);
   if (!tmp)
   {
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyRightPreconditioning:\n"
           << "**ERR**: supplied preconditioner is not a Epetra_Operator!\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
#endif

   // find out a little bit more detailed, just to know....
   ML_NOX::ML_Nox_Preconditioner* tmp1 = dynamic_cast<ML_NOX::ML_Nox_Preconditioner*>(Precptr_);
   if (tmp1 && printlevel_>9 && tmp1->Comm().MyPID()==0)
      cout << "ML (level " << level_ << "): Preconditioner is a ML_NOX::ML_Nox_Preconditioner\n";

   if (tmp1)
   {
     ierr = tmp1->ApplyInverse(in,out);
     if (!ierr) 
        return true;
     else       
        return false;
   }
   
   ML_Epetra::MultiLevelOperator* tmp2 = dynamic_cast<ML_Epetra::MultiLevelOperator*>(Precptr_);
   if (tmp2 && printlevel_>9 && tmp2->Comm().MyPID()==0)
      cout << "ML (level " << level_ << "): Preconditioner is a ML_Epetra::MultiLevelOperator\n";

   if (tmp2)
   {
     ierr = tmp2->ApplyInverse(in,out);
     if (!ierr) 
       return true;
     else       
       return false;
   }
   
   ML_NOX::ML_Nox_ConstrainedMultiLevelOperator* tmp3 = dynamic_cast<ML_NOX::ML_Nox_ConstrainedMultiLevelOperator*>(Precptr_);
   if (tmp3 && printlevel_>9 && tmp3->Comm().MyPID()==0)
      cout << "ML (level " << level_ << "): Preconditioner is a ML_Epetra::MultiLevelOperator\n";

   if (tmp3)
   {
     ierr = tmp3->ApplyInverse(in,out);
     if (!ierr) 
       return true;
     else       
       return false;
   }
   else
   {
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyRightPreconditioning:\n"
           << "**ERR**: unknown type of preconditioning operator\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

}

/*----------------------------------------------------------------------*
 |  computeJacobian (public)                                 m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::computeJacobian(Epetra_Vector& x)
{
   if (ismatrixfree_==false)
   {
      bool status = iJac_.computeJacobian(x);
      return status;
   }
   else
   {
      NOX::EpetraNew::MatrixFree* Jac = 
                              dynamic_cast<NOX::EpetraNew::MatrixFree*>(&J_);
      if (!Jac)
      {
         cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::computeJacobian:\n"
              << "**ERR**: Jacobian is not a NOX::EpetraNew::MatrixFree object\n"  
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      bool status = Jac->computeJacobian(x);
      return status;
   }
}

/*----------------------------------------------------------------------*
 |  createPreconditioner (public)                            m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::createPreconditioner(
                                                  Epetra_Vector& x, 
			                          NOX::Parameter::List& p,
			                          bool recomputeGraph) const
{
   if (level_ == 0)
   {
      // on the fine level, the preconditioner can either be a ML_NOX::ML_Nox_Preconditioner 
      //                                                  or a ML_Epetra::MultiLevelOperator
      ML_NOX::ML_Nox_Preconditioner* Prec1 = dynamic_cast<ML_NOX::ML_Nox_Preconditioner*>(Precptr_);
      if (Prec1)
      {
         bool status = Prec1->computePreconditioner(x);
         return status;
      }
      ML_Epetra::MultiLevelOperator* Prec2 = dynamic_cast<ML_Epetra::MultiLevelOperator*>(Precptr_);
      if (Prec2)
         return true;
      ML_NOX::ML_Nox_ConstrainedMultiLevelOperator* Prec3 = dynamic_cast<ML_NOX::ML_Nox_ConstrainedMultiLevelOperator*>(Precptr_);
      if (Prec3)
         return true;
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::createPreconditioner:\n"
           << "**ERR**: type of preconditioner unknown\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      return false;
   }
   else
   {
      // on a coarse level the preconditioner can only be a ML_Epetra::MultiLevelOperator
      // there is no need to compute it here, it's ready
      ML_Epetra::MultiLevelOperator* Prec2 = dynamic_cast<ML_Epetra::MultiLevelOperator*>(Precptr_);
      if (Prec2)
         return true;
      ML_NOX::ML_Nox_ConstrainedMultiLevelOperator* Prec3 = dynamic_cast<ML_NOX::ML_Nox_ConstrainedMultiLevelOperator*>(Precptr_);
      if (Prec3)
         return true;
      cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::createPreconditioner:\n"
           << "**ERR**: type of preconditioner unknown\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      return false;
   }
}

/*----------------------------------------------------------------------*
 |  destroyPreconditioner (public)                           m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::destroyPreconditioner() const
{
   // do nothing, don't want NOX to destroy my preconditioner, as 
   // I prefer doing it myself (using an offset for recomputing the preconditioner)
   return true;
}

/*----------------------------------------------------------------------*
 |  getJacobianOperator (public)                             m.gee 12/04|
 *----------------------------------------------------------------------*/
const Epetra_Operator& ML_NOX::Ml_Nox_LinearSystem::getJacobianOperator() const
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::getJacobianOperator:\n"
        << "**ERR**: not impl.\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return J_;
}

/*----------------------------------------------------------------------*
 |  getJacobianOperator (public)                             m.gee 12/04|
 *----------------------------------------------------------------------*/
Epetra_Operator& ML_NOX::Ml_Nox_LinearSystem::getJacobianOperator()
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::getJacobianOperator:\n"
        << "**ERR**: not impl.\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return J_;
}

/*----------------------------------------------------------------------*
 |  getGeneratedPrecOperator (public)                        m.gee 12/04|
 *----------------------------------------------------------------------*/
const Epetra_Operator& ML_NOX::Ml_Nox_LinearSystem::getGeneratedPrecOperator() const
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::getGeneratedPrecOperator:\n"
        << "**ERR**: not impl.\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return *Precptr_;
}

/*----------------------------------------------------------------------*
 |  getGeneratedPrecOperator (public)                        m.gee 12/04|
 *----------------------------------------------------------------------*/
Epetra_Operator& ML_NOX::Ml_Nox_LinearSystem::getGeneratedPrecOperator()
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::getGeneratedPrecOperator:\n"
        << "**ERR**: not impl.\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return *Precptr_;
}

/*----------------------------------------------------------------------*
 |  setJacobianOperatorForSolve (public)                     m.gee 12/04|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_LinearSystem::setJacobianOperatorForSolve(const Epetra_Operator& solveJacOp)
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::setJacobianOperatorForSolve:\n"
        << "**ERR**: not impl.\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return;
}

/*----------------------------------------------------------------------*
 |  setPrecOperatorForSolve (public)                         m.gee 08/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_LinearSystem::setPrecOperatorForSolve(const Epetra_Operator& solvePrecOp)
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::setPrecOperatorForSolve:\n"
        << "**ERR**: not impl.\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return;
}

/*----------------------------------------------------------------------*
 |  applyJacobian (public)                                   m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::applyJacobian(const NOX::Epetra::Vector& input, 
		                                NOX::Epetra::Vector& result) const 
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyJacobian:\n"
        << "**ERR**: not supposed to be called???????\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return false;
}

/*----------------------------------------------------------------------*
 |  applyJacobianTranspose (public)                          m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::applyJacobianTranspose(const NOX::Epetra::Vector& input, 
			                                 NOX::Epetra::Vector& result) const 
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyJacobianTranspose:\n"
        << "**ERR**: not supposed to be called???????\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return false;
}

/*----------------------------------------------------------------------*
 |  applyJacobianInverse (public)                            m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Ml_Nox_LinearSystem::applyJacobianInverse(NOX::Parameter::List &params, 
			                               const NOX::Epetra::Vector &input, 
			                               NOX::Epetra::Vector &result) 
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::applyJacobianInverse:\n"
        << "**ERR**: not supposed to be called???????\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return false;
}

/*----------------------------------------------------------------------*
 |  resetScaling (public)                                    m.gee 12/04|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_LinearSystem::resetScaling(NOX::EpetraNew::Scaling& scalingObject)
{
   cout << "**ERR**: ML_Epetra::Ml_Nox_LinearSystem::resetScaling:\n"
        << "**ERR**: not supposed to be called???????\n"  
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   return;
}

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
