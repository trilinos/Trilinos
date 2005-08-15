/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolve.hpp"
#include "Teuchos_dyn_cast.hpp"

#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_PASTIX
#include "Amesos_Pastix.h"
#endif
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_PARDISO
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif

namespace Thyra {

// Constructors/initializers/accessors

AmesosLinearOpWithSolveFactory::AmesosLinearOpWithSolveFactory(
    const Amesos::ESolverType               solverType
    ,const Amesos::ERefactorizationPolicy   refactorizationPolicy
    )
  :solverType_(solverType)
  ,refactorizationPolicy_(refactorizationPolicy)
{}

// Overridden from LinearOpWithSolveFactoryBase

bool AmesosLinearOpWithSolveFactory::isCompatible(
  const LinearOpBase<double> &fwdOp
  ) const
{
  const EpetraLinearOpBase *eFwdOp = NULL;
  if( ! (eFwdOp = dynamic_cast<const EpetraLinearOpBase*>(&fwdOp)) )
    return false;
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*eFwdOp->epetra_op()) )
    return false;
  return true;
}

Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
AmesosLinearOpWithSolveFactory::createOp() const
{
  return Teuchos::rcp(new AmesosLinearOpWithSolve());
}

void AmesosLinearOpWithSolveFactory::initializeOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,LinearOpWithSolveBase<double>                             *Op
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  //
  // Get the derived subclasses with checked dynamic casts
  //
  Teuchos::RefCountPtr<const EpetraLinearOpBase>
    epetraFwdOp = Teuchos::rcp_dynamic_cast<const EpetraLinearOpBase>(fwdOp,true);
  AmesosLinearOpWithSolve
    *epetraOp = &Teuchos::dyn_cast<AmesosLinearOpWithSolve>(*Op);
  //
  // Update the amesos solver
  //
  if(epetraOp->get_amesosSolver() == Teuchos::null) {
    //
    // This LOWS object has not be initialized yet so this is where we setup
    // everything from the ground up.
    //
    // Create the linear problem and set the operator
    Teuchos::RefCountPtr<Epetra_LinearProblem>
      epetraLP = Teuchos::rcp(new Epetra_LinearProblem());
    epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp->epetra_op()));
    // Create the concrete solver
    Teuchos::RefCountPtr<Amesos_BaseSolver>
      amesosSolver;
    switch(solverType()) {
      case Thyra::Amesos::LAPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Lapack(*epetraLP));
        break;
#ifdef HAVE_AMESOS_KLU
      case Thyra::Amesos::KLU :
        amesosSolver = Teuchos::rcp(new Amesos_Klu(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_PASTIX
      case Thyra::Amesos::PASTIX :
        amesosSolver = Teuchos::rcp(new Amesos_Pastix(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_MUMPS
      case Thyra::Amesos::MUMPS :
        amesosSolver = Teuchos::rcp(new Amesos_Mumps(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_SCALAPACK
      case Thyra::Amesos::SCALAPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Scalapack(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_UMFPACK
      case Thyra::Amesos::UMFPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Umfpack(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
      case Thyra::Amesos::SUPERLUDIST :
        amesosSolver = Teuchos::rcp(new Amesos_Superludist(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_SUPERLU
      case Thyra::Amesos::SUPERLU :
        amesosSolver = Teuchos::rcp(new Amesos_Superlu(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_DSCPACK
      case Thyra::Amesos::DSCPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Dscpack(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_PARDISO
      case Thyra::Amesos::PARDISO :
        amesosSolver = Teuchos::rcp(new Amesos_Pardiso(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_TAUCS
      case Thyra::Amesos::TAUCS :
        amesosSolver = Teuchos::rcp(new Amesos_Taucs(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_PARAKLETE
      case Thyra::Amesos::PARAKLETE :
        amesosSolver = Teuchos::rcp(new Amesos_Paraklete(*epetraLP));
        break;
#endif
      default:
        TEST_FOR_EXCEPTION(
          true, std::logic_error
          ,"Error, the solver type ID = " << solverType() << " is invalid!"
          );
    }
    // Set the parameters
    if(paramList_.get()) amesosSolver->SetParameters(*paramList_);
    // Do the initial factorization
    amesosSolver->SymbolicFactorization();
    amesosSolver->NumericFactorization();
    // Initialize the LOWS object and we are done!
    epetraOp->initialize(epetraFwdOp,epetraLP,amesosSolver);
  }
  else {
    //
    // This LOWS object has already be initialized once so we must just reset
    // the matrix and refactor it.
    //
    // Get non-const pointers to the linear problem and the amesos solver.
    // These const-casts are just fine since the epetraOp in non-const.
    Epetra_LinearProblem
      *epetraLP = const_cast<Epetra_LinearProblem*>(&*epetraOp->get_epetraLP());
    Amesos_BaseSolver
      *amesosSolver = &*epetraOp->get_amesosSolver();
    // Reset the matrix object
    epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp->epetra_op()));
    // Reset the parameters
    if(paramList_.get()) amesosSolver->SetParameters(*paramList_);
    // Repivot if asked
    if(refactorizationPolicy()==Amesos::REPIVOT_ON_REFACTORIZATION)
      amesosSolver->SymbolicFactorization();
    amesosSolver->NumericFactorization();
    // Reset the forward operator and we are done!
    epetraOp->reset_epetraFwdOp(epetraFwdOp);
  }
}

void AmesosLinearOpWithSolveFactory::uninitializeOp(
  LinearOpWithSolveBase<double>                       *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *fwdOp
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  AmesosLinearOpWithSolve
    *epetraOp = &Teuchos::dyn_cast<AmesosLinearOpWithSolve>(*Op);
  Teuchos::RefCountPtr<const EpetraLinearOpBase>
    epetraFwdOp = epetraOp->extract_epetraFwdOp(); // Will be null if uninitialized!
  if(fwdOp) *fwdOp = epetraFwdOp;                  // It is fine if the client does not want this
}

} // namespace Thyra
