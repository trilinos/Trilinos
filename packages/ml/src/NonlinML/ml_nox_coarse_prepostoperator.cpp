/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>
#include "ml_nox_coarse_prepostoperator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_CoarsePrePostOperator::Ml_Nox_CoarsePrePostOperator(
                               ML_NOX::Nox_CoarseProblem_Interface& coarseinterface,
                               ML_NOX::Ml_Nox_Fineinterface& fineinterface) :
coarseinterface_(coarseinterface),
fineinterface_(fineinterface)                               
{
  type_ = "Ml_Nox_CoarsePrePostOperator";
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       m.gee 10/05|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_CoarsePrePostOperator::Ml_Nox_CoarsePrePostOperator(const ML_NOX::Ml_Nox_CoarsePrePostOperator& source) :
coarseinterface_(source.coarseinterface_),
fineinterface_(source.fineinterface_)                               
{
  type_ = source.type_;
  return;
}

/*----------------------------------------------------------------------*
 |  clone (public)                                           m.gee 10/05|
 *----------------------------------------------------------------------*/
NOX::Parameter::Arbitrary* ML_NOX::Ml_Nox_CoarsePrePostOperator::clone() const
{
  ML_NOX::Ml_Nox_CoarsePrePostOperator* tmp 
    = new ML_NOX::Ml_Nox_CoarsePrePostOperator(*this);
  return tmp;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_CoarsePrePostOperator::~Ml_Nox_CoarsePrePostOperator()
{ 
   return; 
}

/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPreSolve(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPreSolve called\n";    
#endif
  return;
}
/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPostSolve(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPostSolve called\n";
#endif
  return;
}

/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPreIterate(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPreIterate called\n";

  // get the current solution state out of the nox solver
  const NOX::Epetra::Vector* noxsolnptr =
    dynamic_cast<const NOX::Epetra::Vector*>(&(solver.getSolutionGroup().getX()));
  if( !noxsolnptr )
  {
    cout << "**ERR**: ML_NOX::Ml_Nox_CoarsePrePostOperator::runPreSolve:\n"
         << "**ERR**: cast of const NOX::Epetra::Vector* noxsolnptr failed\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  // get the Epetra_MultiVector out of the nox vector
  // FIXME: is this a view?
  const Epetra_Vector& csol = noxsolnptr->getEpetraVector();
  //cout << csol;
  
  // move that current coarse solution vector to the fine grid
  Epetra_Vector* fsol = coarseinterface_.prolong_this_to_fine(csol);
  //cout << *fsol;
  
  // build a pseudo nox group and a pseudo nox solver to put this solution in
  // build pseudo nox parameters (tell this nox to shut up in the printparams)
  NOX::Parameter::List* nlParams = new NOX::Parameter::List();
  NOX::Parameter::List& printParams = nlParams->sublist("Printing"); 
  printParams.setParameter("Output Information",0);
  
  // build a pseudo group
  NOX::Epetra::Vector nfsol(*fsol,NOX::ShapeCopy,true);
  NOX::EpetraNew::Group group(printParams,fineinterface_,nfsol);  
  
  // build a pseudo convergence test
  NOX::StatusTest::NormF normf(1.0,NOX::StatusTest::NormF::Unscaled);
  
  // build a pseudo nox solver manager
  NOX::Solver::Manager fsolver(group,normf,*nlParams);

  // call the fineinterface that implements the fine prepostoperator
  // in there, the application does whatever it need to do to e.g.
  // enforce constraints
  fineinterface_.runPreIterate(fsolver);
  
  // get the current x out of the fsolver again as 
  const NOX::Epetra::Vector* noxfsolnptr =
    dynamic_cast<const NOX::Epetra::Vector*>(&(fsolver.getSolutionGroup().getX()));
  if( !noxfsolnptr )
  {
    cout << "**ERR**: ML_NOX::Ml_Nox_CoarsePrePostOperator::runPreSolve:\n"
         << "**ERR**: cast of const NOX::Epetra::Vector* noxfsolnptr failed\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  // get the Epetra_MultiVector out of the nox vector
  const Epetra_Vector& fsolconstrained = noxfsolnptr->getEpetraVector();
  
  // move the constrained fine grid current solution up to the current level 
  Epetra_Vector* csolconstrained = coarseinterface_.restrict_fine_to_this(fsolconstrained);

#if 0
  // test what has changed
  for (int i=0; i<csol.MyLength(); ++i)
  {
    double diff = csol[i] - (*csolconstrained)[i];
    if (fabs(diff)>1.0e-08)
      printf("Proc %d: Diff in csol[%d]: %20.10f csol[%d]: %20.10f csolconstrained[%d]: %20.10f\n",
              csol.Comm().MyPID(),i,diff,i,csol[i],i,(*csolconstrained)[i]);
  }
#endif  

  // replace the coarse solution with the constrained one
  const NOX::Abstract::Group& origgroup = solver.getSolutionGroup();
  const NOX::EpetraNew::Group& origgroup2 = dynamic_cast<const NOX::EpetraNew::Group&>(origgroup);
  NOX::EpetraNew::Group& origgroup3 = const_cast<NOX::EpetraNew::Group&>(origgroup2);
  origgroup3.setX(*csolconstrained);

  // tidy up
  delete csolconstrained;
  delete fsol;
  delete nlParams;
#endif  
  return;
}

/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPostIterate(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPostIterate called\n";
#endif
  return;
}


#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
