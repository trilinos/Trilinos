// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

/*! \file BelosEpetraOperator.cpp
    \brief This file provides the implementation for an Epetra_Operator 
     interface so Belos can be integrated into other codes as an abstract operator.
*/

#include "BelosEpetraOperator.h"

using namespace Belos;

//--------------------------------------------------------------
//
// implementation of the EpetraOperator class.
//
// Constructor.
//
EpetraOperator::EpetraOperator( const RCP<LinearProblem<double,Epetra_MultiVector,Epetra_Operator> >& lp,
				const RCP<Teuchos::ParameterList>& plist,
                                bool initSolnVec )
  : lp_(lp), 
    plist_(plist),
    initSolnVec_(initSolnVec)
{
  // Get the solver's name from the parameter list, use block Gmres by default.
  std::string solver = plist_->get("Solver", "BlockGmres");

  // Create a label for this Epetra_Operator.
  std::string solver_name = "Belos " + solver + " Solver";
 
  // Copy std::string to character array.  
  // Not using conversion routine copy() because it's not supported by RW on Janus. (HKT 11/13/2003) 
  Solver.resize(solver_name.length()+1);
  for (int i=0; i<(int)solver_name.length(); i++) {
    Solver[i] = solver_name[i];
  } 
  Solver[solver_name.length()] = 0;

  //
  // Create solver and solve problem.  This is inefficient, an instance of the solver should
  // exist already and just be reset with a new RHS.
  //
  if (solver == "BlockGmres") {
    solver_ = Teuchos::rcp( new BlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator>( lp_, plist_ ) );
  } 
  else if (solver == "PseudoBlockGmres") {
    solver_ = Teuchos::rcp( new PseudoBlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator>( lp_, plist_ ) );
  }
  else if (solver == "BlockCG") {
    solver_ = Teuchos::rcp( new BlockCGSolMgr<double,Epetra_MultiVector,Epetra_Operator>( lp_, plist) );
  }
}

const Epetra_Comm& EpetraOperator::Comm() const 
{ 
  return (lp_->getOperator()->Comm());
}
  
const Epetra_Map& EpetraOperator::OperatorDomainMap() const 
{ 
  return (lp_->getOperator()->OperatorDomainMap());
}

const Epetra_Map& EpetraOperator::OperatorRangeMap() const 
{ 
  return (lp_->getOperator()->OperatorRangeMap());
}

int EpetraOperator::Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
  RCP<const Epetra_MultiVector> vec_X;
  RCP<Epetra_MultiVector> vec_Y;
  vec_X = rcp( &X, false );
  vec_Y = rcp( &Y, false );
  if (initSolnVec_)
    vec_Y->PutScalar( 0.0 );
  lp_->setProblem( vec_Y, vec_X );
  Belos::ReturnType ret = solver_->solve();
  
  if (ret != Converged) 
    return(-1);

  return(0);
}

int EpetraOperator::ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
  RCP<const Epetra_MultiVector> vec_X;
  RCP<Epetra_MultiVector> vec_Y;
  vec_X = rcp( &X, false );
  vec_Y = rcp( &Y, false );
  if (initSolnVec_)
    vec_Y->PutScalar( 0.0 );
  lp_->setProblem( vec_Y, vec_X );
  Belos::ReturnType ret = solver_->solve();
  
  if (ret != Converged) 
    return(-1);

  return(0);
}


