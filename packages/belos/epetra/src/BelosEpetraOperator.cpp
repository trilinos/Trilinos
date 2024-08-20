// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file BelosEpetraOperator.cpp
    \brief This file provides the implementation for an Epetra_Operator 
     interface so Belos can be integrated into other codes as an abstract operator.
*/

#include "BelosEpetraOperator.h"
#include "BelosSolverFactory.hpp"

namespace Belos {

//--------------------------------------------------------------
//
// implementation of the EpetraOperator class.
//
// Constructor.
//
EpetraOperator::EpetraOperator( const Teuchos::RCP<LinearProblem<double,Epetra_MultiVector,Epetra_Operator> >& lp,
				const Teuchos::RCP<Teuchos::ParameterList>& plist,
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
    solver = "Block Gmres";
  } 
  else if (solver == "PseudoBlockGmres") {
    solver = "Pseudo Block Gmres";
  }
  else if (solver == "BlockCG") {
    solver = "Block CG";
  }
  else if (solver == "PseudoBlockCG") {
    solver = "Pseudo Block CG";
  }

  plist_->remove( "Solver" );
  
  Belos::SolverFactory<double,Epetra_MultiVector,Epetra_Operator> factory;
  solver_ = factory.create( solver, plist );
  solver_->setProblem( lp_ );

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
  Teuchos::RCP<const Epetra_MultiVector> vec_X = Teuchos::rcp( &X, false );
  Teuchos::RCP<Epetra_MultiVector> vec_Y = Teuchos::rcp( &Y, false );
  if (initSolnVec_)
  {
    vec_Y->PutScalar( 0.0 );
    lp_->setInitResVec( vec_X );
  }
  lp_->setProblem( vec_Y, vec_X );
  Belos::ReturnType ret = solver_->solve();
  
  if (ret != Converged) 
    return(-1);

  return(0);
}

int EpetraOperator::ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
  Teuchos::RCP<const Epetra_MultiVector> vec_X = Teuchos::rcp( &X, false );
  Teuchos::RCP<Epetra_MultiVector> vec_Y = Teuchos::rcp( &Y, false );
  if (initSolnVec_)
  {
    vec_Y->PutScalar( 0.0 );
    lp_->setInitResVec( vec_X );
  }
  lp_->setProblem( vec_Y, vec_X ); 
  Belos::ReturnType ret = solver_->solve();
  
  if (ret != Converged) 
    return(-1);

  return(0);
}

} // end namespace Belos
