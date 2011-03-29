//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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
  Teuchos::RCP<const Epetra_MultiVector> vec_X;
  Teuchos::RCP<Epetra_MultiVector> vec_Y;
  vec_X = Teuchos::rcp( &X, false );
  vec_Y = Teuchos::rcp( &Y, false );
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
  Teuchos::RCP<const Epetra_MultiVector> vec_X;
  Teuchos::RCP<Epetra_MultiVector> vec_Y;
  vec_X = Teuchos::rcp( &X, false );
  vec_Y = Teuchos::rcp( &Y, false );
  if (initSolnVec_)
    vec_Y->PutScalar( 0.0 );
  lp_->setProblem( vec_Y, vec_X );
  Belos::ReturnType ret = solver_->solve();
  
  if (ret != Converged) 
    return(-1);

  return(0);
}


