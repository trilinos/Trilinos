//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
//@HEADER

#include "AztecOO_ConditionNumber.h"	// class definition

#include "Epetra_Map.h"	        // class data element
#include "Epetra_Vector.h"	// class data element
#include "Epetra_Operator.h"	// class data element
#include "AztecOO.h"            // class data element

// ***********************************************************************
AztecOOConditionNumber::AztecOOConditionNumber() :
  conditionNumber_(0.0),
  domainMap_(NULL),
  rangeMap_(NULL),
  operator_(NULL),
  rhs_(NULL),
  dummy_(NULL),
  solver_(NULL),
  printSolve_(false)
{
  
}

// ***********************************************************************
AztecOOConditionNumber::~AztecOOConditionNumber()
{
  this->freeMemory();
}

// ***********************************************************************
void AztecOOConditionNumber::
initialize(const Epetra_Operator& op,
	   SolverType solverType,
	   int krylovSubspaceSize,
	   bool printSolve)
{
  // cast off const due to aztecoo requiring non-const operator.
  operator_ = const_cast<Epetra_Operator*>(&op);

  // Make sure to clear out any older runs.
  this->freeMemory();

  // Create new vectors and initialize.
  // Call to free memory should zero out map, vectors, and solver.
  domainMap_ = new Epetra_Map(operator_->OperatorDomainMap());
  rangeMap_ = new Epetra_Map(operator_->OperatorRangeMap());
  rhs_ = new Epetra_Vector(*rangeMap_, false);
  rhs_->Random();
  dummy_ = new Epetra_Vector(*domainMap_, true);
  
  // Create solver
  solver_ = new AztecOO(operator_, dummy_, rhs_);
  
  // Set Parameters
  if (solverType == GMRES_)
    solver_->SetAztecOption(AZ_solver, AZ_gmres_condnum);
  else
    solver_->SetAztecOption(AZ_solver, AZ_cg_condnum);

  solver_->SetAztecOption(AZ_kspace, krylovSubspaceSize);

  solver_->SetAztecOption(AZ_precond, AZ_none);

  if (printSolve)
    solver_->SetAztecOption(AZ_output, AZ_last);
  else
    solver_->SetAztecOption(AZ_output, AZ_none);

}

// ***********************************************************************
int AztecOOConditionNumber::computeConditionNumber(int maxIters, double tol)
{
  int status = solver_->Iterate(maxIters, tol);
  const double* aztecStatus = solver_->GetAztecStatus();
  conditionNumber_ = aztecStatus[AZ_condnum];
  
  return status;
}

// ***********************************************************************
double AztecOOConditionNumber::getConditionNumber()
{
  return conditionNumber_;
}

// ***********************************************************************
void AztecOOConditionNumber::freeMemory()
{
  delete solver_;
  solver_ = NULL;
  delete dummy_;
  dummy_ = NULL;
  delete rhs_;
  rhs_ = NULL;
  delete domainMap_;
  domainMap_ = NULL;
  delete rangeMap_;
  rangeMap_ = NULL;
}
