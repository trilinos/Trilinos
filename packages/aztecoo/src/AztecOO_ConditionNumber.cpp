//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
