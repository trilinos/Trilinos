/*
//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
*/

#include "LOCAInterface.H"

LOCAInterface::
LOCAInterface( Teuchos::RCP <ProblemLOCAPrototype> & aProblem ,
    Teuchos::RCP <ContinuationManager> aContinuationManager):
  continuationManager(aContinuationManager),
  problem(aProblem)
{
}

LOCAInterface::
~LOCAInterface()
{
}
      
bool LOCAInterface::
computeF(const Epetra_Vector& x, Epetra_Vector& f, 
	 const NOX::Epetra::Interface::Required::FillType F)
{
  problem->ComputeF(x,f);
  return true;
}
    
bool LOCAInterface::
computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  problem->ComputeJacF(x);
  return true;
}

void LOCAInterface::
setParameters(const LOCA::ParameterVector& params)
{
  // Setting the continuable parameters
  for(int i = 0; i < params.length(); i++ ) {
    problem->SetContinuableParameter(params.getLabel(i), params.getValue(i));
  }

  return;
}

void LOCAInterface::
printSolution (const Epetra_Vector &x, const double conParam)
{

  // Printing a Solution file ****************************************
  if ( continuationManager->GetSolutionFileAttribute() == 
      ContinuationManager::Print)
  {
    // Setting the parameters to be printed in the solution file
    problem->SetSolutionFileParameters(x);

    // Getting the solution parameter list 
    Teuchos::RCP <Teuchos::ParameterList> solutionFileParams = 
      problem->GetSolutionFileParameters();

    // Retrieving the file name from the continuation Manager
    std::string solutionFileName = continuationManager->GetSolutionFileName();

    // Printing a Solution File
    problem->PrintSolutionFile(solutionFileName,x,*solutionFileParams);
  }

  // Printing a continuation step in the continuation file ***********

  // Setting the parameters to be printed in the continuation file
  problem->SetContinuationFileParameters(x);

  // Getting the continuation file parameter list 
  Teuchos::RCP <Teuchos::ParameterList> continuationFileParams = 
    problem->GetContinuationFileParameters();

  // Getting the continuation file name from the Continuation manager
  std::string continuationFileName = continuationManager->GetContinuationFileName();

  // Getting the step id from the continuation manager
  int stepId = continuationManager->GetStepID();

  // Updating the continuation file
  problem->UpdateContinuationFile(continuationFileName,stepId,*continuationFileParams);

  return;
}

bool LOCAInterface::
computeShiftedMatrix (double alpha, double beta, 
                           const Epetra_Vector &x, Epetra_Operator &A)
{

cout << " AGS HACK -- LOCAInterface::computeShiftedMatrix RETURNS JACOBIAN!!! " << std::endl;
  problem->ComputeJacF(x);
  problem->GetJacF()->Scale(alpha);

  // Need to add  beta * I for ODEs of the form:  u_dot = f(u)

  //problem->ComputeShiftedJacobian(alpha,beta);
  return true;
}

void LOCAInterface::setXdot(const Epetra_Vector& xdot, const double time) {
  // current problem does not depend on xdot or t
  t = time;
}
