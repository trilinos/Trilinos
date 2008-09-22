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
    NOX::Epetra::Interface::Required::FillType F)
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
    string solutionFileName = continuationManager->GetSolutionFileName();

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
  string continuationFileName = continuationManager->GetContinuationFileName();

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

cout << " AGS HACK -- LOCAInterface::computeShiftedMatrix RETURNS JACOBIAN!!! " << endl;
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
