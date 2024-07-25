// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PeriodicLinearSystem.H"

// Paracont headers
#include "IOContFileUtils.H"

// Trilinos headers
#include "EpetraExt_MultiVectorOut.h"

PeriodicLinearSystem::
PeriodicLinearSystem(
    const Teuchos::RCP <Epetra_Comm> aComm):
  comm(aComm),
  continuableParams(LOCA::ParameterVector()),
  continuationFileParams(Teuchos::null),
  initialGuess(Teuchos::null),
  jacobian(Teuchos::null),
  myGlobalElements(NULL),
  numMyElements(0),
  vectorMap(Teuchos::null)
{

  // Check if we are running with too many processors
  TEUCHOS_TEST_FOR_EXCEPTION( comm->NumProc() > 3,
      std::logic_error,
      "Run this code with 3 processors at most!");

  // Parameters -------------------------------------------
  //
  // Building the list of continuable parameters
  // (only one parameter: p)
  continuableParams.addParameter("p",1.0);

  // In the continuation file, we want p and the three
  // components of x,
  continuationFileParams =
    Teuchos::rcp (new Teuchos::ParameterList());
  continuationFileParams->set<double>("p",
      continuableParams.getValue("p"));
  continuationFileParams->set<double>("x1",0.0);
  continuationFileParams->set<double>("x2",0.0);
  continuationFileParams->set<double>("x3",0.0);

  // Vector distribution ---------------------------------
  //
  // Vector map for three components
  vectorMap = Teuchos::rcp (new Epetra_Map(3,0,*comm));

  // Number of local elements
  numMyElements = vectorMap->NumMyElements();

  // Quering local to global indexing for vectors
  myGlobalElements = vectorMap->MyGlobalElements();

  // Initial guess ---------------------------------------
  initialGuess =
    Teuchos::rcp (new Epetra_Vector(*vectorMap));
  initialGuess->PutScalar(0.0);

  // Jacobian --------------------------------------------
  // Number of nonzero elements
  int * numNonzeros = new int[3];
  for ( int i = 0; i < numMyElements; i++ )
    if ( myGlobalElements[i] == 0 )
      numNonzeros[i] = 3;
    else
      numNonzeros[i] = 2;

  // Create jacobian
  jacobian = Teuchos::rcp (new Epetra_CrsMatrix(Copy,*vectorMap,numNonzeros));

  // Filling the jacobian
  int * indices = new int[3];
  double * values = new double[3];
  for ( int i = 0; i < numMyElements; i++ )
  {
    switch (myGlobalElements[i]) {
      case 0:
    indices[0] = 0;
    values[0]  = 1.0;
    indices[1] = 1;
    values[1]  = 1.0;
    indices[2] = 2;
    values[2]  = 1.0;
    break;
      case 1:
    indices[0] = 0;
    values[0]  = -2.0;
    indices[1] = 1;
    values[1]  = -2.0;
    break;
      case 2:
    indices[0] = 0;
    values[0]  = 2.0;
    indices[1] = 1;
    values[1]  = 1.0;
    break;
      default:
    throw "Thrown exception in PeriodicLinearSystem: myGlobalElements[i] must be between 0 to 2";
    }
    jacobian->InsertGlobalValues(myGlobalElements[i],numNonzeros[i],values,indices);
  }

  // Optimise storage
  jacobian->FillComplete();
  jacobian->OptimizeStorage();

  savedJacobian = Teuchos::rcp (new Epetra_CrsMatrix(*jacobian));

  // Clean
  delete [] numNonzeros;
  delete [] indices;
  delete [] values;

}

PeriodicLinearSystem::
~PeriodicLinearSystem()
{
}

bool PeriodicLinearSystem::
ComputeF(const Epetra_Vector & x, Epetra_Vector & f)
{
  *jacobian = *savedJacobian;

  // Matrix-vector multiplication
  jacobian->Multiply(false,x,f);

  // Getting the continuable parameter
  double p  = continuableParams.getValue("p");

  // Filling in the forcing terms
  for ( int i = 0; i < numMyElements; i++ )
   if (myGlobalElements[i] == 0 || myGlobalElements[i] == 1)
     f[i] -= p;

  return true;
}

bool PeriodicLinearSystem::
ComputeJacF(const Epetra_Vector & x)
{
  *jacobian = *savedJacobian;
  return true;
}

Teuchos::RCP <Epetra_CrsMatrix> PeriodicLinearSystem::
GetJacF() const
{
  return jacobian;
}

Teuchos::RCP <Epetra_Vector> PeriodicLinearSystem::
GetInitialGuess() const
{
  return initialGuess;
}

LOCA::ParameterVector PeriodicLinearSystem::
GetContinuableParams() const
{
  return continuableParams;
}

bool PeriodicLinearSystem::
SetContinuableParameter(std::string label,double value)
{

  // These are the continuable parameters
  if (label == "p")
    continuableParams.setValue("p",value);
  else
    throw "Thrown exception in SetParameter(): label not known";

  return true;
}

bool PeriodicLinearSystem::
UpdateContinuationFile( const std::string & fileName,
          const int & idStep,
          const Teuchos::ParameterList & continuationFileParams)
{

  // Here we are using the coninuation file utilities
  if (comm->MyPID() == 0)
    UpdateContFile(fileName,idStep,continuationFileParams);

  return true;
}

bool PeriodicLinearSystem::
SetContinuationFileParameters(const Epetra_Vector & x)
{

  // Parameter p
  continuationFileParams->set<double>("p",continuableParams.getValue("p"));

  // Three components of the solution
  for ( int i = 0; i < numMyElements; i++ )
    switch (myGlobalElements[i]) {
      case 0:
       continuationFileParams->set<double>("x1",x[i]);
       break;
      case 1:
       continuationFileParams->set<double>("x2",x[i]);
       break;
      case 2:
       continuationFileParams->set<double>("x3",x[i]);
       break;
      default:
       throw "Thrown exception in PeriodicLinearSystem: myGlobalElements[i] must be between 0 to 2";
    }

  return true;
}

Teuchos::RCP <Teuchos::ParameterList> PeriodicLinearSystem::
GetContinuationFileParameters()
{
  return continuationFileParams;
}

bool PeriodicLinearSystem::
PrintSolutionFile(const std::string & fileName, const Epetra_Vector & x,
    const Teuchos::ParameterList & xParams)
{

  // Here we are using Trilinos output to Matlab file
  EpetraExt::MultiVectorToMatlabFile(fileName.c_str(),x);

  return true;
}
