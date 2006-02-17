//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "GenericEpetraProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(GenericEpetraProblem& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag)
{
  return problem.evaluate(flag, &x, &FVec);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
					Epetra_Operator& Jac)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Prec, &x, 0);
}
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x, 
					      Epetra_Operator& Prec,
                       NOX::Parameter::List* precParams)
{
  // Pass through to let the Problem fill its owned matrix
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0);
}
//-----------------------------------------------------------------------------

// This is a new class that may evantually get moved into NOX.  For now,
// this is simply used as a testbed for driving NOX using Matlab

// Declare and define a static data member
map< NOX::Abstract::Group::ReturnType, string > Matlab_Interface::returnMsg;
  
Matlab_Interface::Matlab_Interface(NOX::Solver::Manager &  noxSolver) :
  solver(noxSolver),
  enginePtr( new EpetraExt::EpetraExt_MatlabEngine(comm) ),
  engine(*enginePtr)
{
  returnMsg[ NOX::Abstract::Group::Ok ]            = "Ok"            ;
  returnMsg[ NOX::Abstract::Group::NotDefined ]    = "NotDefined"    ;
  returnMsg[ NOX::Abstract::Group::BadDependency ] = "BadDependency" ;
  returnMsg[ NOX::Abstract::Group::NotConverged ]  = "NotConverged"  ;
  returnMsg[ NOX::Abstract::Group::Failed ]        = "Failed"        ;

  cout << "matlab started\n";

  // Verify we are using an Epetra concrete implemntation
  const NOX::Epetra::Group * testGroup = &(dynamic_cast<const NOX::Epetra::Group &>(solver.getSolutionGroup()));
  if( NULL == testGroup )
  {
    throw "Matlab_Interface ERROR: NOX solver not using Epetra implementation.";
  }

  groupPtr = const_cast<NOX::Epetra::Group*>(testGroup);

  solnPtr = const_cast<Epetra_Vector*>(&(dynamic_cast<const NOX::Epetra::Vector&>(groupPtr->getX()).getEpetraVector()));

}

//-----------------------------------------------------------------------------
Matlab_Interface::~Matlab_Interface()
{
  delete enginePtr;
}

//-----------------------------------------------------------------------------
void Matlab_Interface::interact()
{
  char s [BUFSIZE] ;
  char matlabBuffer [MATLABBUF];

  int err;
  while(1) 
  {
      // Prompt the user and get a string
      printf(">> ");
      if (fgets(s, BUFSIZE, stdin) == NULL) {
          printf("Bye\n");
          break ;
      }
      printf ("command :%s:\n", s) ;

      // Parse for NOX commands
      std::string inStr(s);
      if( inStr[0] == '#' )
      {
        inStr.erase(0,1);
        // Remove any carriage returns
        inStr.replace( inStr.find('\n'), 1, "");

        doCommand( inStr );
      }
      else
      {
        // Send the command to MATLAB
        // output goes to stdout
        err = engine.EvalString(s, matlabBuffer, MATLABBUF);
        if (err != 0) {
            printf("there was an error: %d", err);
                    err = 0;
        }
        else {
            printf("Matlab Output:\n%s", matlabBuffer);
        }
      }
  }

  return;
}

//-----------------------------------------------------------------------------
void Matlab_Interface::doCommand( std::string & command )
{

  NOX::Abstract::Group::ReturnType returnStatus;

  cout << "NOX: " << command << endl;

  try 
  {
    // Query methods

    if( command.find("isF") != string::npos )
    {
      std::string isValid = (groupPtr->isF() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }

    if( command.find("isJacobian") != string::npos )
    {
      std::string isValid = (groupPtr->isJacobian() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }

    if( command.find("isGradient") != string::npos )
    {
      std::string isValid = (groupPtr->isGradient() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }
      
    if( command.find("isNewton") != string::npos )
    {
      std::string isValid = (groupPtr->isNewton() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }
      
    if( command.find("isNormNewtonSolveResidual") != string::npos )
    {
      std::string isValid = (groupPtr->isNormNewtonSolveResidual() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }
      
    if( command.find("isPreconditioner") != string::npos )
    {
      std::string isValid = (groupPtr->isPreconditioner() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }
      
    if( command.find("isConditionNumber") != string::npos )
    {
      std::string isValid = (groupPtr->isConditionNumber() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return;
    }
      
    if( command.find("showValid") != string::npos )
    {
      std::string isValid = (groupPtr->isF() ? "True" : "False" );
      cout << " isF              --> " << isValid << endl;
      isValid = (groupPtr->isJacobian() ? "True" : "False" );
      cout << " isJacobian       --> " << isValid << endl;
      isValid = (groupPtr->isGradient() ? "True" : "False" );
      cout << " isGradient       --> " << isValid << endl;
      isValid = (groupPtr->isNewton() ? "True" : "False" );
      cout << " isNewton         --> " << isValid << endl;
      isValid = (groupPtr->isPreconditioner() ? "True" : "False" );
      cout << " isPreconditioner --> " << isValid << endl;
      return;
    }

    // Scalar value query methods

    if( command.find("getJacobianConditionNumber") != string::npos )
    {
      cout << groupPtr->getJacobianConditionNumber() << endl;
      return;
    }
      
    if( command.find("getNormF") != string::npos )
    {
      cout << groupPtr->getNormF() << endl;
      return;
    }
      
    // Compute methods

    if( command.find("setX") != string::npos )
    {
      command.replace( command.find("setX"), 5, ""); 
      engine.GetMultiVector( command.c_str(), *solnPtr );
      groupPtr->setX(*solnPtr);
      return;
    }

    if( command.find("computeF") != string::npos )
    {
      returnStatus = groupPtr->computeF();
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return;
    }

    if( command.find("computeJacobian") != string::npos )
    {
      returnStatus = groupPtr->computeJacobian();
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return;
    }

    if( command.find("computeGradient") != string::npos )
    {
      returnStatus = groupPtr->computeGradient();
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return;
    }

    if( command.find("computeNewton") != string::npos )
    {
      const NOX::Parameter::List & const_lsParams = solver.getParameterList().
                                                           sublist("Direction").
                                                           sublist("Newton").
                                                           sublist("Linear Solver");
      NOX::Parameter::List & lsParams = const_cast<NOX::Parameter::List &>(const_lsParams);
      returnStatus = groupPtr->computeNewton(lsParams);
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return;
    }

    // Jacobian operations

    // Get operations

    if( command.find("getX") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getX()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "X" );
      return;
    }

    if( command.find("getF") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getF()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "F" );
      return;
    }

    if( command.find("getGradient") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getGradient()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "Gradient" );
      return;
    }

    if( command.find("getNewton") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getNewton()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "Newton" );
      return;
    }

    if( command.find("getJacobian") != string::npos )
    {
      Epetra_Operator * jacOp = (groupPtr->getLinearSystem()->getJacobianOperator().get());
      Epetra_RowMatrix * rowMatrix = dynamic_cast<Epetra_RowMatrix *>(jacOp);
      if( rowMatrix )
      {
        engine.PutRowMatrix( *rowMatrix, "Jacobian", false );
        return;
      }
      NOX::Epetra::FiniteDifference * fdOp = dynamic_cast<NOX::Epetra::FiniteDifference *>(jacOp);
      if( fdOp )
      {
        engine.PutRowMatrix( fdOp->getUnderlyingMatrix(), "Jacobian", false );
        return; 
      }

      cout << "Could not get a valid matrix." << endl;
      return;
    }
  }
  catch( const char * msg )
  {
    cout << msg << endl;
    return;
  }

  return;
}
