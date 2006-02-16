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

Matlab_Interface::Matlab_Interface(NOX::Solver::Manager & solver) :
  enginePtr( new EpetraExt::EpetraExt_MatlabEngine(comm) ),
  engine(*enginePtr)
{
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

Matlab_Interface::~Matlab_Interface()
{
  delete enginePtr;
}

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

void Matlab_Interface::doCommand( std::string & command )
{
  cout << "NOX: " << command << endl;

  map<string, int> commands;
  commands["getX"] = 1;

  cout << "Trying \"" << command << "\" (" << (command=="getX") << ")" << endl;

  if( commands[command] == 1 )
  {
    cout << "Doing " << command << endl;
    engine.PutMultiVector( *solnPtr, "X" );
    return;
  }
  if( command.find("setX") != string::npos )
  {
    cout << "Doing " << command << endl;
    command.replace( command.find("setX"), 5, ""); 

    engine.GetMultiVector( command.c_str(), *solnPtr );
    cout << "New solution vector :\n" << *solnPtr << endl;
  }

  return;
}
