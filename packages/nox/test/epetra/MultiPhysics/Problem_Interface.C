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
#include "Problem_Manager.H"

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
      if (fgets(s, BUFSIZE, stdin) == NULL) 
      {
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
        if (err != 0) 
        {
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
bool Matlab_Interface::doCommand( std::string & command )
{

  NOX::Abstract::Group::ReturnType returnStatus;

  cout << "NOX: " << command << endl;

  // Give precedence to mapped commands
  if( userMaps.end() != userMaps.find(command) )
    command = userMaps[command];

  try 
  {

    // Convenience methods

    if( command.find("?") != string::npos )
    {
      cout << "\n\tCommand Summary:\n" << endl;

      cout << "\t\t?\t\tShow commands" << endl;
      cout << "\t\t#map [arg1] [arg2] \tMap a command" << endl;
      cout << "\t\t#isF" << endl;
      cout << "\t\t#isJacobian" << endl;
      cout << "\t\t#isGradient" << endl;
      cout << "\t\t#isNewton" << endl;
      cout << "\t\t#isNormNewtonSolveResidual" << endl;
      cout << "\t\t#showValid" << endl;
      cout << "\t\t#isPreconditioner" << endl;
      cout << "\t\t#isConditionNumber" << endl;
      cout << "\t\t#computeF" << endl;
      cout << "\t\t#computeJacobian" << endl;
      cout << "\t\t#computeGradient" << endl;
      cout << "\t\t#computeNewton" << endl;
      cout << "\t\t#setX" << endl;
      cout << "\t\t#getNormF" << endl;

      if( !userMaps.empty() )
      {
        cout << "\n\tUser Defined (Mapped) Command Summary:\n" << endl;

        for( std::map<string, string>::const_iterator iter = userMaps.begin(); iter != userMaps.end(); ++iter )
          cout << "\t\t#" << (*iter).first << " --> " << (*iter).second << endl;
      }

      return true;
    }

    if( command.find("map") != string::npos )
    {
      command.replace( command.find("map "), 4, ""); 
      std::string::size_type loc = command.find(" ");
      if( std::string::npos == loc )
      {
        cout << "Could not get two valid arguments." << endl;
        return false;
      }
      std::string arg1 = command.substr(0, loc);
      command.replace( 0, loc+1, ""); 
      std::string arg2 = command;
      cout << "Mapping \"" << arg1 << "\" to \"" << arg2 << "\"" << endl;
      userMaps[ arg1 ] = arg2;
      return true;
    }

    // Query methods

    if( command.find("isF") != string::npos )
    {
      std::string isValid = (groupPtr->isF() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
    }

    if( command.find("isJacobian") != string::npos )
    {
      std::string isValid = (groupPtr->isJacobian() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
    }

    if( command.find("isGradient") != string::npos )
    {
      std::string isValid = (groupPtr->isGradient() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
    }
      
    if( command.find("isNewton") != string::npos )
    {
      std::string isValid = (groupPtr->isNewton() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
    }
      
    if( command.find("isNormNewtonSolveResidual") != string::npos )
    {
      std::string isValid = (groupPtr->isNormNewtonSolveResidual() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
    }
      
    if( command.find("isPreconditioner") != string::npos )
    {
      std::string isValid = (groupPtr->isPreconditioner() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
    }
      
    if( command.find("isConditionNumber") != string::npos )
    {
      std::string isValid = (groupPtr->isConditionNumber() ? "True" : "False" );
      cout << " --> " << isValid << endl;
      return true;
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
      return true;
    }

    // Scalar value query methods

    if( command.find("getJacobianConditionNumber") != string::npos )
    {
      cout << groupPtr->getJacobianConditionNumber() << endl;
      return true;
    }
      
    if( command.find("getNormF") != string::npos )
    {
      cout << groupPtr->getNormF() << endl;
      return true;
    }
      
    // Compute methods

    if( command.find("setX") != string::npos )
    {
      command.replace( command.find("setX"), 5, ""); 
      engine.GetMultiVector( command.c_str(), *solnPtr );
      groupPtr->setX(*solnPtr);
      return true;
    }

    if( command.find("computeF") != string::npos )
    {
      returnStatus = groupPtr->computeF();
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return true;
    }

    if( command.find("computeJacobian") != string::npos )
    {
      returnStatus = groupPtr->computeJacobian();
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return true;
    }

    if( command.find("computeGradient") != string::npos )
    {
      returnStatus = groupPtr->computeGradient();
      cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
      return true;
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
      return true;
    }

    // Jacobian operations

    // Get operations

    if( command.find("getX") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getX()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "X" );
      return true;
    }

    if( command.find("getF") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getF()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "F" );
      return true;
    }

    if( command.find("getGradient") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getGradient()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "Gradient" );
      return true;
    }

    if( command.find("getNewton") != string::npos )
    {
      const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                         (groupPtr->getNewton()).getEpetraVector());
      engine.PutMultiVector( *tmpVec, "Newton" );
      return true;
    }

    if( command.find("getJacobian") != string::npos )
    {
      Epetra_Operator * jacOp = (groupPtr->getLinearSystem()->getJacobianOperator().get());
      Epetra_RowMatrix * rowMatrix = dynamic_cast<Epetra_RowMatrix *>(jacOp);
      if( rowMatrix )
      {
        engine.PutRowMatrix( *rowMatrix, "Jacobian", false );
        return true;
      }
      NOX::Epetra::FiniteDifference * fdOp = dynamic_cast<NOX::Epetra::FiniteDifference *>(jacOp);
      if( fdOp )
      {
        engine.PutRowMatrix( fdOp->getUnderlyingMatrix(), "Jacobian", false );
        return true; 
      }

      cout << "Could not get a valid matrix." << endl;
      return false;
    }

  }
  catch( const char * msg )
  {
    cout << msg << endl;
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------

Coupling_Matlab_Interface::Coupling_Matlab_Interface(Problem_Manager &  manager_) :
  Matlab_Interface( *(manager_.getCompositeSolver()) ),
  problemManager(manager_)
{
}

//-----------------------------------------------------------------------------
bool Coupling_Matlab_Interface::doCommand( std::string & command )
{
  NOX::Abstract::Group::ReturnType returnStatus;

  cout << "Cpl_NOX: " << command << endl;

  // Give precedence to mapped commands
  if( userMaps.end() != userMaps.find(command) )
    command = userMaps[command];

  try 
  {

    // Give priority to mapping method

    if( command.find("map") != string::npos )
      return Matlab_Interface::doCommand( command );

    // Convenience methods

    if( command.find("?") != string::npos )
    {
      cout << "\n\tCoupling Command Summary:\n" << endl;

      cout << "\t\t#pSummary \t Show coupled problems summary" << endl;
      cout << "\t\t#cJ[i j] \t Compute i,j block Jacobian" << endl;
      cout << "\t\t#cF[i] \t Compute i Residual" << endl;
      cout << "\t\t#getF[i j] \t Get i Residual" << endl;
      cout << "\t\t#getJ[i j] \t Get i,j block Jacobian" << endl;
      cout << "\t\t#getAllF \t Get all Residuals" << endl;
      cout << "\t\t#getAllJ \t Get all block Jacobians" << endl;

      return Matlab_Interface::doCommand( command );
    }

    // Query methods

    if( command.find("pSummary") != string::npos )
    {
      problemManager.outputStatus();
      return true;
    }

    // Compute methods

    if( command.find("cJ[" ) != string::npos )
    {
      command.replace( command.find("cJ["), 3, ""); 
      std::string::size_type loc = command.find(" ");
      if( std::string::npos == loc )
      {
        cout << "Could not get two valid arguments." << endl;
        return false;
      }
      std::string arg1 = command.substr(0, loc);
      command.replace( 0, loc+1, ""); 
      loc = command.find("]");
      std::string arg2 = command.substr(0, loc);
      int probId = atoi( arg1.c_str()) ,
          depID =  atoi( arg2.c_str()) ;

      cout << "Computing Jacobian Block " << probId << "," << depID << endl;

      problemManager.computeBlockJacobian( probId, depID );

      return true;
    }

    if( command.find("cF[" ) != string::npos )
    {
      command.replace( command.find("cR["), 3, ""); 
      std::string::size_type loc = command.find("]");
      if( std::string::npos == loc )
      {
        cout << "Could not get a valid argument." << endl;
        return false;
      }
      std::string arg1 = command.substr(0, loc);
      int probId = atoi( arg1.c_str()) ;

      cout << "Computing Residual for Problem " << probId << endl;

      problemManager.computeGroupF( probId );

      return true;
    }

    // Get methods

    if( command.find("getJ[") != string::npos )
    {
      command.replace( command.find("getJ["), 5, ""); 
      std::string::size_type loc = command.find(" ");
      if( std::string::npos == loc )
      {
        cout << "Could not get two valid arguments." << endl;
        return false;
      }
      std::string arg1 = command.substr(0, loc);
      command.replace( 0, loc+1, ""); 
      loc = command.find("]");
      std::string arg2 = command.substr(0, loc);
      int probId = atoi( arg1.c_str()) ,
          depID =  atoi( arg2.c_str()) ;

      Epetra_RowMatrix * rowMatrix = problemManager.getBlockJacobianMatrix( probId, depID );
      if( rowMatrix )
      {
        std::string name = "BJac_" + arg1 + "_" + arg2;
        engine.PutRowMatrix( *rowMatrix, name.c_str(), false );
        cout << "Stored Block Jacobian (" << probId << "," << depID << ") in \""
             << name << "\"" << endl;
        return true;
      }

      cout << "Could not get a valid matrix." << endl;
      return false;
    }

    if( command.find("getAllJ") != string::npos )
    {
      map<int, GenericEpetraProblem*>::iterator problemIter = problemManager.getProblems().begin();
      map<int, GenericEpetraProblem*>::iterator problemLast = problemManager.getProblems().end();

      // Do diagonal blocks
      for( ; problemLast != problemIter; ++problemIter )
      {

        GenericEpetraProblem & problem = *(*problemIter).second;
        int                    probId  = (*problemIter).first;

        Epetra_RowMatrix * rowMatrix = problemManager.getBlockJacobianMatrix( probId );

        if( rowMatrix )
        {
          ostringstream sval1, sval2;
          sval1 << probId << flush;
          std::string name = "BJac_" + sval1.str() + "_" + sval1.str();
          engine.PutRowMatrix( *rowMatrix, name.c_str(), false );
          cout << "Stored Block Jacobian (" << probId << "," << probId << ") in \""
               << name << "\"" << endl;
        }

        // Do off-diagoanl blocks if appropriate
        if( problemManager.useOffBlocks() ) 
        {
#ifdef HAVE_NOX_EPETRAEXT
          for( unsigned int k = 0; k < problem.getDependentProblems().size(); ++k) 
          {
            int depId = problem.getDependentProblems()[k];

            Epetra_RowMatrix * rowMatrix = problemManager.getBlockJacobianMatrix( probId, depId );

            if( rowMatrix )
            {
              ostringstream sval1, sval2;
              sval1 << probId << flush;
              sval2 << depId  << flush;
              std::string name = "BJac_" + sval1.str() + "_" + sval2.str();
              engine.PutRowMatrix( *rowMatrix, name.c_str(), false );
              cout << "Stored Block Jacobian (" << probId << "," << depId << ") in \""
                   << name << "\"" << endl;
            }
          }
#endif
        }
      }
      return true;
    }

    if( command.find("getAllF") != string::npos )
    {
      map<int, GenericEpetraProblem*>::iterator problemIter = problemManager.getProblems().begin();
      map<int, GenericEpetraProblem*>::iterator problemLast = problemManager.getProblems().end();

      // Do diagonal blocks
      for( ; problemLast != problemIter; ++problemIter )
      {
        GenericEpetraProblem & problem = *(*problemIter).second;
        int                    probId  = (*problemIter).first;

        const Epetra_Vector * resVec = problemManager.getResidual( probId );

        if( resVec )
        {
          ostringstream sval1;
          sval1 << probId << flush;
          std::string name = "F_" + sval1.str();
          engine.PutMultiVector( *resVec, name.c_str() );
          cout << "Stored Residual (" << probId << ") in \"" << name << "\"" << endl;
        }
      }
      return true;
    }

  }
  catch( const char * msg )
  {
    cout << msg << endl;
    return false;
  }

  // If no coupling commands found, fall through to nonlinear solver commands
  return Matlab_Interface::doCommand( command );
}
