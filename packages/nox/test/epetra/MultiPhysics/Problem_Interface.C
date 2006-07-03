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
                       Teuchos::ParameterList* precParams)
{
  // Pass through to let the Problem fill its owned matrix
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0);
}
//-----------------------------------------------------------------------------

// This is a new class that may evantually get moved into NOX.  For now,
// this is simply used as a testbed for driving NOX using Matlab

// Declare and define static data member
map< NOX::Abstract::Group::ReturnType, string > CommandBase::returnMsg;
std::istream * CMD_newMacro::is = &(std::cin);
  
Matlab_Interface::Matlab_Interface(NOX::Solver::Manager &  noxSolver) :
  solver(noxSolver),
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

  // Pack our list of supported commands
  mapHandler = new CMD_map( engine, solver )  ;
  commands.push_back( mapHandler );

  CMD_read  * reader = new CMD_read                         ( engine, solver );
  reader->registerDriver( this );
  commands.push_back( reader );

  CMD_write * writer = new CMD_write                        ( engine, solver );
  writer->registerDriver( this );
  commands.push_back( writer );

  CMD_newMacro * newMacro = new CMD_newMacro                ( engine, solver );
  newMacro->registerDriver( this );
  commands.push_back( newMacro );

  CMD_showStack * showStack = new CMD_showStack             ( engine, solver );
  showStack->registerDriver( this );
  commands.push_back( showStack );

  CMD_clearStack * clearStack = new CMD_clearStack          ( engine, solver );
  clearStack->registerDriver( this );
  commands.push_back( clearStack );

  commands.push_back( new CMD_isF                           ( engine, solver ) );
  commands.push_back( new CMD_isJacobian                    ( engine, solver ) );
  commands.push_back( new CMD_isGradient                    ( engine, solver ) );
  commands.push_back( new CMD_isNewton                      ( engine, solver ) );
  commands.push_back( new CMD_isNormNewtonSolveResidual     ( engine, solver ) );
  commands.push_back( new CMD_isPreconditioner              ( engine, solver ) );
  commands.push_back( new CMD_isConditionNumber             ( engine, solver ) );
  commands.push_back( new CMD_showValid                     ( engine, solver ) );
  commands.push_back( new CMD_getJacobianConditionNumber    ( engine, solver ) );
  commands.push_back( new CMD_getNormF                      ( engine, solver ) );
  commands.push_back( new CMD_setX                          ( engine, solver ) );
  commands.push_back( new CMD_computeF                      ( engine, solver ) );
  commands.push_back( new CMD_computeJacobian               ( engine, solver ) );
  commands.push_back( new CMD_computeGradient               ( engine, solver ) );
  commands.push_back( new CMD_computeNewton                 ( engine, solver ) );
  commands.push_back( new CMD_getX                          ( engine, solver ) );
  commands.push_back( new CMD_getF                          ( engine, solver ) );
  commands.push_back( new CMD_getGradient                   ( engine, solver ) );
  commands.push_back( new CMD_getNewton                     ( engine, solver ) );
  commands.push_back( new CMD_getJacobian                   ( engine, solver ) );
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

  bool status;

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

    std::string inStr(s);
    status = doCommand( inStr );

    if( status )
      commandStack.push_back( inStr );
  }

  return;
}

//-----------------------------------------------------------------------------

bool Matlab_Interface::doCommand( std::string command )
{
  char matlabBuffer [MATLABBUF];

  bool status;

  // Parse for NOX commands
  if( command[0] == '#' )
  {
    command.erase(0,1);
    // Remove any carriage returns
    if( command.find('\n') != string::npos )
      command.replace( command.find('\n'), 1, "");

    status = doNOXCommand( command );
  }
  else
  {
    // Send the command to MATLAB - output goes to stdout
    char * c_comm = const_cast<char *>(command.c_str());
    int err = engine.EvalString(c_comm, matlabBuffer, MATLABBUF);
    if (err != 0) 
    {
      printf("there was an error: %d", err);
      status = false;
    }
    else 
    {
      printf("Matlab Output:\n%s", matlabBuffer);
      status = true;
    }
  }
  
  return status;
}

//-----------------------------------------------------------------------------

bool Matlab_Interface::doNOXCommand( std::string & command )
{

  NOX::Abstract::Group::ReturnType returnStatus;

  cout << "NOX     : " << command << endl;

  // Allow for "#" no-op
  if( '#' == command[0] )
    return true;

  // Give priority to mapped commands

  if( mapHandler->getUserMaps().end() != mapHandler->getUserMaps().find(command) )
    command = mapHandler->getUserMaps()[command];

  try 
  {

    // Convenience methods

    if( command.find("?") != string::npos )
    {
      cout << "\n\tCommand Summary:\n" << endl;

      std::list<CommandBase *>::const_iterator iter = commands.begin() ,
                                           iter_end = commands.end()   ;
      for(  ; iter_end != iter; ++iter )
        if( NULL == (*iter)->getCoupler() )
          (*iter)->describe();

      if( !mapHandler->getUserMaps().empty() )
      {
        cout << "\n\tUser Defined (Mapped) Command Summary:\n" << endl;

        for( std::map<string, string>::const_iterator iter = 
            mapHandler->getUserMaps().begin(); iter != mapHandler->getUserMaps().end(); ++iter )
          cout << "\t\t" << (*iter).first << " --> " << (*iter).second << endl;
      }

      return true;
    }

    std::list<CommandBase *>::const_iterator iter = commands.begin() ,
                                         iter_end = commands.end()   ;
    for(  ; iter_end != iter; ++iter )
      if( string::npos != command.find((*iter)->getCommand()) )
        return (*iter)->doCommand(command);

    cout << "\n" << command << " : Not found." << endl;
    return false;

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

  commands.push_back( new CMD_problemSummary             ( engine, problemManager ) );
  commands.push_back( new CMD_getAllX                    ( engine, problemManager ) );
  commands.push_back( new CMD_setXvec                    ( engine, problemManager ) );
  commands.push_back( new CMD_compJac                    ( engine, problemManager ) );
  commands.push_back( new CMD_compRes                    ( engine, problemManager ) );
  commands.push_back( new CMD_doXfers                    ( engine, problemManager ) );
  commands.push_back( new CMD_getJac                     ( engine, problemManager ) );
  commands.push_back( new CMD_getRes                     ( engine, problemManager ) );
  commands.push_back( new CMD_getAllJac                  ( engine, problemManager ) );
  commands.push_back( new CMD_getAllRes                  ( engine, problemManager ) );
}

//-----------------------------------------------------------------------------
bool Coupling_Matlab_Interface::doCommand( std::string command )
{
  NOX::Abstract::Group::ReturnType returnStatus;

  cout << "Cpl_NOX : " << command << endl;

  // Give precedence to mapped commands

  std::map< std::string, std::string >::const_iterator iter = mapHandler->getUserMaps().begin(),
                                                   iter_end = mapHandler->getUserMaps().end()  ;

  for( ; iter_end != iter; ++iter )
  {
    std::string mKey = (*iter).first;
    if( command.find(mKey) != string::npos )
    {
      command.replace( command.find(mKey), mKey.size(), (*iter).second );
      break;
    }
  }

  try 
  {
    // Convenience methods

    if( command.find("?") != string::npos )
    {
      cout << "\n\tCoupling Command Summary:\n" << endl;

      std::list<CommandBase *>::const_iterator iter = commands.begin() ,
                                           iter_end = commands.end()   ;
      for(  ; iter_end != iter; ++iter )
        if( NULL != (*iter)->getCoupler() )
          (*iter)->describe();

    }

    return Matlab_Interface::doCommand( command );
  }
  catch( const char * msg )
  {
    cout << msg << endl;
    return false;
  }

  // If no coupling commands found, fall through to nonlinear solver commands
  return Matlab_Interface::doCommand( command );
}

//-----------------------------------------------------------------------------
//----------------------  Commands for NOX solver -----------------------------
//-----------------------------------------------------------------------------

CommandBase::CommandBase( EpetraExt::EpetraExt_MatlabEngine & engine_, 
                          NOX::Solver::Manager & solver_ ) :
  engine(engine_),
  solver(solver_),
  problemManager(0),
  driver(0)
{
  initialize();
}

CommandBase::CommandBase( EpetraExt::EpetraExt_MatlabEngine & engine_, 
                          Problem_Manager & problemManager_) :
  engine(engine_),
  solver( *problemManager_.getCompositeSolver() ),
  problemManager(&problemManager_),
  driver(0)
{
  initialize();
}

//-----------------------------------------------------------------------------

void CommandBase::initialize()
{
  // Verify we are using an Epetra concrete implemntation
  const NOX::Epetra::Group * testGroup = &(dynamic_cast<const NOX::Epetra::Group &>(solver.getSolutionGroup()));
  if( NULL == testGroup )
  {
    throw "Matlab_Interface ERROR: NOX solver not using Epetra implementation.";
  }

  groupPtr = const_cast<NOX::Epetra::Group*>(testGroup);

  solnPtr = const_cast<Epetra_Vector*>(&(dynamic_cast<const NOX::Epetra::Vector&>(groupPtr->getX()).getEpetraVector()));

  returnMsg[ NOX::Abstract::Group::Ok ]            = "Ok"            ;
  returnMsg[ NOX::Abstract::Group::NotDefined ]    = "NotDefined"    ;
  returnMsg[ NOX::Abstract::Group::BadDependency ] = "BadDependency" ;
  returnMsg[ NOX::Abstract::Group::NotConverged ]  = "NotConverged"  ;
  returnMsg[ NOX::Abstract::Group::Failed ]        = "Failed"        ;

}

//-----------------------------------------------------------------------------

bool CMD_map::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find("map "), 4, ""); 
  std::string::size_type loc = commandLine.find(" ");
  if( std::string::npos == loc )
  {
    cout << "Could not get two valid arguments." << endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  commandLine.replace( 0, loc+1, ""); 
  std::string arg2 = commandLine;
  cout << "Mapping \"" << arg1 << "\" to \"" << arg2 << "\"" << endl;
  userMaps[ arg1 ] = arg2;
  return false;
}

//-----------------------------------------------------------------------------

void CMD_map::writeMaps( std::ofstream & os )
{

  std::map< std::string, std::string >::const_iterator iter = userMaps.begin() ,
                                                   iter_end = userMaps.end()   ;
  for( ; iter_end != iter; ++iter )
  {
    std::string cLine = "#map " + (*iter).first + " " + (*iter).second;
    os << cLine << endl;
  }

  return;
}

//-----------------------------------------------------------------------------

bool CMD_showStack::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    cout << "ERROR: No valid driver registered with showStack." << endl;
    return false;
  }

  // Show the command stack
  cout << "Command Stack :\n" << endl;

  std::list< std::string >::const_iterator iter = driver->getCommandStack().begin()  ,
                                       iter_end = driver->getCommandStack().end()    ;
  for( ; iter_end != iter; ++iter )
    cout << (*iter) << endl;

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_clearStack::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    cout << "ERROR: No valid driver registered with clearStack." << endl;
    return false;
  }

  // Show the command stack
  cout << "Command Stack Cleared." << endl;

  driver->getCommandStack().clear();

  return false;
}

//-----------------------------------------------------------------------------

bool CMD_isF::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isF() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_isJacobian::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isJacobian() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_isGradient::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isGradient() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_isNewton::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isNewton() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_isNormNewtonSolveResidual::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isNormNewtonSolveResidual() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_isPreconditioner::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isPreconditioner() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_isConditionNumber::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isConditionNumber() ? "True" : "False" );
  cout << " --> " << isValid << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_showValid::doCommand( std::string commandLine )
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

//-----------------------------------------------------------------------------

bool CMD_getJacobianConditionNumber::doCommand( std::string commandLine )
{
  cout << groupPtr->getJacobianConditionNumber() << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getNormF::doCommand( std::string commandLine )
{
  cout << groupPtr->getNormF() << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_setX::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find("setX"), 5, ""); 
  engine.GetMultiVector( commandLine.c_str(), *solnPtr );
  groupPtr->setX(*solnPtr);
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_computeF::doCommand( std::string commandLine )
{
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeF();
  cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_computeJacobian::doCommand( std::string commandLine )
{
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeJacobian();
  cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_computeGradient::doCommand( std::string commandLine )
{
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeGradient();
  cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_computeNewton::doCommand( std::string commandLine )
{
  const Teuchos::ParameterList & const_lsParams = solver.getList().
                                                       sublist("Direction").
                                                       sublist("Newton").
                                                       sublist("Linear Solver");
  Teuchos::ParameterList & lsParams = const_cast<Teuchos::ParameterList &>(const_lsParams);
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeNewton(lsParams);
  cout << "Return Status = " << returnMsg[ returnStatus ] << endl;
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getX::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getX()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "X" );
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getF::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getF()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "F" );
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getGradient::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getGradient()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "Gradient" );
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getNewton::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getNewton()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "Newton" );
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getJacobian::doCommand( std::string commandLine )
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

//-----------------------------------------------------------------------------

bool CMD_read::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    cout << "ERROR: No valid driver registered with reader." << endl;
    return false;
  }

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  ifstream inFile(commandLine.c_str());
  if( !inFile )
  {
    cout << "ERROR: Could not open file \"" << commandLine << "\"" << endl;
    return false;
  }

  CMD_newMacro::is = &inFile;

  char cLine[256];

  while( !inFile.eof() )
  {
    inFile.getline( cLine, 256 );
    std::string sLine(cLine);
    driver->doCommand( sLine );
  }

  CMD_newMacro::is = &(std::cin);

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_write::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    cout << "ERROR: No valid driver registered with reader." << endl;
    return false;
  }

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  ofstream outFile(commandLine.c_str());
  if( !outFile )
  {
    cout << "ERROR: Could not open file \"" << commandLine << "\"" << endl;
    return false;
  }

  // First, we write out all maps
  outFile << "######################" << endl;
  outFile << "####### Maps #########" << endl;
  outFile << "######################" << endl;
  outFile << "##" << endl;
  driver->getMapHandler()->writeMaps( outFile );
  outFile << "##" << endl;
  outFile << "##" << endl;

  // Then we write out all macros
  outFile << "######################" << endl;
  outFile << "###### Macros ########" << endl;
  outFile << "######################" << endl;
  outFile << "##" << endl;
  std::map< std::string, CMD_macro * >::const_iterator miter = driver->getUserMacros().begin() ,
                                                   miter_end = driver->getUserMacros().end()   ;
  for( ; miter_end != miter; ++miter )
    (*miter).second->writeMacro( outFile );

  outFile << "##" << endl;
  outFile << "##" << endl;

  // Now dump the command stack
  outFile << "######################" << endl;
  outFile << "### Command Stack ####" << endl;
  outFile << "######################" << endl;
  outFile << "##" << endl;
  std::list< std::string >::const_iterator siter = driver->getCommandStack().begin()  ,
                                       siter_end = driver->getCommandStack().end()    ;
  for( ; siter_end != siter; ++siter )
    outFile << (*siter) << endl;

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_newMacro::doCommand( std::string commandLine )
{
  if( 0 == driver )
  {
    cout << "ERROR: No valid driver registered with newMacro." << endl;
    return false;
  }

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  std::map< std::string, CMD_macro * > & userMacros = driver->getUserMacros();

  std::map< std::string, CMD_macro * >::iterator iter = userMacros.find(commandLine);
  if( userMacros.end() != iter )
  {
    cout << "Replacing macro \"" << commandLine << "\"" << endl;
    delete (*iter).second;
    driver->removeCommand( (*iter).second );
  }
  else
    cout << "Creating macro \"" << commandLine << "\"" << endl;

  CMD_macro * p_macro = new CMD_macro( engine, solver, commandLine );
  CMD_macro & macro   = *p_macro;
  macro.registerDriver( driver );

  std::string sLine;

  do
  {
    getline( *is, sLine );
    macro.addLineCommand( sLine );
  } 
  while( sLine != "##" );

  userMacros[commandLine] = p_macro;
  driver->addCommand( p_macro );

  return false;
}

//-----------------------------------------------------------------------------

bool CMD_macro::doCommand( std::string commandLine )
{
  if( 0 == driver )
  {
    cout << "ERROR: No valid driver registered with macro." << endl;
    return false;
  }

  // I may be the wrong macro.  Nevertheless, I get and invoke the correct one.

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  CMD_macro * p_macro = driver->getUserMacros()[commandLine];

  if( !p_macro )
  {
    cout << "Macro \"" << commandLine << "\" not found." << endl;
    return false;
  }

  cout << "Doing macro \"" << commandLine << "\"" << endl;
  
  std::list< std::string >::iterator iter = p_macro->macroCommands.begin() ,
                                 iter_end = p_macro->macroCommands.end()   ;

  for( ; iter_end != iter; ++iter )
    driver->doCommand( *iter );

  return true;
}

//-----------------------------------------------------------------------------

void CMD_macro::writeMacro( std::ofstream & os )
{

  std::string cLine = "#newMacro " + macroName;
  os << cLine << endl;

  std::list< std::string >::const_iterator iter = macroCommands.begin() ,
                                       iter_end = macroCommands.end()   ;
  for( ; iter_end != iter; ++iter )
  {
    os << (*iter) << endl;
  }

  return;
}

//-----------------------------------------------------------------------------

void CMD_macro::describe()
{

  cout << "\n\tMacro \"" << macroName << "\" :" << endl;

  std::list< std::string >::const_iterator iter = macroCommands.begin() ,
                                       iter_end = macroCommands.end()   ;
  for( ; iter_end != iter; ++iter )
    cout << "\t\t" << (*iter) << endl;

  return;
}

//-----------------------------------------------------------------------------
//----------------   Matlab Coupling Interface Commands   ---------------------
//-----------------------------------------------------------------------------

bool CMD_problemSummary::doCommand( std::string commandLine )
{
  problemManager->outputStatus(std::cout);
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_compJac::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), ""); 
  std::string::size_type loc = commandLine.find(" ");
  if( std::string::npos == loc )
  {
    cout << "Could not get two valid arguments." << endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  commandLine.replace( 0, loc+1, ""); 
  loc = commandLine.find("]");
  std::string arg2 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ,
      depID =  atoi( arg2.c_str()) ;

  cout << "Computing Jacobian Block " << probId << "," << depID << endl;

  problemManager->computeBlockJacobian( probId, depID );

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_setXvec::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), ""); 
  std::string::size_type loc = commandLine.find("]");
  if( std::string::npos == loc )
  {
    cout << "Could not get a valid argument." << endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ;

  commandLine.replace( 0, loc+2, ""); 

  cout << "Set solution vector for Problem " << probId << " using \"" 
       << commandLine << "\"" << endl;

  // Note we set both vectors, the one one in the problem, and the one in 
  // the corresponding group

  Epetra_Vector * p_soln = problemManager->getProblem(probId).getSolution().get();

  engine.GetMultiVector( commandLine.c_str(), *p_soln );
  problemManager->getGroup(probId).setX(*p_soln);

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_compRes::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), ""); 
  std::string::size_type loc = commandLine.find("]");
  if( std::string::npos == loc )
  {
    cout << "Could not get a valid argument." << endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ;

  cout << "Computing Residual for Problem " << probId << endl;

  problemManager->computeGroupF( probId );

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_doXfers::doCommand( std::string command )
{
  problemManager->syncAllProblems();
  problemManager->setAllGroupX();
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getAllX::doCommand( std::string command )
{
  map<int, GenericEpetraProblem*>::iterator problemIter = problemManager->getProblems().begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = problemManager->getProblems().end();

  // Do diagonal blocks
  for( ; problemLast != problemIter; ++problemIter )
  {
    GenericEpetraProblem & problem = *(*problemIter).second;
    int                    probId  = (*problemIter).first;

    const Epetra_Vector * tempVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                       (problemManager->getGroup(probId).getX()).getEpetraVector());

    if( tempVec )
    {
      ostringstream sval1;
      sval1 << probId << flush;
      std::string name = "X_" + sval1.str();
      engine.PutMultiVector( *tempVec, name.c_str() );
      cout << "Stored Solution (" << probId << ") in \"" << name << "\"" << endl;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getJac::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), ""); 
  std::string::size_type loc = commandLine.find(" ");
  if( std::string::npos == loc )
  {
    cout << "Could not get two valid arguments." << endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  commandLine.replace( 0, loc+1, ""); 
  loc = commandLine.find("]");
  std::string arg2 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ,
      depID =  atoi( arg2.c_str()) ;

  Epetra_RowMatrix * rowMatrix = problemManager->getBlockJacobianMatrix( probId, depID );
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

//-----------------------------------------------------------------------------

bool CMD_getAllJac::doCommand( std::string commandLine )
{
  map<int, GenericEpetraProblem*>::iterator problemIter = problemManager->getProblems().begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = problemManager->getProblems().end();

  // Do diagonal blocks
  for( ; problemLast != problemIter; ++problemIter )
  {

    GenericEpetraProblem & problem = *(*problemIter).second;
    int                    probId  = (*problemIter).first;

    Epetra_RowMatrix * rowMatrix = problemManager->getBlockJacobianMatrix( probId );

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
    if( problemManager->useOffBlocks() ) 
    {
#ifdef HAVE_NOX_EPETRAEXT
      for( unsigned int k = 0; k < problem.getDependentProblems().size(); ++k) 
      {
        int depId = problem.getDependentProblems()[k];

        Epetra_RowMatrix * rowMatrix = problemManager->getBlockJacobianMatrix( probId, depId );

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

//-----------------------------------------------------------------------------

bool CMD_getRes::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), ""); 
  std::string::size_type loc = commandLine.find("]");
  if( std::string::npos == loc )
  {
    cout << "Could not get a valid argument." << endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ;

  const Epetra_Vector * resVec = problemManager->getResidual( probId );

  ostringstream sval1;
  sval1 << probId << flush;
  std::string name = "F_" + sval1.str();
  engine.PutMultiVector( *resVec, name.c_str() );
  cout << "Stored Residual (" << probId << ") in \"" << name << "\"" << endl;

  return true;
}

//-----------------------------------------------------------------------------

bool CMD_getAllRes::doCommand( std::string command )
{
  map<int, GenericEpetraProblem*>::iterator problemIter = problemManager->getProblems().begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = problemManager->getProblems().end();

  // Do diagonal blocks
  for( ; problemLast != problemIter; ++problemIter )
  {
    GenericEpetraProblem & problem = *(*problemIter).second;
    int                    probId  = (*problemIter).first;

    const Epetra_Vector * resVec = problemManager->getResidual( probId );

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
