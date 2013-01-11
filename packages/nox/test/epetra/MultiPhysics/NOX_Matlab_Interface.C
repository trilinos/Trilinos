//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
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
                                                                                
#include "NOX_Matlab_Interface.H"

#ifdef HAVE_MATLAB

#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_RowMatrixOut.h"

// This is a new class that may evantually get moved into NOX.  For now,
// this is simply used as a testbed for driving NOX using Matlab

// Declare and define static data member
// Declare and define static data member
map< NOX::Abstract::Group::ReturnType, std::string > Matlab_Interface::CommandBase::returnMsg;
map< Matlab_Interface::CommandBase::CommandCategory, std::string > Matlab_Interface::CommandBase::categoryDescription;
std::istream * Matlab_Interface::CMD_newMacro::is = &(std::cin);

Matlab_Interface::Matlab_Interface(NOX::Solver::Manager &  noxSolver) :
  enginePtr( new EpetraExt::EpetraExt_MatlabEngine(comm) ),
  engine(*enginePtr),
  solver(noxSolver)
{
  std::cout << "matlab started\n";

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

void 
Matlab_Interface::interact()
{
  char s [BUFSIZE] ;

  bool status;

  while(1) 
  {
    // Prompt the user and get a std::string
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

bool 
Matlab_Interface::doCommand( std::string command )
{
  char matlabBuffer [MATLABBUF];

  bool status;

  // Parse for NOX commands
  if( command[0] == '#' )
  {
    command.erase(0,1);
    // Remove any carriage returns
    if( command.find('\n') != std::string::npos )
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

bool 
Matlab_Interface::doNOXCommand( std::string & command )
{

  std::cout << "NOX     : " << command << std::endl;

  // Allow for "#" no-op
  if( '#' == command[0] )
    return true;

  // Give priority to mapped commands

  if( mapHandler->getUserMaps().end() != mapHandler->getUserMaps().find(command) )
    command = mapHandler->getUserMaps()[command];

  try 
  {
    // Convenience methods

    if( command.find("?") != std::string::npos )
    {
      std::cout << "\n\tCommand Summary:\n" << std::endl;

      std::map< CommandBase::CommandCategory, std::string >::iterator cIter = CommandBase::categoryDescription.begin() ,
                                                                 cIter_end = CommandBase::categoryDescription.end()     ;

      for(  ; cIter_end != cIter; ++cIter )
      {
        std::cout << "\n\t" << (*cIter).second << std::endl << std::endl;

        std::list<CommandBase *>::const_iterator iter = commands.begin() ,
                                             iter_end = commands.end()   ;
        for(  ; iter_end != iter; ++iter )
          if( (*iter)->isCategory((*cIter).first) )
            (*iter)->describe();
      }

      if( !mapHandler->getUserMaps().empty() )
      {
        std::cout << "\n\tUser Defined (Mapped) Command Summary:\n" << std::endl;

        for( std::map<string, std::string>::const_iterator iter = 
            mapHandler->getUserMaps().begin(); iter != mapHandler->getUserMaps().end(); ++iter )
          std::cout << "\t\t" << (*iter).first << " --> " << (*iter).second << std::endl;
      }

      return true;
    }

    std::list<CommandBase *>::const_iterator iter = commands.begin() ,
                                         iter_end = commands.end()   ;
    for(  ; iter_end != iter; ++iter )
      if( std::string::npos != command.find((*iter)->getCommand()) )
        return (*iter)->doCommand(command);

    std::cout << "\n" << command << " : Not found." << std::endl;
    return false;

  }
  catch( const char * msg )
  {
    std::cout << msg << std::endl;
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
//----------------------  Commands for NOX solver -----------------------------
//-----------------------------------------------------------------------------

Matlab_Interface::CommandBase::CommandBase( EpetraExt::EpetraExt_MatlabEngine & engine_, 
                          NOX::Solver::Manager & solver_ ) :
  engine(engine_),
  solver(solver_),
  driver(0)
{
  initialize();
}

//-----------------------------------------------------------------------------

void 
Matlab_Interface::CommandBase::initialize()
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

  categoryDescription[ GROUP_GET      ]             = "Group 'Get' Methods"      ;
  categoryDescription[ GROUP_SET      ]             = "Group 'Set' Methods"      ;
  categoryDescription[ GROUP_COMPUTE  ]             = "Group 'Compute' Methods"  ;
  categoryDescription[ CALLBACK_FUN   ]             = "Callback Functionality"   ;
  categoryDescription[ QUERY          ]             = "Query Methods"            ;
  categoryDescription[ INOUT          ]             = "Input/Output Methods"     ;
  categoryDescription[ GENERIC        ]             = "Generic Methods"          ;

}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_map::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find("map "), 4, ""); 
  std::string::size_type loc = commandLine.find(" ");
  if( std::string::npos == loc )
  {
    std::cout << "Could not get two valid arguments." << std::endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  commandLine.replace( 0, loc+1, ""); 
  std::string arg2 = commandLine;
  std::cout << "Mapping \"" << arg1 << "\" to \"" << arg2 << "\"" << std::endl;
  userMaps[ arg1 ] = arg2;
  return false;
}

//-----------------------------------------------------------------------------

void 
Matlab_Interface::CMD_map::writeMaps( std::ofstream & os )
{

  std::map< std::string, std::string >::const_iterator iter = userMaps.begin() ,
                                                   iter_end = userMaps.end()   ;
  for( ; iter_end != iter; ++iter )
  {
    std::string cLine = "#map " + (*iter).first + " " + (*iter).second;
    os << cLine << std::endl;
  }

  return;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_showStack::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    std::cout << "ERROR: No valid driver registered with showStack." << std::endl;
    return false;
  }

  // Show the command stack
  std::cout << "Command Stack :\n" << std::endl;

  std::list< std::string >::const_iterator iter = driver->getCommandStack().begin()  ,
                                       iter_end = driver->getCommandStack().end()    ;
  for( ; iter_end != iter; ++iter )
    std::cout << (*iter) << std::endl;

  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_clearStack::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    std::cout << "ERROR: No valid driver registered with clearStack." << std::endl;
    return false;
  }

  // Show the command stack
  std::cout << "Command Stack Cleared." << std::endl;

  driver->getCommandStack().clear();

  return false;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isF::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isF() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isJacobian::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isJacobian() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isGradient::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isGradient() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isNewton::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isNewton() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isNormNewtonSolveResidual::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isNormNewtonSolveResidual() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isPreconditioner::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isPreconditioner() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_isConditionNumber::doCommand( std::string commandLine )
{
  std::string isValid = (groupPtr->isConditionNumber() ? "True" : "False" );
  std::cout << " --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_showValid::doCommand( std::string commandLine )
{
  return showValid( groupPtr );
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_showValid::showValid( const NOX::Epetra::Group * p_Grp )
{
  std::string isValid = (p_Grp->isF() ? "True" : "False" );
  std::cout << " isF              --> " << isValid << std::endl;
  isValid = (p_Grp->isJacobian() ? "True" : "False" );
  std::cout << " isJacobian       --> " << isValid << std::endl;
  isValid = (p_Grp->isGradient() ? "True" : "False" );
  std::cout << " isGradient       --> " << isValid << std::endl;
  isValid = (p_Grp->isNewton() ? "True" : "False" );
  std::cout << " isNewton         --> " << isValid << std::endl;
  isValid = (p_Grp->isPreconditioner() ? "True" : "False" );
  std::cout << " isPreconditioner --> " << isValid << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getJacobianConditionNumber::doCommand( std::string commandLine )
{
  std::cout << groupPtr->getJacobianConditionNumber() << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getNormF::doCommand( std::string commandLine )
{
  std::cout << groupPtr->getNormF() << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_setX::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find("setX"), 5, ""); 
  engine.GetMultiVector( commandLine.c_str(), *solnPtr );
  groupPtr->setX(*solnPtr);
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_computeF::doCommand( std::string commandLine )
{
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeF();
  std::cout << "Return Status = " << returnMsg[ returnStatus ] << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_computeJacobian::doCommand( std::string commandLine )
{
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeJacobian();
  std::cout << "Return Status = " << returnMsg[ returnStatus ] << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_computeGradient::doCommand( std::string commandLine )
{
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeGradient();
  std::cout << "Return Status = " << returnMsg[ returnStatus ] << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_computeNewton::doCommand( std::string commandLine )
{
  const Teuchos::ParameterList & const_lsParams = solver.getList().
                                                       sublist("Direction").
                                                       sublist("Newton").
                                                       sublist("Linear Solver");
  Teuchos::ParameterList & lsParams = const_cast<Teuchos::ParameterList &>(const_lsParams);
  NOX::Abstract::Group::ReturnType returnStatus = groupPtr->computeNewton(lsParams);
  std::cout << "Return Status = " << returnMsg[ returnStatus ] << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getX::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getX()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "X" );
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getF::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getF()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "F" );
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getGradient::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getGradient()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "Gradient" );
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getNewton::doCommand( std::string commandLine )
{
  const Epetra_Vector * tmpVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getNewton()).getEpetraVector());
  engine.PutMultiVector( *tmpVec, "Newton" );
  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_getJacobian::doCommand( std::string commandLine )
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

  std::cout << "Could not get a valid matrix." << std::endl;
  return false;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_read::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    std::cout << "ERROR: No valid driver registered with reader." << std::endl;
    return false;
  }

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  ifstream inFile(commandLine.c_str());
  if( !inFile )
  {
    std::cout << "ERROR: Could not open file \"" << commandLine << "\"" << std::endl;
    return false;
  }

  Matlab_Interface::CMD_newMacro::is = &inFile;

  char cLine[256];

  while( !inFile.eof() )
  {
    inFile.getline( cLine, 256 );
    std::string sLine(cLine);
    driver->doCommand( sLine );
  }

  Matlab_Interface::CMD_newMacro::is = &(std::cin);

  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_write::doCommand( std::string commandLine )
{

  if( 0 == driver )
  {
    std::cout << "ERROR: No valid driver registered with reader." << std::endl;
    return false;
  }

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  std::ofstream outFile(commandLine.c_str());
  if( !outFile )
  {
    std::cout << "ERROR: Could not open file \"" << commandLine << "\"" << std::endl;
    return false;
  }

  // First, we write out all maps
  outFile << "######################" << std::endl;
  outFile << "####### Maps #########" << std::endl;
  outFile << "######################" << std::endl;
  outFile << "##" << std::endl;
  driver->getMapHandler()->writeMaps( outFile );
  outFile << "##" << std::endl;
  outFile << "##" << std::endl;

  // Then we write out all macros
  outFile << "######################" << std::endl;
  outFile << "###### Macros ########" << std::endl;
  outFile << "######################" << std::endl;
  outFile << "##" << std::endl;
  std::map< std::string, Matlab_Interface::CMD_macro * >::const_iterator miter = driver->getUserMacros().begin() ,
                                                   miter_end = driver->getUserMacros().end()   ;
  for( ; miter_end != miter; ++miter )
    (*miter).second->writeMacro( outFile );

  outFile << "##" << std::endl;
  outFile << "##" << std::endl;

  // Now dump the command stack
  outFile << "######################" << std::endl;
  outFile << "### Command Stack ####" << std::endl;
  outFile << "######################" << std::endl;
  outFile << "##" << std::endl;
  std::list< std::string >::const_iterator siter = driver->getCommandStack().begin()  ,
                                       siter_end = driver->getCommandStack().end()    ;
  for( ; siter_end != siter; ++siter )
    outFile << (*siter) << std::endl;

  return true;
}

//-----------------------------------------------------------------------------

bool 
Matlab_Interface::CMD_newMacro::doCommand( std::string commandLine )
{
  if( 0 == driver )
  {
    std::cout << "ERROR: No valid driver registered with newMacro." << std::endl;
    return false;
  }

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  std::map< std::string, Matlab_Interface::CMD_macro * > & userMacros = driver->getUserMacros();

  std::map< std::string, Matlab_Interface::CMD_macro * >::iterator iter = userMacros.find(commandLine);
  if( userMacros.end() != iter )
  {
    std::cout << "Replacing macro \"" << commandLine << "\"" << std::endl;
    delete (*iter).second;
    driver->removeCommand( (*iter).second );
  }
  else
    std::cout << "Creating macro \"" << commandLine << "\"" << std::endl;

  Matlab_Interface::CMD_macro * p_macro = new Matlab_Interface::CMD_macro( engine, solver, commandLine );
  Matlab_Interface::CMD_macro & macro   = *p_macro;
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

bool 
Matlab_Interface::CMD_macro::doCommand( std::string commandLine )
{
  if( 0 == driver )
  {
    std::cout << "ERROR: No valid driver registered with macro." << std::endl;
    return false;
  }

  // I may be the wrong macro.  Nevertheless, I get and invoke the correct one.

  commandLine.replace( commandLine.find(command), command.size()+1, ""); 

  Matlab_Interface::CMD_macro * p_macro = driver->getUserMacros()[commandLine];

  if( !p_macro )
  {
    std::cout << "Macro \"" << commandLine << "\" not found." << std::endl;
    return false;
  }

  std::cout << "Doing macro \"" << commandLine << "\"" << std::endl;
  
  std::list< std::string >::iterator iter = p_macro->macroCommands.begin() ,
                                 iter_end = p_macro->macroCommands.end()   ;

  for( ; iter_end != iter; ++iter )
    driver->doCommand( *iter );

  return true;
}

//-----------------------------------------------------------------------------

void 
Matlab_Interface::CMD_macro::writeMacro( std::ofstream & os )
{

  std::string cLine = "#newMacro " + macroName;
  os << cLine << std::endl;

  std::list< std::string >::const_iterator iter = macroCommands.begin() ,
                                       iter_end = macroCommands.end()   ;
  for( ; iter_end != iter; ++iter )
  {
    os << (*iter) << std::endl;
  }

  return;
}

//-----------------------------------------------------------------------------

void 
Matlab_Interface::CMD_macro::describe()
{

  std::cout << "\n\tMacro \"" << macroName << "\" :" << std::endl;

  std::list< std::string >::const_iterator iter = macroCommands.begin() ,
                                       iter_end = macroCommands.end()   ;
  for( ; iter_end != iter; ++iter )
    std::cout << "\t\t" << (*iter) << std::endl;

  return;
}

#endif
