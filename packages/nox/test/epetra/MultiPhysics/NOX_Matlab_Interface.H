// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_MATLAB_INTERFACE_H
#define NOX_MATLAB_INTERFACE_H

#include "NOX_Epetra.H"

#ifdef HAVE_MATLAB

#include "EpetraExt_MatlabEngine.h"

//-----------------------------------------------------------------------------

#include "EpetraExt_MatlabEngine.h"

// the following deal with matlab provided headers:
#include "engine.h"
#include "mex.h"
#undef printf  // matlab has its own printf that we don't want to use

#define BUFSIZE 200
#define MATLABBUF 1024 * 16

//-----------------------------------------------------------------------------

// This is a new class that may evantually get moved into NOX.  For now,
// this is simply used as a testbed for driving NOX using Matlab

class Matlab_Interface
{
public:

  // A base class for Matlab interface commands to NOX
  class CommandBase
  {

   public:

    enum CommandCategory {
                                GROUP_GET       ,
                                GROUP_SET       ,
                                GROUP_COMPUTE   ,
                                CALLBACK_FUN    ,
                                QUERY           ,
                                INOUT           ,
                                GENERIC
                         };

    static map< NOX::Abstract::Group::ReturnType, std::string > returnMsg;

    static map< CommandCategory, std::string > categoryDescription;

    CommandBase(  EpetraExt::EpetraExt_MatlabEngine & engine_,
                  NOX::Solver::Manager & solver_               );

    virtual ~CommandBase() { };

    // methods

    virtual bool doCommand( std::string ) = 0 ;

    virtual const std::string & getCommand()
    { return command; }

    virtual const std::string & getDescription()
    { return description; }

    virtual void describe()
    {
      cout << "\t\t#" << command << "\t" << description << std::endl;
      return;
    }

    virtual bool isCategory( CommandCategory category )
    {
      return ( myCategory == category );
    }

    // Get the NOX solver
    const NOX::Solver::Manager * const getSolver() const
    { return &solver; }

    virtual void registerDriver( Matlab_Interface * driver_ )
    { driver = driver_; }

   protected:

    // methods

    virtual void initialize();

    // data

    std::string command     ;
    std::string description ;

    std::map< std::string, std::string > userMaps;

    EpetraExt::EpetraExt_MatlabEngine & engine;

    NOX::Solver::Manager & solver;

    Matlab_Interface * driver;

    CommandCategory myCategory;

    // Convenient pointers to underlying data in solver

    NOX::Epetra::Group * groupPtr;

    Epetra_Vector * solnPtr;
  };

  // Derived classes that represent individual commands

  class CMD_map : public CommandBase
  {

   public:

    CMD_map( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      command     = "map"              ;
      description = " - Map a command" ;
    }

    virtual bool doCommand( std::string );

    virtual void writeMaps( std::ofstream & );

    std::map< std::string, std::string > & getUserMaps()
    { return userMaps; }

  };

  //-----------------------------------------------------------------------------
  class CMD_showStack : public CommandBase
  {

   public:

    CMD_showStack( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "showStack"              ;
      description = " - Show the command stack" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_clearStack : public CommandBase
  {

   public:

    CMD_clearStack( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "clearStack"              ;
      description = " - Clear the command stack" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isF : public CommandBase
  {

   public:

    CMD_isF( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isF"              ;
      description = " - Query if Residual is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isJacobian : public CommandBase
  {

   public:

    CMD_isJacobian( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isJacobian"              ;
      description = " - Query if Jacobian is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isGradient : public CommandBase
  {

   public:

    CMD_isGradient( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isGradient"              ;
      description = " - Query if Gradient is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isNewton : public CommandBase
  {

   public:

    CMD_isNewton( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isNewton"              ;
      description = " - Query if Newton vector is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isNormNewtonSolveResidual : public CommandBase
  {

   public:

    CMD_isNormNewtonSolveResidual( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isNormNewtonSolveResidual"              ;
      description = " - Query if norm of the Newton solve is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isPreconditioner : public CommandBase
  {

   public:

    CMD_isPreconditioner( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isPreconditioner"              ;
      description = " - Query if Preconditioner is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_isConditionNumber : public CommandBase
  {

   public:

    CMD_isConditionNumber( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "isConditionNumber"              ;
      description = " - Query if Condition Number is current wrt X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_showValid : public CommandBase
  {

   public:

    CMD_showValid( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "showValid"              ;
      description = " - Show status of all \"is***\" query methods" ;
    }

    virtual bool doCommand( std::string );

    // Static member to allow use by others
    static bool showValid( const NOX::Epetra::Group * );
  };

  //-----------------------------------------------------------------------------
  class CMD_getJacobianConditionNumber : public CommandBase
  {

   public:

    CMD_getJacobianConditionNumber( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getJacobianConditionNumber"              ;
      description = " - Report the condition number of the Jacobian matrix" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_getNormF : public CommandBase
  {

   public:

    CMD_getNormF( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getNormF"              ;
      description = " - Report the Residual (2-)norm" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_setX : public CommandBase
  {

   public:

    CMD_setX( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "setX"              ;
      description = " - Set the group solution, X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_computeF : public CommandBase
  {

   public:

    CMD_computeF( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "computeF"              ;
      description = " - Compute the Residual vector" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_computeJacobian : public CommandBase
  {

   public:

    CMD_computeJacobian( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "computeJacobian"              ;
      description = " - Compute the Jacobian" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_computeGradient : public CommandBase
  {

   public:

    CMD_computeGradient( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "computeGradient"              ;
      description = " - Compute the Gradient vector" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_computeNewton : public CommandBase
  {

   public:

    CMD_computeNewton( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "computeNewton"              ;
      description = " - Compute the Newton vector" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_getX : public CommandBase
  {

   public:

    CMD_getX( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getX"              ;
      description = " - Get the solution vector, X" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_getF : public CommandBase
  {

   public:

    CMD_getF( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getF"              ;
      description = " - Get the Residual vector, F" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_getGradient : public CommandBase
  {

   public:

    CMD_getGradient( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getGradient"              ;
      description = " - Get the Gradient vector" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_getNewton : public CommandBase
  {

   public:

    CMD_getNewton( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getNewton"              ;
      description = " - Get the Newton vector" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_getJacobian : public CommandBase
  {

   public:

    CMD_getJacobian( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getJacobian"              ;
      description = " - Get the Jacobian matrix" ;
    }

    virtual bool doCommand( std::string );
  };

  //-----------------------------------------------------------------------------
  class CMD_newMacro : public CommandBase
  {

   public:

    static std::istream * is;

    CMD_newMacro( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      command     = "newMacro"              ;
      description = " - Create a named macro" ;
    }

    virtual bool doCommand( std::string );

  };

  //-----------------------------------------------------------------------------
  class CMD_macro : public CommandBase
  {

   public:

    CMD_macro( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_  ,
             std::string                         name_     ) :
      CommandBase( engine_, solver_ ),
      macroName( name_ )
    {
      command     = "macro"              ;
      description = " - invoke a user macro" ;
    }

    virtual bool doCommand( std::string );

    virtual bool addLineCommand( std::string & com )
    { macroCommands.push_back(com); return true; }

    virtual void describe();

    virtual void writeMacro( std::ofstream & );

    virtual std::string getName()
    { return macroName; }

   protected:

    std::string macroName;

    std::list< std::string > macroCommands;

  };

  //-----------------------------------------------------------------------------
  class CMD_read : public CommandBase
  {

   public:

    CMD_read( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = INOUT;
      command     = "read"              ;
      description = " - read a script file" ;
    }

    virtual bool doCommand( std::string );

  };

  //-----------------------------------------------------------------------------
  class CMD_write : public CommandBase
  {

   public:

    CMD_write( EpetraExt::EpetraExt_MatlabEngine & engine_,
             NOX::Solver::Manager              & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = INOUT;
      command     = "write"              ;
      description = " - write a script file" ;
    }

    virtual bool doCommand( std::string );

  };

  //-----------------------------------------------------------------------------

  Matlab_Interface(NOX::Solver::Manager & solver);

  virtual ~Matlab_Interface();

  // Begin an interactive session
  virtual void interact();

  virtual bool doCommand( std::string );

  virtual bool doNOXCommand( std::string & );

  virtual void addCommand( CommandBase * com )
  { commands.push_back( com ); return; }

  virtual void removeCommand( CommandBase * com )
  { commands.remove( com ); return; }

  CMD_map * getMapHandler()
  { return mapHandler; }

  std::list< std::string > & getCommandStack()
  { return commandStack; }

  std::map< std::string, CMD_macro * > & getUserMacros()
  { return userMacros; }


protected:

  // Methods

  // Data

  // A needed communicator - can be replaced later by initialization argument
  Epetra_SerialComm comm;

  EpetraExt::EpetraExt_MatlabEngine * enginePtr;
  EpetraExt::EpetraExt_MatlabEngine & engine;

  NOX::Solver::Manager & solver;
  NOX::Epetra::Group * groupPtr;
  Epetra_Vector * solnPtr;

  // A user command mapping handler
  CMD_map * mapHandler;

  std::map< std::string, CMD_macro * > userMacros;

  std::list< CommandBase* > commands;

  std::list< std::string > commandStack;
};

#endif
#endif
