// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_MULTIPHYSICS_MATLAB_INTERFACE_H
#define NOX_MULTIPHYSICS_MATLAB_INTERFACE_H

#include "NOX_Matlab_Interface.H" //base class

#ifdef HAVE_MATLAB

//-----------------------------------------------------------------------------

// Forward declarations
class Problem_Manager;

// This is a new class that may evantually get moved into NOX.  For now,
// this is simply used as a testbed for driving NOX using Matlab

class Coupling_Matlab_Interface : public Matlab_Interface
{
public:

  // A base class for Matlab interface commands to NOX
  class CommandBase : public Matlab_Interface::CommandBase
  {

   public:

    CommandBase(  EpetraExt::EpetraExt_MatlabEngine & engine_,
                  Problem_Manager &                         );

    // methods

    virtual bool doCommand( std::string ) = 0 ;


    // Get the Coupling solver
    const Problem_Manager & getCoupler()
    { return problemManager; }

   protected:

    // methods

    virtual ~CommandBase() { };

    // data

    Problem_Manager      & problemManager;

  };

  // Derived classes that represent individual commands

  class CMD_problemSummary : public CommandBase
  {

    public:

    CMD_problemSummary( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = QUERY;
      command     = "problemSummary"              ;
      description = " - (Coupling) Summarize coupled problems" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_showAllValid : public CommandBase
  {

   public:

    CMD_showAllValid( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager                 & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = QUERY;
      command     = "showAllValid"              ;
      description = " - (Coupling) Show valid status of all groups" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_compPreconditioner : public CommandBase
  {

   public:

    CMD_compPreconditioner( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager                 & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "compPreconditioner"              ;
      description = " - (Coupling) Compute the coupling Preconditioner matrix" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_compJac : public CommandBase
  {

   public:

    CMD_compJac( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "compJac["              ;
      description = " - (Coupling) Compute Block Jacobian (i, j)" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_setXvec : public CommandBase
  {

   public:

    CMD_setXvec( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_SET;
      command     = "setSoln["              ;
      description = " - (Coupling) Set a solution vector of problem i" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_compRes : public CommandBase
  {

   public:

    CMD_compRes( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_COMPUTE;
      command     = "compRes["              ;
      description = " - (Coupling) Compute Residual for problem i" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_syncAllGroupX : public CommandBase
  {

   public:

    CMD_syncAllGroupX( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager                 & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_SET;
      command     = "syncGroupsX"              ;
      description = " - (Coupling) Copy pieces of coupling Group solution into each problem's group solution." ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_doXfers : public CommandBase
  {

   public:

    CMD_doXfers( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = CALLBACK_FUN;
      command     = "doXfers"              ;
      description = " - (Coupling) Do all registered transfers of data" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_getAllX : public CommandBase
  {

   public:

    CMD_getAllX( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_GET;
      command     = "getAllX"              ;
      description = " - (Coupling) Get All Problem Solutions" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_getJac : public CommandBase
  {

   public:

    CMD_getJac( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_GET;
      command     = "getJac["              ;
      description = " - (Coupling) Get Block Jacobian [i j]" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_getAllJac : public CommandBase
  {

   public:

    CMD_getAllJac( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_GET;
      command     = "getAllJac"              ;
      description = " - (Coupling) Get All Block Jacobians" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_getRes : public CommandBase
  {

   public:

    CMD_getRes( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_GET;
      command     = "getRes["              ;
      description = " - (Coupling) Get Residual for Problem [i]" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_getAllRes : public CommandBase
  {

   public:

    CMD_getAllRes( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager              & problemManager_ ) :
      CommandBase( engine_, problemManager_ )
    {
      myCategory = GROUP_GET;
      command     = "getAllRes"              ;
      description = " - (Coupling) Get All Problem Residuals" ;
    }

    virtual bool doCommand( std::string );
  };

  class CMD_getPrecMatrix : public CommandBase
  {

   public:

    CMD_getPrecMatrix( EpetraExt::EpetraExt_MatlabEngine & engine_,
             Problem_Manager                 & solver_ ) :
      CommandBase( engine_, solver_ )
    {
      myCategory = GROUP_GET;
      command     = "getPrecMatrix"              ;
      description = " - (Coupling) Get the coupling Preconditioner matrix" ;
    }

    virtual bool doCommand( std::string );
  };

  Coupling_Matlab_Interface(Problem_Manager & manager);

  virtual ~Coupling_Matlab_Interface() { };

protected:

  // Methods

  // Data

  Problem_Manager & problemManager;
};

#endif
#endif // NOX_MULTIPHYSICS_MATLAB_INTERFACE_H
