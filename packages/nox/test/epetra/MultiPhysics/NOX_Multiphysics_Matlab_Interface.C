// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Multiphysics_Matlab_Interface.H"

#ifdef HAVE_MATLAB

#include "Problem_Manager.H"

// This is a new class that may evantually get moved into NOX.  For now,
// this is simply used as a testbed for driving NOX using Matlab

Coupling_Matlab_Interface::Coupling_Matlab_Interface(Problem_Manager &  manager_) :
  Matlab_Interface( *(manager_.getCompositeSolver()) ),
  problemManager(manager_)
{
  std::cout << "coupling matlab started\n";

  commands.push_back( new CMD_problemSummary             ( engine, problemManager ) );
  commands.push_back( new CMD_showAllValid               ( engine, problemManager ) );
  commands.push_back( new CMD_getAllX                    ( engine, problemManager ) );
  commands.push_back( new CMD_setXvec                    ( engine, problemManager ) );
  commands.push_back( new CMD_compJac                    ( engine, problemManager ) );
  commands.push_back( new CMD_compPreconditioner         ( engine, problemManager ) );
  commands.push_back( new CMD_compRes                    ( engine, problemManager ) );
  commands.push_back( new CMD_syncAllGroupX              ( engine, problemManager ) );
  commands.push_back( new CMD_doXfers                    ( engine, problemManager ) );
  commands.push_back( new CMD_getJac                     ( engine, problemManager ) );
  commands.push_back( new CMD_getRes                     ( engine, problemManager ) );
  commands.push_back( new CMD_getAllJac                  ( engine, problemManager ) );
  commands.push_back( new CMD_getAllRes                  ( engine, problemManager ) );
  commands.push_back( new CMD_getPrecMatrix              ( engine, problemManager ) );
}

//-----------------------------------------------------------------------------
//----------------------  Commands for Coupling solver ------------------------
//-----------------------------------------------------------------------------

Coupling_Matlab_Interface::CommandBase::CommandBase( EpetraExt::EpetraExt_MatlabEngine & engine_,
                          Problem_Manager & problemManager_) :
  Matlab_Interface::CommandBase( engine_, *problemManager_.getCompositeSolver() ),
  problemManager(problemManager_)
{
  initialize();
}

//-----------------------------------------------------------------------------
//----------------   Matlab Coupling Interface Commands   ---------------------
//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_problemSummary::doCommand( std::string commandLine )
{
  problemManager.outputStatus(std::cout);
  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_showAllValid::doCommand( std::string commandLine )
{

  std::cout << std::endl;
  std::cout << "\tCoupling Group status " << std::endl;
  std::cout << "\t--------------------- " << std::endl;
  Matlab_Interface::CMD_showValid::showValid( groupPtr );
  std::cout << std::endl;

  for( int probId = 1; probId <= problemManager.getProblemCount(); ++probId )
  {
    std::cout << "\tGroup status for problem \"" << problemManager.getNames()[probId] << "\"" << std::endl;
    std::cout << "\t------------------------ " << std::endl;

    const NOX::Epetra::Group * p_probGrp = &(problemManager.getSolutionGroup(probId));

    Matlab_Interface::CMD_showValid::showValid( p_probGrp );
    std::cout << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_compPreconditioner::doCommand( std::string commandLine )
{
    //Epetra_Operator * dummyOp = NULL;

    std::cout << "Command currently unspupported." << std::endl;
    return false;
    //cout << "Computing coupling matrix preconditioner." << std::endl;
    //return problemManager.computeJacobian( *solnPtr, *dummyOp );
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_compJac::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), "");
  std::string::size_type loc = commandLine.find(" ");
  if( std::string::npos == loc )
  {
    std::cout << "Could not get two valid arguments." << std::endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  commandLine.replace( 0, loc+1, "");
  loc = commandLine.find("]");
  std::string arg2 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ,
      depID =  atoi( arg2.c_str()) ;

  std::cout << "Computing Jacobian Block " << probId << "," << depID << std::endl;

  problemManager.computeBlockJacobian( probId, depID );

  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_setXvec::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), "");
  std::string::size_type loc = commandLine.find("]");
  if( std::string::npos == loc )
  {
    std::cout << "Could not get a valid argument." << std::endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ;

  commandLine.replace( 0, loc+2, "");

  std::cout << "Set solution vector for Problem " << probId << " using \""
       << commandLine << "\"" << std::endl;

  // Note we set both vectors, the one one in the problem, and the one in
  // the corresponding group

  Epetra_Vector * p_soln = problemManager.getProblem(probId).getSolution().get();

  engine.GetMultiVector( commandLine.c_str(), *p_soln );
  const_cast<NOX::Epetra::Group&>(problemManager.getSolutionGroup(probId)).setX(*p_soln);

  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_compRes::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), "");
  std::string::size_type loc = commandLine.find("]");
  if( std::string::npos == loc )
  {
    std::cout << "Could not get a valid argument." << std::endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ;

  std::cout << "Computing Residual for Problem " << probId << std::endl;

  problemManager.computeGroupF( probId );

  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_syncAllGroupX::doCommand( std::string command )
{

  const Epetra_Vector * compositeSoln = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (groupPtr->getX()).getEpetraVector());

  for( int probId = 1; probId <= problemManager.getProblemCount(); ++probId )
  {
    Epetra_Vector * tempSoln = new Epetra_Vector( problemManager.getSolutionVec(probId) );

    problemManager.copyCompositeToVector( *compositeSoln, probId, *tempSoln );

    NOX::Epetra::Group & grp = const_cast<NOX::Epetra::Group&>( problemManager.getSolutionGroup(probId) );
    grp.setX(*tempSoln);

    delete tempSoln; tempSoln = 0;

    std::cout << "Copied composite solution into problem group # " << probId << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_doXfers::doCommand( std::string command )
{
  problemManager.syncAllProblems();
  problemManager.setAllGroupX();
  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_getAllX::doCommand( std::string command )
{
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = problemManager.getProblems().begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = problemManager.getProblems().end();

  // Do diagonal blocks
  for( ; problemLast != problemIter; ++problemIter )
  {
    int probId  = (*problemIter).first;

    const Epetra_Vector * tempVec = &( problemManager.getSolutionVec(probId) );

    if( tempVec )
    {
      ostringstream sval1;
      sval1 << probId << flush;
      std::string name = "X_" + sval1.str();
      engine.PutMultiVector( *tempVec, name.c_str() );
      std::cout << "Stored Solution (" << probId << ") in \"" << name << "\"" << std::endl;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_getJac::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), "");
  std::string::size_type loc = commandLine.find(" ");
  if( std::string::npos == loc )
  {
    std::cout << "Could not get two valid arguments." << std::endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  commandLine.replace( 0, loc+1, "");
  loc = commandLine.find("]");
  std::string arg2 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ,
      depID =  atoi( arg2.c_str()) ;

  Epetra_RowMatrix * rowMatrix = problemManager.getBlockJacobianMatrix( probId, depID ).get();
  if( rowMatrix )
  {
    std::string name = "BJac_" + arg1 + "_" + arg2;
    engine.PutRowMatrix( *rowMatrix, name.c_str(), false );
    std::cout << "Stored Block Jacobian (" << probId << "," << depID << ") in \""
         << name << "\"" << std::endl;
    return true;
  }

  std::cout << "Could not get a valid matrix." << std::endl;
  return false;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_getAllJac::doCommand( std::string commandLine )
{
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = problemManager.getProblems().begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = problemManager.getProblems().end();

  // Do diagonal blocks
  for( ; problemLast != problemIter; ++problemIter )
  {

    GenericEpetraProblem & problem = *(*problemIter).second;
    int                    probId  = (*problemIter).first;

    Epetra_RowMatrix * rowMatrix = problemManager.getBlockJacobianMatrix( probId ).get();

    if( rowMatrix )
    {
      ostringstream sval1, sval2;
      sval1 << probId << flush;
      std::string name = "BJac_" + sval1.str() + "_" + sval1.str();
      engine.PutRowMatrix( *rowMatrix, name.c_str(), false );
      std::cout << "Stored Block Jacobian (" << probId << "," << probId << ") in \""
           << name << "\"" << std::endl;
    }

    // Do off-diagoanl blocks if appropriate
    if( problemManager.useOffBlocks() )
    {
#ifdef HAVE_NOX_EPETRAEXT
      for( unsigned int k = 0; k < problem.getDependentProblems().size(); ++k)
      {
        int depId = problem.getDependentProblems()[k];

        Epetra_RowMatrix * rowMatrix = problemManager.getBlockJacobianMatrix( probId, depId ).get();

        if( rowMatrix )
        {
          ostringstream sval1, sval2;
          sval1 << probId << flush;
          sval2 << depId  << flush;
          std::string name = "BJac_" + sval1.str() + "_" + sval2.str();
          engine.PutRowMatrix( *rowMatrix, name.c_str(), false );
          std::cout << "Stored Block Jacobian (" << probId << "," << depId << ") in \""
               << name << "\"" << std::endl;
        }
      }
#endif
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_getRes::doCommand( std::string commandLine )
{
  commandLine.replace( commandLine.find(command), command.size(), "");
  std::string::size_type loc = commandLine.find("]");
  if( std::string::npos == loc )
  {
    std::cout << "Could not get a valid argument." << std::endl;
    return false;
  }
  std::string arg1 = commandLine.substr(0, loc);
  int probId = atoi( arg1.c_str()) ;

  const Epetra_Vector * resVec = &(problemManager.getResidualVec( probId ));

  ostringstream sval1;
  sval1 << probId << flush;
  std::string name = "F_" + sval1.str();
  engine.PutMultiVector( *resVec, name.c_str() );
  std::cout << "Stored Residual (" << probId << ") in \"" << name << "\"" << std::endl;

  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_getAllRes::doCommand( std::string command )
{
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = problemManager.getProblems().begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = problemManager.getProblems().end();

  // Do diagonal blocks
  for( ; problemLast != problemIter; ++problemIter )
  {
    //GenericEpetraProblem & problem = *(*problemIter).second;
    int                    probId  = (*problemIter).first;

    const Epetra_Vector * resVec = &(problemManager.getResidualVec( probId ));

    if( resVec )
    {
      ostringstream sval1;
      sval1 << probId << flush;
      std::string name = "F_" + sval1.str();
      engine.PutMultiVector( *resVec, name.c_str() );
      std::cout << "Stored Residual (" << probId << ") in \"" << name << "\"" << std::endl;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
Coupling_Matlab_Interface::CMD_getPrecMatrix::doCommand( std::string command )
{

  std::cout << "Command currently unspupported." << std::endl;
  return false;
  //Epetra_RowMatrix & rowMatrix = *(problemManager.get_matrix());
  //engine.PutRowMatrix( rowMatrix, "PrecMatrix", false );
  return true;

}

#endif
