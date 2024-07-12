// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Problem_Manager.H"

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_DebugTools.H"
#include "NOX_Epetra_SchurCoupler.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_ConditionNumber.h"

#include "GenericEpetraProblem.H"
#include "Xfer_Operator.H"
#include "OffBlock_Manager.H"

// Headers needed for Coloring
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif

// Needed for block inverse operators
#include "Ifpack.h"

// Header for Timing info
#include "Epetra_Time.h"

#ifdef HAVE_MATLAB
#include "NOX_Multiphysics_Matlab_Interface.H"
#endif

// Hard-wired switch to turn on/off use of FD Coloring
#define USE_FD
//#undef USE_FD

// Simple macro for turning on debug code
#undef DEBUG_BLOCKGRAPH
#ifdef ENABLE_DEBUG_BLOCKGRAPH
  #define DEBUG_BLOCKGRAPH(a) a
#else
  #define DEBUG_BLOCKGRAPH(a)
#endif

Problem_Manager::Problem_Manager(Epetra_Comm& comm,
                                 bool doOffBlocks_,
                                 int numGlobalElements,
                                 bool useMatlab_ ) :
  GenericEpetraProblem(comm, numGlobalElements),
  problemCount(0),
  doOffBlocks(doOffBlocks_),
  useMatlab(useMatlab_)
{
  // Unset doOffBlocks flag if this build does not include the required
  // EpetraExt library intreface
#ifndef HAVE_NOX_EPETRAEXT
  doOffBlocks = false;
#endif

  // Reset composite number of dofs (overwrites base constructor assignment)
  NumMyNodes = 0;

  // Do all setup in registerComplete after all registrations have been
  // performed
}

//-----------------------------------------------------------------------------

Problem_Manager::~Problem_Manager()
{
}

//-----------------------------------------------------------------------------

void
Problem_Manager::addProblem(GenericEpetraProblem & problem)
{
  Problems[++problemCount] = Teuchos::rcp( &problem, false );
  problem.setId(problemCount);
  problem.setManager(this);

  // Give this problem a name if it doesn't already have one
  if( problem.getName() == "" )
  {
    std::string name = "Problem_";
    char id_str[4];
    (void) sprintf(id_str, "%d",problemCount);
    name += id_str;
    problem.setName(name);
  }

  Names[problemCount] = problem.getName();
  NameLookup[problem.getName()] = problemCount;

  // Keep a running total of dofs for use in constructing composite objects
  NumMyNodes += problem.NumMyNodes;

  return;

}

//-----------------------------------------------------------------------------

GenericEpetraProblem &
Problem_Manager::getProblem( int id_ )
{
  // Get a problem given its unique id
  Teuchos::RCP<GenericEpetraProblem> problem = Problems[id_];

  if( Teuchos::is_null(problem) )
  {
    std::cout << "ERROR: Problem with id --> " << id_ << " not registered with "
         << "Problem_Manager !!" << std::endl;
    outputStatus(std::cout);

    throw "Problem_Manager ERROR";
  }
  else
    return *problem;
}

//-----------------------------------------------------------------------------

GenericEpetraProblem &
Problem_Manager::getProblem( std::string name )
{
  // Get a problem given its name
  map<string, int>::iterator iter = NameLookup.find(name);

  if( iter == NameLookup.end() )
  {
    std::cout << "ERROR: Could not find lookup id for Problem --> " << name
         << std::endl;
    outputStatus(std::cout);

    throw "Problem_Manager ERROR";
  }
  else
    return getProblem( (*iter).second );
}

//-----------------------------------------------------------------------------

const NOX::Epetra::Group &
Problem_Manager::getProblemSolutionGroup(int id)
{
  // Get a group given its unique id
  //if( Teuchos::is_null(Solvers[id]->rcpSolver) )
  if( Teuchos::is_null(Solvers[id]) )
  {
    std::cout << "Could not get a valid NOX::Solver::Manager object for id = " << id << std::endl;
    throw "Problem_Manager ERROR";
  }

  const NOX::Epetra::Group & const_epetraSolnGrp =
    dynamic_cast<const NOX::Epetra::Group&>(Solvers[id]->getSolutionGroup());

  return const_epetraSolnGrp;
}

//-----------------------------------------------------------------------------

const Epetra_Vector &
Problem_Manager::getSolutionVec(int id)
{
  //if( Teuchos::is_null(noxEpetraSolvers[id]->rcpSolver) )
  if( Teuchos::is_null(Solvers[id]) )
  {
    std::cout << "Could not get a valid NOX::Solver::Manager object for id = " << id << std::endl;
    throw "Problem_Manager ERROR";
  }

  const Epetra_Vector & epetraSolnVec = dynamic_cast<const NOX::Epetra::Vector&>
    (getProblemSolutionGroup(id).getX()).getEpetraVector();

  return epetraSolnVec;
}

//-----------------------------------------------------------------------------

NOX::Epetra::Group &
Problem_Manager::getGroup(int id_)
{
  // Get a group given its unique id
  Teuchos::RCP<NOX::Epetra::Group> group = Groups[id_];

  if( Teuchos::is_null(group) )
  {
    std::cout << "ERROR: Could not get Group for Problem with id --> " << id_
         << " !!" << std::endl;
    throw "Problem_Manager ERROR";
  }
  else
    return *(group.get());
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_Vector>
Problem_Manager::getCompositeSoln()
{
  if( !compositeSoln.get() ) {
    std::cout << "ERROR: No valid Composite Solution vector with Problem Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }
  return compositeSoln;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::createDependency( std::string nameA, std::string nameB, bool isInterfacial )
{
  // Create a dependence of Problem A equations on Problem B variables
  int probId_A = (*(NameLookup.find(nameA))).second;
  int probId_B = (*(NameLookup.find(nameB))).second;

  if( !probId_A || !probId_B )
  {
    std::cout << "ERROR: Could not create dependency of \"" << nameA << "\" on \""
         << nameB << "\" !!" << std::endl;
    throw "Problem_Manager ERROR";
  }

  GenericEpetraProblem &probA = *(*(Problems.find(probId_A))).second,
                       &probB = *(*(Problems.find(probId_B))).second;

  createDependency( probA, probB, isInterfacial );

  return;
}

//-----------------------------------------------------------------------------


void
Problem_Manager::createDependency( GenericEpetraProblem & problemA,
                                        GenericEpetraProblem & problemB,
                                        bool isInterfacial )
{
  // Ensure that both problems already exist
  if( !problemA.getId() || !problemB.getId() ) {
    std::cout << "ERROR: Invalid dependency since at least one problem is "
         << "not registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  if( !isInterfacial )
    problemA.addTransferOp(problemB);

  problemA.addProblemDependence(problemB);

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::registerParameters(const Teuchos::RCP<Teuchos::ParameterList>& List)
{
  nlParams = List;

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::registerStatusTest(const Teuchos::RCP<NOX::StatusTest::Combo>& comboTest)
{
  statusTest = comboTest;

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::registerComplete()
{
  if(Problems.empty())
  {
    std::cout << "ERROR: No problems registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  if(Teuchos::is_null(nlParams) || Teuchos::is_null(statusTest))
  {
    std::cout << "ERROR: No nlParams and/or statusTest registered with "
         << "Problem_Manager !!" << std::endl;
    throw "Problem_Manager ERROR";
  }

  // Make sure everything is starting clean
  assert(Groups.empty() && Interfaces.empty() && Solvers.empty()
         && ProblemToCompositeIndices.empty());

  // Iterate over each problem and construct the necessary objects

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator iter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator last = Problems.end();

  iter = Problems.begin();

  int icount = 0; // Problem counter
  int runningProblemNodeCount = 0;
  int runningMaxGlobalId = 0; // Accruing max index used for establishing
                  // each problem's mapping into the composite problem

  // Create an index array for use in constructing a composite Epetra_Map
  int* compositeGlobalNodes = new int[NumMyNodes];

  // Do first pass over problems to allocate needed data.
  // This is needed prior to calling computeF below.
  while( iter != last)
  {
    // Get a convenient reference to the current problem
    GenericEpetraProblem& problem = *(*iter).second;
    int probId = problem.getId();

    // Create dependent vectors for this problem
    problem.createDependentVectors();

    // Give each problem a chance to initialize based on registrations
    problem.initialize();

    // Create index mapping for this problem into the composite problem
    ProblemToCompositeIndices[probId] = Teuchos::rcp( new Epetra_IntVector(*problem.StandardMap) );
    Epetra_IntVector &indices = *(ProblemToCompositeIndices[probId]);

    for (int i=0; i<problem.NumMyNodes; i++)
    {
      int compositeId = runningMaxGlobalId + problem.StandardMap->GID(i);
      compositeGlobalNodes[i + runningProblemNodeCount] = compositeId;
      indices[i] = compositeId;
    }

    runningProblemNodeCount += problem.NumMyNodes;
    runningMaxGlobalId += problem.StandardMap->MaxAllGID() + 1;

    iter++;
  }

  // Do second pass to setup each problem
  iter = Problems.begin();

  // Now that all problems have created vectors to receive fields from
  // problems on which they depend, sync everyone so that the subsequent
  // call to computeF below will use the correct values transferred from
  // each dependent problem.
  syncAllProblems();

  while( iter != last)
  {
    // Get a convenient reference to the current problem
    GenericEpetraProblem& problem = *(*iter).second;
    int probId = problem.getId();

    Interfaces[probId] = Teuchos::rcp(new Problem_Interface(problem));

    NOX::Epetra::Vector nox_soln( problem.getSolution(), NOX::Epetra::Vector::CreateView);

    // Use this for analytic Matrix Fills
#ifndef USE_FD
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInt = Interfaces[probId];
    LinearSystems[probId] = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
      nlParams->sublist("Printing"),
      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      Interfaces[probId],
      jacInt,
      problem.getJacobian(),
      problem.getSolution() ));

    Groups[probId] = Teuchos::rcp( new NOX::Epetra::Group(
      nlParams->sublist("Printing"),
      Interfaces[probId],
      nox_soln,
      LinearSystems[probId] ));

    // OR use this to fill matrices using Finite-Differences with Coloring
#else
#ifdef HAVE_NOX_EPETRAEXT
    // Create a timer for performance
    Epetra_Time fillTime(*Comm);

    EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
      EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN;
    // NOTE: GREEDY causes a core dump????
    int reordering = 0;
    bool useParallel = true;
    bool distance1 = false;
    int verbose = 0;

    TmpMapColorings[probId] = Teuchos::rcp(new EpetraExt::CrsGraph_MapColoring(algType, reordering, distance1, verbose));
    ColorMaps[probId] = Teuchos::rcp(&( (*(TmpMapColorings[probId]))( *problem.getGraph() ) ));
    ColorMapIndexSets[probId] = Teuchos::rcp(new EpetraExt::CrsGraph_MapColoringIndex( *(ColorMaps[probId]) ));
    ColumnsSets[probId] = Teuchos::rcp(& (*(ColorMapIndexSets[probId]))( *problem.getGraph() ));

    Teuchos::RCP<Epetra_CrsGraph> & problemGraph = problem.getGraph();

    if (MyPID == 0)
      printf("\n\tTime to color Jacobian # %d --> %e sec. \n\n",
                  icount++,fillTime.ElapsedTime());
    MatrixOperators[probId] = Teuchos::rcp( new NOX::Epetra::FiniteDifferenceColoring(
            nlParams->sublist("Printing"),
                Interfaces[probId],
        nox_soln,
                problemGraph,
                ColorMaps[probId],
                ColumnsSets[probId],
            useParallel,
            distance1) );
    //MatrixOperators[probId] = Teuchos::rcp(new NOX::Epetra::FiniteDifference(
    //            nlParams->sublist("Printing"),
    //            Interfaces[probId],
    //            nox_soln,
    //            problemGraph ));
    //NOX::Epetra::Interface::Jacobian& jacInt =
    NOX::Epetra::Interface::Jacobian * p_jacInt = dynamic_cast<NOX::Epetra::FiniteDifference*>(MatrixOperators[probId].get());
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInt = Teuchos::rcp(p_jacInt, false);
      //dynamic_cast<NOX::Epetra::Interface::Jacobian&>(*(*(MatrixOperators[probId]);
    LinearSystems[probId] = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
      nlParams->sublist("Printing"),
      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      Interfaces[probId],
      jacInt,
      MatrixOperators[probId],
      *problem.getSolution() ) );

    Groups[probId] = Teuchos::rcp( new NOX::Epetra::Group(
      nlParams->sublist("Printing"),
      Interfaces[probId],
      nox_soln,
      LinearSystems[probId] ));

#else
    if(MyPID==0)
      std::cout << "ERROR: Cannot use EpetraExt with this build !!" << std::endl;
    exit(0);
#endif
#endif

    // Needed to establish initial convergence state
    Groups[probId]->computeF();

    Solvers[probId] =
      NOX::Solver::buildSolver(Groups[probId], statusTest, nlParams);

    ++iter;
  }


  compositeMap = Teuchos::rcp( new Epetra_Map(-1, NumMyNodes, compositeGlobalNodes, 0, *Comm) );

  delete [] compositeGlobalNodes; compositeGlobalNodes = 0;

  compositeSoln = Teuchos::rcp( new Epetra_Vector(*compositeMap) );

  // --- moved from solveMF -------

  // Sync all the problems and get initial convergence state
  syncAllProblems();

  setAllGroupX();

  computeAllF();

  double normSum = getNormSum();
  std::cout << "Initial 2-Norm of composite Problem --> " << normSum;

  // Fill initial composite solution with values from each problem
  copyProblemsToComposite(*(compositeSoln.get()), SOLUTION);

  // Set up a problem interface for the Problem Manager
  Teuchos::RCP<Problem_Interface> interface = Teuchos::rcp(new Problem_Interface(*this));
  // Now create a composite matrix graph needed for preconditioning
  AA = Teuchos::rcp( new Epetra_CrsGraph(Copy, *compositeMap, 0) );
  generateGraph();

/* --------------  Block for Coloring Preconditioner Operator ------

  // NOT YET WORKING
  // This needs more work to deal with the parallel-use coloring capability.

  // We now attempt to use Coloring on the global preconditioning matrix
  // Create the Epetra_RowMatrix using Finite Difference with Coloring
#ifdef HAVE_NOX_EPETRAEXT
  bool verbose = false;
  EpetraExt::CrsGraph_MapColoring tmpMapColoring( verbose );
  Epetra_MapColoring* colorMap = &tmpMapColoring(*AA);
  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  std::vector<Epetra_IntVector>* columns = &colorMapIndex(*AA);
#else
  if(MyPID==0)
    std::cout << "ERROR: Cannot use EpetraExt with this build !!" << std::endl;
  exit(0);
#endif

  // Use this constructor to create the graph numerically as a means of timing
  // the old way of looping without colors :
  //  NOX::Epetra::FiniteDifferenceColoring A(interface, soln,
  //                                          *colorMap, *columns);
  // Or use this as the standard way of using finite differencing with coloring
  // where the application is responsible for creating the matrix graph
  // beforehand, ie as is done in Problem.
  NOX::Epetra::FiniteDifferenceColoring* A =
    new NOX::Epetra::FiniteDifferenceColoring(interface, compositeSoln, *AA,
                                              *colorMap, *columns);
// --------  End of Block for Coloring Preconditioner Operator ------ */


  // Create a preconditioning matrix using the graph just created - this
  // creates a static graph so we can refill the new matirx after
  // FillComplete()  is called.
  A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *AA) );
  A->FillComplete();

  //NOX::Epetra::Interface::Required& reqInt =
  //  dynamic_cast<NOX::Epetra::Interface::Required&>(interface);
  Teuchos::RCP<NOX::Epetra::Vector> compositeNOXSoln =
    Teuchos::rcp( new NOX::Epetra::Vector(*(compositeSoln.get())) );

  Teuchos::ParameterList & printParams = nlParams->sublist("Printing");
  Teuchos::ParameterList & lsParams    = nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver");
  Teuchos::RCP<NOX::Utils> utilsPtr = Teuchos::rcp( new NOX::Utils(printParams) );

  // Setup appropriate operators and their use

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> composite_linearSystem;

  if( 1 ) // default for MultiPhyics tests
  {
    // Use Matrix-Free Jacobian Operator with FDC block preconditioner
    jacOperator = Teuchos::rcp( new NOX::Epetra::MatrixFree( printParams, interface, *compositeSoln) );
    NOX::Epetra::Interface::Jacobian * p_jacInt = dynamic_cast<NOX::Epetra::MatrixFree*>(jacOperator.get());
    jacInterface = Teuchos::rcp(p_jacInt, false);
    precOperator = Teuchos::rcp( A.get(), false );
    precInterface = interface;

    composite_linearSystem = Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO(
        nlParams->sublist("Printing"),
        lsParams,
        jacInterface, jacOperator,
        precInterface, precOperator,
        compositeSoln) );

    dynamic_cast<NOX::Epetra::MatrixFree*>(jacOperator.get())->computeJacobian(*compositeSoln, *jacOperator);
  }

  if( 0 ) // Matrix-Free Jacobian and Broyden preconditioner
  {
    evaluate( NOX::Epetra::Interface::Required::Jac, &(*compositeSoln), NULL );
    jacOperator = Teuchos::rcp( new NOX::Epetra::MatrixFree( printParams, interface, *compositeSoln) );
    NOX::Epetra::Interface::Jacobian * p_jacInt = dynamic_cast<NOX::Epetra::MatrixFree*>(jacOperator.get());
    jacInterface = Teuchos::rcp(p_jacInt, false);

    Teuchos::RCP<NOX::Epetra::BroydenOperator> broydenOp =
      Teuchos::rcp(new NOX::Epetra::BroydenOperator( *nlParams, utilsPtr, *compositeSoln, A, true) );
    precOperator = broydenOp;
    precInterface = broydenOp;

    composite_linearSystem = Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO(
        nlParams->sublist("Printing"),
        lsParams,
        jacInterface, jacOperator,
        precInterface, precOperator,
        compositeSoln) );
  }

  if( 0 ) // Broyden Jacobian and preconditioner
  {
    // Turn on some debugging
    Teuchos::ParameterList & broydenParams =
      nlParams->sublist("Direction").sublist("Newton").sublist("Broyden Op");
    broydenParams.set("Write Broyden Info", true);

    // Create some temp objects needed for construction
    Teuchos::RCP<Epetra_CrsGraph>  blockAA = generateBlockDiagonalGraph( *compositeMap, false );
    blockAA->FillComplete();
    Teuchos::RCP<Epetra_CrsMatrix> blockA = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *blockAA) );
    blockA->FillComplete();
    evaluate( NOX::Epetra::Interface::Required::Jac, &(*compositeSoln), NULL );
    copyProblemJacobiansToComposite( *blockA );

    std::cout << blockA->Graph() << std::endl;
    std::cout << *blockA << std::endl;

    Teuchos::RCP<NOX::Epetra::BroydenOperator> broydenOp =
      Teuchos::rcp(new NOX::Epetra::BroydenOperator( *nlParams, utilsPtr, *compositeSoln, blockA, true) );
    jacOperator = broydenOp;
    jacInterface = broydenOp;

    composite_linearSystem = Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO(
        nlParams->sublist("Printing"),
        lsParams,
        interface,
        jacInterface, jacOperator,
        compositeSoln) );

    // To allow replacement of diagonal block values, register ourself as the
    // needed interface
    //broydenOp->addReplacementInterface( this );
    //broydenOp->removeEntriesFromBroydenUpdate( A->Graph() );
  }

  if( 0 )
  {
    // Use FDC Jacobian and itself as its own preconditioner

    jacOperator = A;
    jacInterface = interface;

    composite_linearSystem = Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO(
        nlParams->sublist("Printing"),
        lsParams,
        interface,
        jacInterface, jacOperator,
        compositeSoln) );
  }

  //lsParams.set("Preconditioning", "None");
  lsParams.set("Preconditioner", "AztecOO");
  compositeGroup = Teuchos::rcp(new NOX::Epetra::Group(nlParams->sublist("Printing"),
                    interface, *compositeNOXSoln,
                    composite_linearSystem));
  compositeGroup->computeF();

  compositeSolver =
    NOX::Solver::buildSolver(compositeGroup, statusTest, nlParams);

#ifdef HAVE_MATLAB
  // Just a test to see if the Matlab engine is working - RWH
  if( useMatlab )
  {
    //Matlab_Interface testMatlab(*compositeSolver);
    Coupling_Matlab_Interface testMatlab(*this);

    testMatlab.interact();
  }
#endif

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::syncAllProblems()
{
  if(Problems.empty())
  {
    std::cout << "ERROR: No problems registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter;
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemBegin = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast  = Problems.end();

  // Give each problem an opportunity to do things before transferring data, eg compute fluxes
  for( problemIter = problemBegin ; problemIter != problemLast; ++problemIter )
    (*problemIter).second->prepare_data_for_transfer();

  // Loop over each problem being managed and invoke its transfer requests
  for( problemIter = problemBegin ; problemIter != problemLast; ++problemIter )
    (*problemIter).second->doTransfer();

  // Give each problem an opportunity to do things after transferring data
  for( problemIter = problemBegin ; problemIter != problemLast; ++problemIter )
    (*problemIter).second->process_transferred_data();

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::setAlldt( double dt )
{
  if(Problems.empty()) {
    std::cout << "ERROR: No problems registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  for( ; problemIter != problemLast; ++problemIter)
    (*problemIter).second->setdt(dt);

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::setGroupX(int probId)
{
  Teuchos::RCP<GenericEpetraProblem> problem = Problems[probId];

  if( Teuchos::is_null(problem) )
  {
    std::cout << "ERROR: Could not get requested Problem to use with group.setX "
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  Teuchos::RCP<NOX::Epetra::Group> grp = Groups[probId];

  if( Teuchos::is_null(grp) ) {
    std::cout << "ERROR: Could not get appropriate group for use in setX !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  grp->setX(*problem->getSolution());

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::setGroupX(int probId, Epetra_Vector & vec)
{
  Teuchos::RCP<GenericEpetraProblem> problem = Problems[probId];

  if( Teuchos::is_null(problem) )
  {
    std::cout << "ERROR: Could not get requested Problem to use with group.setX " << std::endl;
    throw "Problem_Manager ERROR";
  }

  Teuchos::RCP<NOX::Epetra::Group> grp = Groups[probId];

  if( Teuchos::is_null(grp) )
  {
    std::cout << "ERROR: Could not get appropriate group for use in setX !!" << std::endl;
    throw "Problem_Manager ERROR";
  }

  grp->setX( vec );
}

//-----------------------------------------------------------------------------

void
Problem_Manager::setAllGroupX()
{
  if(Problems.empty()) {
    std::cout << "ERROR: No problems registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Loop over each problem being managed and set the corresponding group
  // solution vector (used by NOX) with the problem's (used by application)
  for( ; problemIter != problemLast; problemIter++)
  {
    int probId = (*problemIter).first;
    setGroupX(probId);
  }

  return;
}

//-----------------------------------------------------------------------------

#ifdef HAVE_NOX_EPETRAEXT
void
Problem_Manager::setAllOffBlockGroupX(const Epetra_Vector &inVec)
{
  map<int, std::vector<OffBlock_Manager*> >::iterator offBlockIter = OffBlock_Managers.begin();
  map<int, std::vector<OffBlock_Manager*> >::iterator offBlockLast = OffBlock_Managers.end();

  // Loop over each off-block manager and set the contained groups X-vector
  // with the incoming vector
  for( ; offBlockIter != offBlockLast; ++offBlockIter)
  {
    std::vector<OffBlock_Manager*> managerVec = (*offBlockIter).second;

    // Need to extract block portion of incoming compositeVec
    for( unsigned int i = 0; i < managerVec.size(); ++i )
    {
      // Note that we assign the group soln to be from the problemVar
      int probVarId = managerVec[i]->getProblemVarId();
      copyCompositeToVector(inVec, probVarId, *(managerVec[i]->getRowMapVec()) );
      managerVec[i]->getGroup()->setX(managerVec[i]->getRowMapVec());
    }
  }

  return;
}

#endif

//-----------------------------------------------------------------------------

void
Problem_Manager::resetProblems()
{

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Loop over each problem and copy its solution into its old solution
  for( ; problemIter != problemLast; problemIter++)
  {
    GenericEpetraProblem & problem = *(*problemIter).second;
    problem.reset( *problem.getSolution() );
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::computeAllF()
{
  if(Problems.empty()) {
    std::cout << "ERROR: No problems registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Loop over each problem being managed and invoke the corresponding group's
  // residual evaluation
  for( ; problemIter != problemLast; problemIter++)
  {
    int probId = (*problemIter).first;

    computeGroupF(probId);
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::computeGroupF(int probId)
{
  Teuchos::RCP<NOX::Epetra::Group> grp = Groups[ probId ];

  if( Teuchos::is_null(grp) )
  {
    std::cout << "ERROR: Could not get a group for problem with id --> "
         << probId << std::endl;
    throw "Problem_Manager ERROR";
  }
  grp->computeF();

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::computeAllJacobian()
{
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Do diagoanl blocks
  for( ; problemLast != problemIter; ++problemIter )
    computeBlockJacobian( (*problemIter).first );

  // Do off-diagoanl blocks if appropriate
  if( doOffBlocks )
  {
#ifdef HAVE_NOX_EPETRAEXT
    problemIter = Problems.begin();

    for( ; problemLast != problemIter; ++problemIter )
      for( unsigned int k = 0; k < (*problemIter).second->depProblems.size(); ++k)
        computeBlockJacobian( (*problemIter).first, (*problemIter).second->depProblems[k] );
#endif
  }

  if( 0 )
  {
    // Create a timer for performance
    Epetra_Time fillTime(*Comm);

    // Loop over each problem being managed and invoke its computeJacobian
    // method
    for( ; problemIter != problemLast; problemIter++)
    {
      int probId = (*problemIter).first;
      Teuchos::RCP<NOX::Epetra::Group> grp = Groups[probId];
      if( Teuchos::is_null(grp) )
      {
        std::cout << "ERROR: Could not find valid group for compouteJacobian !!"
             << std::endl;
        throw "Problem_Manager ERROR";
      }

      fillTime.ResetStartTime();

      grp->computeJacobian();

      if (MyPID == 0)
        printf("\n\tTime to fill Jacobian %d --> %e sec. \n\n",
                    probId, fillTime.ElapsedTime());


      if( doOffBlocks )
      {
  #ifdef HAVE_NOX_EPETRAEXT

        fillTime.ResetStartTime();

        std::vector<OffBlock_Manager*> & offBlocksVec = OffBlock_Managers[probId];

        for( unsigned int i = 0; i < offBlocksVec.size(); ++i )
        {
          // Refresh all problem vectors
          copyCompositeToProblems(*(compositeSoln.get()), SOLUTION);

    DEBUG_BLOCKGRAPH( std::cout << "Doing computeJacobian for : " << offBlocksVec[i]->getName() << std::endl;)

          offBlocksVec[i]->getGroup()->computeJacobian();

    DEBUG_BLOCKGRAPH( std::cout << "For block : " << offBlocksVec[i]->getName() << std::endl;)
    DEBUG_BLOCKGRAPH( offBlocksVec[i]->getMatrix().Print(cout);)

          if (MyPID == 0)
            printf("\n\tTime to fill Jacobian %d (%d) --> %e sec. \n\n",
                        probId, i, fillTime.ElapsedTime());
        }
  #endif
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::computeBlockJacobian(int probId, int depId)
{

  Epetra_Time fillTime(*Comm);

  // If this is a self-dependence, simply call the corresponding Group's computeJacobian
  if( (probId == depId) || (-1 == depId) ) // diagonal block contribution
  {
    fillTime.ResetStartTime();

    getGroup(probId).computeJacobian();

    if (MyPID == 0)
      printf("\n\tTime to fill Jacobian %d --> %e sec. \n\n",
                  probId, fillTime.ElapsedTime());
  }
  else // off-diagonal block contribution
  {
    std::vector<OffBlock_Manager*> managerVec = OffBlock_Managers[probId];

    OffBlock_Manager::idToFind = depId;

    std::vector<OffBlock_Manager *>::iterator iter =
      find_if( managerVec.begin(), managerVec.end(), mem_fn( &OffBlock_Manager::isMember) );

    if( managerVec.end() == iter )
    {
      std::string msg =
            "ERROR Problem_Manager::computeBlockJacobian : Problem \""
             + Problems[probId]->getName()
             + "\" does not have a registered dependence on Problem \""
             + Problems[depId]->getName() + "\"\n";
      throw msg;
    }

    fillTime.ResetStartTime();

    // Refresh all problem vectors
    copyCompositeToProblems(*(compositeSoln.get()), SOLUTION);

    (*iter)->getGroup()->computeJacobian();

    if (MyPID == 0)
      printf("\n\tTime to fill Jacobian %d (%d) --> %e sec. \n\n",
                  probId, depId, fillTime.ElapsedTime());
  }

  return;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_CrsMatrix>
Problem_Manager::getBlockJacobianMatrix(int probId, int depId)
{

  Epetra_CrsMatrix * p_problemMatrix = NULL;

  if( (probId == depId) || (-1 == depId) ) // diagonal block
  {
    // Get objects holding our desired matrix
    GenericEpetraProblem & problem = *(Problems[probId]);
    NOX::Epetra::LinearSystemAztecOO & problemLinearSystem = *(LinearSystems[probId]);

    // Use each group's operator test to determine the type of Jacobian
    // operator being used.
    const Epetra_Operator& jacOp = *(problemLinearSystem.getJacobianOperator().get());

    if ( dynamic_cast<const NOX::Epetra::FiniteDifference*>(&jacOp) )
      p_problemMatrix = const_cast<Epetra_CrsMatrix*>(
        &dynamic_cast<const NOX::Epetra::FiniteDifference&>(jacOp).getUnderlyingMatrix());
    else if ( dynamic_cast<const Epetra_CrsMatrix*>(&jacOp) )
      // NOTE: We are getting the matrix from the problem.  This SHOULD be
      // the same matrix wrapped in the group.  A safer alternative would be
      // to get this matrix from the group as above for a more general
      // operator.
      p_problemMatrix = problem.getJacobian().get();

  }
  else // off-diagonal block
  {
    std::vector<OffBlock_Manager*> offVec = OffBlock_Managers[probId];

    OffBlock_Manager::idToFind = depId;

    std::vector<OffBlock_Manager *>::iterator iter =
      find_if( offVec.begin(), offVec.end(), mem_fn( &OffBlock_Manager::isMember) );

    if( offVec.end() == iter )
    {
      std::string msg =
            "ERROR Problem_Manager::getBlockJacobianMatrix : Problem \""
             + Problems[probId]->getName()
             + "\" does not have a registered dependence on Problem \""
             + Problems[depId]->getName() + "\"\n";
      throw msg;
    }

    p_problemMatrix = &((*iter)->getMatrix());
  }

  Teuchos::RCP<Epetra_CrsMatrix> mat = Teuchos::rcp( p_problemMatrix, false );

  return mat;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_Operator>
Problem_Manager::getBlockInverseOperator(int probId)
{

  if( BlockInverseOperators.end() == BlockInverseOperators.find(probId) )
  {
    std::cout << "ERROR: No valid Block Inverse Operator exists for problem # "
         << probId << std::endl;
    throw "Problem_Manager ERROR";
  }

  Teuchos::RCP<Epetra_Operator> pInvOp = BlockInverseOperators[probId];

  return pInvOp;
}

//-----------------------------------------------------------------------------

const Epetra_Vector &
Problem_Manager::getResidualVec(int probId )
{

  const Epetra_Vector & vec = dynamic_cast<const NOX::Epetra::Vector&>
                 (getProblemSolutionGroup(probId).getF()).getEpetraVector();

  return vec;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyGroupCurrentXtoProblemX(int probId)
{
  // Copy Current solution X from NOX solver group into the problem's solution vector

  Teuchos::RCP<GenericEpetraProblem> problem = Problems[probId];

  if( Teuchos::is_null(problem) )
  {
    std::cout << "ERROR: Could not get requested Problem to update with final "
         << "solution" << std::endl;
    throw "Problem_Manager ERROR";
  }

  Teuchos::RCP<NOX::Solver::Generic> solver = Solvers[probId];

  if( Teuchos::is_null(solver) )
  {
    std::cout << "ERROR: Could not get appropriate Solver for use in update !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>
      (finalGroup.getX())).getEpetraVector();

  *problem->getSolution() = finalSolution;

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyAllGroupXtoProblems()
{
  // Copy final solution from NOX solvers into each problem's solution vector

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  for( ; problemIter != problemLast; ++problemIter )
  {
    int probId = (*problemIter).first;
    copyGroupCurrentXtoProblemX(probId);
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::resetCurrentGroupX(int probId)
{
  // Reset each solvers current solution group X vector with itself, thereby resetting all valid flags

  Teuchos::RCP<NOX::Solver::Generic> solver = Solvers[probId];

  if( Teuchos::is_null(solver) )
  {
    std::cout << "ERROR: Could not get appropriate Solver for use in update !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  const NOX::Epetra::Group& currentGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const_cast<NOX::Epetra::Group&>(currentGroup).setX( currentGroup.getX() );

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::resetAllCurrentGroupX()
{
  // Reset each solvers current solution group X vector with itself, thereby resetting all valid flags

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  for( ; problemIter != problemLast; ++problemIter )
  {
    int probId = (*problemIter).first;
    resetCurrentGroupX(probId);
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyCompositeToProblems( const Epetra_Vector& compositeVec, Problem_Manager::VectorType vecType )
{
  // Copy a composite problem vector to each problem's vector

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  Epetra_Vector* problemVec(0);

  // Loop over each problem being managed and copy into the correct problem vector
  for( ; problemIter != problemLast; ++problemIter )
  {
    int probId = (*problemIter).first;
    switch (vecType)
    {
      case SOLUTION :
        problemVec = ((*problemIter).second->getSolution()).get();
    break;

      case GROUP_F :
      default :
        std::cout << "ERROR: vecType not supported for copy FROM composite!!"
             << std::endl;
        throw "Problem_Manager ERROR";
    }

    copyCompositeToVector(compositeVec, probId, *problemVec);
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyProblemsToComposite( Epetra_Vector& compositeVec, Problem_Manager::VectorType vecType )
{
  // Copy vectors from each problem into a composite problem vector

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  const Epetra_Vector* problemVec(0);

  // Loop over each problem being managed and copy from the correct
  // problem vector
  for( ; problemIter != problemLast; ++problemIter )
  {
    int probId = (*problemIter).first;
    switch (vecType) {

      case SOLUTION :
        problemVec = ((*problemIter).second->getSolution()).get();
    break;

      case GROUP_F :
        problemVec = &(dynamic_cast<const NOX::Epetra::Vector&>
               (Groups[probId]->getF()).getEpetraVector());
    break;

      default :
        std::cout << "ERROR: vecType not supported for copy TO composite!!"
             << std::endl;
        throw "Problem_Manager ERROR";
    }

    copyVectorToComposite(compositeVec, probId, *problemVec);
  }

  return;
}

//-----------------------------------------------------------------------------


void
Problem_Manager::copyCompositeToVector( const Epetra_Vector& compositeVec, int id, Epetra_Vector& problemVec )
{
  // Copy part of a composite problem vector to a problem's vector
  Epetra_IntVector& indices = *(ProblemToCompositeIndices[id]);
  // This map is needed to get correct LID of indices
  const Epetra_BlockMap &map = compositeVec.Map();

  for( int i = 0; i < problemVec.MyLength(); ++i )
    problemVec[i] = compositeVec[map.LID(indices[i])];

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyVectorToComposite( Epetra_Vector& compositeVec, int id, const Epetra_Vector& problemVec)
{
  // Copy a vector from a problem into part of a composite problem vector
  Epetra_IntVector& indices = *(ProblemToCompositeIndices[id]);

  for( int i = 0; i < problemVec.MyLength(); ++i )
    compositeVec[indices[i]] = problemVec[i];

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyProblemJacobiansToComposite()
{
  copyProblemJacobiansToComposite(*A);

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::copyProblemJacobiansToComposite( Epetra_CrsMatrix & mat )
{
  // Copy problem Jacobians as block diagonal contributions to
  // composite Jacobian

  Epetra_CrsMatrix & compositeMatrix = mat;

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  int problemMaxNodes = compositeSoln.get()->GlobalLength();

  // Loop over each problem being managed and copy its Jacobian into
  // the composite diagonal blocks
  for( ; problemIter != problemLast; problemIter++)
  {
    int probId = (*problemIter).first;

    // Get the problem, its Jacobian graph and its linear system
    GenericEpetraProblem & problem = *((*problemIter).second);
    Epetra_CrsGraph & problemGraph = *problem.getGraph();
    NOX::Epetra::LinearSystemAztecOO & problemLinearSystem = *(LinearSystems[probId]);

    // Get the indices map for copying data from this problem into
    // the composite problem
    Epetra_IntVector & indices = *(ProblemToCompositeIndices[probId]);

    // Each problem's Jacobian will be determined by type
    Epetra_CrsMatrix * p_problemMatrix(0);

    // Use each group's operator test to determine the type of Jacobian
    // operator being used.
    const Epetra_Operator& jacOp = *(problemLinearSystem.getJacobianOperator().get());

    if ( dynamic_cast<const NOX::Epetra::FiniteDifference*>(&jacOp) )
      p_problemMatrix = const_cast<Epetra_CrsMatrix*>(
        &dynamic_cast<const NOX::Epetra::FiniteDifference&>(jacOp).getUnderlyingMatrix());
    else if ( dynamic_cast<const Epetra_CrsMatrix*>(&jacOp) )
      // NOTE: We are getting the matrix from the problem.  This SHOULD be
      // the same matrix wrapped in the group.  A safer alternative would be
      // to get this matrix from the group as above for a more general
      // operator.
      p_problemMatrix = problem.getJacobian().get();
    else
    {
      if (MyPID==0)
        std::cout << "Jacobian operator for Problem not supported for "
             << "preconditioning Matrix-Free coupling solver." << std::endl;
      throw "Problem_Manager ERROR";
    }

    // Create convenient reference for each Jacobian matrix
    Epetra_CrsMatrix &problemMatrix = *p_problemMatrix;


    // ROGER Block Norm
    if (1)
    {
      std::cout << "Block=" << probId << "  Inf Norm=" << problemMatrix.NormInf()
       << "  One Norm=" << problemMatrix.NormOne() << std::endl;
    }

    // Temporary storage arrays for extracting/inserting matrix row data
    int* columnIndices = new int[problemMaxNodes];
    double* values = new double[problemMaxNodes];

    int problemRow, compositeRow, numCols, numValues;

    for (int i = 0; i<problemMatrix.NumMyRows(); i++)
    {
      problemRow = problemMatrix.Map().GID(i);
      problemGraph.ExtractGlobalRowCopy(problemRow, problemMaxNodes, numCols, columnIndices);
      problemMatrix.ExtractGlobalRowCopy(problemRow, problemMaxNodes, numValues, values);
      if( numCols != numValues )
      {
        if (MyPID==0)
          std::cout << "ERROR: Num Columns != Num Values from problem Matrix !!"
               << std::endl;
        throw "Problem_Manager ERROR";
      }

      // Convert row/column indices to composite problem
      compositeRow = indices[problemRow];

      for (int j = 0; j<numCols; j++)
        columnIndices[j] = indices[columnIndices[j]];

      int ierr = compositeMatrix.ReplaceGlobalValues(compositeRow,
                       numValues, values, columnIndices);
      if( ierr )
      {
        if (MyPID==0)
          std::cout << "ERROR: compositeMatrix.ReplaceGlobalValues(...)" << std::endl;
        throw "Problem_Manager ERROR";
      }
    }
    delete [] values; values = 0;
    delete [] columnIndices; columnIndices = 0;

    // Sync up processors to be safe
    Comm->Barrier();

    // Add off-diagonal FD block contributions if waranted
    if( doOffBlocks )
    {
#ifdef HAVE_NOX_EPETRAEXT
      // Loop over each problem on which this one depends
      for( unsigned int k = 0; k < problem.depProblems.size(); ++k)
      {
        Teuchos::RCP<GenericEpetraProblem> p_depProblem = Problems[problem.depProblems[k]];

        // Copy the off-block jacobian matrices for this
        // problem-problem coupling
        //  !!! --------  THESE COMMENTS ARE NO LONGER VALID --------------!!!
        //  *******************************************************************
        //  *** NOTE: the map used for the off-block graph is the composite Map
        //  *** to allow valid global indexing.  This also allows us to
        //  *** simply copy values directly from the off-block matrices to the
        //  *** composite Jacobian matrix
        //  *******************************************************************
        OffBlock_Manager * p_offBlockMgr = (OffBlock_Managers[probId])[k];

    if( !p_offBlockMgr )
        {
          std::cout << "ERROR: Unable to get OffBlock_Manager for dependence of problem "
               << problem.getName() << " on problem " << p_depProblem->getName()
               << " !!" << std::endl;
          throw "Problem_Manager ERROR";
        }
    Epetra_CrsMatrix & offMatrix = p_offBlockMgr->getMatrix();

        int *        blockColIndices; //  = new int[problemMaxNodes];
        int          numCols = -1;
        double *     values; //  = new double[problemMaxNodes];
        std::vector<int>  compositeColIndices;
        int          compositeRow;

        // Loop over each row and copy into composite matrix
        for (int i = 0; i < offMatrix.NumMyRows(); ++i)
        {
          offMatrix.ExtractMyRowView(i, numCols, values, blockColIndices);

          compositeColIndices.resize(numCols);

          p_offBlockMgr->convertBlockRowIndicesToComposite(1, &i, &compositeRow);
          p_offBlockMgr->convertBlockColIndicesToComposite( numCols, blockColIndices, &(compositeColIndices[0]) );
          compositeMatrix.ReplaceMyValues( compositeRow, numCols, values, &(compositeColIndices[0]) );
        }

        // Sync up processors to be safe
        Comm->Barrier();
      }
#endif
    }
  }

#ifdef DEBUG_PROBLEM_MANAGER
    compositeMatrix.Print(cout);
#endif

  return;
}

//-----------------------------------------------------------------------------

double
Problem_Manager::getNormSum()
{
  // Get each problem's residual norm and sum into a total
  double normSum(0.0);

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  for( ; problemIter != problemLast; problemIter++) {
    int probId = (*problemIter).first;

    //Teuchos::RCP<NOX::Epetra::Group> grp = Groups[ probId ];
    NOX::Epetra::Group grp = const_cast<NOX::Epetra::Group&>(dynamic_cast<const NOX::Epetra::Group&>(Solvers[ probId ]->getSolutionGroup()));

    //if( Teuchos::is_null(grp) )
    //{
    //  std::cout << "ERROR: Could not get appropriate group for use in NormSum !!"
    //       << std::endl;
    //  throw "Problem_Manager ERROR";
    //}

    //double problemNorm = grp->getNormF();
    double problemNorm = grp.getNormF();
    std::cout << "2-Norm of Problem " << probId << " --> " << problemNorm
         << std::endl;
    normSum += problemNorm * problemNorm;
  }
  normSum = sqrt(normSum);

  return normSum;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::solve()
{
  if(Problems.empty())
  {
    std::cout << "ERROR: No problems registered with Problem_Manager !!" << std::endl;
    throw "Problem_Manager ERROR";
  }

  if( Groups.empty() || Interfaces.empty() || Solvers.empty() )
  {
    std::cout << "ERROR: Groups, Interfaces and/or Solvers are emptry !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Sync all the problems and get initial convergence state
  syncAllProblems();
  setAllGroupX();
  computeAllF();

  double normSum = getNormSum();
  std::cout << "Initial 2-Norm of composite Problem --> " << normSum;

  // Now do the decoupled solve
  int iter = 0;
  NOX::StatusTest::StatusType status;

  // RPP: Inner status test
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(15));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteValue =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
  combo->addStatusTest(finiteValue);

  while( iter < 50 && normSum > 1.0e-8 ) // Hard-coded convergence criterion for now.
  {
    ++iter;

    problemIter = Problems.begin();

    // Solve each problem in the order it was registered
    for( ; problemIter != problemLast; problemIter++)
    {
      GenericEpetraProblem & problem = *(*problemIter).second;

      int probId = problem.getId();

      Teuchos::RCP<NOX::Epetra::Group> problemGroup = Groups[probId];
      NOX::Solver::Generic & problemSolver = *(Solvers[probId]);

      // Sync all dependent data with this problem
      problem.doTransfer();
      // Sync the problem solution with its solver group
      setGroupX(probId);
      // Reset the solver for this problem and solve
      //problemSolver.reset(problemGroup, *statusTest, *nlParams);

      problemSolver.reset(problemGroup->getX(), combo);
      status = problemSolver.solve();
      if( status != NOX::StatusTest::Converged )
      {
        if (MyPID==0)
          std::cout << "\nRegistered Problem ## failed to converge !!"  << std::endl;
      }

      // Copy final solution from group into problem solution vector
      copyGroupCurrentXtoProblemX(probId);
    }

    // Determine final residuals for use in testing convergence
    syncAllProblems();
    setAllGroupX();
    computeAllF();

    normSum = getNormSum();
    std::cout << "iter #" << iter << ", 2-Norm of composite Problem --> "
         << normSum << std::endl;
  }

  if (normSum > 1.0e-8)
  {
    std::cout << "Warning: composite problem failed to converge after "
         << iter << " fixed-point iterations." << std::endl;

    return false;
  }

  std::cout << "\nDecoupled solution required --> " << iter << " iterations.\n"
       << std::endl;

  return true;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::solveMF()
{
  if(Problems.empty())
  {
    std::cout << "ERROR: No problems registered with Problem_Manager !!"
         << std::endl;
    throw "Problem_Manager ERROR";
  }

  compositeGroup->setX( compositeGroup->getX() );

  compositeSolver->reset(compositeGroup->getX(), statusTest);

  NOX::StatusTest::StatusType status = compositeSolver->solve();
  if( status != NOX::StatusTest::Converged )
    if (MyPID==0)
      std::cout << "\nMatrix-Free coupled Problem failed to converge !!"  << std::endl;

  // Update all problem's solutions with final solutions from solvers
  copyAllGroupXtoProblems();

  return true;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::solveSchurBased()
{
  if( 2 != Problems.size() )
  {
    std::cout << "ERROR: Schur-based coupling currently requires exactly two problems!" << std::endl;
    throw "Problem_Manager ERROR";
  }

  Teuchos::RCP<Epetra_Vector> resA = Teuchos::rcp( new Epetra_Vector( getResidualVec(1)) );
  Teuchos::RCP<Epetra_Vector> resB = Teuchos::rcp( new Epetra_Vector( getResidualVec(2)) );

  NOX::Epetra::SchurCoupler schurCoupler(*this);
  NOX::Epetra::SchurOp & schurCplOp = *schurCoupler.getSchurOperator(0, 0);
  //NOX::Epetra::SchurOp schurCplOp( 1, 2, *this );

  Teuchos::ParameterList& lsParams = nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver");

  setGroupX(1);
  setGroupX(2);
  getGroup(1).computeF();
  getGroup(2).computeF();
  getGroup(1).computeJacobian();
  getGroup(2).computeJacobian();
  createBlockInverseOperator(1, lsParams);
  createBlockInverseOperator(2, lsParams);

  //schurCplOp.Apply( *tmpA1, *tmpA2 );

  Epetra_Vector solutionA(getSolutionVec(1));
  solutionA.PutScalar(0.0);
  std::cout << "Solution A: \n" << solutionA << std::endl;

  std::cout << "Original RHS :\n" << *resA << std::endl;
  schurCplOp.modifyRHS( resA, resB );
  std::cout << "Modified RHS :\n" << *resA << std::endl;

  Epetra_LinearProblem * linear_problem = new Epetra_LinearProblem(&schurCplOp, &solutionA, &(*resA));
  AztecOO * aztecoo_solver = new AztecOO(*linear_problem);

  aztecoo_solver->SetAztecOption(AZ_precond, AZ_none);

  // Test use of a diagonal block preconditioner
  if( 1 )
  {
    double conditionNumberEstimate = -1.0;
    Teuchos::RCP<Epetra_CrsMatrix> pMatrix = getBlockJacobianMatrix(1);
    aztecoo_solver->SetPrecMatrix(pMatrix.get());
    aztecoo_solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    aztecoo_solver->SetAztecOption(AZ_overlap, 0);
    aztecoo_solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
    aztecoo_solver->SetAztecOption(AZ_graph_fill, 0);
    aztecoo_solver->ConstructPreconditioner(conditionNumberEstimate);
  }

  std::cout << "... Solving.\n";
  aztecoo_solver->Iterate(100, 1.0e-8);
  std::cout << "... Solved to a tolerance of " << aztecoo_solver->TrueResidual() << " in " << aztecoo_solver->NumIters() << " iterations.\n";

  std::cout << "\n----- Solution for update vector ---- :\n" << solutionA << std::endl;

  *resA = solutionA;

  solutionA.Update(1.0, getSolutionVec(1), -1.0);

  std::cout << "\n----- Solution vector ---- :\n" << solutionA << std::endl;

  // Now back substitue for other problem solution
  applyBlockAction( 2, 1, *resA, *resB);
  resB->Update(-1.0, getResidualVec(2), 1.0);
  getBlockInverseOperator(2)->ApplyInverse(*resB, *resB);

  std::cout << "\n----- Second Solution for update vector ---- :\n" << *resB << std::endl;

  resB->Update(1.0, getSolutionVec(2), 1.0);

  std::cout << "\n----- Second Solution vector ---- :\n" << *resB << std::endl;

  exit(0);
  // Update all problem's solutions with final solutions from solvers
  copyAllGroupXtoProblems();

  return true;
}

//-----------------------------------------------------------------------------


// These methods are needed to allow inheritance from GenericEpetraProblem base

bool
Problem_Manager::evaluate(
              NOX::Epetra::Interface::Required::FillType flag,
              const Epetra_Vector *solnVector,
              Epetra_Vector *rhsVector)
{
  //Determine what to fill (F or Jacobian)
  bool fillF      = false;
  bool fillMatrix = false;

  if (rhsVector != 0)
    fillF = true;
  else
    fillMatrix = true;

  // Copy incoming vector from NOX solver into our composite solution Vector
  *(compositeSoln.get()) = *solnVector;

  copyCompositeToProblems(*solnVector, SOLUTION);

  // If used, give each off-block FDC manager a copy of the current total
  // solution vector
#ifdef HAVE_NOX_EPETRAEXT
  if( fillMatrix && doOffBlocks )
    setAllOffBlockGroupX(*solnVector);
#endif

  // Do transfers from problem solution vectors into problem dependent vectors
  syncAllProblems();

  // Set each problem group Xvec with its problem solution vector
  setAllGroupX();

  if (fillF)
  {
    computeAllF();
    copyProblemsToComposite(*rhsVector, GROUP_F);

  }

  if (fillMatrix)
  {
    A->PutScalar(0.0);

    computeAllJacobian();

    copyProblemJacobiansToComposite();

    A->FillComplete();

#ifdef DEBUG_PROBLEM_MANAGER
    A->Print(cout);
#endif
  }

  return true;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_CrsGraph>
Problem_Manager::generateBlockDiagonalGraph( const Epetra_Map & rowMap, bool includeProblemOffContribs )
{
  Teuchos::RCP<Epetra_CrsGraph> graphPtr =
    Teuchos::rcp( new Epetra_CrsGraph(Copy, rowMap, 0) );

  // First construct a graph for each problem's self-dependence
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Loop over each problem being managed and ascertain its graph as well
  // as its graph from its dependencies
  for( ; problemIter != problemLast; ++problemIter )
  {
    GenericEpetraProblem & problem = *(*problemIter).second;
    int probId = problem.getId();

    Epetra_CrsGraph & problemGraph = *problem.getGraph();

    // Use max potential number of nonzero columns to dimension index array
    int problemMaxNodes = problemGraph.Map().MaxAllGID();

    // Get the indices map for copying data from this problem into
    // the composite problem
    Epetra_IntVector& problemIndices = *(ProblemToCompositeIndices[probId]);

    int problemRow, compositeRow, numCols;

    // First fill composite graph for each problem's self-dependence.  This
    // corresponds to diagonal blocks.
    int* columnIndices = new int[problemMaxNodes];

    for (int i = 0; i < problemGraph.NumMyRows(); ++i)
    {
      problemRow = problemGraph.Map().GID(i);
      problemGraph.ExtractGlobalRowCopy(problemRow, problemMaxNodes,
                           numCols, columnIndices);

      // Convert row/column indices to composite problem
      compositeRow = problemIndices[problemRow];
      for (int j = 0; j<numCols; j++)
        columnIndices[j] = problemIndices[columnIndices[j]];
      int ierr = graphPtr->InsertGlobalIndices(compositeRow, numCols, columnIndices);
      if( ierr )
      {
        if (MyPID==0)
          std::cout << "ERROR: graphPtr->InsertGlobalIndices(...)" << std::endl;
        throw "Problem_Manager ERROR";
      }
    }
    delete [] columnIndices; columnIndices = 0;

    // Now allow the problem to explicitly create known dependencies on other
    // problems
    if( includeProblemOffContribs )
    {
      for( unsigned int k = 0; k < problem.depProblems.size(); ++k)
      {
        // Get the needed objects for the depend problem
        GenericEpetraProblem & dependProblem = *(Problems[problem.depProblems[k]]);
        int dependId                         =  dependProblem.getId();
        Epetra_IntVector & dependIndices     = *(ProblemToCompositeIndices[dependId]);

        if( problem.suppliesOffBlockStructure(k) )
        {
          map<int, std::vector<int> > offBlockIndices;
          problem.getOffBlockIndices( offBlockIndices );

          map<int, std::vector<int> >::iterator indIter     = offBlockIndices.begin(),
                                           indIter_end = offBlockIndices.end()   ;

          for( ; indIter != indIter_end; ++indIter )
          {
            int compositeRow = problemIndices[(*indIter).first];
            std::vector<int> & colIndices = (*indIter).second;

            // Convert column indices to composite values
            for( unsigned int cols = 0; cols < colIndices.size(); ++cols )
              colIndices[cols] = dependIndices[ colIndices[cols] ];

            graphPtr->InsertGlobalIndices( compositeRow, colIndices.size(), &colIndices[0] );
          }
        }
      }
    }
  }

  return graphPtr;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_CrsGraph>
Problem_Manager::generateOffDiagonalBlockGraph( const Epetra_Map & rowMap )
{
  Teuchos::RCP<Epetra_CrsGraph> graphPtr =
    Teuchos::rcp( new Epetra_CrsGraph(Copy, rowMap, 0) );

  // Two things are achieved here: 1) The composite Graph is augmented to
  // accommodate off-diagonal blocks and 2) these blocks are packaged as
  // individual NOX::EpetreNew::Group's owned by the manager.

#ifdef HAVE_NOX_EPETRAEXT

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  problemIter = Problems.begin();

  // Loop over each problem being managed and ascertain its graph as well
  // as its graph from its dependencies
  for( ; problemIter != problemLast; problemIter++)
  {
    GenericEpetraProblem & problem = *(*problemIter).second;
    int probId = problem.getId();

    Epetra_CrsGraph & problemGraph = *problem.getGraph();

    // Get the indices map for copying data from this problem into
    // the composite problem
    Epetra_IntVector& problemIndices = *(ProblemToCompositeIndices[probId]);

    // Create containers for the off-block objects
    std::vector<OffBlock_Manager*> OffBlock_ManagersVec;

    int problemMaxNodes = problemGraph.Map().NumGlobalElements();

    int problemRow, compositeRow, numCols, numDepCols;

    // Loop over each problem on which this one depends
    for( unsigned int k = 0; k < problem.depProblems.size(); ++k)
    {
      // Create the off-block graph to be constructed for this
      // problem-problem coupling
      // NOTE: the map used for the off-block graph is the composite Map
      // to allow valid global indexing
      Epetra_CrsGraph* offGraphPtr = new Epetra_CrsGraph(Copy, *compositeMap, 0);

      Epetra_CrsGraph &offGraph = *offGraphPtr;

      // Get the needed objects for the depend problem
      GenericEpetraProblem & dependProblem = *(Problems[problem.depProblems[k]]);

      int dependId = dependProblem.getId();

      XferOp * xferOpPtr = problem.xferOperators[dependId];

      if( !xferOpPtr )
      {
        std::cout << "ERROR: Unable to get Xfer_Operator for dependence of "
             << "problem \"" << problem.getName() << "\" on problem "
             << "\"" << dependProblem.getName() << "\" !!" << std::endl;
        throw "Problem_Manager ERROR";
      }
      XferOp &xferOp = *xferOpPtr;
      multimap<int,int>& depNodesMap = xferOp.getDependentNodesMap();

      // Get the indices map for copying data from the dependent problem into
      // the composite problem
      Epetra_IntVector& dependIndices = *(ProblemToCompositeIndices[dependId]);

      // Dimension nonzero columns index array with upper bound which is
      // the previous definition * 2 since each node in problem could
      // depend at most on 2 nodes in dependProblem
      int* columnIndices = new int[problemMaxNodes];
      int maxDepNodes = 2 * problemGraph.Map().MaxAllGID();
      int* dependentColIndices = new int[maxDepNodes];

      // We must loop over each dependent node of problem and then determine
      // the dependence of each on the nodes of dependProblem

      // Loop over each row in problem and ascertain all dependencies on
      // dependProblem as determined by the xferOp map
      for (int i = 0; i<problemGraph.NumMyRows(); i++)
      {
        problemRow = problemGraph.Map().GID(i);

        problemGraph.ExtractGlobalRowCopy(problemRow, problemMaxNodes,
                             numCols, columnIndices);

        // Convert row/column indices to composite problem
        compositeRow = problemIndices[problemRow];
        numDepCols = 0;
        for (int j = 0; j<numCols; j++)
        {
          pair< multimap<int, int>::iterator,
                multimap<int, int>::iterator > rangeN
              = depNodesMap.equal_range(columnIndices[j]);
          multimap<int, int>::iterator iterN;
          for( iterN = rangeN.first; iterN != rangeN.second; iterN++)
            dependentColIndices[numDepCols++] = dependIndices[(*iterN).second];
        }
        graphPtr->InsertGlobalIndices(compositeRow, numDepCols, dependentColIndices);
        offGraph.InsertGlobalIndices(compositeRow, numDepCols, dependentColIndices);
      }
      delete [] columnIndices; columnIndices = 0;
      delete [] dependentColIndices; dependentColIndices = 0;

      offGraph.FillComplete();
#ifdef DEBUG_PROBLEM_MANAGER
      offGraph.Print(cout);
#endif
      // A new graph is created within the OffBlock_Manager; so we delete our temporary
      OffBlock_ManagersVec.push_back( new OffBlock_Manager(*this, offGraph,
                                                      probId, dependId) );
      delete offGraphPtr;
    }

    OffBlock_Managers[probId] = OffBlock_ManagersVec;
  }
#endif

  return graphPtr;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_CrsGraph>
Problem_Manager::composeGraphs( const Epetra_CrsGraph & graph1, const Epetra_CrsGraph & graph2 )
{
  Teuchos::RCP<Epetra_CrsGraph> graphPtr =
    Teuchos::rcp( new Epetra_CrsGraph(Copy, graph1.RowMap(), 0) );

  if( !graph1.Filled() || !graph2.Filled() )
  {
    std::cout << "ERROR: Problem_Manager::composeGraphs(...) One or both incoming graphs have not "
         << "been transformed to local index space by calling FillComplete()." << std::endl;
    throw "Problem_Manager ERROR";
  }

  if( !graph1.RowMap().SameAs( graph2.RowMap() ) )
  {
    std::cout << "ERROR: Problem_Manager::composeGraphs(...) Incoming graphs do not have "
         << "the same row maps." << std::endl;
    throw "Problem_Manager ERROR";
  }

  int * indices1   ;
  int * indices2   ;
  int * newIndices ;
  int   numInd1    ;
  int   numInd2    ;
  int   numNewInds ;
  list<int> list1  ;
  list<int> list2  ;


  for( int row = 0; row < graph1.NumMyRows(); ++row )
  {
    list1.clear();
    list2.clear();

    graph1.ExtractMyRowView( row, numInd1, indices1 );
    graph2.ExtractMyRowView( row, numInd2, indices2 );

    for( int col = 0; col < numInd1; ++col )
      list1.push_back(graph1.ColMap().GID(indices1[col]));

    for( int col = 0; col < numInd2; ++col )
      list2.push_back(graph2.ColMap().GID(indices2[col]));

    list1.sort();
    list2.sort();

    list1.merge(list2);

    numNewInds = list1.size();
    newIndices = new int[ numNewInds ];

    list<int>::iterator iter = list1.begin();
    for( int col = 0; col < numNewInds; )
      newIndices[col++] = *(iter++);

    graphPtr->InsertGlobalIndices( graph1.RowMap().GID(row), numNewInds, newIndices );

    delete [] newIndices;
  }

  return graphPtr;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::generateGraph()
{

  // First construct a graph for each problem's self-dependence
  Teuchos::RCP<Epetra_CrsGraph> diagBlkGraph = generateBlockDiagonalGraph( *compositeMap, false );
  diagBlkGraph->FillComplete();

  // Next create inter-problem block graph contributions if desired;
  //   default is false
  if( doOffBlocks )
  {
    Teuchos::RCP<Epetra_CrsGraph> offBlkGraph = generateOffDiagonalBlockGraph( *compositeMap );
    offBlkGraph->FillComplete();

    AA = composeGraphs( *diagBlkGraph, *offBlkGraph );
  }
  else
  {
    AA = diagBlkGraph;
  }

  AA->FillComplete();

  return;
}

//-----------------------------------------------------------------------------

string
Problem_Manager::createIOname(GenericEpetraProblem & problem, int timeStep)
{

  int probId = problem.getId();

  char file_name[25];
  (void) sprintf(file_name, "output.%03d.%03d_%05d", probId, MyPID, timeStep);

  std::string name(file_name);

  return name;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::outputSolutions( const std::string outputDir, int timeStep )
{

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = Problems.end();

  // Loop over each problem being managed and write its solution vector
  // to a file.
  for( ; problemIter != problemLast; problemIter++)
  {
    GenericEpetraProblem & problem = *(*problemIter).second;

    Epetra_Vector& xMesh = problem.getMesh();
    Epetra_Vector& problemSoln = *problem.getSolution();

    std::string baseName = createIOname( problem, timeStep );
    std::string fileName = outputDir + baseName;

#ifdef HAVE_NOX_EPETRAEXT
    NOX::Epetra::DebugTools::writeVector( fileName, problemSoln, NOX::Epetra::DebugTools::MATRIX_MARKET );
    fileName = outputDir + baseName + "_mesh";
    NOX::Epetra::DebugTools::writeVector( fileName, xMesh, NOX::Epetra::DebugTools::MATRIX_MARKET );
#else
    int numNodes = xMesh.Map().NumMyElements();
    FILE *ifp;
    ifp = fopen(fileName.c_str(), "w");
    for (int i = 0; i < numNodes; i++)
      fprintf(ifp, "%d  %E  %E \n", i, xMesh[i], problemSoln[i]);
    fclose(ifp);
#endif
  }

  return;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::outputStatus( std::ostream & os )
{

  map<int, Teuchos::RCP<GenericEpetraProblem> >::const_iterator problemIter = Problems.begin();
  map<int, Teuchos::RCP<GenericEpetraProblem> >::const_iterator problemLast = Problems.end();

  os << std::endl << std::endl << "\t\t********************************"   << std::endl;
  os                 << "\t\t*******  PROBLEM SUMMARY  ******"   << std::endl;
  os                 << "\t\t********************************"   << std::endl;
  os << std::endl << std::endl << "\t\tOff-diagonal contributions are "
                       << ( doOffBlocks ? "Enabled" : "Disabled" ) << std::endl;
  os << std::endl << std::endl;

  // Loop over each problem being managed and output its dependencies
  for( ; problemIter != problemLast; problemIter++)
  {
    GenericEpetraProblem & problem = *(*problemIter).second;

    os << "\tProblem \"" << problem.getName() << "\" (" << (*problemIter).first
         << ")\t Depends on:" << std::endl;

    for( unsigned int j = 0; j < problem.depProblems.size(); ++j )
    {
      GenericEpetraProblem & depProblem = *(Problems[ problem.depProblems[j] ]);

      os << "\t\t-------------> \t\t\"" << depProblem.getName()
           << "\" (" << depProblem.getId() << ")" << std::endl;
    }
    os << std::endl;

    // Allow problems to provide additional info if desired
    for( unsigned int j = 0; j < problem.depProblems.size(); ++j )
    {
      //GenericEpetraProblem & depProblem = *(Problems[ problem.depProblems[j] ]);
      problem.outputStatus(os);
    }
    os << std::endl;
  }

  return;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::exchangeAllData()
{
  std::cout << "Problem_Manager::exchangeAllData() .... called." << std::endl;

  // Preceding this call, the solution vector for each problem has been placed
  // in each solver's group X vector.  We need to copy this into each problem's
  // corresponding solution vector and then fire off the appropriate transfers
  // into auxiliary data vectors.

  copyAllGroupXtoProblems();

  syncAllProblems();

  return true;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::exchangeDataTo(int solverId)
{
  std::cout << "Problem_Manager::exchangeDataTo( " << solverId << " ) .... called." << std::endl;

  // Note: the incoming solverId reflects the order in which this problem occurs
  // in the container of solvers passed to the coupling solver and is 0-based.
  // We need to convert this to 1-based nd account for any permutations of problem
  // ordering.

  solverId++;

  // Preceding this call, the solution vector for each problem has been placed
  // in each solver's group X vector.  We need to copy this into each problem's
  // corresponding solution vector and then fire off the appropriate transfers
  // into auxiliary data vectors.

  copyAllGroupXtoProblems();

  Problems[solverId]->doTransfer();

  return false;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_CrsMatrix>
Problem_Manager::getReplacementValuesMatrix( const Epetra_Vector & x, FILL_TYPE )
{
  std::cout << "Problem_Manager::getReplacementValuesMatrix(...) called." << std::endl;

  bool ok = evaluate( NOX::Epetra::Interface::Required::Jac, &x, NULL );

  if( !ok )
  {
    std::cout << "ERROR: Problem_Manager::getReplacementValuesMatrix call to evaluate failed."
         << std::endl;

    throw "Problem_Manager ERROR";
  }

  return A;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::applyBlockAction( int probId, int depId, const Epetra_Vector & x, Epetra_Vector & result )
{
  // We currently rely on an existing MatrixFree jacobian operator for this action
  Teuchos::RCP<NOX::Epetra::MatrixFree> mfOp = Teuchos::rcp_dynamic_cast<NOX::Epetra::MatrixFree>( jacOperator );
  if( Teuchos::is_null(mfOp) )
  {
    std::cout << "ERROR: Problem_Manager::applyBlockAction : jacOperator is not of type NOX::Epetra::MatrixFree." << std::endl;
    throw "Problem_Manager ERROR";
  }

  if( !x     .Map().SameAs(getSolutionVec(depId).Map())  ||
      !result.Map().SameAs(getSolutionVec(probId).Map())   )
  {
    std::cout << "ERROR: Problem_Manager::applyBlockAction : vector sizes not compatible with domain and/or range of block dimensions." << std::endl;
    throw "Problem_Manager ERROR";
  }
  Teuchos::RCP<Epetra_Vector> tmpXcomposite = Teuchos::rcp( new Epetra_Vector(*compositeSoln) );
  Teuchos::RCP<Epetra_Vector> tmpYcomposite = Teuchos::rcp( new Epetra_Vector(*compositeSoln) );

  tmpXcomposite->PutScalar(0.0);
  copyVectorToComposite( *tmpXcomposite, depId, x );

  mfOp->Apply( *tmpXcomposite, *tmpYcomposite );

  copyCompositeToVector( *tmpYcomposite, probId, result );

  return true;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::applyBlockInverseAction( int probId, int depId, const Epetra_Vector & x, Epetra_Vector & result )
{
  // We simply invoke the action of a previously created inverse
  getBlockInverseOperator(probId)->ApplyInverse(x, result);

  return true;
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::hasExplicitOperator( int rowBlock, int colBlock )
{
  return ( rowBlock == colBlock );
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Epetra_Operator>
Problem_Manager::getExplicitOperator( int rowBlock, int colBlock )
{
  return getBlockJacobianMatrix(rowBlock);
}

//-----------------------------------------------------------------------------

bool
Problem_Manager::createBlockInverseOperator( int probId, Teuchos::ParameterList & pList )
{
  // For now, we will hard-code use of Ifpack.  Note that direct inverses via Amesos can be
  // used by propoer specification in the Teuchos::ParameterList.

  Teuchos::RCP<Ifpack_Preconditioner> inverseOperator;

  // Get any relevant parameters for the Inverting package using sublists
  Teuchos::ParameterList& teuchosParams = pList.sublist("Ifpack");

  Ifpack Factory;

  Teuchos::RCP<Epetra_CrsMatrix> pMatrix = getBlockJacobianMatrix(probId);

  if( NULL == pMatrix.get() )
  {
    std::cout << "ERROR: Problem_Manager::createBlockInverseOperator : Cannot get an Epetra_CrsMatrix for diagonal block # " << probId << std::endl;
    throw "Problem_Manager ERROR";
  }

  inverseOperator = Teuchos::rcp(Factory.Create( pList.get("Ifpack Preconditioner", "ILU"),
    &(*pMatrix),
    pList.get("Overlap", 0) ));
  inverseOperator->SetParameters(teuchosParams);
  inverseOperator->Initialize();
  inverseOperator->Compute();

  BlockInverseOperators[probId] = inverseOperator;

  return true;
}

//-----------------------------------------------------------------------------
