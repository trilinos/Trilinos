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
                                                                                
#include "NOX.H"
#include "NOX_Epetra.H"
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

#include "Problem_Manager.H"
#include "GenericEpetraProblem.H"
#include "Xfer_Operator.H"
#include "OffBlock_Manager.H"

// Headers needed for Coloring
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif

// Header for Timing info
#include "Epetra_Time.h"

// Hard-wired switch to turn on/off use of FD Coloring
#define USE_FD
//#undef USE_FD

Problem_Manager::Problem_Manager(Epetra_Comm& comm, 
                                 bool doOffBlocks_,
                                 int numGlobalElements) :
  GenericEpetraProblem(comm, numGlobalElements),
  problemCount(0),
  doOffBlocks(doOffBlocks_),
  compositeMap(0),
  compositeSoln(0),
  nlParams(0),
  statusTest(0)
{
  // Unset doOffBlocks flag if this build does not include the required 
  // EpetraExt library intreface
#ifndef HAVE_NOX_EPETRAEXT
  doOffBlocks = false;
#endif

  // Create a problem interface to the manager
  compositeProblemInterface = new Problem_Interface(*this);

  // Reset composite number of dofs (overwrites base constructor assignment)
  NumMyNodes = 0;

  // Do all setup in registerComplete after all registrations have been
  // performed
}

Problem_Manager::~Problem_Manager()
{
  delete AA; AA = 0;
  delete A; A = 0;

  // Iterate over each problem and destroy/free the necessary objects

  map<int, GenericEpetraProblem*>::iterator iter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator last = Problems.end();

  map<int, NOX::EpetraNew::Group*>::iterator GroupsIter = Groups.begin();   
  map<int, Problem_Interface*>::iterator InterfacesIter = Interfaces.begin();
  map<int, NOX::Solver::Manager*>::iterator SolversIter = Solvers.begin();

#ifdef HAVE_NOX_EPETRAEXT
#ifdef USE_FD
  map<int, EpetraExt::CrsGraph_MapColoring*>::iterator 
	  TmpMapColoringsIter = TmpMapColorings.begin();
  map<int, Epetra_MapColoring*>::iterator 
	  ColorMapsIter = ColorMaps.begin();
  map<int, EpetraExt::CrsGraph_MapColoringIndex*>::iterator 
	  ColorMapIndexSetsIter = ColorMapIndexSets.begin();
  map<int, vector<Epetra_IntVector>*>::iterator 
	  ColumnsSetsIter = ColumnsSets.begin();
  map<int, Epetra_Operator*>::iterator 
	  MatrixOperatorsIter = MatrixOperators.begin();
#endif
#endif

  while( iter != last)
  {
    delete (SolversIter++)->second;
    delete (GroupsIter++)->second;
#ifdef HAVE_NOX_EPETRAEXT
#ifdef USE_FD
    delete (MatrixOperatorsIter++)->second;
    delete (TmpMapColoringsIter++)->second;
    delete (ColorMapsIter++)->second;
    delete (ColorMapIndexSetsIter++)->second;
    //delete *ColumnsSetsIter++;
#endif
#endif
    iter++; // Problems are owned by the app driver (Example.C)
  }
}

void Problem_Manager::addProblem(GenericEpetraProblem& problem)
{
  Problems.insert(pair<int, GenericEpetraProblem*>(++problemCount, &problem));
  problem.setId(problemCount);
  problem.setManager(this);

  // Give this problem a name if it doesn't already have one
  if( problem.getName() == "" ) {
    string name = "Problem_";
    char id_str[4];
    (void) sprintf(id_str, "%d",problemCount);
    name += id_str;
    problem.setName(name);
  }

  Names.insert( pair<int, string> (problemCount, problem.getName()) );
  NameLookup.insert( pair<string, int> (problem.getName(), problemCount) );

  // Keep a running total of dofs for use in constructing composite objects
  NumMyNodes += problem.NumMyNodes;

}

GenericEpetraProblem& Problem_Manager::getProblem(int id_)
{
  // Get a problem given its unique id
  GenericEpetraProblem* problem = Problems.find(id_)->second;
  if( !problem ) {
    cout << "ERROR: Problem with id --> " << id_ << " not registered with "
         << "Problem_Manager !!" << endl;
    outputStatus();
    throw "Problem_Manager ERROR";
  }
  else
    return *problem;
}

GenericEpetraProblem& Problem_Manager::getProblem(string name)
{
  // Get a problem given its name
  map<string, int>::iterator iter = NameLookup.find(name);
  if( iter == NameLookup.end() ) {
    cout << "ERROR: Could not find lookup id for Problem --> " << name 
         << endl;
    outputStatus();
    throw "Problem_Manager ERROR";
  }
  else
    return getProblem( iter->second );
}

NOX::EpetraNew::Group& Problem_Manager::getGroup(int id_)
{
  // Get a group given its unique id
  NOX::EpetraNew::Group* group = Groups.find(id_)->second;
  if( !group ) {
    cout << "ERROR: Could not get Group for Problem with id --> " << id_ 
         << " !!" << endl;
    throw "Problem_Manager ERROR";
  }
  else
    return *group;
}

Epetra_Vector& Problem_Manager::getCompositeSoln()
{
  if( !compositeSoln ) {
    cout << "ERROR: No valid Composite Solution vector with Problem Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }
  return *compositeSoln;
}

void Problem_Manager::createDependency(string nameA, string nameB)
{
  // Create a dependence of Problem A equations on Problem B variables
  int probId_A = NameLookup.find(nameA)->second;
  int probId_B = NameLookup.find(nameB)->second;
  if( !probId_A || !probId_B ) {
    cout << "ERROR: Could not create dependency of \"" << nameA << "\" on \""
         << nameB << "\" !!" << endl;
    throw "Problem_Manager ERROR";
  }
 
  GenericEpetraProblem &probA = *Problems.find(probId_A)->second,
                       &probB = *Problems.find(probId_B)->second;

  createDependency(probA, probB);
}


void Problem_Manager::createDependency(GenericEpetraProblem& problemA,
                                       GenericEpetraProblem& problemB)
{
  // Ensure that both problems already exist
  if( !problemA.getId() || !problemB.getId() ) {
    cout << "ERROR: Invalid dependency since at least one problem is "
         << "not registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  problemA.addTransferOp(problemB);
  problemA.addProblemDependence(problemB);

}

void Problem_Manager::registerParameters(NOX::Parameter::List& List)
{
  nlParams = &List;
}

void Problem_Manager::registerStatusTest(NOX::StatusTest::Combo& comboTest)
{
  statusTest = &comboTest;
}

void Problem_Manager::registerComplete()
{
  if(Problems.empty())
  {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  if(nlParams == 0 || statusTest == 0)
  {
    cout << "ERROR: No nlParams and/or statusTest registered with "
         << "Problem_Manager !!" << endl;
    throw "Problem_Manager ERROR";
  }

  // Iterate over each problem and construct the necessary objects

  map<int, GenericEpetraProblem*>::iterator iter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator last = Problems.end();

  // Make sure everything is starting clean
  assert(Groups.empty() && Interfaces.empty() && Solvers.empty()
         && ProblemToCompositeIndices.empty());

  int icount = 0; // Problem counter
  int runningProblemNodeCount = 0;
  int runningMaxGlobalId = 0; // Accruing max index used for establishing
                  // each problem's mapping into the composite problem
  //
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

    // Create index mapping for this problem into the composite problem
    ProblemToCompositeIndices.insert( pair<int, Epetra_IntVector*>
      (probId, new Epetra_IntVector(*problem.StandardMap)));
    Epetra_IntVector &indices = *ProblemToCompositeIndices.find(probId)->second;

    for (int i=0; i<problem.NumMyNodes; i++) {
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

  while( iter != last)
  {
    // Get a convenient reference to the current problem
    GenericEpetraProblem& problem = *(*iter).second;
    int probId = problem.getId();

    Interfaces.insert( pair<int, Problem_Interface*>
		       (probId, new Problem_Interface(problem)));
    NOX::EpetraNew::Interface::Required& reqInt = 
      dynamic_cast<NOX::EpetraNew::Interface::Required&>
         (*Interfaces.find(probId)->second);

    NOX::Epetra::Vector nox_soln( problem.getSolution() );

    // Use this for analytic Matrix Fills
#ifndef USE_FD
    NOX::EpetraNew::Interface::Jacobian& jacInt = 
      dynamic_cast<NOX::EpetraNew::Interface::Jacobian&>
         (*Interfaces.find(probId)->second);
    LinearSystems.insert( pair<int, NOX::EpetraNew::LinearSystemAztecOO*>
     ( probId, new NOX::EpetraNew::LinearSystemAztecOO(
      nlParams->sublist("Printing"),
      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      reqInt,
      jacInt,
      problem.getJacobian(),
      problem.getSolution())) );

    Groups.insert( pair<int, NOX::EpetraNew::Group*>
     ( probId, new NOX::EpetraNew::Group(
      nlParams->sublist("Printing"),
      reqInt,
      nox_soln,
      *LinearSystems.find(probId)->second)) );

    // OR use this to fill matrices using Finite-Differences with Coloring
#else
#ifdef HAVE_NOX_EPETRAEXT
    // Create a timer for performance
    Epetra_Time fillTime(*Comm);

    bool verbose = false;
    EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
      EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY;
    TmpMapColorings.insert( pair<int, EpetraExt::CrsGraph_MapColoring*>
      (probId, new EpetraExt::CrsGraph_MapColoring(algType, verbose)));
    ColorMaps.insert( pair<int, Epetra_MapColoring*>
      (probId, &((*TmpMapColorings.find(probId)->second)(problem.getGraph()))) );
    ColorMapIndexSets.insert( pair<int, EpetraExt::CrsGraph_MapColoringIndex*>
      (probId, new EpetraExt::CrsGraph_MapColoringIndex(*ColorMaps.find(probId)->second)) );
    ColumnsSets.insert( pair<int, vector<Epetra_IntVector>* >
      (probId, &(*ColorMapIndexSets.find(probId)->second)(problem.getGraph())) );

    if (MyPID == 0)
      printf("\n\tTime to color Jacobian # %d --> %e sec. \n\n",
                  icount++,fillTime.ElapsedTime());
    MatrixOperators.insert( pair<int, NOX::EpetraNew::FiniteDifferenceColoring*>
      (probId, new NOX::EpetraNew::FiniteDifferenceColoring(*Interfaces.find(probId)->second, 
        problem.getSolution(), problem.getGraph(), *ColorMaps.find(probId)->second, 
        *ColumnsSets.find(probId)->second)) );
    NOX::EpetraNew::Interface::Jacobian& jacInt = 
      dynamic_cast<NOX::EpetraNew::Interface::Jacobian&>(*MatrixOperators.find(probId)->second);
    LinearSystems.insert( pair<int, NOX::EpetraNew::LinearSystemAztecOO*>
      (probId, new NOX::EpetraNew::LinearSystemAztecOO(
      nlParams->sublist("Printing"),
      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      reqInt,
      jacInt,
      *MatrixOperators.find(probId)->second,
      problem.getSolution())) );

    Groups.insert( pair<int, NOX::EpetraNew::Group*>
      (probId, new NOX::EpetraNew::Group(
      nlParams->sublist("Printing"),
      reqInt,
      nox_soln,
      *LinearSystems.find(probId)->second)) );
#else
    if(MyPID==0)
      cout << "ERROR: Cannot use EpetraExt with this build !!" << endl;
    exit(0);
#endif
#endif

    // Needed to establish initial convergence state
    Groups.find(probId)->second->computeF(); 
   
    Solvers.insert( pair<int, NOX::Solver::Manager*>
      (probId, new NOX::Solver::Manager(*Groups.find(probId)->second, 
					*statusTest, *nlParams)) );
    iter++;
  }


  compositeMap = new Epetra_Map(-1, NumMyNodes, compositeGlobalNodes, 0, *Comm);
  delete [] compositeGlobalNodes; compositeGlobalNodes = 0;

  compositeSoln = new Epetra_Vector(*compositeMap);

  return;

}

void Problem_Manager::syncAllProblems()
{
  if(Problems.empty()) {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }
  
  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Loop over each problem being managed and invoke its transfer requests
  for( ; problemIter != problemLast; problemIter++)
    (*problemIter).second->doTransfer();
}

void Problem_Manager::setGroupX(int probId)
{
  GenericEpetraProblem *problem = Problems.find(probId)->second;
  if( !problem ) {
    cout << "ERROR: Could not get requested Problem to use with group.setX "
         << endl;
    throw "Problem_Manager ERROR";
  }

  NOX::EpetraNew::Group *grp = Groups.find(probId)->second;
  if( !grp ) {
    cout << "ERROR: Could not get appropriate group for use in setX !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  grp->setX(problem->getSolution());
}

void Problem_Manager::setAllGroupX()
{
  if(Problems.empty()) {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Loop over each problem being managed and set the corresponding group
  // solution vector (used by NOX) with the problem's (used by application)
  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;

    setGroupX(probId);
  }
}

#ifdef HAVE_NOX_EPETRAEXT
void Problem_Manager::setAllOffBlockGroupX(const Epetra_Vector &inVec)
{
  map<int, vector<OffBlock_Manager*> >::iterator offBlockIter = 
                                                   OffBlock_Managers.begin();
  map<int, vector<OffBlock_Manager*> >::iterator offBlockLast = 
                                                   OffBlock_Managers.end();

  // Loop over each off-block manager and set the contained groups X-vector
  // with the incoming vector
  for( ; offBlockIter != offBlockLast; offBlockIter++) {

    vector<OffBlock_Manager*> managerVec = offBlockIter->second;

    for( int i = 0; i<managerVec.size(); i++ )
      managerVec[i]->getGroup().setX(inVec);
  }
}
#endif

void Problem_Manager::resetProblems()
{ 

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Loop over each problem and copy its solution into its old solution
  for( ; problemIter != problemLast; problemIter++) {
    GenericEpetraProblem& problem = *(*problemIter).second;
    problem.reset( problem.getSolution() );
  }
}

void Problem_Manager::computeAllF()
{
  if(Problems.empty()) {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Loop over each problem being managed and invoke the corresponding group's
  // residual evaluation
  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;

    computeGroupF(probId);
  }
}

void Problem_Manager::computeGroupF(int probId)
{
  NOX::EpetraNew::Group *grp = Groups.find(probId)->second;
  if( !grp ) {
    cout << "ERROR: Could not get a group for problem with id --> "
         << probId << endl;
    throw "Problem_Manager ERROR";
  }
  grp->computeF();
}

void Problem_Manager::computeAllJacobian()
{
  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Create a timer for performance
  Epetra_Time fillTime(*Comm);

  // Loop over each problem being managed and invoke its computeJacobian
  // method
  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;
    NOX::EpetraNew::Group *grp = Groups.find(probId)->second;
    if( !grp ) {
      cout << "ERROR: Could not find valid group for compouteJacobian !!"
           << endl;
      throw "Problem_Manager ERROR";
    }

    fillTime.ResetStartTime();

    grp->computeJacobian();

    if (MyPID == 0)
      printf("\n\tTime to fill Jacobian %d --> %e sec. \n\n",
                  probId, fillTime.ElapsedTime());


    if( doOffBlocks ) {
#ifdef HAVE_NOX_EPETRAEXT

      fillTime.ResetStartTime();
  
      vector<OffBlock_Manager*> &offBlocksVec =
        OffBlock_Managers.find(probId)->second;

      for( int i = 0; i<offBlocksVec.size(); i++ ) {
        
        offBlocksVec[i]->getGroup().computeJacobian();
  
        if (MyPID == 0)
          printf("\n\tTime to fill Jacobian %d (%d) --> %e sec. \n\n",
                      probId, i, fillTime.ElapsedTime());
      }
#endif
    }
  }
}

void Problem_Manager::updateWithFinalSolution(int probId)
{
  // Copy final solution from NOX solver into the problem's solution vector

  GenericEpetraProblem *problem = Problems.find(probId)->second;
  if( !problem ) {
    cout << "ERROR: Could not get requested Problem to update with final "
         << "solution" << endl;
    throw "Problem_Manager ERROR";
  }

  NOX::Solver::Manager* solver = Solvers.find(probId)->second;
  if( !solver ) {
    cout << "ERROR: Could not get appropriate Solver for use in update !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  const NOX::EpetraNew::Group& finalGroup =
    dynamic_cast<const NOX::EpetraNew::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>
      (finalGroup.getX())).getEpetraVector();
  
  problem->getSolution() = finalSolution;
}

void Problem_Manager::updateAllWithFinalSolution()
{
  // Copy final solution from NOX solvers into each problem's solution vector

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;

    updateWithFinalSolution(probId);
  }
}

void Problem_Manager::copyCompositeToProblems(
                         const Epetra_Vector& compositeVec, 
                         Problem_Manager::VectorType vecType)
{
  // Copy a composite problem vector to each problem's vector

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  Epetra_Vector* problemVec(0);

  // Loop over each problem being managed and copy into the correct 
  // problem vector
  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;
    switch (vecType) {

      case SOLUTION :
        problemVec = &(problemIter->second->getSolution());
	break;

      case GROUP_F :
      default :
        cout << "ERROR: vecType not supported for copy FROM composite!!" 
             << endl;
        throw "Problem_Manager ERROR";
    }

    copyCompositeToVector(compositeVec, probId, *problemVec); 
  }
}

void Problem_Manager::copyProblemsToComposite(
                         Epetra_Vector& compositeVec,
                         Problem_Manager::VectorType vecType)
{
  // Copy vectors from each problem into a composite problem vector

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  const Epetra_Vector* problemVec(0);

  // Loop over each problem being managed and copy from the correct 
  // problem vector
  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;
    switch (vecType) {

      case SOLUTION :
        problemVec = &(problemIter->second->getSolution());
	break;

      case GROUP_F :
        problemVec = &(dynamic_cast<const NOX::Epetra::Vector&>
                     (Groups.find(probId)->second->getF()).getEpetraVector());
	break;

      default :
        cout << "ERROR: vecType not supported for copy TO composite!!" 
             << endl;
        throw "Problem_Manager ERROR";
    }

    copyVectorToComposite(compositeVec, probId, *problemVec); 
  }
}


void Problem_Manager::copyCompositeToVector(
                      const Epetra_Vector& compositeVec, int id,
                      Epetra_Vector& problemVec)
{
  // Copy part of a composite problem vector to a problem's vector
  Epetra_IntVector& indices = *ProblemToCompositeIndices.find(id)->second;
  // This map is needed to get correct LID of indices
  const Epetra_BlockMap &map = compositeVec.Map();

  for (int i=0; i<problemVec.MyLength(); i++)
    problemVec[i] = compositeVec[map.LID(indices[i])];
}

void Problem_Manager::copyVectorToComposite(
                      Epetra_Vector& compositeVec, int id,
                      const Epetra_Vector& problemVec)
{
  // Copy a vector from a problem into part of a composite problem vector
  Epetra_IntVector& indices = *ProblemToCompositeIndices.find(id)->second;
  for (int i=0; i<problemVec.MyLength(); i++)
    compositeVec[indices[i]] = problemVec[i];
}

void Problem_Manager::copyProblemJacobiansToComposite()
{
  // Copy problem Jacobians as block diagonal contributions to 
  // composite Jacobian

  Epetra_CrsMatrix &compositeMatrix = dynamic_cast<Epetra_CrsMatrix&>(*A);

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  int problemMaxNodes = compositeSoln->GlobalLength();

  // Loop over each problem being managed and copy its Jacobian into 
  // the composite diagonal blocks
  for( ; problemIter != problemLast; problemIter++) {

    int probId = problemIter->first;

    // Get the problem, its Jacobian graph and its linear system
    GenericEpetraProblem &problem = *(problemIter->second);
    Epetra_CrsGraph &problemGraph = problem.getGraph();
    NOX::EpetraNew::LinearSystemAztecOO &problemLinearSystem = 
            *LinearSystems.find(probId)->second;

    // Get the indices map for copying data from this problem into 
    // the composite problem
    Epetra_IntVector& indices = 
            *ProblemToCompositeIndices.find(probId)->second;


    // Each problem's Jacobian will be determined by type 
    Epetra_CrsMatrix *problemMatrixPtr(0);

    // Use each group's operator test to determine the type of Jacobian
    // operator being used.
    const Epetra_Operator& jacOp = problemLinearSystem.getJacobianOperator();

    if ( dynamic_cast<const NOX::EpetraNew::FiniteDifference*>(&jacOp) )
      problemMatrixPtr = const_cast<Epetra_CrsMatrix*>(
        &dynamic_cast<const NOX::EpetraNew::FiniteDifference&>
          (jacOp).getUnderlyingMatrix());
    else if ( dynamic_cast<const Epetra_CrsMatrix*>(&jacOp) )
      // NOTE: We are getting the matrix from the problem.  This SHOULD be
      // the same matrix wrapped in the group.  A safer alternative would be
      // to get this matrix from the group as above for a more general 
      // operator.
      problemMatrixPtr = &(problem.getJacobian());
    else
    {
      if (MyPID==0)
        cout << "Jacobian operator for Problem not supported for "
             << "preconditioning Matrix-Free coupling solver." << endl;
      throw "Problem_Manager ERROR";
    }

    // Create convenient reference for each Jacobian matrix
    Epetra_CrsMatrix &problemMatrix = *problemMatrixPtr;

    // Temporary storage arrays for extracting/inserting matrix row data
    int* columnIndices = new int[problemMaxNodes];
    double* values = new double[problemMaxNodes];

    int problemRow, compositeRow, numCols, numValues;
    for (int i = 0; i<problemMatrix.NumMyRows(); i++)
    { 
      problemRow = problemMatrix.Map().GID(i);
      problemGraph.ExtractGlobalRowCopy(problemRow, problemMaxNodes, 
                           numCols, columnIndices);
      problemMatrix.ExtractGlobalRowCopy(problemRow, problemMaxNodes,
                           numValues, values);
      if( numCols != numValues ) {
        if (MyPID==0)
          cout << "ERROR: Num Columns != Num Values from problem Matrix !!"
               << endl;
        throw "Problem_Manager ERROR";
      }
      // Convert row/column indices to composite problem
      compositeRow = indices[problemRow];
      for (int j = 0; j<numCols; j++)
        columnIndices[j] = indices[columnIndices[j]];
      int ierr = compositeMatrix.ReplaceGlobalValues(compositeRow, 
                       numValues, values, columnIndices);
    }
    delete [] values; values = 0;
    delete [] columnIndices; columnIndices = 0;

    // Sync up processors to be safe
    Comm->Barrier();

    // Add off-diagonal FDC block contributions if waranted
    if( doOffBlocks ) {
#ifdef HAVE_NOX_EPETRAEXT
      // Loop over each problem on which this one depends
      for( int k = 0; k<problem.depProblems.size(); k++) {

        // Copy the off-block jacobian matrices for this 
        // problem-problem coupling
        // NOTE: the map used for the off-block graph is the composite Map
        // to allow valid global indexing.  This also allows us to
        // simply copy values directly from the off-block matrices to the
        // composite Jacobian matrix
	Epetra_CrsMatrix *offMatrixPtr = 
          &( (OffBlock_Managers.find(probId)->second)[k]->getMatrix() );
	if( !offMatrixPtr ) {
          cout << "ERROR: Unable to get FDColoring underlying matrix for "
               << "dependence of problem " << probId << " on problem "
               << problem.depProblems[k] << " !!" << endl;
          throw "Problem_Manager ERROR";
        }
        Epetra_CrsMatrix &offMatrix = *offMatrixPtr;

        int* columnIndices = new int[problemMaxNodes];
        double* values = new double[problemMaxNodes];

        // Loop over each row and copy into composite matrix
        for (int i = 0; i<offMatrix.NumMyRows(); i++) {

          problemRow = offMatrix.Map().GID(i);

          offMatrix.ExtractGlobalRowCopy(problemRow, problemMaxNodes,
                           numValues, values, columnIndices);
          compositeMatrix.ReplaceGlobalValues(problemRow, 
                       numValues, values, columnIndices);
        }
        delete [] columnIndices; columnIndices = 0;
        delete [] values; values = 0;

        // Sync up processors to be safe
        Comm->Barrier();
      }
#endif
    }
  }
#ifdef DEBUG
    compositeMatrix.Print(cout);
#endif
}

double Problem_Manager::getNormSum()
{
  // Get each problem's residual norm and sum into a total
  double normSum(0.0);

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  for( ; problemIter != problemLast; problemIter++) {
    int probId = problemIter->first;

    NOX::EpetraNew::Group *grp = Groups.find(probId)->second;
    if( !grp ) {
      cout << "ERROR: Could not get appropriate group for use in NormSum !!"
           << endl;
      throw "Problem_Manager ERROR";
    }

    double problemNorm = grp->getNormF();
    cout << "2-Norm of Problem " << probId << " --> " << problemNorm
         << endl;
    normSum += problemNorm * problemNorm;
  }
  normSum = sqrt(normSum);

  return normSum;
}

bool Problem_Manager::solve()
{
  if(Problems.empty())
  {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  if( Groups.empty() || Interfaces.empty() || Solvers.empty() )
  {
    cout << "ERROR: Groups, Interfaces and/or Solvers are emptry !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Sync all the problems and get initial convergence state
  syncAllProblems();
  setAllGroupX();
  computeAllF();

  double normSum = getNormSum();
  cout << "Initial 2-Norm of composite Problem --> " << normSum;

  // Now do the decoupled solve
  int iter = 0;
  NOX::StatusTest::StatusType status;

//  while( normSum > 2.e-0 ) // Hard-coded convergence criterion for now.
  int nlIter = 0; 
  while( nlIter < 5 ) // Hard-coded convergence criterion for now.
  {
    iter++;

    problemIter = Problems.begin();

    // Solve each problem in the order it was registered
    for( ; problemIter != problemLast; problemIter++) {

      GenericEpetraProblem& problem = *(*problemIter).second;
      int probId = problem.getId();

      NOX::EpetraNew::Group &problemGroup = *Groups.find(probId)->second;
      NOX::Solver::Manager &problemSolver = *Solvers.find(probId)->second;
    
      // Sync all dependent data with this problem 
      problem.doTransfer();
      // Sync the problem solution with its solver group
      setGroupX(probId);
      // Reset the solver for this problem and solve
      problemSolver.reset(problemGroup, *statusTest, *nlParams);
      status = problemSolver.solve();
      if( status != NOX::StatusTest::Converged )
      { 
        if (MyPID==0)
          cout << "\nRegistered Problem ## failed to converge !!"  << endl;
        //exit(0);
      }

      updateWithFinalSolution(probId);
    }

    // Determine final residuals for use in testing convergence
    syncAllProblems();
    setAllGroupX();
    computeAllF();

    normSum = getNormSum();
    cout << "iter #" << iter << ", 2-Norm of composite Problem --> " 
         << normSum << endl;

    nlIter++;
  }
  
  cout << "\nDecoupled solution required --> " << iter << " iterations.\n" 
       << endl;

  return true;
}

bool Problem_Manager::solveMF()
{
  if(Problems.empty())
  {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  // Sync all the problems and get initial convergence state
  syncAllProblems();

  setAllGroupX();

  computeAllF();

  double normSum = getNormSum();
  cout << "Initial 2-Norm of composite Problem --> " << normSum;

  // Fill initial composite solution with values from each problem
  copyProblemsToComposite(*compositeSoln, SOLUTION);

  // Set up a problem interface for the Problem Manager
  Problem_Interface interface(*this);

  // Now create a composite matrix graph needed for preconditioning
  AA = new Epetra_CrsGraph(Copy, *compositeMap, 0);
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
  vector<Epetra_IntVector>* columns = &colorMapIndex(*AA);
#else
  if(MyPID==0)
    cout << "ERROR: Cannot use EpetraExt with this build !!" << endl;
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
  // TransformToLocal()  is called.  
  A = new Epetra_CrsMatrix(Copy, *AA); 
  A->TransformToLocal();

  NOX::EpetraNew::Interface::Required& reqInt = 
    dynamic_cast<NOX::EpetraNew::Interface::Required&>(interface);
  NOX::Epetra::Vector nox_soln(*compositeSoln);

  // Create the Matrix-Free Jacobian Operator
  NOX::EpetraNew::MatrixFree Jac(interface, *compositeSoln);
  NOX::EpetraNew::Interface::Jacobian& jacInt = Jac;

  NOX::EpetraNew::Interface::Preconditioner& precInt = 
    dynamic_cast<NOX::EpetraNew::Interface::Preconditioner&>(interface);

  NOX::Parameter::List& lsParams = 
    nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver");

  NOX::EpetraNew::LinearSystemAztecOO composite_linearSystem(
    nlParams->sublist("Printing"),
    lsParams,
    jacInt, Jac,
    precInt, *A,
    *compositeSoln);

  //lsParams.setParameter("Preconditioning", "None");
  lsParams.setParameter("Preconditioner", "AztecOO");
  NOX::EpetraNew::Group grp(nlParams->sublist("Printing"), 
    interface, nox_soln, composite_linearSystem);
  grp.computeF();

  NOX::Solver::Manager solver(grp, *statusTest, *nlParams);
  NOX::StatusTest::StatusType status = solver.solve();
  if( status != NOX::StatusTest::Converged )
  { 
    if (MyPID==0)
      cout << "\nMatrix-Free coupled Problem failed to converge !!"  << endl;
    exit(0);
  }

  // Update all problem's solutions with final solutions from solvers
  updateAllWithFinalSolution();

  // Cleanup locally allocated memory
  delete A; A = 0;
  delete AA; AA = 0;

  return true;
}


// These methods are needed to allow inheritance from GenericEpetraProblem base

bool Problem_Manager::evaluate(
              NOX::EpetraNew::Interface::Required::FillType flag,
              const Epetra_Vector *solnVector,
              Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix)
{
  //Determine what to fill (F or Jacobian)
  bool fillF = false;
  bool fillMatrix = false;
  if (rhsVector != 0) {
    fillF = true;
  }
  else {
    fillMatrix = true;
  }

  // Note that incoming matrix is no longer used.  Instead, the problem
  // should own the matrix to be filled into.

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

  if (fillF) {
    computeAllF();
    copyProblemsToComposite(*rhsVector, GROUP_F);
  }

  if (fillMatrix) {
    
    Epetra_CrsMatrix* Matrix = dynamic_cast<Epetra_CrsMatrix*>(A);
    Matrix->PutScalar(0.0);
    
    computeAllJacobian();

    copyProblemJacobiansToComposite();

    Matrix->TransformToLocal();

#ifdef DEBUG
    Matrix->Print(cout);
#endif
  }

  return true;
}

void Problem_Manager::generateGraph()
{ 

  // First construct a graph for each problem's self-dependence
  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Loop over each problem being managed and ascertain its graph as well
  // as its graph from its dependencies
  for( ; problemIter != problemLast; problemIter++) {

    GenericEpetraProblem& problem = *(*problemIter).second;
    int probId = problem.getId();

    Epetra_CrsGraph &problemGraph = problem.getGraph();

    // Use max potential number of nonzero columns to dimension index array
    int problemMaxNodes = problemGraph.Map().MaxAllGID();

    // Get the indices map for copying data from this problem into 
    // the composite problem
    Epetra_IntVector& problemIndices = 
            *ProblemToCompositeIndices.find(probId)->second;

    int problemRow, compositeRow, numCols;

    // First fill composite graph for each problem's self-dependence.  This
    // corresponds to diagonal blocks.
    int* columnIndices = new int[problemMaxNodes];

    for (int i = 0; i<problemGraph.NumMyRows(); i++)
    { 
      problemRow = problemGraph.Map().GID(i);
      problemGraph.ExtractGlobalRowCopy(problemRow, problemMaxNodes, 
                           numCols, columnIndices);

      // Convert row/column indices to composite problem
      compositeRow = problemIndices[problemRow];
      for (int j = 0; j<numCols; j++)
        columnIndices[j] = problemIndices[columnIndices[j]];
      int ierr = AA->InsertGlobalIndices(compositeRow, numCols, columnIndices);
    }
    delete [] columnIndices; columnIndices = 0;
  }

  // Next create inter-problem block graph contributions if desired;
  // default is false

  // Two things are achieved here: 1) The composite Graph is augmented to
  // accommodate off-diagonal blocks and 2) these blocks are packaged as
  // individual NOX::EpetreNew::Group's owneed by the manager.

  if( doOffBlocks ) {
#ifdef HAVE_NOX_EPETRAEXT

    problemIter = Problems.begin();

    // Loop over each problem being managed and ascertain its graph as well
    // as its graph from its dependencies
    for( ; problemIter != problemLast; problemIter++) {
  
      GenericEpetraProblem& problem = *(*problemIter).second;
      int probId = problem.getId();
  
      Epetra_CrsGraph &problemGraph = problem.getGraph();

      // Get the indices map for copying data from this problem into 
      // the composite problem
      Epetra_IntVector& problemIndices = 
              *ProblemToCompositeIndices.find(probId)->second;
  
      // Create containers for the off-block objects
//      vector<Epetra_CrsGraph*> offGraphs;
      vector<OffBlock_Manager*> OffBlock_ManagersVec;

      int problemMaxNodes = problemGraph.Map().NumGlobalElements();

      int problemRow, compositeRow, numCols, numDepCols;
  
      // Loop over each problem on which this one depends
      for( int k = 0; k<problem.depProblems.size(); k++) {

        // Create the off-block graph to be constructed for this 
        // problem-problem coupling
        // NOTE: the map used for the off-block graph is the composite Map
        // to allow valid global indexing
        Epetra_CrsGraph* offGraphPtr = 
          new Epetra_CrsGraph(Copy, *compositeMap, 0);
	Epetra_CrsGraph &offGraph = *offGraphPtr;

//        offGraphs.push_back(new Epetra_CrsGraph(Copy, *compositeMap, 0));
//	Epetra_CrsGraph &offGraph = *offGraphs.back();

        // Get the needed objects for the depend problem
        GenericEpetraProblem &dependProblem = 
          *(Problems.find(problem.depProblems[k])->second);
        int dependId = dependProblem.getId();
        XferOp *xferOpPtr = problem.xferOperators.find(dependId)->second;
	if( !xferOpPtr ) {
          cout << "ERROR: Unable to get Xfer_Operator for dependence of "
               << "problem \"" << problem.getName() << "\" on problem "
               << "\"" << dependProblem.getName() << "\" !!" << endl;
          throw "Problem_Manager ERROR";
        }
        XferOp &xferOp = *xferOpPtr;
        multimap<int,int>& depNodesMap = xferOp.getDependentNodesMap();

        // Get the indices map for copying data from the dependent problem into 
        // the composite problem
        Epetra_IntVector& dependIndices = 
                *ProblemToCompositeIndices.find(dependId)->second;

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
        for (int i = 0; i<problemGraph.NumMyRows(); i++) {

          problemRow = problemGraph.Map().GID(i);

          problemGraph.ExtractGlobalRowCopy(problemRow, problemMaxNodes, 
                               numCols, columnIndices);
  
          // Convert row/column indices to composite problem
          compositeRow = problemIndices[problemRow];
          numDepCols = 0;
          for (int j = 0; j<numCols; j++) {
            int numDepNodes = depNodesMap.count(columnIndices[j]);
            pair< multimap<int, int>::iterator,
                  multimap<int, int>::iterator > rangeN
                = depNodesMap.equal_range(columnIndices[j]);
            multimap<int, int>::iterator iterN;
            for( iterN = rangeN.first; iterN != rangeN.second; iterN++)
              dependentColIndices[numDepCols++] = dependIndices[iterN->second];
          }
          AA->InsertGlobalIndices(compositeRow, numDepCols, dependentColIndices);
          offGraph.InsertGlobalIndices(compositeRow, numDepCols, dependentColIndices);
        }
        delete [] columnIndices; columnIndices = 0;
        delete [] dependentColIndices; dependentColIndices = 0;

	offGraph.TransformToLocal();
	offGraph.SortIndices();
	offGraph.RemoveRedundantIndices();
#ifdef DEBUG
	offGraph.Print(cout);
#endif
        OffBlock_ManagersVec.push_back( new OffBlock_Manager(*this, offGraph,
                                                        probId, dependId) );
      }
//      Off_Graphs.insert( pair<int, vector<Epetra_CrsGraph*> >
//                         (probId, offGraphs) );
//      createFDCobjects(probId);
      OffBlock_Managers.insert( pair<int, vector<OffBlock_Manager*> >
                         (probId, OffBlock_ManagersVec) );
    }
#endif
  } // end doOffBlocks

  AA->TransformToLocal();
  AA->SortIndices();
  AA->RemoveRedundantIndices();

#ifdef DEBUG
  AA->Print(cout);
#endif

  return;
}

void Problem_Manager::outputSolutions(int timeStep)
{ 

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // Loop over each problem being managed and write its solution vector
  // to a file.
  for( ; problemIter != problemLast; problemIter++) {

    GenericEpetraProblem& problem = *(*problemIter).second;
    int probId = problem.getId();

    Epetra_Vector& xMesh = problem.getMesh();
    Epetra_Vector& problemSoln = problem.getSolution();
    int NumMyNodes = xMesh.Map().NumMyElements();
    int NumGlobalNodes = xMesh.Map().NumGlobalElements();

    char file_name[25];
    FILE *ifp;
    (void) sprintf(file_name, "output.%03d.%03d_%05d",probId,MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (int i=0; i<NumMyNodes; i++)
      fprintf(ifp, "%d  %E  %E \n", i, xMesh[i], problemSoln[i]);
    fclose(ifp);
  }
}

void Problem_Manager::outputStatus()
{ 

  map<int, GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  map<int, GenericEpetraProblem*>::iterator problemLast = Problems.end();

  map<int, GenericEpetraProblem*>::iterator dependIter;

  cout << endl << endl << "\t\t********************************" << endl;
  cout                 << "\t\t*******  Problem Summary  ******" << endl;
  cout                 << "\t\t********************************" << endl;
  cout << endl << endl;

  // Loop over each problem being managed and write its solution vector
  // to a file.
  for( ; problemIter != problemLast; problemIter++) {

    GenericEpetraProblem& problem = *(*problemIter).second;
    cout << "\tProblem \"" << problem.getName() << "\" (" << problemIter->first
         << ")\t Depends on:" << endl;
    
    for( int j = 0; j<problem.depProblems.size(); j++ ) {
      dependIter = Problems.find( problem.depProblems[j] );
      cout << "\t\t-------------> \t\t\"" << dependIter->second->getName() 
           << "\"" << endl;
    }
    cout << endl;
  }
}
