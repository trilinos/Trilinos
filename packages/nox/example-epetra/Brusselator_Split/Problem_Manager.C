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
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

#include "Problem_Manager.H"
#include "GenericEpetraProblem.H"

// Headers needed for Coloring
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif

// Header for Timing info
#include "Epetra_Time.h"

Problem_Manager::Problem_Manager(Epetra_Comm& comm, 
                                 int numGlobalElements) :
  GenericEpetraProblem(comm, numGlobalElements),
  nlParams(0),
  statusTest(0)
{
}

Problem_Manager::~Problem_Manager()
{
  delete AA; AA = 0;
  delete A; A = 0;

  // Iterate over each problem and destroy/free the necessary objects

  vector<GenericEpetraProblem*>::iterator iter = Problems.begin();
  vector<GenericEpetraProblem*>::iterator last = Problems.end();

  vector<NOX::Epetra::Group*>::iterator GroupsIter = Groups.begin();   
  vector<Problem_Interface*>::iterator InterfacesIter = Interfaces.begin();
  vector<NOX::Solver::Manager*>::iterator SolversIter = Solvers.begin();

#ifdef HAVE_NOX_EPETRAEXT
  vector<EpetraExt::CrsGraph_MapColoring*>::iterator TmpMapColoringsIter = TmpMapColorings.begin();
  vector<Epetra_MapColoring*>::iterator ColorMapsIter = ColorMaps.begin();
  vector<EpetraExt::CrsGraph_MapColoringIndex*>::iterator ColorMapIndexSetsIter = ColorMapIndexSets.begin();
  vector<vector<Epetra_IntVector>*>::iterator ColumnsSetsIter = ColumnsSets.begin();
  vector<Epetra_Operator*>::iterator MatrixOperatorsIter = MatrixOperators.begin();
#endif

  while( iter != last)
  {
    delete *SolversIter++;
    delete *GroupsIter++;
#ifdef HAVE_NOX_EPETRAEXT
    delete *MatrixOperatorsIter++;
    delete *TmpMapColoringsIter++;
    delete *ColorMapsIter++;
    delete *ColorMapIndexSetsIter++;
    //delete *ColumnsSetsIter++;
#endif
    iter++; // Problems are owned by the app driver (Example.C)
  }
}

void Problem_Manager::addProblem(GenericEpetraProblem& problem)
{
  Problems.push_back(&problem);
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

  vector<GenericEpetraProblem*>::iterator iter = Problems.begin();
  vector<GenericEpetraProblem*>::iterator last = Problems.end();

  // Make sure everything is starting clean
  assert(Groups.empty() && Interfaces.empty() && Solvers.empty());

  int icount = 0; // Problem counter

  while( iter != last)
  {
    Interfaces.push_back(new Problem_Interface(**iter));

    // Use this for analytic Matrix Fills
//    Groups.push_back(new NOX::Epetra::Group(nlParams->sublist("Printing"),
//      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
//      *Interfaces.back(), (*iter)->getSolution(), (*iter)->getJacobian()));

    // OR use this to fill matrices using Finite-Differences with Coloring
#ifdef HAVE_NOX_EPETRAEXT
    // Create a timer for performance
    Epetra_Time fillTime(*Comm);

    bool verbose = false;
    EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
      EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY;
    TmpMapColorings.push_back(new 
      EpetraExt::CrsGraph_MapColoring(algType, verbose));
    ColorMaps.push_back(&((*TmpMapColorings.back())((*iter)->getGraph())));
    ColorMapIndexSets.push_back(new 
      EpetraExt::CrsGraph_MapColoringIndex(*ColorMaps.back()));
    ColumnsSets.push_back(&(*ColorMapIndexSets.back())((*iter)->getGraph()));

    if (MyPID == 0)
      printf("\n\tTime to color Jacobian # %d --> %e sec. \n\n",
                  icount++,fillTime.ElapsedTime());
    MatrixOperators.push_back(new
      NOX::Epetra::FiniteDifferenceColoring(*Interfaces.back(), 
        (*iter)->getSolution(), (*iter)->getGraph(), *ColorMaps.back(), 
        *ColumnsSets.back()));

    Groups.push_back(new NOX::Epetra::Group(nlParams->sublist("Printing"),
      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      *Interfaces.back(), (*iter)->getSolution(), *MatrixOperators.back()));
#else
    if(MyPID==0)
      cout << "ERROR: Cannot use EpetraExt with this build !!" << endl;
    exit(0);
#endif

    Groups.back()->computeF(); // Needed to establish convergence state
   
    Solvers.push_back(new NOX::Solver::Manager(*Groups.back(), *statusTest,
                                               *nlParams));
    iter++;
  }

  return;

}

bool Problem_Manager::solve()
{
  if(Problems.empty())
  {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  assert( !Groups.empty() );
  assert( !Interfaces.empty() );
  assert( !Solvers.empty() );

  vector<GenericEpetraProblem*>::iterator problemIter = Problems.begin();
  vector<GenericEpetraProblem*>::iterator problemLast = Problems.end();

  // These iterators would be needed in general, but we later specialize
  // for the case of a 2-problem system.
  vector<NOX::Epetra::Group*>::iterator     groupIter = Groups.begin();
  vector<Problem_Interface*>::iterator  interfaceIter = Interfaces.begin();
  vector<NOX::Solver::Manager*>::iterator  solverIter = Solvers.begin();

  // This is set up for more than 2 problems, but for now, we deal explicitly
  // with just 2.

  GenericEpetraProblem &problemA = *Problems[0],
                       &problemB = *Problems[1];

  NOX::Epetra::Group &grpA = *Groups[0],
                     &grpB = *Groups[1];

  NOX::Solver::Manager &solverA = *Solvers[0],
                       &solverB = *Solvers[1];

  // Sync the two problems and get initial convergence state
  problemA.setAuxillarySolution(problemB.getSolution());
  problemB.setAuxillarySolution(problemA.getSolution());
  grpA.setX(problemA.getSolution());
  grpB.setX(problemB.getSolution());
  grpA.computeF();
  grpB.computeF();
  cout << "Initial 2-Norms of F (A, B) --> " << grpA.getNormF() << " ,  "
       << grpB.getNormF() << endl;
  double normSum = grpA.getNormF() + grpB.getNormF();

  // Now do the decoupled solve
  int iter = 0;
  NOX::StatusTest::StatusType status;

  while( normSum > 1.e-5 ) // Hard-coded convergence criterion for now.
  {
    iter++;

    solverA.reset(grpA, *statusTest, *nlParams);
    status = solverA.solve();
    if( status != NOX::StatusTest::Converged )
    { 
      if (MyPID==0)
        cout << "\nRegistered Problem A failed to converge !!"  << endl;
      exit(0);
    }

    // Extract and use final solution
    const NOX::Epetra::Group& finalGroupA =
      dynamic_cast<const NOX::Epetra::Group&>(solverA.getSolutionGroup());
    const Epetra_Vector& finalSolutionA =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroupA.getX())).getEpetraVector();

    problemB.setAuxillarySolution(finalSolutionA);
    solverB.reset(grpB, *statusTest, *nlParams);
    status = solverB.solve();
    if( status != NOX::StatusTest::Converged )
    { 
      if (MyPID==0)
        cout << "\nRegistered Problem B failed to converge !!"  << endl;
      exit(0);
    }

    // Extract and use final solution
    const NOX::Epetra::Group& finalGroupB =
      dynamic_cast<const NOX::Epetra::Group&>(solverB.getSolutionGroup());
    const Epetra_Vector& finalSolutionB =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroupB.getX())).getEpetraVector();
  
    problemA.setAuxillarySolution(finalSolutionB);
    grpA.setX(finalSolutionA);
    grpA.computeF();
    problemB.setAuxillarySolution(finalSolutionA);
    grpB.setX(finalSolutionB);
    grpB.computeF();

    cout << "Decoupled iteration #" << iter << " : 2-Norms of F (A, B) --> " 
         << grpA.getNormF() << " ,  " << grpB.getNormF() << endl;
    normSum = grpA.getNormF() + grpB.getNormF();
  }
  
  cout << "\nDecoupled solution required --> " << iter << " iterations.\n" 
       << endl;

  // Extract and use final solutions
  const NOX::Epetra::Group& finalGroupA =
    dynamic_cast<const NOX::Epetra::Group&>(solverA.getSolutionGroup());
  const NOX::Epetra::Group& finalGroupB =
    dynamic_cast<const NOX::Epetra::Group&>(solverB.getSolutionGroup());
  const Epetra_Vector& finalSolutionA =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroupA.getX())).getEpetraVector();
  const Epetra_Vector& finalSolutionB =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroupB.getX())).getEpetraVector();
  
  // Put final solutions back into problems
  problemA.setSolution(finalSolutionA);
  problemB.setSolution(finalSolutionB);

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

  // This is set up for more than 2 problems, but for now, we deal explicitly
  // with just 2.

  GenericEpetraProblem &problemA = *Problems[0],
                       &problemB = *Problems[1];

  NOX::Epetra::Group   &grpA = *Groups[0],
                       &grpB = *Groups[1];

  NOX::Solver::Manager &solverA = *Solvers[0],
                       &solverB = *Solvers[1];

  // Sync the two problems and get initial convergence state
  problemA.setAuxillarySolution(problemB.getSolution());
  problemB.setAuxillarySolution(problemA.getSolution());
  grpA.setX(problemA.getSolution());
  grpB.setX(problemB.getSolution());
  grpA.computeF();
  grpB.computeF();
  cout << "Initial 2-Norms of F (A, B) --> " << grpA.getNormF() << " ,  "
       << grpB.getNormF() << endl;

  // Construct a composite Epetra_Map
  NumMyNodes = problemA.NumMyNodes + problemB.NumMyNodes;
  int* myGlobalNodes = new int[NumMyNodes];
  for (int i=0; i<problemA.NumMyNodes; i++)
    myGlobalNodes[i] = problemA.StandardMap->GID(i);
  int maxAllAGID = problemA.StandardMap->MaxAllGID();
  for (int i=0; i<problemB.NumMyNodes; i++)
    myGlobalNodes[i+problemA.NumMyNodes] = maxAllAGID + 1 +
                                           problemB.StandardMap->GID(i);
  Epetra_Map compositeMap(-1, NumMyNodes, myGlobalNodes, 0, *Comm);
  delete [] myGlobalNodes; myGlobalNodes = 0;

  Epetra_Vector compositeSoln(compositeMap);

  // Fill initial composite solution with values from each problem
  int solnAlength = problemA.getSolution().MyLength();
  for (int i=0; i<solnAlength; i++)
    compositeSoln[i] = problemA.getSolution()[i];
  for (int i=0; i<problemB.getSolution().MyLength(); i++)
    compositeSoln[i+solnAlength] = problemB.getSolution()[i];

  // Set up a problem interface for the Problem Manager
  Problem_Interface interface(*this);

  // Now create a composite matrix graph needed for preconditioning
  AA = new Epetra_CrsGraph(Copy, compositeMap, 0);
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

  // Create the Matrix-Free Jacobian Operator
  NOX::Epetra::MatrixFree Jac(interface, compositeSoln);

  NOX::Parameter::List& lsParams = 
    nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver");
  //lsParams.setParameter("Preconditioning", "None");
  lsParams.setParameter("Preconditioning", "AztecOO: User RowMatrix");
  NOX::Epetra::Group grp(nlParams->sublist("Printing"), lsParams,
    interface, compositeSoln, Jac, *A);
  grp.computeF();

  NOX::Solver::Manager solver(grp, *statusTest, *nlParams);
  NOX::StatusTest::StatusType status = solver.solve();
  if( status != NOX::StatusTest::Converged )
  { 
    if (MyPID==0)
      cout << "\nMatrix-Free coupled Problem failed to converge !!"  << endl;
    exit(0);
  }

  // Extract and use final solutions
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
  
  for (int i=0; i<solnAlength; i++)
    problemA.getSolution()[i] = finalSolution[i];
  for (int i=0; i<problemB.getSolution().MyLength(); i++)
    problemB.getSolution()[i] = finalSolution[i+solnAlength];

  // Cleanup locally allocated memory
  delete A; A = 0;
  delete AA; AA = 0;

  return true;
}


// These methods are needed to allow inheritance from GenericEpetraProblem base

bool Problem_Manager::evaluate(FillType f, const Epetra_Vector *solnVector,
              Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix,
              NOX::Epetra::Interface::FillType fill)
{

  GenericEpetraProblem &problemA = *Problems[0],
                       &problemB = *Problems[1];

  NOX::Epetra::Group   &grpA = *Groups[0],
                       &grpB = *Groups[1];

  Epetra_Vector solnA(*problemA.StandardMap);
  Epetra_Vector solnB(*problemB.StandardMap);

  for (int i=0; i<solnA.MyLength(); i++)
    solnA[i] = (*solnVector)[i];
  for (int i=0; i<solnB.MyLength(); i++)
    solnB[i] = (*solnVector)[i+solnA.MyLength()];

  // Pass solutions and compute residuals
  problemA.setAuxillarySolution(solnB);
  problemB.setAuxillarySolution(solnA);
  grpA.setX(solnA);
  grpB.setX(solnB);

  if (f == F_ONLY || f == ALL)
  {
    grpA.computeF();
    grpB.computeF();

    for (int i=0; i<solnA.MyLength(); i++)
      (*rhsVector)[i] = dynamic_cast<const NOX::Epetra::Vector&>(grpA.getF()).
                        getEpetraVector()[i];
    for (int i=0; i<solnB.MyLength(); i++)
      (*rhsVector)[i+solnA.MyLength()] = 
                        dynamic_cast<const NOX::Epetra::Vector&>(grpB.getF()).
                        getEpetraVector()[i];
  }

  if (f == MATRIX_ONLY || f == ALL)
  {
    
    Epetra_CrsMatrix* Matrix = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    Matrix->PutScalar(0.0);
    
    // Create a timer for performance
    Epetra_Time fillTime(*Comm);

    grpA.computeJacobian();
    if (MyPID == 0)
      printf("\n\tTime to fill Jacobian A --> %e sec. \n\n",
                  fillTime.ElapsedTime());
    fillTime.ResetStartTime();
    grpB.computeJacobian();
    if (MyPID == 0)
      printf("\n\tTime to fill Jacobian B --> %e sec. \n\n",
                  fillTime.ElapsedTime());

    Epetra_CrsGraph &graphA = (*Problems[0]).getGraph(),
                    &graphB = (*Problems[1]).getGraph();

    Epetra_CrsMatrix *matrixAPtr,
                     *matrixBPtr;

    // Use each group's operator test to determine the type of Jacobian
    // operator being used by each.
    NOX::Epetra::Group::OperatorType opTypeA = grpA.getOperatorType(
                                     grpA.getSharedJacobian().getOperator());
    NOX::Epetra::Group::OperatorType opTypeB = grpB.getOperatorType(
                                     grpB.getSharedJacobian().getOperator());

    if (opTypeA == NOX::Epetra::Group::NoxFiniteDifferenceRowMatrix )
      matrixAPtr = const_cast<Epetra_CrsMatrix*>(
        &dynamic_cast<const NOX::Epetra::FiniteDifference&>
        (grpA.getSharedJacobian().getOperator()).getUnderlyingMatrix());
    else if (opTypeA == NOX::Epetra::Group::EpetraRowMatrix )
      // NOTE: We are getting the matrix from the problem.  This SHOULD be
      // the same matrix wrapped in the group.  A safer alternative would be
      // to get this matrix from the group as above for a more general 
      // operator.
      matrixAPtr = &(*Problems[0]).getJacobian();
    else
    {
      if (MyPID==0)
        cout << "Jacobian operator for Problem A not supported for "
             << "preconditioning Matrix-Free coupling solver." << endl;
      exit(0);
    }

    if (opTypeB == NOX::Epetra::Group::NoxFiniteDifferenceRowMatrix )
      matrixBPtr = const_cast<Epetra_CrsMatrix*>(
        &dynamic_cast<const NOX::Epetra::FiniteDifference&>
        (grpB.getSharedJacobian().getOperator()).getUnderlyingMatrix());
    else if (opTypeB == NOX::Epetra::Group::EpetraRowMatrix )
      // NOTE: We are getting the matrix from the problem.  This SHOULD be
      // the same matrix wrapped in the group.  A safer alternative would be
      // to get this matrix from the group as above for a more general 
      // operator.
      matrixBPtr = &(*Problems[1]).getJacobian();
    else
    {
      if (MyPID==0)
        cout << "Jacobian operator for Problem B not supported for "
             << "preconditioning Matrix-Free coupling solver." << endl;
      exit(0);
    }

    // Create convenient references for each Jacobian operator
    Epetra_CrsMatrix &matrixA = *matrixAPtr,
                     &matrixB = *matrixBPtr;

    int maxAllAGID = matrixA.Map().MaxAllGID();
    int* indices = new int[maxAllAGID];
    double* values = new double[maxAllAGID];
    int i, row, numCols, numValues;
    for (i=0; i<matrixA.NumMyRows(); i++)
    { 
      row = matrixA.Map().GID(i);
      graphA.ExtractGlobalRowCopy(row, maxAllAGID, numCols, indices);
      matrixA.ExtractGlobalRowCopy(row, maxAllAGID, numValues, values);
      //printf("MatrixA, row --> %d\tnumValues --> %d\n",row,numValues);
      //for (int j=0; j<numValues; j++)
      //  cout << "\t[" << indices[j] << "]  " << values[j];
      //cout << endl;
      int ierr = Matrix->ReplaceGlobalValues(row, numValues, values, indices);
      //printf("\nAfter insertion, ierr --> %d\n\n",ierr);
    }
    delete [] values; values = 0;
    delete [] indices; indices = 0;

    int maxAllBGID = graphB.Map().MaxAllGID();
    indices = new int[maxAllBGID];
    values = new double[maxAllBGID];
    for (i=0; i<matrixB.NumMyRows(); i++)
    { 
      row = matrixB.Map().GID(i);
      graphB.ExtractGlobalRowCopy(row, maxAllBGID, numCols, indices);
      matrixB.ExtractGlobalRowCopy(row, maxAllBGID, numValues, values);
      for (int j=0; j<numCols; j++)
        indices[j] += maxAllAGID + 1;
      row += maxAllAGID + 1;
      Matrix->ReplaceGlobalValues(row, numValues, values, indices);
    }
    delete [] values; values = 0;
    delete [] indices; indices = 0;
  
    // Sync up processors to be safe
    Comm->Barrier();
  
    Matrix->TransformToLocal();

    //matrixA.Print(cout);
    //matrixB.Print(cout);
    //Matrix->Print(cout);
  }

  return true;
}

void Problem_Manager::generateGraph()
{ 

  // Here again, a general capability has been specialized to 2 problems
  
  Epetra_CrsGraph &graphA = (*Problems[0]).getGraph(),
                  &graphB = (*Problems[1]).getGraph();

  int maxAllAGID = graphA.Map().MaxAllGID();
  int* indices = new int[maxAllAGID];
  int i, row, numCols;
  for (i=0; i<graphA.NumMyRows(); i++)
  {
    row = graphA.Map().GID(i);
    graphA.ExtractGlobalRowCopy(row, maxAllAGID, numCols, indices);
    AA->InsertGlobalIndices(row, numCols, indices);
  }
  delete [] indices; indices = 0;

  int maxAllBGID = graphB.Map().MaxAllGID();
  indices = new int[maxAllBGID];
  for (i=0; i<graphB.NumMyRows(); i++)
  {
    row = graphB.Map().GID(i);
    graphA.ExtractGlobalRowCopy(row, maxAllBGID, numCols, indices);
    for (int j=0; j<numCols; j++)
      indices[j] += maxAllAGID + 1;
    row += maxAllAGID + 1;
    AA->InsertGlobalIndices(row, numCols, indices);
  }
  delete [] indices; indices = 0;

  AA->TransformToLocal();
  AA->SortIndices();
  AA->RemoveRedundantIndices();

  //AA->Print(cout);

  return;
}

