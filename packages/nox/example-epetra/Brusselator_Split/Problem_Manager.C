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

Problem_Manager::Problem_Manager(Epetra_Comm& comm, 
                                 int numGlobalElements) :
  GenericEpetraProblem(comm, numGlobalElements),
  nlParams(0),
  statusTest(0)
{
}

Problem_Manager::~Problem_Manager()
{
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

  while( iter != last)
  {
    Interfaces.push_back(new Problem_Interface(**iter));

    Groups.push_back(new NOX::Epetra::Group(nlParams->sublist("Printing"),
      nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver"),
      *Interfaces.back(), (*iter)->getSolution(), (*iter)->getJacobian()));
    Groups.back()->computeF();
   
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

  while( normSum > 1.e-5 ) // Hard-coded convergence criterion for now.
  {
    iter++;

    solverA.reset(grpA, *statusTest, *nlParams);
    solverA.solve();

    // Extract and use final solution
    const NOX::Epetra::Group& finalGroupA =
      dynamic_cast<const NOX::Epetra::Group&>(solverA.getSolutionGroup());
    const Epetra_Vector& finalSolutionA =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroupA.getX())).getEpetraVector();

    problemB.setAuxillarySolution(finalSolutionA);
    solverB.reset(grpB, *statusTest, *nlParams);
    solverB.solve();

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

  // Now create a composite matrix graph needed for preconditioning
  AA = new Epetra_CrsGraph(Copy, compositeMap, 0);
  generateGraph();

  // Create a preconditioning matrix using the graph just created - this 
  // creates a static graph so we can refill the new matirx after
  // TransformToLocal()  is called.  
  A = new Epetra_CrsMatrix(Copy, *AA); 
  A->TransformToLocal();

  Problem_Interface interface(*this);

  NOX::Epetra::MatrixFree Jac(interface, compositeSoln);
  Epetra_CrsMatrix& Prec = *A;

  NOX::Parameter::List& lsParams = 
    nlParams->sublist("Direction").sublist("Newton").sublist("Linear Solver");
//  lsParams.setParameter("Preconditioning", "None");
  lsParams.setParameter("Preconditioning", "AztecOO: User RowMatrix");
  NOX::Epetra::Group grp(nlParams->sublist("Printing"), lsParams,
    interface, compositeSoln, Jac, Prec);
  grp.computeF();

  NOX::Solver::Manager solver(grp, *statusTest, *nlParams);
  NOX::StatusTest::StatusType status = solver.solve();

  if (status != NOX::StatusTest::Converged)
    if (MyPID==0)
      cout << "Nonlinear solver failed to converge!" << endl;

  // Extract and use final solutions
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
  
  for (int i=0; i<solnAlength; i++)
    problemA.getSolution()[i] = finalSolution[i];
  for (int i=0; i<problemB.getSolution().MyLength(); i++)
    problemB.getSolution()[i] = finalSolution[i+solnAlength];

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
    
    grpA.computeJacobian();
    grpB.computeJacobian();

    Epetra_CrsGraph &graphA = (*Problems[0]).getGraph(),
                    &graphB = (*Problems[1]).getGraph();

    Epetra_CrsMatrix &matrixA = (*Problems[0]).getJacobian(),
                     &matrixB = (*Problems[1]).getJacobian();
  
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

