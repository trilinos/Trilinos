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

Problem_Manager::Problem_Manager() :
  nlParams(0),
  statusTest(0)
{
}

Problem_Manager::Problem_Manager(NOX::Parameter::List& List) :
  statusTest(0)
{
  nlParams = &List;
}

Problem_Manager::Problem_Manager(NOX::Parameter::List& List, 
                                 GenericEpetraProblem& problem) :
  statusTest(0)
{
  nlParams = &List;
  addProblem(problem);
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

  vector<NOX::Epetra::Group*>::iterator     groupIter = Groups.begin();
  vector<Problem_Interface*>::iterator  interfaceIter = Interfaces.begin();
  vector<NOX::Solver::Manager*>::iterator  solverIter = Solvers.begin();

  while( problemIter != problemLast )
  {
    (*solverIter)->solve();
    cout << "\n\n\n\t\t------- Finished solving a problem !! -------" << endl;
    problemIter++;
    groupIter++;
    interfaceIter++;
    solverIter++;
    //cout << (*iter++)->getSolution() << endl << endl;
  }

  // Setup the problem and solve
  return true;
}

bool Problem_Manager::reset()
{
  // Iterate over the problem and invoke their reset methods
  return true;
}
