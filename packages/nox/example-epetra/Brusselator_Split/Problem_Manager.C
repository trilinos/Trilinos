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

Problem_Manager::Problem_Manager()
{
}

Problem_Manager::Problem_Manager(NOX::Parameter::List& List)
{
  nlParams = &List;
}

Problem_Manager::Problem_Manager(NOX::Parameter::List& List, 
                                 GenericEpetraProblem& problem)
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

bool Problem_Manager::solve()
{
  if(Problems.empty())
  {
    cout << "ERROR: No problems registered with Problem_Manager !!"
         << endl;
    throw "Problem_Manager ERROR";
  }

  vector<GenericEpetraProblem*>::iterator iter = Problems.begin();
  vector<GenericEpetraProblem*>::iterator last = Problems.end();

  while( iter != last)
    cout << (*iter++)->getSolution() << endl << endl;

  // Setup the problem and solve
  return true;
}

bool Problem_Manager::reset()
{
  // Iterate over the problem and invoke their reset methods
  return true;
}
