#ifndef _NOX_PROBLEM_MANAGER_H
#define _NOX_PROBLEM_MANAGER_H

#include "NOX.H"
#include "NOX_Epetra.H"

#include "Problem_Interface.H"
#include "GenericEpetraProblem.H"

//Forward declarations
//class GenericEpetraProblem;

#include <vector>

class Problem_Manager : GenericEpetraProblem { 

public:


  // Constructor

  // Constructor
  Problem_Manager(Epetra_Comm& Comm, int NumGlobalElements = 0);

  // Destructor
  virtual ~Problem_Manager();

  // Add a problem to the Manager
  void addProblem(GenericEpetraProblem& problem);

  // Register NLSolver parameter list
  void registerParameters(NOX::Parameter::List& List);

  // Register Convergence Test(s)
  void registerStatusTest(NOX::StatusTest::Combo& comboTest);

  // Signal completion of problem, parameter, and convergence test
  // regitrations and trigger the problem setup
  void registerComplete();

  // Solve problems 
  bool solve();

  // Solve problems using matrix-free coupling
  bool solveMF();


  // These allow inheritance from GenericEpetraProblem base class

  bool evaluate(FillType f, const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix,
               NOX::Epetra::Interface::FillType fill=NOX::Epetra::Interface::F);

private:

  vector<GenericEpetraProblem*> Problems;   // Problem container
  vector<NOX::Epetra::Group*> Groups;       // Groups container
  vector<Problem_Interface*> Interfaces;    // Interfaces container
  vector<NOX::Solver::Manager*> Solvers;    // Solvers container

  NOX::Parameter::List* nlParams;          // NLSolver parameter list
  NOX::StatusTest::Combo* statusTest;      // Convergence criterion (ia)

};
#endif
