#ifndef _NOX_PROBLEM_MANAGER_H
#define _NOX_PROBLEM_MANAGER_H

#include "NOX.H"
#include "NOX_Epetra.H"

#include "Problem_Interface.H"

//Forward declarations
class GenericEpetraProblem;

#include <vector>

class Problem_Manager { 

public:

  // Constructor
  Problem_Manager();

  // Constructor
  Problem_Manager(NOX::Parameter::List& List);

  // Constructor
  Problem_Manager(NOX::Parameter::List& List, GenericEpetraProblem&);

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

  // Reset problems
  bool reset();

private:

  vector<GenericEpetraProblem*> Problems;  // Problem container
  vector<NOX::Epetra::Group*> Groups;       // Groups container
  vector<Problem_Interface*> Interfaces;    // Interfaces container
  vector<NOX::Solver::Manager*> Solvers;    // Solvers container

  NOX::Parameter::List* nlParams;          // NLSolver parameter list
  NOX::StatusTest::Combo* statusTest;      // Convergence criterion (ia)

};
#endif
