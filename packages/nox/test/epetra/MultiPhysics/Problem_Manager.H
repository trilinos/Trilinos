//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#ifndef NOX_PROBLEM_MANAGER_H
#define NOX_PROBLEM_MANAGER_H

#include "NOX.H"
#include "NOX_Multiphysics_DataExchange_Interface.H"
#include "NOX_Epetra.H"

#include "Problem_Interface.H"
#include "GenericEpetraProblem.H"
#include "Teuchos_RefCountPtr.hpp"

#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
class OffBlock_Manager;
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif  

#include <vector>

class Problem_Manager : public GenericEpetraProblem ,
                        public NOX::Multiphysics::DataExchange::Interface
{

public:

  enum VectorType { SOLUTION, GROUP_F };

  // Constructor
  Problem_Manager(Epetra_Comm& Comm, bool doOffBlocks = false,
                  int NumGlobalElements = 0, bool useMatlab = false);

  // Destructor
  virtual ~Problem_Manager();

  // Add a problem to the Manager
  void addProblem(GenericEpetraProblem& problem);

  // Get the number of registered problems
  int getProblemCount() { return problemCount; }

  // Get the problem names container
  map<int, std::string> & getNames() { return Names; }

  // Return my Problems container
  map<int, Teuchos::RefCountPtr<GenericEpetraProblem> > & getProblems() 
  { return Problems; }

  // Get a problem given its unique id
  GenericEpetraProblem & getProblem(int id);

  // Get a problem given its name
  GenericEpetraProblem & getProblem(string name);

  // Get a group given its unique id
  const NOX::Epetra::Group & getSolutionGroup(int id);

  // Get the solution vector from a problem group
  const Epetra_Vector & getSolutionVec(int id);

  // Get an individual residual vector given an integer problem id
  const Epetra_Vector & getResidualVec(int i);

  // Get a group given its unique id
  NOX::Epetra::Group& getSolutionGroup(int id);

  // Get the solution (Epetra) vector from a problem (Solv_NOX_Epetra)
  const Epetra_Vector & getSolutionVec(int id);

  // Get a group given its unique id
  NOX::Epetra::Group& getGroup(int id);

  // Get the compositeSoln vector
  Teuchos::RefCountPtr<Epetra_Vector> getCompositeSoln();

  // Get the container of solver managers
  map<int, Teuchos::RefCountPtr<NOX::Solver::Manager> > & getSolvers()
  {  return Solvers; }

  // Get the composite NOX solver
  Teuchos::RefCountPtr<NOX::Solver::Manager> & getCompositeSolver()
  {  return compositeSolver; }

  // Create a dependence of Problem A equations on Problem B variables
  void createDependency(string nameA, string nameB, bool interfacialCoupling = false);

  // Create a dependence of Problem A equations on Problem B variables
  void createDependency(GenericEpetraProblem& problemA, GenericEpetraProblem& problemB, 
                        bool interfacialCoupling = false);

  // Register NLSolver parameter list
  void registerParameters(const Teuchos::RefCountPtr<Teuchos::ParameterList>& List);

  // Register Convergence Test(s)
  void registerStatusTest(const Teuchos::RefCountPtr<NOX::StatusTest::Combo>& comboTest);

  // Signal completion of problem, parameter, and convergence test
  // regitrations and trigger the problem setup
  void registerComplete();

  // Have each problem receive dependent data from every other problem 
  // on which it depends
  void syncAllProblems();

  // Set all problems with the same time step size
  void setAlldt( double dt);

  // Set the solver group solution using the current solution contained 
  // in the problem
  void setGroupX(int id);

  // Set the solver group solution using a provided vector
  void setGroupX(int id, Epetra_Vector &);

  // Set all groups with the current solution contained in each problem
  void setAllGroupX();

#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
  // Set all off-block groups with the current composite-size solution
  void setAllOffBlockGroupX(const Epetra_Vector&);
#endif

  // Reset all problems by copying the current solution into the old solution
  // Needed for time-dependent problems
  void resetProblems();

  // Compute residual (F) in solution groups for all registered problems
  void computeAllF();

  // Compute residual (F) in solution group specified by id
  void computeGroupF(int id);

  // Compute Jacobian in solution groups for all registered problems
  void computeAllJacobian();

  // Compute an individual block Jacobian given an integer pairing of problem and
  // the problem variables on which it depends
  void computeBlockJacobian(int i, int j = -1);

  // Get an individual block Jacobian given an integer pairing of problem and
  // the problem variables on which it depends
  Teuchos::RefCountPtr<Epetra_CrsMatrix> getBlockJacobianMatrix(int i, int j = -1);

  // Compute sum of norms of each problem residual
  double getNormSum();

  // Copy final solution from NOX solver into the problem's solution vector
  void copyGroupCurrentXtoProblemX(int id);

  // Copy final solution from NOX solvers into each problem's solution vector
  void copyAllGroupXtoProblems();

  // Reset specified solvers current solution group X vector with itself, thereby resetting all valid flags
  void resetCurrentGroupX(int probId);

  // Reset each solvers current solution group X vector with itself, thereby resetting all valid flags
  void resetAllCurrentGroupX();

  // Copy a composite problem vector to each problem's vector
  void copyCompositeToProblems(const Epetra_Vector& compositeVec, VectorType);

  // Copy part of a composite problem vector to a problem's vector
  void copyCompositeToVector(const Epetra_Vector& compositeVec, int id,
                             Epetra_Vector& problemVec);

  // Copy vectors from each problem into a composite problem vector
  void copyProblemsToComposite(Epetra_Vector& compositeVec, VectorType);

  // Copy a vector from a problem into part of a composite problem vector
  void copyVectorToComposite(Epetra_Vector& compositeVec, int id,
                             const Epetra_Vector& problemVec);

  // Copy problem Jacobians as block diagonal contributions to 
  // composite Jacobian
  void copyProblemJacobiansToComposite();

  // Solve problem using fixed-point iteration
  bool solve();

  // Solve problems using matrix-free coupling
  bool solveMF();

  // Output summary of copled problems
  void outputStatus( ostream & os );

  // Solve problems using matrix-free coupling
  void outputSolutions(int timeStep = 0);

  // Query whether or not offBlocks ar used
  bool useOffBlocks()
  { return doOffBlocks; }

  // Return my name
  string getName() const 
  { return GenericEpetraProblem::getName(); };

  // Return the problem name queries
  string getProblemName(int id) 
  { return Names[id]; };

  // These allow inheritance from GenericEpetraProblem base class

  bool evaluate(NOX::Epetra::Interface::Required::FillType type,
               const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector);

  // Inherited methods to allow use of library fixed-point
  virtual bool exchangeAllData();

  virtual bool exchangeDataTo(int solverId);

private:

  // Creates a composite graph from existing graphs of registered problems
  void generateGraph();

public:
  // Made these public for now
  Teuchos::RefCountPtr<Teuchos::ParameterList> nlParams; // NLSolver parameter list
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> statusTest;       // Convergence criterion (ia)

private:

  int problemCount;                         // Number of registered Problems
  bool doOffBlocks;	// Flag to specify whether or not to compute
			// off-diagonal composite Jacobian blocks via FDC
  bool useMatlab;	// Flag to specify whether or not to use a Matlab interactive session

  Teuchos::RefCountPtr<Epetra_Map> compositeMap;	// Composite map for all registered problems
  Teuchos::RefCountPtr<Epetra_Vector> compositeSoln;	// Composite map for all registered problems
  map<int, Teuchos::RefCountPtr<GenericEpetraProblem> > Problems; // Problem container
  map<int, string> Names;                   // Problem names
  map<string, int> NameLookup;              // Needed to lookup problems by name
  map<int, Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> > LinearSystems;
  					    // Linear Systems container
  map<int, Teuchos::RefCountPtr<NOX::Epetra::Group> >Groups;     // Groups container
  map<int, Teuchos::RefCountPtr<Problem_Interface>  >Interfaces; // Interfaces container
  map<int, Teuchos::RefCountPtr<NOX::Solver::Manager> > Solvers; // Solvers container
  map<int, Teuchos::RefCountPtr<Epetra_IntVector> > ProblemToCompositeIndices;
    // Essentially a mapping from a problems indices to the composite indices

  Teuchos::RefCountPtr<Epetra_Operator> jacOperator   ; // The composite Jacobian operator
  Teuchos::RefCountPtr<Epetra_Operator> precOperator  ; // The composite Preconditioner operator
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> jacInterface;
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> precInterface;
  Teuchos::RefCountPtr<NOX::Epetra::Group> compositeGroup  ; // The composite solution group
  Teuchos::RefCountPtr<NOX::Solver::Manager> compositeSolver ; // The composite solver

#ifdef HAVE_NOX_EPETRAEXT
  // Additional containers to allow each problem to use FD Coloring
  map<int, Teuchos::RefCountPtr<EpetraExt::CrsGraph_MapColoring> > TmpMapColorings;
  map<int, Teuchos::RefCountPtr<Epetra_MapColoring> > ColorMaps;
  map<int, Teuchos::RefCountPtr<EpetraExt::CrsGraph_MapColoringIndex> > ColorMapIndexSets;
  map<int, Teuchos::RefCountPtr<vector<Epetra_IntVector> > > ColumnsSets;
  map<int, Teuchos::RefCountPtr<Epetra_Operator> > MatrixOperators;

  // These containers allow off-diagonal block contributions representing
  // inter-problem coupling
  map<int, vector<OffBlock_Manager*> > OffBlock_Managers;    
#endif


};
#endif
