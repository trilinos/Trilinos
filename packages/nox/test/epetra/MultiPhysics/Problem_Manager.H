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
                                                                                
#ifndef _NOX_PROBLEM_MANAGER_H
#define _NOX_PROBLEM_MANAGER_H

#include "NOX.H"
#include "NOX_Epetra.H"

#include "Problem_Interface.H"
#include "GenericEpetraProblem.H"

#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
class OffBlock_Manager;
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif  

#include <vector>

class Problem_Manager : GenericEpetraProblem { 

public:

  enum VectorType { SOLUTION, GROUP_F };

  // Constructor
  Problem_Manager(Epetra_Comm& Comm, bool doOffBlocks = false,
                  int NumGlobalElements = 0);

  // Destructor
  virtual ~Problem_Manager();

  // Add a problem to the Manager
  void addProblem(GenericEpetraProblem& problem);

  // Get a problem given its unique id
  GenericEpetraProblem& getProblem(int id);

  // Get a group given its unique id
  NOX::EpetraNew::Group& getGroup(int id);

  // Get the compositeSoln vector
  Epetra_Vector& getCompositeSoln();

  // Create a dependence of Problem A equations on Problem B variables
  void createDependency(string nameA, string nameB);

  // Create a dependence of Problem A equations on Problem B variables
  void createDependency(GenericEpetraProblem& problemA,
                        GenericEpetraProblem& problemB);

  // Register NLSolver parameter list
  void registerParameters(NOX::Parameter::List& List);

  // Register Convergence Test(s)
  void registerStatusTest(NOX::StatusTest::Combo& comboTest);

  // Signal completion of problem, parameter, and convergence test
  // regitrations and trigger the problem setup
  void registerComplete();

  // Have each problem receive auxillary data from every other problem 
  // on which it depends
  bool syncAllProblems();

  // Set the solver group solution using the current solution contained 
  // in the problem
  bool setGroupX(int id);

  // Set all groups with the current solution contained in each problem
  bool setAllGroupX();

#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
  // Set all off-block groups with the current composite-size solution
  bool setAllOffBlockGroupX(const Epetra_Vector&);
#endif

  // Reset all problems by copying the current solution into the old solution
  // Needed for time-dependent problems
  void resetProblems();

  // Compute residual (F) in solution groups for all registered problems
  bool computeAllF();

  // Compute residual (F) in solution group specified by id
  bool computeGroupF(int id);

  // Compute Jacobian in solution groups for all registered problems
  bool computeAllJacobian();

  // Compute sum of norms of each problem residual
  double getNormSum();

  // Copy final solution from NOX solver into the problem's solution vector
  void updateWithFinalSolution(int id);

  // Copy final solution from NOX solvers into each problem's solution vector
  void updateAllWithFinalSolution();

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

  // Solve problems using matrix-free coupling
  void outputSolutions(int timeStep = 0);

  // Print status of Problem Manager, problems and dependencies
  void outputStatus();

  // Return the problem name queries
  string getName(int id) { return Names.find(id)->second; };


  // These allow inheritance from GenericEpetraProblem base class

  bool evaluate(NOX::EpetraNew::Interface::Required::FillType type,
               const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix);

private:

  // Creates a composite graph from existing graphs of registered problems
  void generateGraph();

public:
  // Made these public for now
  NOX::Parameter::List* nlParams;           // NLSolver parameter list
  NOX::StatusTest::Combo* statusTest;       // Convergence criterion (ia)

private:

  int problemCount;                         // Number of registered Problems
  bool doOffBlocks;	// Flag to specify whether or not to compute
			// off-diagonal composite Jacobian blocks via FDC

  Epetra_Map* compositeMap;	// Composite map for all registered problems
  Epetra_Vector* compositeSoln;	// Composite map for all registered problems

  map<int, GenericEpetraProblem*> Problems; // Problem container
  map<int, string> Names;                   // Problem names
  map<string, int> NameLookup;              // Needed to lookup problems by name
  map<int, NOX::EpetraNew::LinearSystemAztecOO*> LinearSystems;
  					    // Linear Systems container
  map<int, NOX::EpetraNew::Group*> Groups;    // Groups container
  map<int, Problem_Interface*> Interfaces;    // Interfaces container
  map<int, NOX::Solver::Manager*> Solvers;    // Solvers container
  map<int, Epetra_IntVector*> ProblemToCompositeIndices;
    // Essentially a mapping from a problems indices to the composite indices

  Problem_Interface *compositeProblemInterface; // A problem interface to 
  // the composite problem manager.  Needed for off-diagonal FDC blocks

#ifdef HAVE_NOX_EPETRAEXT
  // Additional containers to allow each problem to use FD Coloring
  map<int, EpetraExt::CrsGraph_MapColoring*> TmpMapColorings;
  map<int, Epetra_MapColoring*> ColorMaps;
  map<int, EpetraExt::CrsGraph_MapColoringIndex*> ColorMapIndexSets;
  map<int, vector<Epetra_IntVector>*> ColumnsSets;
  map<int, Epetra_Operator*> MatrixOperators;

  // These containers allow off-diagonal block contributions representing
  // inter-problem coupling
  map<int, vector<OffBlock_Manager*> > OffBlock_Managers;    
#endif


};
#endif
