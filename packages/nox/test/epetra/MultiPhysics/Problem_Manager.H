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
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif  

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

  bool evaluate(NOX::EpetraNew::Interface::Required::FillType type,
               const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix);

private:

  // Creates a composite graph from existing graphs of registered problems
  void generateGraph();

private:

  vector<GenericEpetraProblem*> Problems;   // Problem container
  vector<NOX::EpetraNew::LinearSystemAztecOO*> LinearSystems;
  					    // Groups container
  vector<NOX::EpetraNew::Group*> Groups;    // Groups container
  vector<Problem_Interface*> Interfaces;    // Interfaces container
  vector<NOX::Solver::Manager*> Solvers;    // Solvers container

  NOX::Parameter::List* nlParams;           // NLSolver parameter list
  NOX::StatusTest::Combo* statusTest;       // Convergence criterion (ia)

#ifdef HAVE_NOX_EPETRAEXT
  // Additional containers to allow each problem to use FD Coloring
  vector<EpetraExt::CrsGraph_MapColoring*> TmpMapColorings;
  vector<Epetra_MapColoring*> ColorMaps;
  vector<EpetraExt::CrsGraph_MapColoringIndex*> ColorMapIndexSets;
  vector<vector<Epetra_IntVector>*> ColumnsSets;
  vector<Epetra_Operator*> MatrixOperators;
#endif


};
#endif
