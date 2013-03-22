//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include "NOX_Epetra_SchurInterface.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_BroydenOperator.H"

#include "Problem_Interface.H"
#include "GenericEpetraProblem.H"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
class OffBlock_Manager;
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif  

#include <vector>

enum CouplingSolveMethod { LOOSE_HARDCODED,
                           LOOSE_LIBRARY  ,
                           MATRIX_FREE    ,
                           SCHUR_BASED      };

const CouplingSolveMethod SolveMethodValues[] = { LOOSE_HARDCODED, LOOSE_LIBRARY, MATRIX_FREE, SCHUR_BASED };

static const char * SolveMethodNames[] = { "loose-hardcoded" ,
                                           "loose-library"   ,
                                           "matrix-free"     ,
                                           "schur-based"       };

class Problem_Manager : public GenericEpetraProblem                               ,
                        public NOX::Multiphysics::DataExchange::Interface         ,
                        public NOX::Epetra::SchurInterface                        ,
                        public NOX::Epetra::BroydenOperator::ReplacementInterface
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
  map<int, Teuchos::RCP<GenericEpetraProblem> > & getProblems() 
  { return Problems; }

  // Get a problem given its unique id
  GenericEpetraProblem & getProblem(int id);

  // Get a problem given its name
  GenericEpetraProblem & getProblem(std::string name);

  // Get a group given its unique id
  const NOX::Epetra::Group & getProblemSolutionGroup(int id);

  // Get the solution vector from a problem group
  const Epetra_Vector & getSolutionVec(int id);

  // Get an individual residual vector given an integer problem id
  const Epetra_Vector & getResidualVec(int i);

  // Get a group given its unique id
  NOX::Epetra::Group& getGroup(int id);

  // Get the compositeSoln vector
  Teuchos::RCP<Epetra_Vector> getCompositeSoln();

  // Get the compositeSoln Map
  Teuchos::RCP<Epetra_Map> getCompositeMap()
    { return compositeMap; }

  // Get the container of solver managers
  map<int, Teuchos::RCP<NOX::Solver::Generic> > & getSolvers()
  {  return Solvers; }

  // Get the composite NOX solver
  Teuchos::RCP<NOX::Solver::Generic> & getCompositeSolver()
  {  return compositeSolver; }

  // Get mapping of problem indices to composite
  map<int, Teuchos::RCP<Epetra_IntVector> > & getProblemToCompositeIndices()
  {  return ProblemToCompositeIndices; }

  // Get jacobian operator
  Teuchos::RCP<Epetra_Operator> & getJacobianOperator()
  {  return jacOperator; }

  // Create a dependence of Problem A equations on Problem B variables
  void createDependency(std::string nameA, std::string nameB, bool interfacialCoupling = false);

  // Create a dependence of Problem A equations on Problem B variables
  void createDependency(GenericEpetraProblem& problemA, GenericEpetraProblem& problemB, 
                        bool interfacialCoupling = false);

  // Register NLSolver parameter list
  void registerParameters(const Teuchos::RCP<Teuchos::ParameterList>& List);

  // Register Convergence Test(s)
  void registerStatusTest(const Teuchos::RCP<NOX::StatusTest::Combo>& comboTest);

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
  Teuchos::RCP<Epetra_CrsMatrix> getBlockJacobianMatrix(int i, int j = -1);

  // Get an individual block Jacobian inverse operator 
  Teuchos::RCP<Epetra_Operator> getBlockInverseOperator(int i);

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

  // Copy problem Jacobians as block diagonal contributions to 
  // specified composite matrix
  void copyProblemJacobiansToComposite( Epetra_CrsMatrix & mat );

  // Solve problem using fixed-point iteration
  bool solve();

  // Solve problems using matrix-free coupling
  bool solveMF();

  // Solve problems using schur-based coupling
  bool solveSchurBased();

  // Output summary of copled problems
  void outputStatus( std::ostream & os );

  // Create a filename for given problem, time-step, etc.
  std::string createIOname(GenericEpetraProblem &, int timeStep = 0);

  // Solve problems using matrix-free coupling
  void outputSolutions( const std::string outputDir, int timeStep = 0 );

  // Query whether or not offBlocks ar used
  bool useOffBlocks()
  { return doOffBlocks; }

  // Return my name
  std::string getName() const 
  { return GenericEpetraProblem::getName(); };

  // Return the problem name queries
  std::string getProblemName(int id) 
  { return Names[id]; };

  // These allow inheritance from GenericEpetraProblem base class

  bool evaluate(NOX::Epetra::Interface::Required::FillType type,
               const Epetra_Vector *solnVector,
               Epetra_Vector *rhsVector);

  // Inherited methods to allow use of library fixed-point
  virtual bool exchangeAllData();

  virtual bool exchangeDataTo(int solverId);

  // Inherited method to allow matrix replacement with BroydenOperator
  virtual Teuchos::RCP<const Epetra_CrsMatrix> 
  getReplacementValuesMatrix( const Epetra_Vector & x, NOX::Epetra::BroydenOperator::ReplacementInterface::FILL_TYPE );

  // method to generate only a diagonal block composite graph
  virtual Teuchos::RCP<Epetra_CrsGraph> 
  generateBlockDiagonalGraph( const Epetra_Map &, bool includeProblemOffContribs = false );

  // method to generate only off-diagonal block composite graph
  virtual Teuchos::RCP<Epetra_CrsGraph> 
  generateOffDiagonalBlockGraph( const Epetra_Map & );

  // method to combine (compose) two composite graphs into a new one
  virtual Teuchos::RCP<Epetra_CrsGraph> 
  composeGraphs( const Epetra_CrsGraph &, const Epetra_CrsGraph & );

  // Inherited from NOX::Epetra::SchurInterface

  // method to compute the action of a block on a vector using a directional derivative
  virtual bool applyBlockAction( int probId, int depId, const Epetra_Vector &, Epetra_Vector & );

  // method to compute the action of a block inverse on a vector
  virtual bool applyBlockInverseAction( int probId, int depId, const Epetra_Vector &, Epetra_Vector & );

  // Accessor for the communicator
  virtual const Epetra_Comm & getComm()
    { return *Comm; }

  //! Query availability of explicit operator
  virtual bool hasExplicitOperator( int rowBlock, int colBlock );

  //! Access explicit operator
  virtual Teuchos::RCP<Epetra_Operator> getExplicitOperator( int rowBlock, int colBlock );

  // method to create an inverse or its approximation to a diagonal block matrix
  virtual bool createBlockInverseOperator( int probId, Teuchos::ParameterList & );

protected:

  // Creates a composite graph from existing graphs of registered problems
  void generateGraph();

public:
  // Made these public for now
  Teuchos::RCP<Teuchos::ParameterList> nlParams; // NLSolver parameter list
  Teuchos::RCP<NOX::StatusTest::Combo> statusTest;       // Convergence criterion (ia)

private:

  int problemCount;                         // Number of registered Problems
  bool doOffBlocks;	// Flag to specify whether or not to compute
			// off-diagonal composite Jacobian blocks via FDC
  bool useMatlab;	// Flag to specify whether or not to use a Matlab interactive session

  Teuchos::RCP<Epetra_Map> compositeMap;	// Composite map for all registered problems
  Teuchos::RCP<Epetra_Vector> compositeSoln;	// Composite map for all registered problems
  map<int, Teuchos::RCP<GenericEpetraProblem> > Problems; // Problem container
  map<int, std::string> Names;                   // Problem names
  map<string, int> NameLookup;              // Needed to lookup problems by name
  map<int, Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> > LinearSystems;
  					    // Linear Systems container
  map<int, Teuchos::RCP<NOX::Epetra::Group> >Groups;     // Groups container
  map<int, Teuchos::RCP<Problem_Interface>  >Interfaces; // Interfaces container
  map<int, Teuchos::RCP<NOX::Solver::Generic> > Solvers; // Solvers container
  map<int, Teuchos::RCP<Epetra_IntVector> > ProblemToCompositeIndices;
    // Essentially a mapping from a problems indices to the composite indices

  Teuchos::RCP<Epetra_Operator> jacOperator   ; // The composite Jacobian operator
  Teuchos::RCP<Epetra_Operator> precOperator  ; // The composite Preconditioner operator
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterface;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> precInterface;
  Teuchos::RCP<NOX::Epetra::Group> compositeGroup  ; // The composite solution group
  Teuchos::RCP<NOX::Solver::Generic> compositeSolver ; // The composite solver

#ifdef HAVE_NOX_EPETRAEXT
  // Additional containers to allow each problem to use FD Coloring
  map<int, Teuchos::RCP<EpetraExt::CrsGraph_MapColoring> > TmpMapColorings;
  map<int, Teuchos::RCP<Epetra_MapColoring> > ColorMaps;
  map<int, Teuchos::RCP<EpetraExt::CrsGraph_MapColoringIndex> > ColorMapIndexSets;
  map<int, Teuchos::RCP<vector<Epetra_IntVector> > > ColumnsSets;
  map<int, Teuchos::RCP<Epetra_Operator> > MatrixOperators;

  // These containers allow off-diagonal block contributions representing
  // inter-problem coupling
  map<int, std::vector<OffBlock_Manager*> > OffBlock_Managers;    
#endif

  map<int, Teuchos::RCP<Epetra_Operator> > BlockInverseOperators; // inverses (or approximations) to diagonal blocks

};
#endif
