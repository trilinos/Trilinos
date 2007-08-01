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
                                                                                
//  A pure virtual problem class to be used (via encapsulation) in
//  NOX_Epetra_Interface.H

#ifndef _NOX_GENERIC_PROBLEM_H
#define _NOX_GENERIC_PROBLEM_H

#include "NOX.H"
#include "NOX_Epetra.H"

// This is needed to get the FillType enum values from NOX_Epetra_Interface
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

//#include "Problem_Manager.H"
class Problem_Manager;

// Forward Declarations
class Epetra_Comm;
class Epetra_Map; 
class Epetra_Vector;
class Epetra_Import;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;

class XferOp;

// Flag to tell the evaluate routine what objects to fill
enum FillType {F_ONLY, MATRIX_ONLY, ALL}; 

// Finite Element Problem Class
class GenericEpetraProblem { 

  friend class Problem_Manager;

public:

  // Constructor
  GenericEpetraProblem(const Epetra_Comm& Comm, int NumGlobalElements = 0,
                       string name = "");

  // Destructor
  virtual ~GenericEpetraProblem();

  // Pure virtual methods 

  // Evaluates the function (RHS) and/or the Jacobian using the solution 
  // values in solnVector.
  virtual bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
                           const Epetra_Vector *solnVector, 
                           Epetra_Vector *rhsVector ) = 0;

  // Virtual methods with defaults provided

  // Allow a reset capability, eg for time dependent or parameter
  // continuation problems.
  virtual void reset() { /* default no-op */ }; 

  virtual void reset(const Epetra_Vector& x) { /* default no-op */ };

  // Set the unique int id for this problem
  virtual void setId(const int id_) { myId = id_; }; 

  // Query the unique int id for this problem
  virtual int getId() const { return myId; }; 

  // Query the name for this problem
  virtual string getName() const { return myName; }; 

  // Set the name for this problem
  virtual void setName(string name_) { myName = name_; }; 

  // Set the problem Manager with which this problem interacts 
  // (ie for transfers)
  virtual void setManager(Problem_Manager* manager_) 
    { myManager =  manager_; }; 

  // Add a problem to the list of that this one depends on
  virtual void addProblemDependence(const GenericEpetraProblem& problemB); 

  // Add a transfer operator to get fields from another problem
  virtual void addTransferOp(const GenericEpetraProblem& problemB); 

  // Signal whether the problem directly provides off-block structural info, 
  // ie global composite indices
  virtual bool suppliesOffBlockStructure()
  { return false; }

  // Signal whether the problem directly provides off-block structural info for a
  // particular dependent problem
  virtual bool suppliesOffBlockStructure(int depID)
  { return false; }

  // Default no-op method that allows a child to provide off-block structural info,
  virtual void getOffBlockIndices( map<int, vector<int> > & )
  { return; }

  // Do pre-transfer stuff, eg compute fluxes
  virtual void prepare_data_for_transfer()
  { /* default no-op */ }; 

  // Do transfers from each problem to this
  virtual void doTransfer(); 

  // Do post-transfer stuff
  virtual void process_transferred_data()
  { /* default no-op */ }; 

  // Set the solution vector (initialSolution)
  virtual void setSolution(const Epetra_Vector&); 

  // Create storage for dependent data received from other problems via
  // transfer operations
  virtual void createDependentVectors(); 

  virtual void initialize() { return; } //  default no-op

  // Create Standard and Overlap (ghost) maps
  virtual void createMaps(); 

  // Method needed to satisfy Problem_Interface.  
  virtual bool computePrecMatrix(const Epetra_Vector& solnVector, 
                               Epetra_RowMatrix& matrix);

  // Method needed to satisfy Problem_Interface.  
  virtual bool computePreconditioner(const Epetra_Vector& solnVector, 
                                   Epetra_Operator& precOperator);

  // Return a reference to the mesh for this problem
  virtual Epetra_Vector& getMesh();

  // Return a reference to the Epetra_Vector with the initial guess
  // that is generated by this class.
  virtual Teuchos::RCP<Epetra_Vector> getSolution();
  
  // Return a reference to the Epetra_CrsMatrix if one is used
  virtual Teuchos::RCP<Epetra_CrsMatrix> getJacobian();

  // Return a reference to the Epetra_CrsGraph if one is valid
  virtual Teuchos::RCP<Epetra_CrsGraph> & getGraph()
  { return AA; }

  // Output results 
  virtual void outputResults(const NOX::Solver::Generic& solver, 
			     Teuchos::ParameterList& printParams);

  // Output problem description
  virtual void outputStatus( ostream & os ) 
  { return; }

  // Output current solution state 
  virtual void outputSolutionStatus( ostream & os ) const;

  // Provide for transient inherited problems
  virtual void setdt( double dt );

  // Provide for transient inherited problems
  virtual double getdt() const;

  // Accessor for container of dependent problems
  virtual vector<int> & getDependentProblems()
  { return depProblems; }

protected:

  // inserts the global column indices into the Graph
  virtual void generateGraph() { };

protected:

  int                                           myId            ; // My unique integer id 
  string                                        myName          ; 
  Problem_Manager *                             myManager       ; // Problem Manager with which this problem interacts
  Epetra_Map *                                  OverlapMap      ;
  Epetra_Import *                               Importer        ;
  Teuchos::RCP<Epetra_Vector>           xptr            ; // Mesh data
  Teuchos::RCP<Epetra_Vector>           initialSolution ;
  Teuchos::RCP<Epetra_CrsGraph>         AA              ;
  Teuchos::RCP<Epetra_CrsMatrix>        A               ;
  vector<int>                                   depProblems     ; // List of problems that this one depends on
  map<string, int>                              nameToMyIndex   ;
  map<int, Epetra_Vector*>                      depSolutions    ; // Solutions from other Problems needed by this Problem
  map<int, XferOp*>                             xferOperators   ; // Container of transfer operators

public: // for now

  const Epetra_Comm *                           Comm            ;
  Epetra_Map *                                  StandardMap     ;
  int                                           MyPID           ; // Process number
  int                                           NumProc         ; // Total number of processes

  int                                           NumMyNodes      ; // Number of nodes (dofs) owned by this process for this problem
  int                                           NumGlobalNodes  ; // Total Number of nodes (dofs) for this problem

};
#endif
