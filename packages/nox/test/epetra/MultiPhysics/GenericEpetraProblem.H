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
                                                                                
//  A pure virtual problem class to be used (via encapsulation) in
//  NOX_Epetra_Interface.H

#ifndef _NOX_GENERIC_PROBLEM_H
#define _NOX_GENERIC_PROBLEM_H

#include "NOX.H"
#include "NOX_Epetra.H"

// This is needed to get the FillType enum values from NOX_Epetra_Interface
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"

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

// Flag to tell how to order the unknowns
enum StoreUnknowns {Interleaved, Staggered}; // Not currently used

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
  virtual bool evaluate( NOX::EpetraNew::Interface::Required::FillType fill,
                           const Epetra_Vector *solnVector, 
                           Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix)
                           = 0;

  // Virtual methods with defaults provided

  // Allow a reset capability, eg for time dependent or parameter
  // continuation problems.
  virtual void reset() {}; 

  virtual void reset(const Epetra_Vector& x) {};

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

  // Do transfers from each problem to this
  virtual void doTransfer(); 

  // Set the solution vector (initialSolution)
  virtual void setSolution(const Epetra_Vector&); 

  // Create storage for auxillary data received from other problems via
  // transfer operations
  virtual void createAuxillaryVectors(); 

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
  virtual Epetra_Vector& getSolution();
  
  // Return a reference to the Epetra_CrsMatrix if one is used
  virtual Epetra_CrsMatrix& getJacobian();

  // Return a reference to the Epetra_CrsGraph if one is valid
  virtual Epetra_CrsGraph& getGraph();

  // Output results 
  virtual void outputResults(NOX::Solver::Manager& solver, 
                             NOX::Parameter::List& printParams);

protected:

  // inserts the global column indices into the Graph
  virtual void generateGraph() {};

protected:

  int myId;                     // Unique integer id for this Problem
  string myName;		// String identifier for convenience
  Problem_Manager* myManager;   // Problem Manager with which this problem interacts
  //Epetra_Map *StandardMap;
  Epetra_Map *OverlapMap;
  Epetra_Import *Importer;
  Epetra_Vector *xptr;		// Spatial locations of mesh points
  Epetra_Vector *initialSolution;
  Epetra_CrsGraph *AA;
  Epetra_CrsMatrix *A;
  vector<int> auxProblems;      // List of problems that this one depends on
  map<int, Epetra_Vector*> auxSolutions; // Auxillary data; solutions from 
                                  // other Problems needed by this Problem
  map<int, XferOp*> xferOperators;   // Container of transfer operators

public: // for now

  const Epetra_Comm *Comm;
  Epetra_Map *StandardMap;
  int MyPID;              // Process number
  int NumProc;            // Total number of processes

  int NumMyNodes;         // Number of nodes (dofs) owned by this process
                          // for this problem
  int NumGlobalNodes;     // Total Number of nodes (dofs) for this problem
};
#endif
