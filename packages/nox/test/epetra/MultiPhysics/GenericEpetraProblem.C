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
#include "Xfer_Operator.H"

GenericEpetraProblem::GenericEpetraProblem(const Epetra_Comm& comm, 
                                           int numGlobalNodes,
                                           string name_) :
  Comm(&comm),
  NumGlobalNodes(numGlobalNodes),
  myId(0),
  myName(name_),
  StandardMap(0),
  OverlapMap(0),
  Importer(0),
  xptr(0),
  initialSolution(0),
//  auxSolution(0),
  AA(0),  
  A(0) 
{

  // Commonly used variables
  int i;
  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Construct data layout
  createMaps();

  // Construction and initialization of mesh and solution vectors are
  // left to the derived problem class.
}


void GenericEpetraProblem::createMaps()
{

  if (NumGlobalNodes > 0)
  { 
    // Construct a Source Map that puts approximately the same 
    // number of equations on each processor 
    
    // Begin by distributing nodes fairly equally
    StandardMap = new Epetra_Map(NumGlobalNodes, 0, *Comm);
  
    // Get the number of nodes owned by this processor
    NumMyNodes = StandardMap->NumMyElements();
  
    // Construct an overlap node map for the finite element fill
    // For single processor jobs, the overlap and standard maps are the same
    if (NumProc == 1) { 
      OverlapMap = new Epetra_Map(*StandardMap);
    } else {
      
      int OverlapNumMyNodes;
      int OverlapMinMyNodeGID;
      OverlapNumMyNodes = NumMyNodes + 2;
      if ((MyPID == 0) || (MyPID == NumProc - 1))
        OverlapNumMyNodes --;
  
      if (MyPID==0)
        OverlapMinMyNodeGID = StandardMap->MinMyGID();
      else
        OverlapMinMyNodeGID = StandardMap->MinMyGID() - 1;
      
      int* OverlapMyGlobalNodes = new int[OverlapNumMyNodes];
      
      for (int i = 0; i < OverlapNumMyNodes; i ++)
        OverlapMyGlobalNodes[i] = OverlapMinMyNodeGID + i;
      
      OverlapMap = new Epetra_Map(-1, OverlapNumMyNodes,
                              OverlapMyGlobalNodes, 0, *Comm);
  
      delete [] OverlapMyGlobalNodes;
  
    } // End Overlap node map construction ********************************
  
    Importer = new Epetra_Import(*OverlapMap, *StandardMap);

#ifdef DEBUG
    // Output to check progress so far
    printf("NumMyNodes, NumGlobalNodes --> %d\t%d\n",NumMyNodes, 
		                                     NumGlobalNodes);
    cout << *StandardMap << endl;
    cout << *OverlapMap << endl;
    Importer->Print(cout);
#endif
  }

  return;
}


// Destructor
GenericEpetraProblem::~GenericEpetraProblem()
{
  delete A; A = 0;
  delete AA; AA = 0;
  delete initialSolution; initialSolution = 0;
  //  Need to fix this !!!  RHooper
//  delete auxSolution; auxSolution = 0;
  delete Importer; Importer = 0;
  delete OverlapMap; OverlapMap = 0;
  delete StandardMap; StandardMap = 0;
}

void GenericEpetraProblem::outputResults(NOX::Solver::Manager& solver, 
                   NOX::Parameter::List& printParams)
{
  // Output the parameter list
  NOX::Utils utils(printParams);
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getParameterList().print(cout);
    cout << endl;
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>
        (finalGroup.getX())).getEpetraVector();

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = finalSolution.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",finalSolution.Map().Comm().MyPID());
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp,"%d  %E\n",finalSolution.Map().MinMyGID()+i,finalSolution[i]);
  fclose(ifp);
}

Epetra_Vector& GenericEpetraProblem::getMesh()
{
  assert( xptr != 0 ); // Mesh vector had better exist
  return *xptr;
}

Epetra_Vector& GenericEpetraProblem::getSolution()
{
  assert( initialSolution != 0 ); // Solution vector had better exist
  return *initialSolution;
}

Epetra_CrsGraph& GenericEpetraProblem::getGraph()
{
  if(AA)
    return *AA;
  else
  {
    cout << "ERROR: No valid Matrix Graph exists for this problem !!" << endl;
    return *AA;
  }
}

Epetra_CrsMatrix& GenericEpetraProblem::getJacobian()
{
  if(A)
    return *A;
  else
  {
    cout << "ERROR: No valid Jacobian exists for this problem !!" << endl;
    return *A;
  }
}

void GenericEpetraProblem::setSolution(const Epetra_Vector& data)
{
  // Ensure that the derived problem class created the solution vector
  if(!initialSolution)
  {
    cout << "ERROR: No solution vector exists for this problem !!" << endl;
    throw "GenericEpetraProblem ERROR";
  }

  *initialSolution = data;
}

void GenericEpetraProblem::createAuxillaryVectors()
{
  // Create the auxillary vectors needed to receive data from other problems
  if( !initialSolution ) {
    cout << "ERROR: Cannot create auxillary data without an existing solution "
         << "vector for this problem !!" << endl;
    throw "GenericEpetraProblem ERROR";
  }
  for( int i = 0; i<auxProblems.size(); i++ ) {
#ifdef DEBUG
    cout << "For problem : " << myId << "  Creating AuxVec for problem : "
         << auxProblems[i] << endl;
#endif
    auxSolutions.insert( pair<int, Epetra_Vector*>(auxProblems[i],
			    new Epetra_Vector(*initialSolution)) );
  }
}

void GenericEpetraProblem::addProblemDependence(
		const GenericEpetraProblem& problemB)
{
  // Add a problem to the list of those this one depends on
  auxProblems.push_back(problemB.getId());
}

void GenericEpetraProblem::addTransferOp(const GenericEpetraProblem& problemB)
{
  // Add a transfer operator to get fields from another problem
  xferOperators.insert(pair<int, XferOp*>(problemB.getId(), 
			  new XferOp(*this, problemB)));
}

void GenericEpetraProblem::doTransfer()
{
  // Do transfers from each dependent problem to this one
  for( int i = 0; i<auxProblems.size(); i++) {
    int auxId = auxProblems[i];
    XferOp* xfer = xferOperators.find(auxId)->second;

    if( !xfer ) {
      cout << "ERROR: doTransfer: No valid transfer operator !!" << endl;
      throw "GenericEpetraProblem ERROR";
    }

    else {
      // NOTE that we are transferring (by default) to/from each problem's
      // solution vector which may not be the same as in each respective
      // group.
      Epetra_Vector& fromVec = myManager->getProblem(auxId).getSolution();
      Epetra_Vector& toVec = *(auxSolutions.find(auxId)->second);
      xfer->transferField(toVec, fromVec); 
    }  
  }
}

bool GenericEpetraProblem::computePrecMatrix(const Epetra_Vector& solnVector,
                               Epetra_RowMatrix& matrix)
{
  cout << "WARNING: computePrecMatrix not implemented for this problem !!"
       << endl;

  return false;
}
    
bool GenericEpetraProblem::computePreconditioner(
                                   const Epetra_Vector& solnVector,
                                   Epetra_Operator& precOperator)
{
  cout << "WARNING: computePreconditioner not implemented for this problem !!"
       << endl;

  return false;
}

