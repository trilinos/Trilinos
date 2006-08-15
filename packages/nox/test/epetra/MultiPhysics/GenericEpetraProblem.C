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
  AA(0)
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

//-----------------------------------------------------------------------------

// Destructor
GenericEpetraProblem::~GenericEpetraProblem()
{
  delete AA; AA = 0;
  delete Importer; Importer = 0;
  delete OverlapMap; OverlapMap = 0;
  delete StandardMap; StandardMap = 0;
  
  for( map<int, Epetra_Vector*>::iterator iter = depSolutions.begin();
         iter != depSolutions.end(); delete (*iter).second, ++iter );
  for( map<int, XferOp*>::iterator iter = xferOperators.begin();
         iter != xferOperators.end(); delete (*iter).second, ++iter );

}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::createMaps()
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

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::outputResults(NOX::Solver::Manager& solver, 
                   Teuchos::ParameterList& printParams)
{
  // Output the parameter list
  NOX::Utils utils(printParams);
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getList().print(cout);
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

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::outputSolutionStatus( ostream & os ) const
{
  // Output all solution vectors
  os << "\nProblem Status for " << myName << endl;
  os << "  ----------------------------------" << endl;
  os << *initialSolution << endl;

  for( int i = 0; i < depProblems.size(); ++i) {
    int depId = depProblems[i];
    os << "\tDependent Status for " << myManager->getProblem(depId).getName() << endl;
    os << "\t----------------------------------" << endl;
    Epetra_Vector & depVec = *( (*(depSolutions.find(depId))).second );
    os << depVec << endl;
  }
}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::setdt( double dt )
{
  cout << "No-op : Implement time dependence in inherited problem !!" << endl;
}

//-----------------------------------------------------------------------------

double 
GenericEpetraProblem::getdt() const
{
  cout << "No-op : Implement time dependence in inherited problem !!" << endl;
  return 0.0;
}

//-----------------------------------------------------------------------------

Epetra_Vector & 
GenericEpetraProblem::getMesh()
{
  assert( xptr != 0 ); // Mesh vector had better exist
  return *xptr;
}

//-----------------------------------------------------------------------------

Teuchos::RefCountPtr<Epetra_Vector> 
GenericEpetraProblem::getSolution()
{
  return initialSolution;
}

//-----------------------------------------------------------------------------

Epetra_CrsGraph & 
GenericEpetraProblem::getGraph()
{
  if(AA)
    return *AA;
  else
  {
    cout << "ERROR: No valid Matrix Graph exists for this problem !!" << endl;
    return *AA;
  }
}

//-----------------------------------------------------------------------------

Teuchos::RefCountPtr<Epetra_CrsMatrix> 
GenericEpetraProblem::getJacobian()
{
  if(A.get())
    return A;
  else
  {
    cout << "ERROR: No valid Jacobian exists for this problem !!" << endl;
    return A;
  }
}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::setSolution(const Epetra_Vector& data)
{
  // Ensure that the derived problem class created the solution vector
  if(!initialSolution.get())
  {
    cout << "ERROR: No solution vector exists for this problem !!" << endl;
    throw "GenericEpetraProblem ERROR";
  }

  (*initialSolution.get()) = data;
}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::createDependentVectors()
{
  // Create the dependent vectors needed to receive data from other problems
  if( !initialSolution.get() ) 
  {
    cout << "ERROR: Cannot create dependent data without an existing solution "
         << "vector for this problem !!" << endl;
    throw "GenericEpetraProblem ERROR";
  }
  for( int i = 0; i<depProblems.size(); i++ ) 
  {
#ifdef DEBUG
    cout << "For problem : " << myId << "  Creating DepVec for problem : "
         << depProblems[i] << endl;
#endif
    depSolutions[ depProblems[i] ] = new Epetra_Vector(*initialSolution);
  }
}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::addProblemDependence( const GenericEpetraProblem& problemB )
{
  // Add a problem to the list of those this one depends on
  depProblems.push_back(problemB.getId());

  // Add to the Name-to-My-Index lookup map
  nameToMyIndex.insert( pair<string, int>(problemB.getName(),
                          depProblems.size() - 1) );
}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::addTransferOp( const GenericEpetraProblem& problemB )
{
  // Add a transfer operator to get fields from another problem
  xferOperators[problemB.getId()] = new XferOp(*this, problemB);
}

//-----------------------------------------------------------------------------

void 
GenericEpetraProblem::doTransfer()
{
  // Do transfers from each dependent problem to this one
  for( int i = 0; i < depProblems.size(); ++i)  
  {
    int depId = depProblems[i];
    XferOp* xfer = (*xferOperators.find(depId)).second;

    if( !xfer ) 
    {
      cout << "ERROR: doTransfer: No valid transfer operator !!" << endl;
      throw "GenericEpetraProblem ERROR";
    }
    else 
    {
      // NOTE that we are transferring (by default) to/from each problem's
      // solution vector which may be different data than in each respective
      // group.
      Epetra_Vector & fromVec = *(myManager->getProblem(depId).getSolution());
      Epetra_Vector & toVec = *( (*(depSolutions.find(depId))).second );
      xfer->transferField( toVec, fromVec ); 
    }  
  }
}

//-----------------------------------------------------------------------------

bool 
GenericEpetraProblem::computePrecMatrix( const Epetra_Vector& solnVector, Epetra_RowMatrix& matrix )
{
  cout << "WARNING: computePrecMatrix not implemented for this problem !!"
       << endl;

  return false;
}
    
//-----------------------------------------------------------------------------

bool 
GenericEpetraProblem::computePreconditioner( const Epetra_Vector& solnVector, Epetra_Operator& precOperator )
{
  cout << "WARNING: computePreconditioner not implemented for this problem !!"
       << endl;

  return false;
}

//-----------------------------------------------------------------------------

