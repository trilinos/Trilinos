#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.H"

#include "FiniteElementProblem.H"

// Constructor - creates the Epetra objects (maps and vectors) 
FiniteElementProblem::FiniteElementProblem(int numGlobalElements, Epetra_Comm& comm) :
  Comm(&comm),
  NumGlobalElements(numGlobalElements)
{
  // Set the default parameter values 
  double factor = 1000.0;
  double leftBCValue = 0.0;
  double rightBCValue = 1.0;

  // Commonly used variables
  int i;
  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering
  StandardMap = new Epetra_Map(NumGlobalElements, 0, *Comm);

  // Get the number of elements owned by this processor
  NumMyElements = StandardMap->NumMyElements();

  // Construct an overlaped map for the finite element fill *************
  // For single processor jobs, the overlap and standard map are the same
  if (NumProc == 1) {
    OverlapMap = new Epetra_Map(*StandardMap);
  } else {

    int OverlapNumMyElements;
    int OverlapMinMyGID;
    OverlapNumMyElements = NumMyElements + 2;
    if ((MyPID == 0) || (MyPID == NumProc - 1)) 
      OverlapNumMyElements --;
    
    if (MyPID==0) 
      OverlapMinMyGID = StandardMap->MinMyGID();
    else 
      OverlapMinMyGID = StandardMap->MinMyGID() - 1;
    
    int* OverlapMyGlobalElements = new int[OverlapNumMyElements];
    
    for (i = 0; i < OverlapNumMyElements; i ++) 
      OverlapMyGlobalElements[i] = OverlapMinMyGID + i;
    
    OverlapMap = new Epetra_Map(-1, OverlapNumMyElements, 
			    OverlapMyGlobalElements, 0, *Comm);

    delete [] OverlapMyGlobalElements;

  } // End Overlap map construction *************************************

  // Construct Linear Objects  
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  initialSolution = new Epetra_Vector(*StandardMap);
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 5);

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  generateGraph(*AA);

  // Create a second matrix using graph of first matrix - this creates a 
  // static graph so we can refill the new matirx after TransformToLocal()
  // is called.
  A = new Epetra_CrsMatrix (Copy, *AA);
  A->TransformToLocal();
}

// Destructor
FiniteElementProblem::~FiniteElementProblem()
{
  delete AA;
  delete A;
  delete initialSolution;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

// Matrix and Residual Fills
bool FiniteElementProblem::evaluate(FillType f, 
				    const Epetra_Vector* soln, 
				    Epetra_Vector* tmp_rhs, 
				    Epetra_RowMatrix* tmp_matrix)
{
  flag = f;

  // Set the incoming linear objects
  if (flag == RHS_ONLY) {
    rhs = tmp_rhs;
  } else if (flag == MATRIX_ONLY) {
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else if (flag == ALL) { 
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else {
    cout << "ERROR: FiniteElementProblem::fillMatrix() - FillType flag is broken" << endl;
    throw;
  }

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector x(*OverlapMap);

  // Export Solution to Overlap vector
  u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int i,j,ierr;
  int OverlapNumMyElements = OverlapMap->NumMyElements();

  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;

  int row, column;
  double jac;
  double xx[2];
  double uu[2];
  Basis basis;

  // Create the nodal coordinates
  double Length=1.0;
  double dx=Length/((double) NumGlobalElements-1);
  for (i=0; i < OverlapNumMyElements; i++) {
    x[i]=dx*((double) OverlapMinMyGID+i);
  }
  
  // Zero out the objects that will be filled
  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == RHS_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
    
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes 
      xx[0]=x[ne];
      xx[1]=x[ne+1];
      uu[0]=u[ne];
      uu[1]=u[ne+1];
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu);
	            
      // Loop over Nodes in Element
      for (i=0; i< 2; i++) {
	row=OverlapMap->GID(ne+i);
	//printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
	//     MyPID, row, ne+i,StandardMap.MyGID(row));
	if (StandardMap->MyGID(row)) {
	  if ((flag == RHS_ONLY)    || (flag == ALL)) {
	    (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
	      +basis.wt*basis.dx
	      *((1.0/(basis.dx*basis.dx))*basis.duu*
		basis.dphide[i]+factor*basis.uu*basis.uu*basis.uu
		*basis.phi[i]);
	  }
	}
	// Loop over Trial Functions
	if ((flag == MATRIX_ONLY) || (flag == ALL)) {
	  for(j=0;j < 2; j++) {
	    if (StandardMap->MyGID(row)) {
	      column=OverlapMap->GID(ne+j);
	      jac=basis.wt*basis.dx*((1.0/(basis.dx*basis.dx))*
				     basis.dphide[j]*basis.dphide[i]
				     +3.0*factor*basis.uu*basis.uu*
				     basis.phi[j]*basis.phi[i]);  
	      ierr=A->SumIntoGlobalValues(row, 1, &jac, &column);
	    }
	  }
	}
      }
    }
  } 

  // Insert Boundary Conditions and modify Jacobian and function (RHS)
  // U(0)=1
  if (MyPID == 0) {
    if ((flag == RHS_ONLY)    || (flag == ALL))
      (*rhs)[0] = (*soln)[0] - leftBCValue;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      column=0;
      jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }
  
  // U(Npts) = bcValue
  if (MyPID == (NumProc - 1)) {
    if ((flag == RHS_ONLY)    || (flag == ALL))
      (*rhs)[NumMyElements-1] = (*soln)[NumMyElements-1] - rightBCValue;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      column=NumMyElements-1;
      jac=1.0;
      A->ReplaceGlobalValues((NumMyElements-1), 1, &jac, &column);
      column=NumMyElements-2;
      jac=0.0;
      A->ReplaceGlobalValues((NumMyElements-1), 1, &jac, &column);
    }
  }
  
  // Sync up processors to be safe
  Comm->Barrier();
 
  A->TransformToLocal();

  return true;
}

Epetra_Vector& FiniteElementProblem::getSolution()
{
  return *initialSolution;
}
  
Epetra_CrsMatrix& FiniteElementProblem::getJacobian()
{
  return *A;
}

Epetra_CrsGraph& FiniteElementProblem::generateGraph(Epetra_CrsGraph& AA)
{
  
  // Declare required variables
  int i,j;
  int row, column;
  int OverlapNumMyElements = OverlapMap->NumMyElements();
  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;
  
  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
          
    // Loop over Nodes in Element
    for (i=0; i< 2; i++) {
      row=OverlapMap->GID(ne+i);
      
      // Loop over Trial Functions
      for(j=0;j < 2; j++) {
	
	// If this row is owned by current processor, add the index
	if (StandardMap->MyGID(row)) {
	  column=OverlapMap->GID(ne+j);
	  AA.InsertGlobalIndices(row, 1, &column);
	}
      } 	
    }
  }
  AA.TransformToLocal();
  AA.SortIndices();
  AA.RemoveRedundantIndices();
  return AA;
}

bool FiniteElementProblem::setParameter(string name, double value)
{
  if (name == "Nonlinear Factor")
    factor = value;
  else if (name == "Left BC")
    leftBCValue = value;
  else if (name == "Right BC")
    rightBCValue = value;
  else {
    cout << "ERROR: FiniteElementProblem::setParameter() - the string "
	 << name << " is not a valid parameter name!"
	 << endl;
    //return false;
    throw "NOX Error";
  }

  return true;
} 
