#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.H"

#include "Brusselator.H"

// Constructor - creates the Epetra objects (maps and vectors) 
Brusselator::Brusselator(int numGlobalNodes, Epetra_Comm& comm) :
  xmin(0.0),
  xmax(1.0),
  dt(1.0e-01),
  Comm(&comm),
  NumGlobalNodes(numGlobalNodes)
{

  // Commonly used variables
  int i;
  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Here we assume a 2-species Brusselator model, ie 2 dofs per node
  NumSpecies = 2;
  NumGlobalUnknowns = NumSpecies * NumGlobalNodes;

  // Construct a Source Map that puts approximately the same 
  // number of equations on each processor 

  // Begin by distributing nodes fairly equally
  StandardNodeMap = new Epetra_Map(NumGlobalNodes, 0, *Comm);

  // Get the number of nodes owned by this processor
  NumMyNodes = StandardNodeMap->NumMyElements();

  // Construct an overlap node map for the finite element fill
  // For single processor jobs, the overlap and standard maps are the same
  if (NumProc == 1) {
    OverlapNodeMap = new Epetra_Map(*StandardNodeMap);
  } else {

    int OverlapNumMyNodes;
    int OverlapMinMyNodeGID;
    OverlapNumMyNodes = NumMyNodes + 2;
    if ((MyPID == 0) || (MyPID == NumProc - 1)) 
      OverlapNumMyNodes --;
    
    if (MyPID==0) 
      OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();
    else 
      OverlapMinMyNodeGID = StandardNodeMap->MinMyGID() - 1;
    
    int* OverlapMyGlobalNodes = new int[OverlapNumMyNodes];
    
    for (i = 0; i < OverlapNumMyNodes; i ++) 
      OverlapMyGlobalNodes[i] = OverlapMinMyNodeGID + i;
    
    OverlapNodeMap = new Epetra_Map(-1, OverlapNumMyNodes, 
			    OverlapMyGlobalNodes, 0, *Comm);

    delete [] OverlapMyGlobalNodes;

  } // End Overlap node map construction ********************************

  // Now create the unknowns maps from the node maps
  NumMyUnknowns = NumSpecies * NumMyNodes;
  int* StandardMyGlobalUnknowns = new int[NumMyUnknowns];
  for (int k=0; k<NumSpecies; k++) 
    for (i=0; i<NumMyNodes; i++)

      // For now, we employ an interleave of unknowns

      StandardMyGlobalUnknowns[ NumSpecies * i + k ] = 
               NumSpecies * StandardNodeMap->GID(i) + k;

  StandardMap = new Epetra_Map(-1, NumMyUnknowns, StandardMyGlobalUnknowns,
                                 0, *Comm);
  delete [] StandardMyGlobalUnknowns;

  assert(StandardMap->NumGlobalElements() == NumGlobalUnknowns);

  if (NumProc == 1) {
    OverlapMap = new Epetra_Map(*StandardMap);
  } 
  else {
    int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();
    int OverlapNumMyUnknowns = NumSpecies * OverlapNumMyNodes;
    int* OverlapMyGlobalUnknowns = new int[OverlapNumMyUnknowns];
    for (int k=0; k<NumSpecies; k++) 
      for (i=0; i<OverlapNumMyNodes; i++)
        OverlapMyGlobalUnknowns[ NumSpecies * i + k ] = 
                 NumSpecies * OverlapNodeMap->GID(i) + k;

    OverlapMap = new Epetra_Map(-1, OverlapNumMyUnknowns, 
                                OverlapMyGlobalUnknowns, 0, *Comm);
    delete [] OverlapMyGlobalUnknowns;

  } // End Overlap unknowns map construction ***************************


#ifdef DEBUG
  // Output to check progress so far
  printf("NumSpecies, NumMyNodes, NumGlobalNodes, NumMyUnknowns, NumGlobalUnknowns --> %d\t%d\t%d\t%d\n",NumSpecies, NumMyNodes, NumGlobalNodes, NumMyUnknowns, NumGlobalUnknowns);

  cout << *StandardNodeMap << endl;
  cout << *OverlapNodeMap << endl;
  cout << *StandardMap << endl;
  cout << *OverlapMap << endl;
#endif

  // Construct Linear Objects  
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  nodeImporter = new Epetra_Import(*OverlapNodeMap, *StandardNodeMap);
  initialSolution = new Epetra_Vector(*StandardMap);
  oldSolution = new Epetra_Vector(*StandardMap);
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 0);

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  generateGraph(*AA);

#ifdef DEBUG
  AA->Print(cout);
  Importer->Print(cout);
#endif

  // Create a second matrix using graph of first matrix - this creates a 
  // static graph so we can refill the new matirx after TransformToLocal()
  // is called.
  A = new Epetra_CrsMatrix (Copy, *AA);
  A->TransformToLocal();

  // Create the nodal coordinates
  xptr = new Epetra_Vector(*StandardNodeMap);
  double Length= xmax - xmin;
  dx=Length/((double) NumGlobalNodes-1);
  for (i=0; i < NumMyNodes; i++) {
    (*xptr)[i]=xmin + dx*((double) StandardNodeMap->MinMyGID()+i);
  }

  initializeSoln();

}

// Destructor
Brusselator::~Brusselator()
{
  delete A;
  delete AA;
  delete xptr;
  delete initialSolution;
  delete oldSolution;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

// Reset function
void Brusselator::reset(const Epetra_Vector& x)
{
  *oldSolution = x;
}

// Set initialSolution to desired initial condition
void Brusselator::initializeSoln()
{
  Epetra_Vector& soln = *initialSolution;
  Epetra_Vector& x = *xptr;

  // Here we do a sinusoidal perturbation of the unstable
  // steady state.

  double pi = 4.*atan(1.0);

  for (int i=0; i<x.MyLength(); i++) {
    soln[2*i] = 0.6 + 1.e-1*sin(1.0*pi*x[i]);
    soln[2*i+1] = 10.0/3.0 + 1.e-1*sin(1.0*pi*x[i]);
  }
  
  *oldSolution = soln;
} 

// Matrix and Residual Fills
bool Brusselator::evaluate(FillType f, 
				    const Epetra_Vector* soln, 
				    Epetra_Vector* tmp_rhs, 
				    Epetra_RowMatrix* tmp_matrix)
{
  flag = f;

  // Set the incoming linear objects
  if (flag == F_ONLY) {
    rhs = tmp_rhs;
  } else if (flag == MATRIX_ONLY) {
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else if (flag == ALL) { 
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else {
    cout << "ERROR: Brusselator::fillMatrix() - FillType flag is broken" << endl;
    throw;
  }

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  Epetra_Vector xvec(*OverlapNodeMap);

  // Export Solution to Overlap vector
  u.Import(*soln, *Importer, Insert);
  uold.Import(*oldSolution, *Importer, Insert);
  xvec.Import(*xptr, *nodeImporter, Insert);

  // Declare required variables
  int i,j,ierr;
  int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();
  int OverlapNumMyUnknowns = OverlapMap->NumMyElements();

  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardNodeMap->MinMyGID()-1;

  int row1, row2, column1, column2;
  double factor1=1000.0;
  double factor2=1000.0;
  double Dcoeff1 = 0.025;
  double Dcoeff2 = 0.025;
  double alpha = 0.6;
  double beta = 2.0;
  double jac11, jac12, jac21, jac22;
  double xx[2];
  double uu[2*NumSpecies];
  double uuold[2*NumSpecies];
  Basis basis(NumSpecies);

  
  // Zero out the objects that will be filled
  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == F_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {
    
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes 
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      for (int k=0; k<NumSpecies; k++) {
        uu[NumSpecies * 0 + k] = u[NumSpecies * ne + k];
        uu[NumSpecies * 1 + k] = u[NumSpecies * (ne+1) + k];
        uuold[NumSpecies * 0 + k] = uold[NumSpecies * ne + k];
        uuold[NumSpecies * 1 + k] = uold[NumSpecies * (ne+1) + k];
      }
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu, uuold);
	            
      // Loop over Nodes in Element
      for (i=0; i< 2; i++) {
	row1=OverlapMap->GID(NumSpecies * (ne+i));
	row2=OverlapMap->GID(NumSpecies * (ne+i) + 1);
	//printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
	//     MyPID, row, ne+i,StandardMap.MyGID(row));
	if (StandardMap->MyGID(row1)) {
	  if ((flag == F_ONLY)    || (flag == ALL)) {
	    (*rhs)[StandardMap->LID(OverlapMap->GID(NumSpecies*(ne+i)))]+=
	      +basis.wt*basis.dx
	      *((basis.uu[0] - basis.uuold[0])/dt * basis.phi[i] 
              +(1.0/(basis.dx*basis.dx))*Dcoeff1*basis.duu[0]*basis.dphide[i]
              + basis.phi[i] * ( -alpha + (beta+1.0)*basis.uu[0]
                - basis.uu[0]*basis.uu[0]*basis.uu[1]) );
	    (*rhs)[StandardMap->LID(OverlapMap->GID(NumSpecies*(ne+i)+1))]+=
	      +basis.wt*basis.dx
	      *((basis.uu[1] - basis.uuold[1])/dt * basis.phi[i] 
	      +(1.0/(basis.dx*basis.dx))*Dcoeff2*basis.duu[1]*basis.dphide[i]
              + basis.phi[i] * ( -beta*basis.uu[0]
                 + basis.uu[0]*basis.uu[0]*basis.uu[1]) );
	  }
	}
	// Loop over Trial Functions
	if ((flag == MATRIX_ONLY) || (flag == ALL)) {
	  for(j=0;j < 2; j++) {
	    if (StandardMap->MyGID(row1)) {
	      column1=OverlapMap->GID(NumSpecies * (ne+j));
	      column2=OverlapMap->GID(NumSpecies * (ne+j) + 1);
	      jac11=basis.wt*basis.dx*(
                      basis.phi[j]/dt*basis.phi[i] 
                      +(1.0/(basis.dx*basis.dx))*Dcoeff1*basis.dphide[j]*
                                                        basis.dphide[i]
                      + basis.phi[i] * ( (beta+1.0)*basis.phi[j]
                      - 2.0*basis.uu[0]*basis.phi[j]*basis.uu[1]) );  
	      jac12=basis.wt*basis.dx*(
                      basis.phi[i] * ( -basis.uu[0]*basis.uu[0]*basis.phi[j]) );
	      jac21=basis.wt*basis.dx*(
                      basis.phi[i] * ( -beta*basis.phi[j] 
                      + 2.0*basis.uu[0]*basis.phi[j]*basis.uu[1]) );  
	      jac22=basis.wt*basis.dx*(
                      basis.phi[j]/dt*basis.phi[i] 
                      +(1.0/(basis.dx*basis.dx))*Dcoeff2*basis.dphide[j]*
                                                        basis.dphide[i]
                      + basis.phi[i] * basis.uu[0]*basis.uu[0]*basis.phi[j] );
	      ierr=A->SumIntoGlobalValues(row1, 1, &jac11, &column1);
              column1++;
	      ierr=A->SumIntoGlobalValues(row1, 1, &jac12, &column1);
	      ierr=A->SumIntoGlobalValues(row2, 1, &jac22, &column2);
              column2--;
	      ierr=A->SumIntoGlobalValues(row2, 1, &jac21, &column2);
	    }
	  }
	}
      }
    }
  } 

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(0)=1
  if (MyPID==0) {
    if ((flag == F_ONLY)    || (flag == ALL)) {
      (*rhs)[0]= (*soln)[0] - 0.6;
      (*rhs)[1]= (*soln)[1] - 10.0/3.0;
    }
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      int column=0;
      double jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(1, 1, &jac, &column);
      jac=0.0;
      column=0;
      A->ReplaceGlobalValues(1, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      A->ReplaceGlobalValues(1, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      A->ReplaceGlobalValues(1, 1, &jac, &column);
    }
  }
  // U(1)=1
  if ( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 ) {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    if ((flag == F_ONLY)    || (flag == ALL)) {
      (*rhs)[lastDof - 1] = (*soln)[lastDof - 1] - 0.6;
      (*rhs)[lastDof] = (*soln)[lastDof] - 10.0/3.0;
    }
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      int row=StandardMap->MaxAllGID() - 1;
      int column = row;
      double jac = 1.0;
      A->ReplaceGlobalValues(row++, 1, &jac, &column);
      column++;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      jac=0.0;
      row = column - 1;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      column--;
      A->ReplaceGlobalValues(row+1, 1, &jac, &column);
      column--;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      A->ReplaceGlobalValues(row+1, 1, &jac, &column);
      column--;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      A->ReplaceGlobalValues(row+1, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();
 
  A->TransformToLocal();

#ifdef DEBUG
  A->Print(cout);
  cin.get();
#endif

/* 
  if( (flag == F_ONLY) || (flag == ALL) )
    cout << "For residual fill :" << endl << *rhs << endl;

  if( (flag == MATRIX_ONLY) || (flag == ALL) ) {
    cout << "For jacobian fill :" << endl;
    A->Print(cout);
  }

  cin.get();
*/


  return true;
}

Epetra_Vector& Brusselator::getSolution()
{
  return *initialSolution;
}
  
Epetra_CrsMatrix& Brusselator::getJacobian()
{
  return *A;
}

Epetra_Vector& Brusselator::getMesh()
{ 
  return *xptr;
}

Epetra_Vector& Brusselator::getOldSoln()
{
  return *oldSolution;
} 
  
double Brusselator::getdt()
{
  return dt;
}

Epetra_CrsGraph& Brusselator::generateGraph(Epetra_CrsGraph& AA)
{
  
  // Declare required variables
  int i,j;
  int row, column;
  int OverlapNumMyNodes = OverlapNodeMap->NumMyElements();
  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardNodeMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardNodeMap->MinMyGID()-1;
  
  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {
          
    // Loop over Nodes in Element
    for (i=0; i<2; i++) {

      // If this node is owned by current processor, add indices
      if (StandardNodeMap->MyGID(OverlapNodeMap->GID(ne+i))) {

        // Loop over unknowns in Node
        for (int k=0; k<NumSpecies; k++) {
          row=OverlapMap->GID( NumSpecies*(ne+i) + k); // Interleave scheme

          // Loop over supporting nodes
          for(j=0; j<2; j++) {
	
            // Loop over unknowns at supporting nodes
            for (int m=0; m<NumSpecies; m++) {
	      column=OverlapMap->GID( NumSpecies*(ne+j) + m);
              //printf("\t\tWould like to insert -> (%d, %d)\n",row,column);
	      AA.InsertGlobalIndices(row, 1, &column);
            }
	  }
        } 	
      }
    }
  }
  AA.TransformToLocal();
  AA.SortIndices();
  AA.RemoveRedundantIndices();
  return AA;
}
