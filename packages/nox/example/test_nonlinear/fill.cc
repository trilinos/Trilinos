#include<iostream.h>
#include<math.h>
#include "fill.h"
#include "basis.h"

void Fill::registerFillObjects(
		      Epetra_Map &STANDARDMAP, 
		      Epetra_Map &OVERLAPMAP,
		      Epetra_Import & IMPORTER,
		      Epetra_Comm &COMM) 
{
  StandardMap=&STANDARDMAP; 
  OverlapMap=&OVERLAPMAP;
  Importer=&IMPORTER;
  Comm=&COMM;
}

// Matrix and Residual Fills
void Fill::fillMatrix(Epetra_Vector *tmp_soln, Epetra_Vector *tmp_rhs, 
		      Epetra_RowMatrix *tmp_matrix) 
{

  // Decide what kind of fill call this will be
  if (tmp_matrix==NULL) {
    flag = RHS_ONLY;
    soln = tmp_soln;
    rhs = tmp_rhs;
    //if (Comm->MyPID()==0) printf("Starting RESIDUAL ONLY Fill....");
  } else if (tmp_rhs==NULL) {
    flag = MATRIX_ONLY;
    soln = tmp_soln;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
    //if (Comm->MyPID()==0) printf("Starting MATRIX ONLY Fill....");
  } else { 
    flag = ALL;
    soln = tmp_soln;
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
    //if (Comm->MyPID()==0) printf("Starting Matrix and RESIDUAL Fill....");
  }

  // Create Linear Objects for Fill
  Epetra_Vector& u = *new Epetra_Vector(*OverlapMap);
  Epetra_Vector& x = *new Epetra_Vector(*OverlapMap);

  int i,j,ierr;  
  int MyPID = Comm->MyPID();
  int NumProc = Comm->NumProc();
  int NumMyElements = StandardMap->NumMyElements();
  int NumGlobalElements = StandardMap->NumGlobalElements();
  int OverlapNumMyElements = OverlapMap->NumMyElements();
  int * StandardMyGlobalElements = new int[NumMyElements];
  StandardMap->MyGlobalElements(StandardMyGlobalElements);

  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;

  // Apply Essential BC to initial solution vector
  if (MyPID==0) {
    (*soln)[0]=1.0;
  }

  // Export Solution to Overlap vector
  Epetra_Import & Importer2 = *new Epetra_Import(*OverlapMap, *StandardMap);
  assert(u.Import(*soln, Importer2, Insert)==0);
  delete &Importer2;

  /*
  cout << "Soln Vector" << soln[0] << endl;
  Comm->Barrier();
  cout << "Overlapped Vector" << (&u)[0] << endl;
  */
   
  int row;
  double eta;
  //double factor=20.0;
  double factor=1000.0;
  double *xx = new double[2];
  double *uu = new double[2];
  double jac;
  int column;
  int *indicies = new int[3];
  double *RowValues = new double[OverlapMap->NumGlobalElements()];
  Basis basis;

  // Create the nodal position variables
  double Length=1.0;
  double dx=Length/((double) NumGlobalElements-1);

  for (i=0; i < OverlapNumMyElements; i++) {
    x[i]=dx*((double) OverlapMinMyGID+i);
  }
  
  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == RHS_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
    
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      xx[0]=x[ne];
      xx[1]=x[ne+1];
      uu[0]=u[ne];
      uu[1]=u[ne+1];
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
		basis.dphide[i]+factor*basis.uu*basis.uu*basis.phi[i]);
	  }
	}
	// Loop over Trial Functions
	if ((flag == MATRIX_ONLY) || (flag == ALL)) {
	for(j=0;j < 2; j++) {
	  column=OverlapMap->GID(ne+j);
	  jac=basis.wt*basis.dx*((1.0/(basis.dx*basis.dx))*basis.dphide[j]*
	      basis.dphide[i]+2.0*factor*basis.uu*basis.phi[j]*basis.phi[i]);  
	  ierr=A->SumIntoGlobalValues(row, 1, &jac, &column);
	  if (ierr!=0) {	    
	    //printf("SumInto failed at (%d,%d)!!\n",row,column);
	    ierr=A->InsertGlobalValues(row, 1, &jac, &column);
	 //if (ierr==0) printf("Insert SUCCEEDED at (%d,%d)!!\n",row,column);
	    //else printf("Insert FAILED at (%d,%d)!!\n",row,column);
	  } //else if (ierr==0) 
	    //printf("SumInto SUCCEEDED at (%d,%d)!!\n",row,column);
	}
	}
      }
    }
  } 

  // Insert Boundary Conditions
  // U(0)=1.0
  if (MyPID==0) {
    (*soln)[0]=1.0;
    if ((flag == RHS_ONLY)    || (flag == ALL)) (*rhs)[0]=0.0;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      column=0;
      jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }

  Comm->Barrier();
 
  //Epetra_Import & Importer3 = *new Epetra_Import(StandardMap, OverlapMap);
  //assert(rhs.Import(rhs1, Importer3, Add)==0);
  //assert(A.Import(A1, Importer3, Add)==0);
  //delete &Importer3;

  A->TransformToLocal();
  /*
  // Print RHS   
  printf("Proc=%d  #Elements=%d\n",MyPID,NumMyElements);
  for (i=0; i< NumMyElements; i++) {
    printf("Proc=%d, GID=%d, rhs=%e\n",MyPID ,
	   StandardMap->GID(i),(*rhs)[i]);
  }
  */
  /*
  // Print Matrix
  int StandardNumMyRows = A->NumMyRows();
  int StandardNumEntries; 
  int * StandardIndices; 
  double * StandardValues;
  for (i=0; i< StandardNumMyRows; i++) {
    A->ExtractMyRowView(i, StandardNumEntries, 
		       StandardValues, StandardIndices);
    for (j=0; j < StandardNumEntries; j++) {
      printf("MyPID=%d, J[%d,%d]=%e\n",MyPID,StandardMap->GID(i),
	     OverlapMap->GID(StandardIndices[j]),StandardValues[j]);
    }
  }
  */

  //if (MyPID==0) printf("Completed!!!\n");

  delete &u;
  delete &x;
  
  delete [] xx;
  delete [] uu;
  delete [] indicies;
  delete RowValues;

  return ;
}
