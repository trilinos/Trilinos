#include<iostream.h>
#include<math.h>
#include "fill.h"

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
void Fill::fillMatrix(const Epetra_Vector *tmp_soln, Epetra_Vector *tmp_rhs, 
		      Epetra_RowMatrix *tmp_matrix) 
{

  // Decide what kind of fill call this will be
  if (tmp_matrix==NULL) {
    flag = RHS_ONLY;
    soln = const_cast<Epetra_Vector*>(tmp_soln);
//    soln = tmp_soln;
    rhs = tmp_rhs;
    //if (Comm->MyPID()==0) printf("Starting RESIDUAL ONLY Fill....");
  } else if (tmp_rhs==NULL) {
    flag = MATRIX_ONLY;
    soln = const_cast<Epetra_Vector*>(tmp_soln);
//    soln = tmp_soln;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
    //if (Comm->MyPID()==0) printf("Starting MATRIX ONLY Fill....");
  } else { 
    flag = ALL;
    soln = const_cast<Epetra_Vector*>(tmp_soln);
//    soln = tmp_soln;
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
    //if (Comm->MyPID()==0) printf("Starting Matrix and RESIDUAL Fill....");
  }


  // Create Linear Objects for Fill
  Epetra_Vector& u = *new Epetra_Vector(*OverlapMap);

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

  // Export Solution to Overlap vector
  Epetra_Import & Importer2 = *new Epetra_Import(*OverlapMap, *StandardMap);
  assert(u.Import(*soln, Importer2, Insert)==0);
  delete &Importer2;

/*
  printf("Inside fill.cc Proc=%d  #Elements=%d\n",MyPID,NumMyElements);
  for (i=0; i< NumMyElements; i++)    {
    printf("Proc=%d, GID=%d, soln=%e\n",MyPID ,
      StandardMap->GID(i),(*soln)[i]); }
*/

  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == RHS_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);


  if((flag==RHS_ONLY) || (flag==ALL))
  {
    if (MyPID==0) 
      (*rhs)[0]=(u[0]*u[0] + u[1]*u[1] - 2.);
    else 
      (*rhs)[0]=(exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.);
  }

  if((flag==MATRIX_ONLY) || (flag==ALL))
  {
    if (MyPID==0) 
    {
      int row = 0; 
      int column = 0;
      double jac = 2.*u[0];
      ierr=A->ReplaceGlobalValues(row, 1, &jac, &column);
      if (ierr!=0) ierr=A->InsertGlobalValues(row, 1, &jac, &column);
      column++;
      jac = 2.*u[1];
      ierr=A->ReplaceGlobalValues(row, 1, &jac, &column);
      if (ierr!=0) ierr=A->InsertGlobalValues(row, 1, &jac, &column);
    }
    else 
    {
      int row = 1; 
      int column = 0;
      double jac = exp(u[0]-1.);
      ierr=A->ReplaceGlobalValues(row, 1, &jac, &column);
      if (ierr!=0) ierr=A->InsertGlobalValues(row, 1, &jac, &column);
      column++;
      jac = 3.*u[1]*u[1];
      ierr=A->ReplaceGlobalValues(row, 1, &jac, &column);
      if (ierr!=0) ierr=A->InsertGlobalValues(row, 1, &jac, &column);
    }
  }

  Comm->Barrier();
 
  A->TransformToLocal();
  
/*
  // Print RHS   
  if((flag==RHS_ONLY) || (flag==ALL))
  {
    printf("Proc=%d  #Elements=%d\n",MyPID,NumMyElements);
    for (i=0; i< NumMyElements; i++) {
      printf("Proc=%d, GID=%d, rhs=%e\n",MyPID ,
        StandardMap->GID(i),(*rhs)[i]);
    }
  }
 

  // Print Matrix
  if((flag==MATRIX_ONLY) || (flag==ALL))
  {
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
  }

  Comm->Barrier();

  if (MyPID==0) 
  {
    printf("Completed!!!\n");
    getchar();
  }
*/

  delete &u;
  

  return ;
}
