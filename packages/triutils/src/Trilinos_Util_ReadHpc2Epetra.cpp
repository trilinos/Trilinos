#include <stdlib.h>
#include <stdio.h>
#include "Trilinos_Util.h"
#include "iohb.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

void Trilinos_Util_ReadHpc2Epetra(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {


  FILE *in_file ;

  int l;
  int * lp = &l;
  double v;
  double * vp = &v;
#ifdef DEBUG
  bool debug = true;
#else
  bool debug = false;
#endif

  int size = comm.NumProc();
  int rank = comm.MyPID();
  printf("Reading matrix info from %s...\n",data_file);
  
  in_file = fopen( data_file, "r");
  if (in_file == NULL)
    {
      printf("Error: Cannot open file: %s\n",data_file);
      exit(1);
    }
  int numGlobalEquations, total_nnz;
  fscanf(in_file,"%d",&numGlobalEquations);
  fscanf(in_file,"%d",&total_nnz);
  map = new Epetra_Map(numGlobalEquations, 0, comm); // Create map with uniform distribution
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix

  x = new Epetra_Vector(*map);
  b = new Epetra_Vector(*map);
  xexact = new Epetra_Vector(*map);
  int numMyEquations = map->NumMyPoints();

  // Allocate arrays that are of length numMyEquations

  // Find max nnz per row for this processor

  int max_nnz = 0;

  for (int i=0; i<numGlobalEquations; i++) {
      fscanf(in_file, "%d",lp); /* row #, nnz in row */
      if (map->MyGID(i)) max_nnz = EPETRA_MAX(max_nnz,l);
    }


  // Allocate arrays that are of length local_nnz
  double * list_of_vals = new double[max_nnz];
  int *list_of_inds = new int   [max_nnz];

  for (int i=0; i<numGlobalEquations; i++)
    {
      int cur_nnz;
      fscanf(in_file, "%d",&cur_nnz);
      if (map->MyGID(i)) // See if nnz for row should be added
	{
	  if (debug) cout << "Process "<<rank
			  <<" of "<<size<<" getting row "<<i<<endl;
	  int nnz_kept = 0;
	  for (int j=0; j<cur_nnz; j++) 
	    {
	      fscanf(in_file, "%lf %d",vp,lp);
	      if (v!=0.0) {
		list_of_vals[nnz_kept] = v;
		list_of_inds[nnz_kept] = l;
		nnz_kept++;
	      }
	    }
	  A->InsertGlobalValues(i, nnz_kept, list_of_vals, list_of_inds);
	}
      else
	for (int j=0; j<cur_nnz; j++) fscanf(in_file, "%lf %d",vp,lp); // otherwise read and discard
    }

  double xt, bt, xxt;
  for (int i=0; i<numGlobalEquations; i++) 
    {
      if (map->MyGID(i)) // See if entry should be added
	{
	  if (debug) cout << "Process "<<rank<<" of "
                       <<size<<" getting RHS "<<i<<endl;
	  fscanf(in_file, "%lf %lf %lf",&xt, &bt, &xxt);
	  int cur_local_row = map->LID(i);
	  (*x)[cur_local_row] = xt;
	  (*b)[cur_local_row] = bt;
	  (*xexact)[cur_local_row] = xxt;
	}
      else
	fscanf(in_file, "%lf %lf %lf",vp, vp, vp); // or thrown away
    }

  fclose(in_file);

  
  if (debug)
    cout << "Process "<<rank<<" of "<<size<<" has "<<numMyEquations
	 << " rows. Min global row "<< map->MinMyGID()
	 <<" Max global row "<< map->MaxMyGID() <<endl
	 <<" and "<<A->NumMyNonzeros()<<" nonzeros."<<endl;

  A->TransformToLocal();
  

  Epetra_Vector bcomp(*map);

  A->Multiply(false, *xexact, bcomp);
  double residual;
  bcomp.Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of computed b = " << residual << endl;
  b->Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of given b    = " << residual << endl;
  bcomp.Update(-1.0, *b, 1.0);
  bcomp.Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of difference between computed b and given b for xexact = " << residual << endl;
  
  delete [] list_of_vals;
  delete []list_of_inds;

  return;
}
