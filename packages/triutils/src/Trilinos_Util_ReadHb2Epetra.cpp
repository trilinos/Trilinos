#include <stdlib.h>
#include <stdio.h>
#include "Trilinos_Util.h"
#include "iohb.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

void Trilinos_Util_ReadHb2Epetra(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {
  FILE *in_file ;
  int numGlobalEquations, N_columns, n_entries, Nrhs;
  char Title[73], Key[9], Rhstype[4];
  char Type[4] = "XXX";
  char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
  int Ptrcrd, Indcrd, Valcrd, Rhscrd;
  
  int * bindx, * pntr, * indx1, * pntr1;
  double * val, * val1, * hbx, * hbxexact, * hbb;

  hbb = 0; hbxexact = 0; hbb = 0;

  if(comm.MyPID() == 0)  { 
      in_file = fopen( data_file, "r");
      if (in_file == NULL)
	{
	  printf("Error: Cannot open file: %s\n",data_file);
	  exit(1);
	}

      /* Get information about the array stored in the file specified in the  */
      /* argument list:                                                       */

      printf("Reading matrix info from %s...\n",data_file);
      
      in_file = fopen( data_file, "r");
      if (in_file == NULL)
	{
	  printf("Error: Cannot open file: %s\n",data_file);
	  exit(1);
	}

      readHB_header(in_file, Title, Key, Type, &numGlobalEquations, &N_columns, 
		    &n_entries, &Nrhs,
		    Ptrfmt, Indfmt, Valfmt, Rhsfmt,
		    &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
      fclose(in_file);

      if (Nrhs < 0 ) Nrhs = 0;

      printf("***************************************************************\n");
      printf("Matrix in file %s is %d x %d, \n",data_file, numGlobalEquations, N_columns);
      printf("with %d nonzeros with type %3s;\n", n_entries, Type);
      printf("***************************************************************\n");
      printf("Title: %72s\n",Title);
      printf("***************************************************************\n");
      /*Nrhs = 0; */
      printf("%d right-hand-side(s) available.\n",Nrhs);

      if (Type[0] != 'R') perror("Can only handle real valued matrices");
      int isym = 0;
      if (Type[1] == 'S') 
	{
	  printf("Converting symmetric matrix to nonsymmetric storage\n");
	  n_entries = 2*n_entries - N_columns;
	  isym = 1;
	}
      if (Type[2] != 'A') perror("Can only handle assembled matrices");
      if (N_columns != numGlobalEquations) perror("Matrix dimensions must be the same");
      
      /* Read the matrix information, generating the associated storage arrays  */
      printf("Reading the matrix from %s...\n",data_file);

      /* Allocate space.  Note that we add extra storage in case of zero
	 diagonals.  This is necessary for conversion to MSR format. */

      pntr   = (int    *) calloc(N_columns+1,sizeof(int)) ;
      bindx = (int    *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val   = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      readHB_mat_double(data_file, pntr, bindx, val);

      /* Translate integer arrays to zero base */
      for (int i = 0; i <= numGlobalEquations; i++) pntr[i]--;
      {for (int i = 0; i <= n_entries; i++) bindx[i]--;}

      /* If a rhs is specified in the file, read one, 
	 generating the associate storage */
      if (Nrhs > 0 && Rhstype[2] =='X')
	{
	  printf("Reading right-hand-side vector(s) from %s...\n",data_file);
	  hbb = (double *) calloc(N_columns,sizeof(double));
	  readHB_aux_double(data_file, 'F', hbb);
	  printf("Reading exact solution  vector(s) from %s...\n",data_file);
	  hbxexact = (double *) calloc(N_columns,sizeof(double));
	      readHB_aux_double(data_file, 'X', hbxexact);

	}
      else
	{
	  
	  /* Set Xexact to a random vector */

	  printf("Setting  random exact solution  vector\n");
	  hbxexact = (double *) calloc(N_columns,sizeof(double));
	  
	  for (int i=0;i<numGlobalEquations;i++)	 hbxexact[i] = 
                                              ((double)
                                               rand())/((double) RAND_MAX);
	  
	  /* Compute b to match xexact */
	  
	 hbb = (double   *) calloc(N_columns,sizeof(double)) ;
	  if (hbb == NULL) perror("Error: Not enough space to create rhs");
 

      Trilinos_Util_scscmv (isym, N_columns, N_columns, val, bindx, pntr, hbxexact, hbb);
	}

      /* Compute residual using CSC format */

      double res = Trilinos_Util_scscres(isym, numGlobalEquations, numGlobalEquations, val, bindx, pntr, 
		    hbxexact, hbb);
      printf(
	      "The residual using CSC format and exact solution is %12.4g\n",
	      res);

      
      /* Set initial guess to zero */
      
      hbx = (double   *) calloc(numGlobalEquations,sizeof(double)) ;
      
      if (hbx == NULL) 
	perror("Error: Not enough space to create guess");
      
      
      /* Set RHS to a random vector, initial guess to zero */
      {for (int i=0;i<numGlobalEquations;i++) hbx[i] = 0.0;}
      
      
      /* Allocate temporary space */
      
      pntr1 = (int   *) calloc(N_columns+1,sizeof(int)) ;
      indx1 = (int   *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val1 = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      
      /* Convert in the following way:
	 - CSC to CSR 
	 - CSR to MSR
      */
      Trilinos_Util_csrcsc(numGlobalEquations, numGlobalEquations, 0, 0, val, bindx, pntr, val1, indx1, pntr1);
      
      if (Type[1] == 'S') 
	{
	  int *indu, *iwk;
	  int ierr;
	  indu = new int[N_columns];
	  iwk = new int[N_columns+1];
	  ierr = Trilinos_Util_ssrcsr(3, 1, N_columns, val1, indx1, pntr1, n_entries,
	  		  val1, indx1, pntr1, indu, iwk);
	  delete [] indu;
	  delete [] iwk;
	  if (ierr !=0 ) 
	    {
	    printf(" Error in converting from symmetric form\n  IERR = %d\n",ierr);
	    abort();
	    }
	}
  }
  comm.Broadcast(&numGlobalEquations, 1, 0);
  int nlocal = 0;
  if (comm.MyPID()==0) nlocal = numGlobalEquations;
  map = new Epetra_Map(numGlobalEquations, nlocal, 0, comm); // Create map with all elements on PE 0
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0
  if (comm.MyPID()==0)
    for (int i=0; i<numGlobalEquations; i++)
      A->InsertGlobalValues(i, pntr1[i+1]-pntr1[i], val1+pntr1[i], indx1+pntr1[i]);
  A->TransformToLocal();
  
  x = new Epetra_Vector(Copy, *map, hbx);
  b = new Epetra_Vector(Copy, *map, hbb);
  xexact = new Epetra_Vector(Copy, *map, hbxexact);

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
  
  /* Release unneeded space */
  if (comm.MyPID()==0) {
    if (hbb!=0) free((void *) hbb);
    if (hbx!=0) free((void *) hbx);
    if (hbxexact!=0) free((void *) hbxexact);
    free((void *) val);
    free((void *) bindx);
    free((void *) val1);
    free((void *) indx1);
    free((void *) pntr1);
    free((void *) pntr);
  }
  return;
}
