//
// file distrib_msr_matrix.cpp
//
// Modified version of Trilinos_Util_distrib_msr_matrix.cpp
// This version only distributes the matrix from an HB
// file. It ignores any right-hand side, initial guess,
// and exact solution information that may be in the file.

#include <stdlib.h>
#include <stdio.h>
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Util.h"

void distrib_msr_matrix(const Epetra_Comm & Comm, int *N_global,
	      int *n_nonzeros, int *N_update, int **update,
	      double **val, int **bindx)
#undef DEBUG 

{
  int i, n_global_nonzeros;
  int j, row, have_xexact = 0 ;
  int *bindx1;
  double *val1;

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  printf("Processor %d of %d entering distrib_matrix.\n",
	 MyPID,NumProc) ;

  /*************** Distribute global matrix to all processors ************/


  if(NumProc  > 1)
    { 

      Comm.Broadcast( N_global, 1, 0);
      Comm.Broadcast( n_nonzeros, 1, 0);

      if(MyPID != 0)
	{
	  (*bindx) = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
	  (*val) = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
	}

      Comm.Broadcast( (*bindx), (*n_nonzeros+1),  0);
      Comm.Broadcast( (*val),  (*n_nonzeros+1),   0);

      printf("Processor %d of %d done with matrix broadcast.\n",
	     MyPID,NumProc) ;
 
    }

  /********************** Generate update map  *************************/


  Epetra_Map map(*N_global, 0, Comm);
  *N_update = map.NumMyElements();
  (*update) = (int *) calloc(*N_update,sizeof(int)) ;
  map.MyGlobalElements(*update);
  
  printf("Processor %d of %d has %d rows of %d total rows.\n",
	 MyPID,NumProc,*N_update,(*N_global)) ;

  /*************** Construct local matrix from global matrix ************/

  /* The local matrix is a copy of the rows assigned to this processor.  
     It is stored in MSR format and still has global indices
  */

  if(NumProc  > 1)
    { 
      n_global_nonzeros = *n_nonzeros;

      *n_nonzeros = *N_update;
      
      for (i=0; i<*N_update; i++)
	*n_nonzeros += (*bindx)[(*update)[i]+1] - (*bindx)[(*update)[i]];

      printf("Processor %d of %d has %d nonzeros of %d total nonzeros.\n",
	     MyPID,NumProc,
	     *n_nonzeros,n_global_nonzeros) ;

      /* Allocate memory for local matrix */

      bindx1 = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
      val1 = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
     
      bindx1[0] = *N_update+1;
      
      for (i=0; i<*N_update; i++)
	{
	  row = (*update)[i];
	  val1[i] = (*val)[row];
	  bindx1[i+1] = bindx1[i];

#ifdef DEBUG	  
	  printf("Proc %d of %d: Global row = %d: Local row = %d: bindx = %d: val = %12.4g \n",
		 MyPID,NumProc, row, i, bindx1[i], val1[i]) ;
#endif

	  for (j = (*bindx)[row]; j < (*bindx)[row+1]; j++)
	    {
	      val1[  bindx1 [i+1] ] = (*val)[j];
	      bindx1[bindx1 [i+1] ] = (*bindx)[j];
	      bindx1[i+1] ++;
	    }
	}

      printf("Processor %d of %d done with extracting local operators.\n",
	     MyPID,NumProc) ;

  
      /* Release memory for global matrix, rhs and solution */
      
      free ((void *) (*val));
      free ((void *) (*bindx));
      //free ((void *) (*b));
      //free ((void *) (*x));
      //if (have_xexact) free((void *) *xexact);

      /* Return local matrix through same pointers. */
      
      *val = val1;
      *bindx = bindx1;

    }
     
  
  printf("Processor %d of %d leaving distrib_matrix.\n",
	 MyPID,NumProc) ;
  
  /* end distrib_matrix */
}
