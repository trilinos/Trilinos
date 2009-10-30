// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

void Trilinos_Util_distrib_msr_matrix(const Epetra_Comm & Comm, int *N_global,
	      int *n_nonzeros, int *N_update, int **update,
	      double **val, int **bindx,
	      double **x, double **b, double **xexact)
#undef DEBUG 

{
  int i, n_global_nonzeros;
  int j, row, have_xexact = 0 ;
  int *bindx1;
  double *val1, *b1, *x1, *xexact1=0;

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  printf("Processor %d of %d entering distrib_matrix.\n",
	 MyPID,NumProc) ;

  /*************** Distribute global matrix to all processors ************/

  if(MyPID == 0)
    {
      if ((*xexact) != NULL) have_xexact = 1;
      printf("%s", "Broadcasting exact solution\n");
    }

  if(NumProc  > 1)
    { 

      Comm.Broadcast( N_global, 1, 0);
      Comm.Broadcast( n_nonzeros, 1, 0);
      Comm.Broadcast( &have_xexact, 1, 0);

      if(MyPID != 0)
	{
	  (*bindx) = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
	  (*val) = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
	}

      Comm.Broadcast( (*bindx), (*n_nonzeros+1),  0);
      Comm.Broadcast( (*val),  (*n_nonzeros+1),   0);

      printf("Processor %d of %d done with matrix broadcast.\n",
	     MyPID,NumProc) ;
 
      /* Set rhs and initialize guess */
      if(MyPID != 0)
	{
	  (*b) = (double *) calloc(*N_global,sizeof(double)) ;
	  (*x) = (double *) calloc(*N_global,sizeof(double)) ;
	  if (have_xexact)
	  (*xexact) =   (double *) calloc(*N_global,sizeof(double)) ;
	}

      Comm.Broadcast( (*x), (*N_global), 0);
      Comm.Broadcast( (*b), (*N_global), 0);
      if (have_xexact)
	Comm.Broadcast( (*xexact), (*N_global), 0);
      printf("Processor %d of %d done with rhs/guess broadcast.\n",
	     MyPID,NumProc) ;

    }

  /********************** Generate update map  *************************/

  // read_update(N_update, update, proc_config, (*N_global), 1, linear) ;

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

#ifdef DEBUG
      { double sum1 = 0.0;
      for (i=0;i<(*N_global); i++) sum1 += (*b)[i];

      printf("Processor %d of %d has sum of b = %12.4g.\n",
	     MyPID,NumProc,sum1) ;
      }
#endif /* DEBUG */

      /* Allocate memory for local matrix */

      bindx1 = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
      val1 = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
      b1 =   (double *) calloc(*N_update,sizeof(double)) ;
      x1 =   (double *) calloc(*N_update,sizeof(double)) ;
      if (have_xexact)
      xexact1 =   (double *) calloc(*N_update,sizeof(double)) ;
     
      bindx1[0] = *N_update+1;
      
      for (i=0; i<*N_update; i++)
	{
	  row = (*update)[i];
	  b1[i] = (*b)[row];
	  x1[i] = (*x)[row];
	  if (have_xexact) xexact1[i] = (*xexact)[row];
	  val1[i] = (*val)[row];
	  bindx1[i+1] = bindx1[i];

#ifdef DEBUG	  
	  printf("Proc %d of %d: Global row = %d: Local row = %d: b = %12.4g: x = %12.4g: bindx = %d: val = %12.4g \n",
		 MyPID,NumProc, row, i, b1[i], x1[i], bindx1[i], val1[i]) ;
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

      if (have_xexact)
	{
	  printf(
     "The residual using MSR format and exact solution on processor %d is %12.4g\n",
	      MyPID,
	      Trilinos_Util_smsrres (*N_update, (*N_global), val1, bindx1, xexact1, (*xexact), b1));
	}
  
      /* Release memory for global matrix, rhs and solution */
      
      free ((void *) (*val));
      free ((void *) (*bindx));
      free ((void *) (*b));
      free ((void *) (*x));
      if (have_xexact) free((void *) *xexact);

      /* Return local matrix through same pointers. */
      
      *val = val1;
      *bindx = bindx1;
      *b = b1;
      *x = x1;
      if (have_xexact) *xexact = xexact1;

    }
      if (have_xexact && NumProc  == 1)
	{
	  printf(
     "The residual using MSR format and exact solution on processor %d is %12.4g\n",
	      MyPID,
	      Trilinos_Util_smsrres (*N_update, (*N_global), (*val), 
              (*bindx), (*xexact), (*xexact), (*b)));
	}

  
  printf("Processor %d of %d leaving distrib_matrix.\n",
	 MyPID,NumProc) ;
  
  /* end distrib_matrix */
}
//
//---------------------------------------------------------------------------
// Alternate version of Trilinos_Util_distrib_msr_matrix.cpp
// This version only distributes the matrix from an HB
// file. It ignores any right-hand side, initial guess,
// and exact solution information that may be in the file.
//---------------------------------------------------------------------------
//

void Trilinos_Util_distrib_msr_matrix(const Epetra_Comm & Comm, int *N_global,
	      int *n_nonzeros, int *N_update, int **update,
	      double **val, int **bindx)
#undef DEBUG 

{
  int i, n_global_nonzeros;
  int j, row;
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

      /* Return local matrix through same pointers. */
      
      *val = val1;
      *bindx = bindx1;

    }
     
  
  printf("Processor %d of %d leaving distrib_matrix.\n",
	 MyPID,NumProc) ;
  
  /* end distrib_matrix */
}
