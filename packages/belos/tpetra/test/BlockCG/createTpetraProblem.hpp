// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
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
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace Belos {

template <class Scalar>
bool createEpetraProblem(
    Teuchos::RCP<Tpetra::Platform<int> >           &platform,
    std::string                                    &filename,
    Teuchos::RCP<Tpetra::Map<int> >                &rowMap,
    Teuchos::RCP<Tpetra::CrsMatrix<int,Scalar> >   &A,
    Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > &B,
    Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > &X
);

template <class Scalar>
void distrib_matrix(const Teuchos::Comm<int> &comm, int *N_global,
                    int *n_nonzeros, int *N_update, int **update,
                    Scalar **val, int **bindx,
                    Scalar **x, Scalar **b, Scalar **xexact);

template <class Scalar>
bool createEpetraProblem(
    Teuchos::RCP<Tpetra::Platform<int> >           &platform,
    std::string                                    &filename,
    Teuchos::RCP<Tpetra::Map<int> >                &rowMap,
    Teuchos::RCP<Tpetra::CrsMatrix<int,Scalar> >   &A,
    Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > &B,
    Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > &X
    )
{
  using Teuchos::RCP; 
  using Teuchos::rcp;
  RCP<Tpetra::Comm<int> > comm = platform->createComm();
  const int myImageID = comm->getRank();
  const int numImages = comm->getSize();
  //
  int i;
  int n_nonzeros, N_update;
  int *bindx=0, *update=0, *col_inds=0;
  double *val=0, *row_vals=0;
  double *xguess=0, *b=0, *xexact=0;
  char have_exact = 0;

  //
  // **********************************************************************
  // ******************Set up the problem to be solved*********************
  // **********************************************************************
  //
  int NumGlobalElements;  // total # of rows in matrix
  //
  // *****Read in matrix from HB file******
  //
  Trilinos_Util_read_hb(const_cast<char *>(filename.c_str()), myImageID, &NumGlobalElements, &n_nonzeros,
      &val, &bindx, &xguess, &b, &xexact);

  // 
  // *****Distribute data among processors*****
  //
  distrib_matrix(platform, &NumGlobalElements, &n_nonzeros, &N_update, 
      &update, &val, &bindx, &xguess, &b, &xexact);
  //
  // *****Construct the matrix*****
  //
  int NumMyElements = N_update; // # local rows of matrix on processor
  //
  // Create an integer std::vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  //
  int * NumNz = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++) {
    NumNz[i] = bindx[i+1] - bindx[i] + 1;
  }
  //
  RCP<Tpetra::Map<int> > tpetraMap = rcp(new Tpetra::Map<int>(NumGlobalElements, Teuchos::arrayView(update,NumMyElements), 0, platform) );
  rowMap = tpetraMap;
  //
  // Create a Tpetra::CrsMatrix
  //
  A = rcp(new Tpetra::CrsMatrix(*tpetraMap));
  //
  // Add rows one-at-a-time
  //
  int NumEntries;
  for (i=0; i<NumMyElements; i++) {
    row_vals = val + bindx[i];
    col_inds = bindx + bindx[i];
    NumEntries = bindx[i+1] - bindx[i];
    // insert off-diags
    A->submitEntries(update[i], Teuchos::arrayView(col_inds,NumEntries), Teuchos::arrayView(row_vals,NumEntries) );
    // insert diag
    A->submitEntry(update[i],update[i],val[i]);
  }
  //
  // Finish up
  //
  A->fillComplete();
  //
  // Construct the right-hand side and solution multivectors.
  //
  if (B) {
    B = rcp(new Tpetra::MultiVector<int,Scalar>(*tpetraMap,Teuchos::arrayView(b,
    *B = rcp(new Epetra_MultiVector(::Copy, *epetraMap, b, NumMyElements, 1 ));
    Teuchos::set_extra_data( epetraMap, "B::Map", B );
  }
  if (X) {
    *X = rcp(new Epetra_MultiVector(*epetraMap, 1 ));
    Teuchos::set_extra_data( epetraMap, "X::Map", X );
  }
  //
  // Create workspace
  //
  Teuchos::set_default_workspace_store(
      Teuchos::rcp(new Teuchos::WorkspaceStoreInitializeable(static_cast<size_t>(2e+6)))
      );
  //
  // Free up memory
  //
  delete [] NumNz;
  free(update);
  free(val);
  free(bindx);
  if (xexact) free(xexact);
  if (xguess) free(xguess);
  if (b) free(b);
}


template <class Scalar>
void distrib_matrix(Teuchos::RCP<Tpetra::Platform<int> > &platform, int *N_global,
                    int *n_nonzeros, int *N_update, int **update,
                    Scalar **val, int **bindx,
                    Scalar **x, Scalar **b, Scalar **xexact)
{
  int i, n_global_nonzeros;
  int j, row, have_xexact = 0 ;
  int *bindx1;
  double *val1, *b1, *x1, *xexact1;

  Teuchos::RCP<Teuchos::Comm<int> > comm = platform->createComm();

  int MyPID = Teuchos::rank(*comm);
  int NumProc = Teuchos::size(*comm);

  printf("Processor %d of %d entering distrib_matrix.\n", MyPID,NumProc) ;

  /*************** Distribute global matrix to all processors ************/

  if (MyPID == 0) {
    if ((*xexact) != NULL) have_xexact = 1;
    printf("%s", "Broadcasting exact solution\n");
  }

  if (NumProc  > 1) {
    Teuchos::broadcast(*comm,0,&N_global);
    Teuchos::broadcast(*comm,0,&n_nonzeros);
    Teuchos::broadcast(*comm,0,&have_xexact);

    if (MyPID != 0) {
      (*bindx) = (int    *) calloc(*n_nonzeros+1,sizeof(int)) ;
      (*val)   = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
    }

    Teuchos::broadcast(*comm,0,*n_nonzeros+1,*bindx);
    Teuchos::broadcast(*comm,0,*n_nonzeros+1,*val);

    printf("Processor %d of %d done with matrix broadcast.\n",
        MyPID,NumProc) ;

    /* Set rhs and initialize guess */
    if(MyPID != 0) {
      (*b) = (double *) calloc(*N_global,sizeof(double)) ;
      (*x) = (double *) calloc(*N_global,sizeof(double)) ;
      if (have_xexact) {
        (*xexact) =   (double *) calloc(*N_global,sizeof(double)) ;
      }
    }

    Teuchos::broadcast(*comm,0,*N_global,*x);
    Teuchos::broadcast(*comm,0,*N_global,*b);
    if (have_xexact) {
      Teuchos::broadcast(*comm,0,*N_global,*xexact);
    }
    printf("Processor %d of %d done with rhs/guess broadcast.\n", MyPID,NumProc) ;
  }

  /********************** Generate update map  *************************/

  Tpetra::Map<int> map(*N_global,0,*platform);
  *N_update = map.getNumMyEntries();
  (*update) = (int *) calloc(*N_update,sizeof(int)) ;
  {
    Teuchos::ArrayView<int> av = map.getMyGlobalEntries();
    std::copy(av.begin(),av.end(),*update);
  }

  printf("Processor %d of %d has %d rows of %d total rows.\n",
      MyPID,NumProc,*N_update,(*N_global)) ;

  /*************** Construct local matrix from global matrix ************/

  /* The local matrix is a copy of the rows assigned to this processor.  
     It is stored in MSR format and still has global indices
   */

  if (NumProc  > 1)
  { 
    n_global_nonzeros = *n_nonzeros;
    *n_nonzeros = *N_update;
    for (i=0; i<*N_update; i++) {
      *n_nonzeros += (*bindx)[(*update)[i]+1] - (*bindx)[(*update)[i]];
    }
    printf("Processor %d of %d has %d nonzeros of %d total nonzeros.\n", MyPID,NumProc, *n_nonzeros,n_global_nonzeros) ;
#ifdef DEBUG
    { 
      double sum1 = 0.0;
      for (i=0;i<(*N_global); i++) sum1 += (*b)[i];
      printf("Processor %d of %d has sum of b = %12.4g.\n", MyPID,NumProc,sum1) ;
    }
#endif /* DEBUG */

    /* Allocate memory for local matrix */
    bindx1 = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
    val1 = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
    b1 =   (double *) calloc(*N_update,sizeof(double)) ;
    x1 =   (double *) calloc(*N_update,sizeof(double)) ;
    if (have_xexact) {
      xexact1 =   (double *) calloc(*N_update,sizeof(double)) ;
    }
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

    printf("Processor %d of %d done with extracting local operators.\n", MyPID,NumProc) ;

    if (have_xexact) {
      printf("The residual using MSR format and exact solution on processor %d is %12.4g\n",
             MyPID,Trilinos_Util_smsrres (*N_update, (*N_global), val1, bindx1, xexact1, (*xexact), b1));
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

  if (have_xexact && NumProc  == 1) {
    printf("The residual using MSR format and exact solution on processor %d is %12.4g\n",
           MyPID, Trilinos_Util_smsrres (*N_update, (*N_global), (*val), (*bindx), (*xexact), (*xexact), (*b)));
  }

  printf("Processor %d of %d leaving distrib_matrix.\n", MyPID,NumProc) ;
}

} // namespace Belos
