#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Trilinos_Util_CountTriples.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

int Trilinos_Util_ReadTriples2Epetra( char *data_file,
				      bool symmetric, 
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact ) {
  FILE *in_file ;
  int N_rows, nnz ; 

  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  vector<int> non_zeros;   // Number of non-zeros in each row
  Trilinos_Util_CountTriples( data_file, symmetric, non_zeros, N_rows, nnz, comm ) ; 

  vector<int> ptrs(N_rows+1) ; // Pointers into inds and vals for the start of each row
  vector<int> inds(nnz);     //  Column Indices
  vector<double> vals(nnz);  //  Matrix values

  if(comm.MyPID() == 0)  { 
    //  ptrs, inds and vals together constitute a compressed row storage of the matrix.


    in_file = fopen( data_file, "r");
    assert (in_file != NULL) ;  // Checked in T_U_CountTriples() 

    ptrs[0] = 0 ; 
    for( int i=0; i< N_rows; i++ ) { 
      ptrs[i+1] = ptrs[i] + non_zeros[i]; 
    }

    vector<int> iptrs = ptrs ; //  Current pointers into inds and vals for each row

    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int i, j; 
      double val ; 
      sscanf( buffer, "%d %d %lg", &i, &j, &val ) ; 
      if ( i == j && i%100 == 0 ) cout << " i,j = " << i << " val = " << val << endl ; 
      int iptr = iptrs[i-1] ; 
      iptrs[i-1]++ ;
      vals[iptr] = val ; 
      inds[iptr] = j-1 ; 
      //
      //  If this is a symmetric matrix, we need to enter the entry 
      //  for the other triangular half
      //
      if (symmetric && i != j ) {
	iptr = iptrs[j-1] ; 
	iptrs[j-1]++;
	vals[iptr] = val ; 
	inds[iptr] = i-1 ; 
      }
    } 
    fclose(in_file);
  }

  int nlocal = 0;
  if (comm.MyPID()==0) nlocal = N_rows;
  map = new Epetra_Map(N_rows, nlocal, 0, comm); // Create map with all elements on PE 0
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0
  if (comm.MyPID()==0)
    for (int i=0; i<N_rows; i++)
      A->InsertGlobalValues(i, ptrs[i+1]-ptrs[i], &vals[ptrs[i]], &inds[ptrs[i]]);
  A->TransformToLocal();

  vector<double> hbx(N_rows);
  vector<double> hbb(N_rows);
  vector<double> hbxexact(N_rows);

  
  x = new Epetra_Vector(Copy, *map, &hbx[0]);
  b = new Epetra_Vector(Copy, *map, &hbb[0]);
  xexact = new Epetra_Vector(Copy, *map, &hbxexact[0]);

  EPETRA_CHK_ERR( x->PutScalar( 0.0 ) );
  EPETRA_CHK_ERR( xexact->Random( ) ) ; 
  EPETRA_CHK_ERR( A->Multiply( false, *xexact, *b ) ); 

  return 0;
}
