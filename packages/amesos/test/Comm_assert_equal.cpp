#include "Comm_assert_equal.h"

bool Comm_assert_equal( const Epetra_Comm *Comm, int value ) {

  int max_value, min_value, diff; 

  Comm->MaxAll( &value, &max_value, 1) ; 
  Comm->MinAll( &value, &min_value, 1) ; 
  diff = max_value - min_value ;
  //  The broadcast can be omitted on homogeneous architectures
  Comm->Broadcast( &diff, 1, 0 ) ; 
  return( diff == 0 ) ; 
}

bool Comm_assert_equal( const Epetra_Comm *Comm, double value ) {

  double max_value, min_value, diff; 

  Comm->MaxAll( &value, &max_value, 1) ; 
  Comm->MinAll( &value, &min_value, 1) ; 
  diff = max_value - min_value ;
  //  The broadcast can be omitted on homogeneous architectures
  Comm->Broadcast( &diff, 1, 0 ) ; 
  return( diff == 0 ) ; 
}

