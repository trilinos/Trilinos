#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "FEpetra_Map.h"
#include "FEpetra_Vector.h"


int main( int argc, char* argv[] )
{
  
  /*
   * Data declarations (how old-school is this!)
   */

  int numGlobalElements, numGlobalElements_rtn;
  MapID mapID;
  VectorID xID, bID;
  double bnorm, xnorm;

  /*
   * Executable code
   */
  
  /* Create a map */
  numGlobalElements = 4;
  mapID = FEpetra_Map_Create(numGlobalElements);

  numGlobalElements_rtn = FEpetra_Map_NumGlobalElements(mapID);
  printf( "NumGlobalElements = %d\n", numGlobalElements_rtn );
  assert( numGlobalElements == numGlobalElements_rtn );
  
  /* Create vectors */
  xID = FEpetra_Vector_Create(mapID);
  bID = FEpetra_Vector_Create(mapID);

  /* Do some vector operations */
  FEpetra_Vector_Random(bID);
  FEpetra_Vector_Update(xID,2.0,bID,0.0); /* x = 2*b */

  bnorm = FEpetra_Vector_Norm2(bID);
  xnorm = FEpetra_Vector_Norm2(xID);

  printf( "2 norm of x = %f\n", xnorm );
  printf( "2 norm of b = %f\n", bnorm );

  /* Clean up memory (in reverse order)! */
  FEpetra_Vector_Destroy(bID);
  FEpetra_Vector_Destroy(xID);
  FEpetra_Map_Destroy(mapID);

  /* This should throw an exception and print an error message! */
  /* FEpetra_Map_NumGlobalElements(mapID); */

  return 0;

}
