#include <stdlib.h>
#include <stdio.h>
#include "Trilinos_Util.h"

void Trilinos_Util_dusds_vbr(SPBLASMAT *A)
 
/*  Destroy handle for a VBR matrix.  Build any auxiliary data structures
    that might be helpful for performance.
*/

{
  int * ncolvec;
  double * buffer;
  
  ncolvec = A->ncolvec;
  buffer  = A->buffer;
  free ((void *) ncolvec);
  free ((void *) buffer);
  
/* end Trilinos_Util_dusds_vbr */
}
