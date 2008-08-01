/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_include.h"
#include "ml_qr_fix.h"
#include <stdio.h>

/* -mb: make this global var for now - decide where to put in ML */
static ML_qr_fix *xCDeadNodDof = NULL;

int ML_qr_fix_Create(const int nCoarseNod)
{
   int    i, nBytes;
   xCDeadNodDof = (ML_qr_fix*) ML_allocate(sizeof(ML_qr_fix));
   xCDeadNodDof->level         = 0;
   xCDeadNodDof->numDeadNodDof = 0;
   xCDeadNodDof->nDeadNodDof   = nCoarseNod;
   nBytes = (nCoarseNod+1)*sizeof(ML_QR_FIX_TYPE *);
   xCDeadNodDof->xDeadNodDof   = (ML_QR_FIX_TYPE *)ML_allocate(nBytes);
   
   for (i=0; i < nCoarseNod; i++) (xCDeadNodDof->xDeadNodDof)[i] = 0;

   return(0);
}

int ML_qr_fix_Destroy(void)
{
  if (xCDeadNodDof == NULL)
    return(0);

  /* loop over all fields and release them */
  if (xCDeadNodDof->xDeadNodDof) ML_free(xCDeadNodDof->xDeadNodDof);

  /* free the structure itself */

  ML_free(xCDeadNodDof);
  xCDeadNodDof = NULL;

  return(0);
}

int ML_qr_fix_Print(ML_qr_fix* ptr)
{
    int i, cnt;

    if (ptr == NULL)
        return(-1);

    printf("level = %d nodes containing dead dofs:\n", ptr->level);
    cnt = 0;
    for (i=0; i < ptr->nDeadNodDof; i++) {
        if ((ptr->xDeadNodDof)[i]) {
            printf("No. %8d node index %8d\n", ++cnt, i);
        }    
    }

    return(0);
}

int ML_qr_fix_NumDeadNodDof(void) { 
    if (xCDeadNodDof == NULL) return 0;
    return xCDeadNodDof->numDeadNodDof; 
}

ML_QR_FIX_TYPE ML_qr_fix_getDeadNod(const int nodInx)
{
     if (xCDeadNodDof == NULL) return 0;
     return (xCDeadNodDof->xDeadNodDof)[nodInx];
}

void ML_qr_fix_setNumDeadNod(int num)
{
     if (xCDeadNodDof == NULL) return;
     xCDeadNodDof->numDeadNodDof = num;
}

void ML_qr_fix_setDeadNod( 
         const int      nodInx, 
         ML_QR_FIX_TYPE val
     )
{
     if (xCDeadNodDof == NULL) return;
     (xCDeadNodDof->xDeadNodDof)[nodInx] = val;
}

int  ML_qr_fix_Bitsize(void)
{
     if (xCDeadNodDof == NULL) return 0;
     return (int)(8*sizeof((xCDeadNodDof->xDeadNodDof)[0]));
}
