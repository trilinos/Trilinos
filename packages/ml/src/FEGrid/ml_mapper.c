/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Function for manipulating the ML_Mapper structure                    */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_mapper.h"

/* ******************************************************************** */
/* Constructor                                                          */
/* -------------------------------------------------------------------- */

int ML_Mapper_Create(ML_Mapper **mapper)
{
   ML_Mapper *ml_mapper;

   ML_memory_alloc( (void**) mapper, sizeof(ML_Mapper), "MA1" );
   ml_mapper = (*mapper);
   ml_mapper->ML_id = ML_ID_MAPPER;
   ml_mapper->data = NULL; 
   ml_mapper->inlength = 0; 
   ml_mapper->outlength = 0; 
   ml_mapper->map = NULL; 
   return 0;
}

/* ******************************************************************** */
/* Initialize                                                           */
/* -------------------------------------------------------------------- */

int ML_Mapper_Init(ML_Mapper *ml_mapper)
{
   ml_mapper->ML_id = ML_ID_MAPPER;
   ml_mapper->data = NULL; 
   ml_mapper->inlength = 0; 
   ml_mapper->outlength = 0; 
   ml_mapper->map = NULL; 
   return 0;
}

/* ******************************************************************** */
/* Destructor                                                           */
/* -------------------------------------------------------------------- */

int ML_Mapper_Destroy(ML_Mapper **mapper)
{
   ML_Mapper *ml_mapper;

   ml_mapper = (*mapper);
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_Destroy : wrong object.\n");
      exit(1);
   }
   ml_mapper->ML_id = -1;
   ml_mapper->data = NULL; 
   ml_mapper->inlength = 0;
   ml_mapper->outlength = 0;
   ml_mapper->map = NULL;
   ML_memory_free((void**) mapper );
   (*mapper) = NULL;
   return 0;
}

/* ******************************************************************** */
/* Destructor (only destroy internal data structures)                   */
/* -------------------------------------------------------------------- */

int ML_Mapper_Clean(ML_Mapper *ml_mapper)
{
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_Clean : wrong object.\n");
      exit(1);
   }
   ml_mapper->ML_id = -1;
   ml_mapper->data = NULL; 
   ml_mapper->inlength = 0;
   ml_mapper->outlength = 0;
   ml_mapper->map = NULL;
   return 0;
}

/* ******************************************************************** */
/* Check the existence of mapper function                               */
/* -------------------------------------------------------------------- */

int ML_Mapper_Check(ML_Mapper *ml_mapper)
{
#ifdef DEBUG
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_Set error : wrong object.\n");
      exit(1);
   }
#endif
   if ( ml_mapper->map != NULL ) return 1;
   else                          return 0;
} 

/* ******************************************************************** */
/* Set the mapper function                                              */
/* -------------------------------------------------------------------- */

int ML_Mapper_SetFunc(ML_Mapper *ml_mapper, int inlen, int outlen,
                  int (*func)(void*,double*,double*))
{
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_SetFunc error : wrong object.\n");
      exit(1);
   }
   ml_mapper->inlength = inlen;
   ml_mapper->outlength = outlen;
   ml_mapper->map = func;
   return 0;
}

/* ******************************************************************** */
/* Set the mapper data                                                  */
/* -------------------------------------------------------------------- */

int ML_Mapper_GetLength(ML_Mapper *ml_mapper, int *inlen, int *outlen)
{
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_SetData error : wrong object.\n");
      exit(1);
   }
   (*inlen) = ml_mapper->inlength;
   (*outlen) = ml_mapper->outlength;
   return 0;
}

/* ******************************************************************** */
/* Set the mapper data                                                  */
/* -------------------------------------------------------------------- */

int ML_Mapper_SetData(ML_Mapper *ml_mapper, void *data)
{
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_SetData error : wrong object.\n");
      exit(1);
   }
   ml_mapper->data = data;
   return 0;
}

/* ******************************************************************** */
/* apply mapping function                                               */
/* -------------------------------------------------------------------- */

int ML_Mapper_Apply(ML_Mapper *ml_mapper, double *invec, double *outvec)
{
   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_Set error : wrong object.\n");
      exit(1);
   }
   if ( ml_mapper->map != NULL ) {
      ml_mapper->map(ml_mapper->data, invec, outvec);
      return 0;
   } else { 
      return -1;
   }
}




#ifdef WKC
/* WKC
   EXTENSION FOR EPETRA CLASSES
*/
/* ******************************************************************** */
/* apply mapping function                                               */
/* -------------------------------------------------------------------- */
/* NOT CALLED! */

int ML_Mapper_Apply(ML_Mapper *ml_mapper, Epetra_MultiVector &ep_invec, 
                    Epetra_MultiVector &ep_outvec)
{
   double ** pp_invec;
   ep_invec.ExtractView ( &pp_invec );
   double ** pp_outvec;
   ep_outvec.ExtractView ( &pp_outvec );
   int iRetVal = -1;

   for ( int KK = 0 ; KK != ep_invec.NumVectors() ; KK++ ) {
      double *invec = ep_invec[KK];
      double *outvec = ep_outvec[KK];

   if ( ml_mapper->ML_id != ML_ID_MAPPER ) {
      printf("ML_Mapper_Set error : wrong object.\n");
      exit(1);
   }
   if ( ml_mapper->map != NULL ) {
      ml_mapper->map(ml_mapper->data, invec, outvec);
      iRetVal = 0;
   }
   }

   return iRetVal;
}
#endif
