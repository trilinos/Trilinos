/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for manipulating DVector objects                           */ 
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : August, 1998                                         */
/* ******************************************************************** */

#include "ml_vec.h"

/* ******************************************************************** */
/* initialize the Vector structure                                      */
/* -------------------------------------------------------------------- */

int ML_DVector_Create(ML_DVector **vec, ML_Comm *com)
{
   ML_DVector *vec2;

   ML_memory_alloc( (void **) vec, sizeof(ML_DVector), "dv1" );
   vec2                     = (*vec);
   vec2->ML_id              = ML_ID_VEC;
   vec2->SetOrLoad          = 0;
   vec2->VecLength          = 0;
   vec2->VecData            = NULL;
   vec2->comm               = com;
   return 0;
}

/* ******************************************************************** */
/* initialize the Vector structure                                      */
/* -------------------------------------------------------------------- */

int ML_DVector_Init(ML_DVector *vec)
{
   vec->ML_id     = ML_ID_VEC;
   vec->SetOrLoad = 0;
   vec->VecLength = 0;
   vec->VecData   = NULL;
   return 0;
}

/* ******************************************************************** */
/* deallocate the Vector structure                                      */
/* -------------------------------------------------------------------- */

int ML_DVector_Destroy( ML_DVector **vec ) 
{
   ML_DVector *vec2;

   vec2 = (*vec);
   if ( vec2->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_Destroy : wrong object. \n");
      exit(1);
   }
   if ( vec2->VecData != NULL && vec2->SetOrLoad == 2 ) {
      ML_memory_free( (void **) &(vec2->VecData) );
   }
   ML_memory_free ( (void**) vec );
   return 0;
}

/* ******************************************************************** */
/* deallocate the Vector structure                                      */
/* -------------------------------------------------------------------- */

int ML_DVector_Clean( ML_DVector *vec ) 
{
   if ( vec->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_DestroyData : wrong object. \n");
      exit(1);
   }
   if ( vec->VecData != NULL && vec->SetOrLoad == 2 ) {
      ML_memory_free( (void **) &(vec->VecData) );
   }
   return 0;
}

/* ******************************************************************** */
/* load data into the Vector structure                                  */
/* -------------------------------------------------------------------- */

int ML_DVector_LoadData( ML_DVector *vec, int n, double *data )
{
   int i, nbytes;

   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

   if ( vec->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_LoadData : wrong object (%d).\n",vec->ML_id);
      exit(1);
   }
   if ( n < 0 ) 
   {
      printf("ML_DVector_LoadData : length < 0. \n");
      exit(1);
   }

   /* ----------------------------------------------------- */
   /* deallocate previous vector                            */
   /* ----------------------------------------------------- */

   if ( vec->VecData != NULL && vec->SetOrLoad == 2 ) 
   {
      ML_memory_free( (void **) &(vec->VecData) );
   }

   /* ----------------------------------------------------- */
   /* allocate new vector and load                          */
   /* ----------------------------------------------------- */

   nbytes = n * sizeof(double);
   ML_memory_alloc( (void **) &(vec->VecData), nbytes, "dv2" );
   for ( i = 0; i < n; i ++ ) vec->VecData[i] = data[i];
   vec->VecLength = n;
   vec->SetOrLoad = 2;
   return 0;
}

/* ******************************************************************** */
/* load data pointer into the Vector structure                          */
/* -------------------------------------------------------------------- */

int  ML_DVector_SetData( ML_DVector *vec, int n, double *data )
{
   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

   if ( vec->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_SetData : wrong object. \n");
      exit(1);
   }
   if ( n < 0 ) 
   {
      printf("ML_DVector_SetData : length < 0. \n");
      exit(1);
   }

   /* ----------------------------------------------------- */
   /* deallocate previous vector                            */
   /* ----------------------------------------------------- */

   if ( vec->VecData != NULL && vec->SetOrLoad == 2 )
      ML_memory_free( (void **) &(vec->VecData) );

   /* ----------------------------------------------------- */
   /* set data pointer                                      */
   /* ----------------------------------------------------- */

   vec->VecData = data;
   vec->VecLength = n;
   vec->SetOrLoad = 1;
   return 0;
}

/* ******************************************************************** */
/* return length of vector                                              */
/* -------------------------------------------------------------------- */

int  ML_DVector_GetLength( ML_DVector *vec )
{
   if ( vec->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_GetLength : wrong object. \n");
      exit(1);
   }
   return vec->VecLength;
}

/* ******************************************************************** */
/* get data from the Vector structure                                   */
/* -------------------------------------------------------------------- */

int  ML_DVector_GetData( ML_DVector *vec, int *n, double *data )
{
   int i;

   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

   if ( vec->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_GetData : wrong object. \n");
      exit(1);
   }
   if ( vec->VecData == NULL ) 
   {
      printf("ML_DVector_GetData : no data. \n");
      exit(1);
   }

   /* ----------------------------------------------------- */
   /* return data                                           */
   /* ----------------------------------------------------- */

   (*n) = vec->VecLength;
   for ( i = 0; i < (*n); i++ ) data[i] = vec->VecData[i];

   return 0;
}

/* ******************************************************************** */
/* get data pointer from the Vector structure                           */
/* -------------------------------------------------------------------- */

int  ML_DVector_GetDataPtr( ML_DVector *vec, double **data )
{
   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

   if ( vec->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_GetDataPtr : wrong object (%d). \n", vec->ML_id);
      exit(1);
   }

   /* ----------------------------------------------------- */
   /* return data pointer                                   */
   /* ----------------------------------------------------- */

   (*data) = vec->VecData;

   return ( vec->VecLength );
}

/* ******************************************************************** */
/* check whether the object is a ML_DVector                             */
/* -------------------------------------------------------------------- */

int  ML_DVector_Check( ML_DVector *vec )
{
   if ( vec->ML_id != ML_ID_VEC ) 
      printf("ML_DVector_Check : wrong ID %ld (%d).\n",(long int) vec,vec->ML_id);
   else
      printf("ML_DVector_Check : passed (%ld) \n", (long int) vec );
   return 0;
}

/* ******************************************************************** */
/* vector copy                                                          */
/* -------------------------------------------------------------------- */

int ML_DVector_Copy( ML_DVector *src, ML_DVector *dest )
{
   int    i, leng;
   double *srcVec, *destVec;

   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

#ifdef ML_VECTOR_DEBUG 
   if ( src->ML_id != ML_ID_VEC || dest->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_Copy : wrong object (%d). \n", src->ML_id);
      exit(1);
   }
   if ( src->VecLength != dest->VecLength )
   {
      printf("ML_DVector_Copy error : different lengths\n");
      exit(1);
   }
#endif

   leng = src->VecLength;
   srcVec = src->VecData;
   destVec = dest->VecData;
   for ( i = 0; i < leng; i++ ) destVec[i] = srcVec[i];
   return 0;
}

/* ******************************************************************** */
/* vector scale                                                         */
/* -------------------------------------------------------------------- */

int ML_DVector_Scale( double alpha, ML_DVector *src )
{
   int    i, leng;
   double *srcVec;

   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

#ifdef ML_VECTOR_DEBUG 
   if ( src->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_Scale : wrong object (%d). \n", src->ML_id);
      exit(1);
   }
#endif

   leng = src->VecLength;
   srcVec = src->VecData;
   for ( i = 0; i < leng; i++ ) srcVec[i] *= alpha;
   return 0;
}

/* ******************************************************************** */
/* vector Ax+y                                                          */
/* -------------------------------------------------------------------- */

int ML_DVector_Axpy( double alpha, ML_DVector *src, ML_DVector *dest )
{
   int    i, leng;
   double *srcVec, *destVec;

   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

#ifdef ML_VECTOR_DEBUG 
   if ( src->ML_id != ML_ID_VEC || dest->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_Axpy : wrong object (%d). \n", src->ML_id);
      exit(1);
   }
   if ( src->VecLength != dest->VecLength )
   {
      printf("ML_DVector_Axpy error : different lengths\n");
      exit(1);
   }
#endif

   leng = src->VecLength;
   srcVec = src->VecData;
   destVec = dest->VecData;
   for ( i = 0; i < leng; i++ ) destVec[i] += (alpha * srcVec[i]);
   return 0;
}

/* ******************************************************************** */
/* vector Ay+x                                                          */
/* -------------------------------------------------------------------- */

int ML_DVector_Aypx( double alpha, ML_DVector *src, ML_DVector *dest )
{
   int    i, leng;
   double *srcVec, *destVec;

   /* ----------------------------------------------------- */
   /* error checking                                        */
   /* ----------------------------------------------------- */

#ifdef ML_VECTOR_DEBUG 
   if ( src->ML_id != ML_ID_VEC || dest->ML_id != ML_ID_VEC ) 
   {
      printf("ML_DVector_Aypx : wrong object (%d). \n", src->ML_id);
      exit(1);
   }
   if ( src->VecLength != dest->VecLength )
   {
      printf("ML_DVector_Aypx error : different lengths\n");
      exit(1);
   }
#endif

   leng = src->VecLength;
   srcVec = src->VecData;
   destVec = dest->VecData;
   for (i = 0; i < leng; i++) destVec[i] = srcVec[i] + alpha * destVec[i];
   return 0;
}

int ML_DVector_Print(int length, double *data, char *label, ML_Comm *comm)
{
   int i;
   char filename[128];
   FILE *fid;

   if (comm->ML_nprocs == 1)
      sprintf(filename,"%s.serial",label);
   else
      sprintf(filename,"%s.%d",label,comm->ML_mypid);
   printf("Writing matrix to file %s...\n",filename);
   fid = fopen(filename,"w");

   for (i=0; i<length; i++)
   {
      fprintf(fid,"%d      %20.15e\n",i,data[i]);
   }
   fclose(fid);
}
