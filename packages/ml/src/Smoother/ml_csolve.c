/* *********************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact person, */
/* and disclaimer.                                                         */
/* *********************************************************************** */

/* *********************************************************************** */
/* Functions for the ML_CSolve structure                                   */
/* *********************************************************************** */
/* Author        : Charles Tong (LLNL)                                     */
/* Date          : March, 1999                                             */
/* *********************************************************************** */

#include "ml_csolve.h"
#include <string.h>

/* *********************************************************************** */
/* references to external functions                                        */
/* ----------------------------------------------------------------------- */

extern int ML_Clean_CSolveSuperLU( void *, ML_CSolveFunc *); 
extern int SuperLU_Solve(void *,int ,double *,int ,double *);

/* *********************************************************************** */
/* Constructor                                                             */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Create(ML_CSolve **cs)
{
   ML_CSolve *ml_cs;

   ML_memory_alloc((void**) cs, sizeof(ML_CSolve), "CS1" );
   ml_cs = (*cs);
   ml_cs->ML_id = ML_ID_CSOLVE; 
   ml_cs->my_level = NULL;
   ml_cs->ntimes = 0;
   ml_cs->tol = 0;
   ml_cs->data = NULL;
   ML_memory_alloc((void**) &(ml_cs->func),sizeof(ML_CSolveFunc),"CF1" );
   ml_cs->func->ML_id = ML_EMPTY;
   ml_cs->func->internal = NULL;
   ml_cs->func->external = NULL;
   ml_cs->build_time = 0.0;
   ml_cs->apply_time = 0.0;
   ml_cs->label      = NULL;
   ml_cs->data_destroy = NULL;
   return 0;
} 

/* *********************************************************************** */
/* Initialize                                                              */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Init(ML_CSolve *ml_cs)
{
   ml_cs->ML_id = ML_ID_CSOLVE; 
   ml_cs->my_level = NULL;
   ml_cs->ntimes = 0;
   ml_cs->tol = 0;
   ml_cs->data = NULL;
   ML_memory_alloc((void**) &(ml_cs->func),sizeof(ML_CSolveFunc),"CF2" );
   ml_cs->func->ML_id = ML_EMPTY;
   ml_cs->func->internal = NULL;
   ml_cs->func->external = NULL;
   ml_cs->build_time = 0.0;
   ml_cs->apply_time = 0.0;
   ml_cs->label      = NULL;
   ml_cs->data_destroy = NULL;
   return 0;
} 

/* *********************************************************************** */
/* Destructor                                                              */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Destroy(ML_CSolve **cs)
{
   ML_CSolve_Clean(*cs);
   ML_memory_free( (void**) cs );
   (*cs) = NULL; 
   return 0;
}

/* *********************************************************************** */
/* partial destructor                                                      */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Clean(ML_CSolve *ml_cs)
{
#ifdef ML_DETAILED_TIMING
   double t1;
#endif

   if ( ml_cs->ML_id != ML_ID_CSOLVE ) 
   {
      printf("ML_CSolve_Clean error : wrong object.\n");
      exit(-1);
   }
#ifdef ML_DETAILED_TIMING
   if (ml_cs->label != NULL) {
      t1 = ML_gsum_double(ml_cs->build_time, global_comm);
      t1 = t1/((double) global_comm->ML_nprocs);
      if ( (global_comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Build time for %s\t= %e\n",ml_cs->label,t1);
   }

   if  (ml_cs->label != NULL) {
      t1 = ML_gsum_double(ml_cs->apply_time, global_comm);
      t1 = t1/((double) global_comm->ML_nprocs);
      if ( (global_comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s\t= %e\n",ml_cs->label,t1);
    }
#endif

   ml_cs->ML_id = -1; 
   ml_cs->my_level = NULL;
   ml_cs->ntimes = 0;
   ml_cs->tol = 0;
   if (ml_cs->data_destroy != NULL)
      ml_cs->data_destroy( ml_cs->data );

   if ( ml_cs->func->internal == SuperLU_Solve ) {
       ML_Clean_CSolveSuperLU( ml_cs->data, ml_cs->func );
   }
   ML_memory_free( (void**) &(ml_cs->func) );
   ml_cs->data = NULL;
   ml_cs->func = NULL;
   ml_cs->data_destroy = NULL;
   if (ml_cs->label != NULL) { free(ml_cs->label); ml_cs->label = NULL; }
   return 0;
}

/* *********************************************************************** */
/* Check to see if csolve function exists                                  */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Check(ML_CSolve *ml_cs)
{
   if ( ml_cs->ML_id != ML_ID_CSOLVE )
   {
      printf("ML_CSolve_Check : wrong object.\n");
      exit(1);
   }
   if ( ml_cs->func->ML_id == ML_EMPTY ) return 0;
   else                                  return 1;
}

/* *********************************************************************** */
/* set the 1Level parameter                                                */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Set_1Level(ML_CSolve *ml_cs, ML_1Level *mylevel)
{
   if ( ml_cs->ML_id != ML_ID_CSOLVE )
   {
      printf("ML_CSolve_Set_1Level error : wrong object.\n");
      exit(-1);
   }
   ml_cs->my_level = mylevel;
   return 0;
}

/* *********************************************************************** */
/* perform solve                                                           */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Apply(ML_CSolve *csolve, int inlen, double din[], 
                    int outlen, double dout[])
{
#if defined(ML_TIMING) || defined(ML_DETAILED_TIMING)
   double t0;
   t0 = GetClock();
#endif
   if (csolve->func->ML_id == ML_EMPTY) 
      pr_error("ML_CSolve_Apply error : coarse solver not defined\n");

   if (csolve->func->ML_id == ML_EXTERNAL)
        csolve->func->external(csolve->data, inlen, din, outlen, dout);
   else csolve->func->internal(csolve->data, inlen, din, outlen, dout);
#if defined(ML_TIMING) || defined(ML_DETAILED_TIMING)
   csolve->apply_time += (GetClock() - t0);
#endif
   return 0;
}

/* ******************************************************************** */
/* give a label to this object                                          */
/* ******************************************************************** */

int ML_CSolve_Set_Label( ML_CSolve *csolve, char *label)
{
  int size;

   size = strlen(label) + 1;
   csolve->label = (char *) malloc(size*sizeof(char));
   if (csolve->label == NULL) pr_error("Not enough space in ML_CSolve_Set_Label\n");
   strncpy(csolve->label,label,size);
   return(1);
}

