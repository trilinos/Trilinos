/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/*       User Interface Functions                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include <stdlib.h>
#include <assert.h>
#include "ml_struct.h"

/* ************************************************************************* *
/* ------------------------------------------------------------------------- */
/* generate the sparse approximate inverse smoother */
/* ------------------------------------------------------------------------- */

#define ML_AMESOS
#ifdef ML_AMESOS
#include "ml_amesos_wrap.h"

int ML_Smoother_Amesos(void *sm,int inlen,double x[],int outlen,
                        double rhs[])
{
   ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
   void *Amesos_Handle = smooth_ptr->smoother->data;

   ML_Amesos_Solve( Amesos_Handle, x, rhs ) ;

   return 0;
}


void ML_Smoother_Clean_Amesos(void *Amesos_Handle)
{

   ML_Amesos_Destroy(Amesos_Handle);
}
#endif


int ML_Gen_Smoother_Amesos(ML *ml, int nl /* parameter list */)
{
#ifdef ML_AMESOS
   int            (*fun1)(void *, int, double *, int, double *);

   void *Amesos_Handle ;
   int status;

#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif

   if (nl < 0) {
      printf("ML_Gen_Smoother_Amesos: cannot set smoother on level %d\n",nl );
      return 1;
   }
	
   {

     fun1 = ML_Smoother_Amesos;

     printf(" ML_AMESOS      ml = %lx nl=%d ml[].Amat = %lx\n", ml , nl, ml[nl].Amat );
     status = ML_Amesos_Gen( ml, nl, &Amesos_Handle) ; 
     assert( status == 0 ) ; 

     status = ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL,
			      (void *) Amesos_Handle, fun1, NULL, 1, 0.0,NULL);
     assert( status == 0 ) ; 
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Amesos;

#ifdef ML_TIMING
         ml->post_smoother[nl].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[nl].build_time;
#endif

   /* note: in free, post and pre are the same */
   }

   return(status);
#else
   printf("Amesos not linked\n");
   ML_avoid_unused_param((void *) ml);
   ML_avoid_unused_param((void *) &nl);
   return(1);
#endif
}

/* ************************************************************************* */
/* Sparse approximate inverse smoother                                       */
/* ------------------------------------------------------------------------- */

