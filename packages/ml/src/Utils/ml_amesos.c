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
#include "ml_include.h"
#include "ml_struct.h"
#include "ml_amesos_wrap.h"

/* ************************************************************************* */
/* ------------------------------------------------------------------------- */
/* generate the Amesos smoother                                              */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Amesos(ML_Smoother *sm,int inlen,double x[],int outlen,
                        double rhs[])
{

#ifdef HAVE_ML_AMESOS
  ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
  void *Amesos_Handle = smooth_ptr->smoother->data;

  ML_Amesos_Solve( Amesos_Handle, x, rhs ) ;
#else
  fprintf( stderr,
	   "*ML*ERR* you should configure with --with-ml_amesos\n"
	   "*ML*ERR* to use Amesos as smoother\n"
	   "*ML*ERR* (file=%s, line=%d)\n",
	   __FILE__,
	   __LINE__ );
  exit( EXIT_FAILURE );
#endif
  return 0;

} /* ML_Smoother_Amesos */


void ML_Smoother_Clean_Amesos(void *Amesos_Handle)
{

#ifdef HAVE_ML_AMESOS
  ML_Amesos_Destroy(Amesos_Handle);
#else
  fprintf( stderr,
	   "*ML*ERR* you should configure with --with-ml_amesos\n"
	   "*ML*ERR* to use Amesos as smoother\n"
	   "*ERR*ML* (file=%s, line=%d)\n",
	   __FILE__,
	   __LINE__ );
  exit( EXIT_FAILURE );
#endif
  return;
  
} /* ML_Smoother_Clean_Amesos */


int ML_Gen_Smoother_Amesos(ML *ml, int nl, int AmesosSolver,
			   int MaxProcs)
{
#ifdef HAVE_ML_AMESOS
   int            (*fun1)(ML_Smoother *, int, double *, int, double *);

   void *Amesos_Handle ;
   int status;
   char str[80];

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

     status = ML_Amesos_Gen( ml, nl, AmesosSolver, MaxProcs, &Amesos_Handle) ; 
     assert( status == 0 ) ; 

     switch( AmesosSolver ) {
     case ML_AMESOS_KLU:
       sprintf( str, "Amesos_KLU_%d", nl );
       break;
     case ML_AMESOS_UMFPACK:
       sprintf( str, "Amesos_UMFPACK_%d", nl );
       break;
     case ML_AMESOS_SUPERLUDIST:
       sprintf( str, "Amesos_SUPERLUDIST_%d", nl );
       break;
     case ML_AMESOS_MUMPS:
       sprintf( str, "Amesos_MUMPS_%d", nl );
       break;
     default:
      sprintf( str, "Amesos_%d", nl );
      break;

     } 

     status = ML_Smoother_Set(&(ml->post_smoother[nl]), 
			      (void *) Amesos_Handle, fun1, 1, 0.0,str);
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
