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
#include "ml_config.h"
#include "ml_include.h"
#include "ml_struct.h"
#include "ml_ifpack_wrap.h"

/* ------------------------------------------------------------------------- */
/* generate the Ifpack smoother                                              */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
		       double rhs[])
{

#ifdef HAVE_ML_IFPACK
  ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
  void *Ifpack_Handle = smooth_ptr->smoother->data;

  ML_Ifpack_Solve( Ifpack_Handle, x, rhs ) ;
#else
  fprintf( stderr,
	   "*ML*ERR* you should configure with --with-ml_ifpack\n"
	   "*ML*ERR* to use Ifpack as smoother\n"
	   "*ML*ERR* (file=%s, line=%d)\n",
	   __FILE__,
	   __LINE__ );
  exit( EXIT_FAILURE );
#endif
  return 0;

} /* ML_Smoother_Ifpack */


void ML_Smoother_Clean_Ifpack(void *Ifpack_Handle)
{

#ifdef HAVE_ML_IFPACK
  ML_Ifpack_Destroy(Ifpack_Handle);
#else
  fprintf( stderr,
	   "*ML*ERR* you should configure with --with-ml_ifpack\n"
	   "*ML*ERR* to use Ifpack as smoother\n"
	   "*ML*ERR* (file=%s, line=%d)\n",
	   __FILE__,
	   __LINE__ );
  exit( EXIT_FAILURE );
#endif
  return;
  
} /* ML_Smoother_Clean_Ifpack */


int ML_Gen_Smoother_Ifpack(ML *ml, int nl, int choice, int * options, double * params)
{

  /* da mettere POST or PRE or BOTH */
  
#ifdef HAVE_ML_IFPACK
   int            (*fun1)(ML_Smoother *, int, double *, int, double *);

   void *Ifpack_Handle ;
   int status;

#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif

   if (nl < 0) {
      printf("ML_Gen_Smoother_Ifpack: cannot set smoother on level %d\n",nl );
      return 1;
   }
	
   {

     fun1 = ML_Smoother_Ifpack;

     status = ML_Ifpack_Gen( ml, nl, choice, options, params, &Ifpack_Handle) ; 
     assert( status == 0 ) ; 

     status = ML_Smoother_Set(&(ml->post_smoother[nl]), 
			      (void *) Ifpack_Handle, fun1, 1, 0.0,NULL);
     assert( status == 0 ) ; 
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;

#ifdef ML_TIMING
     ml->post_smoother[nl].build_time = GetClock() - t0;
     ml->timing->total_build_time   += ml->post_smoother[nl].build_time;
#endif

     /* note: in free, post and pre are the same */
   }

   return(status);
#else
   printf("Ifpack not linked\n");
   ML_avoid_unused_param((void *) ml);
   ML_avoid_unused_param((void *) &nl);
   return(1);
#endif
}
