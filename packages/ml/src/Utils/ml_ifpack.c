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

  ML_Ifpack_Solve(Ifpack_Handle, x, rhs);
#else
  fprintf( stderr,
	   "*ML*ERR* you should configure ML with --enable-ifpack\n"
	   "*ML*ERR* to use Ifpack smoothers\n"
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
	   "*ML*ERR* you should configure ML with --enable-ifpack\n"
	   "*ML*ERR* to use Ifpack smoothers\n"
	   "*ML*ERR* (file=%s, line=%d)\n",
	   __FILE__,
	   __LINE__ );
  exit( EXIT_FAILURE );
#endif
  return;
  
} /* ML_Smoother_Clean_Ifpack */


/* for choices in options and params, see file
 * ml_ifpack_wrap.cpp 
 * */
int ML_Gen_Smoother_Ifpack(ML *ml, int nl, int pre_or_post,
			   int * options, double * params)
{

  /* da mettere POST or PRE or BOTH */
  
#ifdef HAVE_ML_IFPACK
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, status = 1;
   char str[80];
   void *Ifpack_Handle ;

   if (start_level < 0) {
      printf("ML_Gen_Smoother_Jacobi: cannot set smoother on level %d\n",
	     start_level);
      return 1;
   }

   fun = ML_Smoother_Ifpack;

   /* Creates IFPACK objects */

   status = ML_Ifpack_Gen(ml, nl, options, params, &Ifpack_Handle) ; 
   assert (status == 0); 

   /* Sets function pointers */

   if (pre_or_post == ML_PRESMOOTHER) {
     sprintf(str,"IFPACK_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]), (void*)Ifpack_Handle,
			      fun, 1, 0.0, str);
     ml->pre_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
     sprintf(str,"IFPACK_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]), 
			      (void*)Ifpack_Handle, fun, 1, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
   }
   else if (pre_or_post == ML_BOTH) {
     sprintf(str,"IFPACK_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]),
			      (void*)Ifpack_Handle,
			      fun, 1,  0.0, str);
     sprintf(str,"IFPACK_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]),
			      (void*)Ifpack_Handle, fun, 1, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
   }
   else 
     return(pr_error("ML_Gen_Smoother_Jacobi: unknown pre_or_post choice\n"));

   return(status);
#else
   printf("*ML*ERR* IFPACK not linked\n");
   ML_avoid_unused_param((void *) ml);
   ML_avoid_unused_param((void *) &nl);
   return(1);
#endif
}

int ML_Ifpack_Defaults(int options[], double params[])
{
  options[ML_IFPACK_TYPE] = ML_IFPACK_AMESOS;
  options[ML_IFPACK_OVERLAP] = 0;
  options[ML_IFPACK_LOCAL_PARTS] = 1;
  options[ML_IFPACK_SWEEPS] = 1;
  options[ML_IFPACK_BLOCK_OVERLAP] = 0;
  options[ML_IFPACK_LEVEL_OF_FILL] = 0;

  params[ML_IFPACK_DAMPING_FACTOR] = 1.0;
  return(0);
}
