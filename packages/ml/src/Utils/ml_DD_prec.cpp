#include "ml_common.h"
#include "ml_include.h"

extern "C" {

double ML_DD_OneLevel(ML_1Level *curr, double *sol, double *rhs,
			int approx_all_zeros, ML_Comm *comm,
			    int res_norm_or_not, ML *ml)
{

  ML_Smoother * post = curr->post_smoother;
  ML_Operator * Amat = curr->Amat;
  int lengf = Amat->outvec_leng;

  for ( int i = 0; i < lengf; i++ ) sol[i] = 0.0;
  
  ML_Smoother_Apply(post, lengf, sol, lengf, rhs, approx_all_zeros);

  return 0.0;

} /* ML_DD_OneLevel */

double ML_DD_Additive(ML_1Level *curr, double *sol, double *rhs,
			  int approx_all_zeros, ML_Comm *comm,
			  int res_norm_or_not, ML *ml)
{
   
  ML_Operator * Amat = curr->Amat;
  ML_Operator * Rmat = curr->Rmat;
  ML_Smoother * post      = curr->post_smoother;
  int lengf = Amat->outvec_leng;
  int lengc = Rmat->outvec_leng;

  double * sols = new double[lengf];
  double * rhs2 = new double[lengc];
  double * sol2 = new double[lengc];

  for ( int i = 0; i < lengf; i++ ) sols[i] = 0.0, sol[i] = 0.0;
  for ( int i = 0; i < lengc; i++ ) sol2[i] = 0.0, rhs2[i] = 0.0;

  ML_Smoother_Apply(post, lengf, sol, lengf, rhs, approx_all_zeros);

  ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, rhs, lengc, rhs2);

  ML_Smoother_Apply(Rmat->to->post_smoother, lengc, sol2, lengc, rhs2, ML_NONZERO);

  ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat, lengc, sol2, lengf, sols);

  for ( int i = 0; i < lengf; i++ ) sol[i] += sols[i];

  delete [] sols;
  delete [] rhs2;
  delete [] sol2;
  
  return 0.0;
}


double ML_DD_Hybrid(ML_1Level *curr, double *sol, double *rhs,
		    int approx_all_zeros, ML_Comm *comm,
		    int res_norm_or_not, ML *ml)
{

  ML_Operator *Amat, *Rmat;
  ML_Smoother  *post;
  //  ML_Smoother  *pre;
  //  ML_CSolve   *csolve;
  
  Amat     = curr->Amat;
  Rmat     = curr->Rmat;
  //  pre      = curr->pre_smoother;
  post     = curr->post_smoother;
  // csolve   = curr->csolve;
  int lengf    = Amat->outvec_leng;
  int lengc    = Rmat->outvec_leng;

  double * alpha1 = new double[lengf];
  double * alpha2 = new double[lengf];
  double * tmp_c  = new double[lengc];
  double * tmp2_c = new double[lengc];

  for ( int i = 0; i < lengf ; i++ ) alpha1[i] = 0.0, alpha2[i] = 0.0, sol[i] = 0.0;
  for ( int i = 0; i < lengc ; i++ ) tmp_c[i]  = 0.0, tmp2_c[i] = 0.0;

  // first step : rhs --> alpha1
  ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, rhs, lengc, tmp_c);
  ML_Smoother_Apply(Rmat->to->post_smoother, lengc, tmp2_c, lengc, tmp_c, ML_NONZERO);
  ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat, lengc, tmp2_c, lengf, alpha1);

  // second step
  ML_Operator_ApplyAndResetBdryPts(Amat, lengf, alpha1, lengc, sol);
  for ( int i = 0; i < lengf; i++ ) sol[i] = rhs[i] - sol[i];

  // sol --> alpha2
  ML_Smoother_Apply(post, lengf, alpha2, lengf, sol, approx_all_zeros);
  
  // third step

  for ( int i = 0; i < lengf ; i++ ) alpha1[i] += alpha2[i];
  for ( int i = 0; i < lengf ; i++ ) alpha2[i] = 0.0, sol[i] = 0.0;
  
  ML_Operator_ApplyAndResetBdryPts(Amat, lengf, alpha1, lengc, alpha2);
  
  for ( int i = 0; i < lengf; i++ ) alpha2[i] = rhs[i] - alpha2[i] ;
  //  alpha2 --> sol
  ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, alpha2, lengc, tmp_c);
  ML_Smoother_Apply(Rmat->to->post_smoother, lengc, tmp2_c, lengc, tmp_c, ML_NONZERO);
  ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat, lengc, tmp2_c, lengf, sol);
  
  // compose solution
  for ( int i = 0; i < lengf; i++ ) sol[i] += alpha1[i];

  delete [] alpha1;
  delete [] alpha2;
  delete [] tmp_c;
  delete [] tmp2_c;
  
  return 0.0;
}


double ML_DD_Hybrid_2(ML_1Level *curr, double *sol, double *rhs,
		      int approx_all_zeros, ML_Comm *comm,
		      int res_norm_or_not, ML *ml)
{
  ML_Operator *Amat, *Rmat;
  ML_Smoother *pre,  *post;
  // ML_CSolve   *csolve;
  
  Amat     = curr->Amat;
  Rmat     = curr->Rmat;
  pre      = curr->pre_smoother;
  post     = curr->post_smoother;
  // csolve   = curr->csolve;
  int lengf    = Amat->outvec_leng;
  int lengc    = Rmat->outvec_leng;

  double * alpha1 = new double[lengf];
  double * alpha2 = new double[lengf];
  double * tmp_c  = new double[lengc];
  double * tmp2_c = new double[lengc];

  for ( int i = 0; i < lengf ; i++ ) alpha1[i] = 0.0, alpha2[i] = 0.0, sol[i] = 0.0;
  for ( int i = 0; i < lengc ; i++ ) tmp_c[i]  = 0.0, tmp2_c[i] = 0.0;

  // first step  
  ML_Smoother_Apply(pre, lengf, alpha1, lengf, rhs, approx_all_zeros);

  // second step
  ML_Operator_ApplyAndResetBdryPts(Amat, lengf, alpha1, lengc, sol);
  for ( int i = 0; i < lengf; i++ ) sol[i] = rhs[i] - sol[i];
  
  ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, sol, lengc, tmp_c);
  ML_Smoother_Apply(Rmat->to->post_smoother, lengc, tmp2_c, lengc, tmp_c, ML_NONZERO);

  ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat, lengc, tmp2_c, lengf, alpha2);
  
  // third step

  for ( int i = 0; i < lengf ; i++ ) alpha1[i] += alpha2[i];
  for ( int i = 0; i < lengf ; i++ ) alpha2[i] = 0.0, sol[i] = 0.0;
  
  ML_Operator_ApplyAndResetBdryPts(Amat, lengf, alpha1, lengc, alpha2);
  
  for ( int i = 0; i < lengf; i++ ) alpha2[i] = rhs[i] - alpha2[i] ;
  ML_Smoother_Apply(post, lengf, sol, lengf, alpha2, approx_all_zeros);
  
  // compose solution
  for ( int i = 0; i < lengf; i++ ) sol[i] += alpha1[i];

  delete [] alpha1;
  delete [] alpha2;
  delete [] tmp_c;
  delete [] tmp2_c;

  return 0.0;
}

} // extern "C"
