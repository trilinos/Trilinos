/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/************************************************************************/
/*          Utilities for Trilinos/ML users                             */
/*----------------------------------------------------------------------*/
/* Authors : Mike Heroux (SNL)                                          */
/*           Jonathan Hu  (SNL)                                         */
/*           Ray Tuminaro (SNL)                                         */
/*           Marzio Sala (SNL)                                          */
/************************************************************************/


#include "ml_common.h"
#include "ml_epetra_preconditioner.h"

extern "C" {

extern double ML_DD_OneLevel(ML_1Level *curr, double *sol, double *rhs,
			     int approx_all_zeros, ML_Comm *comm,
			     int res_norm_or_not, ML *ml);
extern double ML_DD_Additive(ML_1Level *curr, double *sol, double *rhs,
			     int approx_all_zeros, ML_Comm *comm,
			     int res_norm_or_not, ML *ml);
extern double ML_DD_Hybrid(ML_1Level *curr, double *sol, double *rhs,
			   int approx_all_zeros, ML_Comm *comm,
			   int res_norm_or_not, ML *ml);
extern double ML_DD_Hybrid_2(ML_1Level *curr, double *sol, double *rhs,
			     int approx_all_zeros, ML_Comm *comm,
			     int res_norm_or_not, ML *ml);
}

extern "C" {

double ML_DD_OneLevel(ML_1Level *curr, double *sol, double *rhs,
			    int approx_all_zeros, ML_Comm *comm,
			    int res_norm_or_not, ML *ml)
{
  ML_Smoother * post = curr->post_smoother;
  ML_Operator * Amat = curr->Amat;
  int lengf = Amat->outvec_leng;
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
  ML_Smoother *pre,  *post;
  ML_CSolve   *csolve;
  
  double * sols;
   
  Amat     = curr->Amat;
  Rmat     = curr->Rmat;
  pre      = curr->pre_smoother;
  post     = curr->post_smoother;
  csolve   = curr->csolve;
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
  ML_CSolve   *csolve;
  
  double * sols;
   
  Amat     = curr->Amat;
  Rmat     = curr->Rmat;
  pre      = curr->pre_smoother;
  post     = curr->post_smoother;
  csolve   = curr->csolve;
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

} /* extern "C" */





#if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"
#include "ml_agg_ParMETIS.h"

void Epetra_ML_Preconditioner::Destroy_ML_Preconditioner()
{
  
  if( agg_ != 0 ) {
    ML_Aggregate_Destroy(&agg_); agg_ = 0;
  }
  if( ml_ != 0 ) {
    ML_Destroy(&ml_);
    ml_ = 0;
  }
  
}

Epetra_ML_Preconditioner::Epetra_ML_Preconditioner( const Epetra_RowMatrix & RowMatrix,
						    ParameterList & List) :
  Comm_(RowMatrix.Comm()),
  DomainMap_(RowMatrix.OperatorDomainMap()),
  RangeMap_(RowMatrix.OperatorRangeMap()),
  List_(List),
  MaxLevels_(20)
{

  Label_ = new char[80];
  
  for( int i=0 ; i<MaxLevels_ ; ++i )  {
    AZ_defaults(SmootherOptions_[i],SmootherParams_[i]);
    SmootherOptions_[i][AZ_precond] = AZ_dom_decomp;
    SmootherOptions_[i][AZ_subdomain_solve] = AZ_lu;
  }
  

#ifdef HAVE_MPI
  const Epetra_MpiComm * MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm_);
  AZ_set_proc_config(ProcConfig_,MpiComm->Comm());
#else
  AZ_set_proc_config(ProcConfig_,AZ_NOT_MPI);
#endif
  
  NumLevels_ = List_.getParameter("max num levels",10);  

  int OutputLevel = List_.getParameter("output level", 10);  
  ML_Set_PrintLevel(OutputLevel);

  ML_Create(&ml_,NumLevels_);

  int NumMyRows, osize;

  NumMyRows = RowMatrix.NumMyRows();
  int N_ghost = RowMatrix.NumMyCols() - NumMyRows;

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  ML_Init_Amatrix(ml_,0,NumMyRows, NumMyRows, (void *) &RowMatrix);
  ML_Set_Amatrix_Getrow(ml_, 0, Epetra_ML_getrow,
			Epetra_ML_comm_wrapper, NumMyRows+N_ghost);

  ML_Set_Amatrix_Matvec(ml_, 0, Epetra_ML_matvec);

  ML_Aggregate_Create(&agg_);
  
  /* ********************************************************************** */
  /* set null space                                                         */
  /* ********************************************************************** */

  int NumPDEEqns = List_.getParameter("num pde eqns", 1);
  ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns,NumPDEEqns,NULL,NumMyRows);

  /* ********************************************************************** */
  /* pick up coarsening strategy. METIS and ParMETIS requires additional    */
  /* lines, as we have to set the number of aggregates                      */
  /* ********************************************************************** */

  char parameter[80];
  int value = -777;
  
  for( int level=0 ; level<NumLevels_-1 ; ++level ) {  

    sprintf(parameter,"aggregation scheme, level %d",level);
    string CoarsenScheme = List_.getParameter(parameter,"METIS");

    if( CoarsenScheme == "METIS" )
      ML_Aggregate_Set_CoarsenSchemeLevel_METIS(level,NumLevels_,agg_);
    else if( CoarsenScheme == "ParMETIS" ) 
      ML_Aggregate_Set_CoarsenSchemeLevel_ParMETIS(level,NumLevels_,agg_);
    else if( CoarsenScheme == "MIS" ) 
      ML_Aggregate_Set_CoarsenSchemeLevel_MIS(level,NumLevels_,agg_);
    else if(  CoarsenScheme == "Uncoupled" ) 
      ML_Aggregate_Set_CoarsenSchemeLevel_Uncoupled(level,NumLevels_,agg_);
    else if(  CoarsenScheme == "Coupled" ) 
      ML_Aggregate_Set_CoarsenSchemeLevel_Coupled(level,NumLevels_,agg_);
    else {
      if( Comm_.MyPID() == 0 ) {
	cout << "ML_Preconditioner:ERROR: specified options ("
	     << CoarsenScheme << ") not valid." << endl;
	cout <<"ML_Preconditioner:ERROR: Using METIS instead" << endl;
      }
      ML_Aggregate_Set_CoarsenSchemeLevel_METIS(level,NumLevels_,agg_);
    } 

    bool isSet = false;
    
    sprintf(parameter,"number of global aggregates, level %d", level);
    if( List_.isParameter(parameter) == true ) {
      value = -777;
      value = List_.getParameter(parameter,value);
      if( value != -777 ) {
	ML_Aggregate_Set_GlobalNumber(ml_,agg_,level,value );
	isSet = true;
      }
    }

    sprintf(parameter,"number of local aggregates, level %d", level);
    if( List_.isParameter(parameter) == true ) {
      value = -777;
      value = List_.getParameter(parameter,value);
      if( value != -777 ) {
	ML_Aggregate_Set_GlobalNumber(ml_,agg_,level,value );
	isSet = true;
      }
    }

    sprintf(parameter,"number of nodes per aggregate, level %d", level);
    if( List_.isParameter(parameter) == true ) {
      value = -777;
      value = List_.getParameter(parameter,value);
      if( value != -777 ) {
	ML_Aggregate_Set_NodesPerAggr(ml_,agg_,level,value );
	isSet = true;
      }
    }

    if( isSet == false ) {
      // put some default value
      ML_Aggregate_Set_NodesPerAggr(ml_,agg_,level,1);
    }

  } /* for */

  /* ********************************************************************** */
  /* minor settings                                                         */
  /* ********************************************************************** */

  double DampingFactor = 0.01;
  if( List_.isParameter("damping factor") == true ) 
    DampingFactor = List_.getParameter("damping factor", DampingFactor);
  ML_Aggregate_Set_DampingFactor( agg_, DampingFactor ); 
  
  int MaxCoarseSize = 512;
  if( List_.isParameter("max coarse size") == true ) 
    MaxCoarseSize = List_.getParameter("max coarse size", MaxCoarseSize);
  ML_Aggregate_Set_MaxCoarseSize(agg_, MaxCoarseSize );

  double Threshold = 0.0;
  if( List_.isParameter("threshold") == true ) 
    Threshold = List_.getParameter("threshold", Threshold);
  ML_Aggregate_Set_Threshold(agg_,Threshold);
  
  int ReqAggrePerProc = 128;
  if( List_.isParameter("req aggregates per process") == true)
    ReqAggrePerProc = List_.getParameter("req aggregates per process", ReqAggrePerProc);
  ML_Aggregate_Set_ReqLocalCoarseSize( ml_, agg_, -1, ReqAggrePerProc);

  if( 5 < ML_Get_PrintLevel() && ProcConfig_[AZ_node] == 0 ) {
    cout << "Damping Factor = " << DampingFactor << endl;
    cout << "Threshold for aggregation =" << Threshold << endl;
    cout << "Max coarse size = " << MaxCoarseSize << endl;
    cout << "Req aggregates per process ParMETIS) " << ReqAggrePerProc << endl;
  }

  /************************************************************************/
  /* Build hierarchy using smoothed aggregation.                          */
  /* NOTE: the first level is 0. This means that I have Nevels, the last  */
  /* one begin Nlevels-1. This also means that I have Nlevels-2 smoothers */
  /* plus the "smoother" on the coarse grid (handled independentely).     */
  /*----------------------------------------------------------------------*/

  NumLevels_ = ML_Gen_MGHierarchy_UsingAggregation(ml_, 0, ML_INCREASING, agg_);

  int num_smoother_steps = 1;
  double omega = 1.0;
  int pre_or_post;
  int * smoother_options =  NULL;
  double * smoother_params  = NULL;
  string PreOrPostSmoother = "post";
  string Smoother = "Jacobi";
  bool AztecSmootherAsASolver = false;
  int aztec_its;
  
  for( int level=0 ; level<NumLevels_-1 ; ++level ) {

    sprintf(parameter,"num smoother steps, level %d",level );
    num_smoother_steps = List_.getParameter(parameter,num_smoother_steps);

    sprintf(parameter,"omega, level %d",level );
    omega = List_.getParameter(parameter,omega);
    
    sprintf(parameter,"pre or post smoother, level %d",level );
    PreOrPostSmoother = List_.getParameter(parameter, PreOrPostSmoother);
    
    if( PreOrPostSmoother == "post" ) {
      pre_or_post = ML_POSTSMOOTHER;
    } else if( PreOrPostSmoother == "pre" ) {
      pre_or_post = ML_PRESMOOTHER;
    } else {
      pre_or_post = ML_BOTH;
    }

    sprintf(parameter,"smoother, level %d", level);
    Smoother = List_.getParameter(parameter,Smoother);
    
    if( Smoother == "Jacobi" ) 
      ML_Gen_Smoother_Jacobi(ml_, level, pre_or_post,
			     num_smoother_steps, omega);
    else if( Smoother == "Gauss-Seidel" ) 
      ML_Gen_Smoother_GaussSeidel(ml_, level, pre_or_post,
				  num_smoother_steps, omega);
    else if( Smoother == "symmetric Gauss-Seidel" ) 
      ML_Gen_Smoother_SymGaussSeidel(ml_, level, pre_or_post,
				     num_smoother_steps, omega);
    else if( Smoother == "Block Gauss-Seidel" )
      ML_Gen_Smoother_BlockGaussSeidel(ml_, level, pre_or_post,
				       num_smoother_steps, omega, NumPDEEqns);
    else if( Smoother == "MLS" )
      ML_Gen_Smoother_MLS(ml_, level, pre_or_post, 30.,
			  num_smoother_steps);
    else if( Smoother == "aztec" ) {
      sprintf(parameter,"aztec smoother options, level %d", level);
      smoother_options = List_.getParameter(parameter, smoother_options);
      sprintf(parameter,"aztec smoother params, level %d", level);
      smoother_params = List_.getParameter(parameter, smoother_params);

      sprintf(parameter,"aztec smoother as solver, level %d", level);
      AztecSmootherAsASolver = List_.getParameter(parameter,AztecSmootherAsASolver);
      if( AztecSmootherAsASolver == false ) aztec_its = AZ_ONLY_PRECONDITIONER;
      else                                  aztec_its = num_smoother_steps;
      
      if( smoother_options == NULL || smoother_params == NULL ) {
	smoother_options = SmootherOptions_[level];
	smoother_params = SmootherParams_[level];
      }
      ML_Gen_SmootherAztec(ml_, level, smoother_options, smoother_params,
			   ProcConfig_, SmootherStatus_,
			   aztec_its, pre_or_post, NULL);
    } else if( Smoother == "ifpack" ) {
      // get ifpack options from list ??? pass list ???
      ML_Gen_Smoother_Ifpack(ml_, level, pre_or_post, NULL, NULL);
    }
    
  } /* for */
  
  /* ********************************************************************** */
  /* solution of the coarse problem                                         */
  /* ********************************************************************** */

  string CoarseSolution = List_.getParameter("coarse solver", "Amesos_KLU");
  int NumSmootherSteps = List_.getParameter("coarse num smoother steps", 1);
  double Omega = List_.getParameter("coarse damping factor", 0.67);

  int MaxProcs = List_.getParameter("max processes for coarse solution", -1);
  
  if( CoarseSolution == "Jacobi" ) 
    ML_Gen_Smoother_Jacobi(ml_, NumLevels_-1, ML_POSTSMOOTHER,
			   NumSmootherSteps, Omega);
  else if( CoarseSolution == "Gauss-Seidel" ) 
    ML_Gen_Smoother_GaussSeidel(ml_, NumLevels_-1, ML_POSTSMOOTHER,
			        NumSmootherSteps, Omega);
  else if( CoarseSolution == "SuperLU" ) 
    ML_Gen_CoarseSolverSuperLU( ml_, NumLevels_-1);
  else if( CoarseSolution == "Amesos_KLU" )
    ML_Gen_Smoother_Amesos(ml_, NumLevels_-1, ML_AMESOS_KLU, MaxProcs);
  else if( CoarseSolution == "Amesos_UMFPACK" )
    ML_Gen_Smoother_Amesos(ml_, NumLevels_-1, ML_AMESOS_UMFPACK, MaxProcs);
  else if(  CoarseSolution == "Amesos_Superludist" )
    ML_Gen_Smoother_Amesos(ml_, NumLevels_-1, ML_AMESOS_SUPERLUDIST, MaxProcs);
  else if( CoarseSolution == "Amesos_MUMPS" )
    ML_Gen_Smoother_Amesos(ml_, NumLevels_-1, ML_AMESOS_MUMPS, MaxProcs);
  else
    ML_Gen_Smoother_Amesos(ml_, NumLevels_-1, ML_AMESOS_KLU, MaxProcs);

  ML_Gen_Solver(ml_, ML_MGV, 0, NumLevels_-1);
  
  ownership_ = false;

  CreateLabel();
  
  /* ********************************************************************** */
  /* Specific   preconditioners                                             */
  /* NOTE: the two-level DD preconditioners are kind of experimental!       */
  /* I suppose that the user knows what he/she is doing.... No real checks  */
  /* are performed. Note also that the coarse solver is somehow supposed to */
  /* be implemented as a post smoother (FIXME)                              */
  /* ********************************************************************** */

  string str = List_.getParameter("prec type","MGV");
  if( str == "fine level postsmoothing" ) {
    
    sprintf(Label_, "1-level postsmoothing only");
    ml_->ML_scheme = ML_ONE_LEVEL_DD;

  } else if( str == "two-level additive DD" ) {

    sprintf(Label_, "two-level additive DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_ADD;
    if( NumLevels_ != 2 ) {
      if( Comm_.MyPID() == 0 ) {
	cerr << "Error : You asked for `two-level additive DD' but you don't have" << endl
	     << "Error : exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level hybrid DD") {

    sprintf(Label_, "two-level hybrid DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID;
    if( NumLevels_ != 2 ) {
      if( Comm_.MyPID() == 0 ) {
	cerr << "Error : You asked for `two-level hybrid DD' but you don't have" << endl
	     << "Error : exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level hybrid DD (2)") {

    sprintf(Label_, "two-level hybrid DD (2)");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID_2;
    if( NumLevels_ != 2 ) {
      if( Comm_.MyPID() == 0 ) {
	cerr << "Error : You asked for `two-level hybrid DD (2)' but you don't have" << endl
	     << "Error : exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else {

    if( Comm_.MyPID() == 0 ) {
      cerr << "Error : `prec type' has an incorrect value. It should be" << endl
	   << "Error : <1-level postsmoothing only> / <two-level additive DD>" << endl
	   << "Error : <two-level hybrid DD> / <two-level hybrid DD (2)>" << endl;
    }
    exit( EXIT_FAILURE );
    
  }

}

int Epetra_ML_Preconditioner::CreateLabel()
{
  
  int i = ml_->ML_finest_level;
  char finest[80];
  char coarsest[80];
  char label[80];
  finest[0] = '\0';
  coarsest[0] = '\0';
  
  if (ml_->pre_smoother[i].ML_id != ML_EMPTY) 
    sprintf(finest, "%s", ml_->pre_smoother[i].label);
  if (ml_->post_smoother[i].ML_id != ML_EMPTY) 
    sprintf(finest, "%s/%s", finest,ml_->post_smoother[i].label);
  
  if (i != ml_->ML_coarsest_level) {
    i = ml_->ML_coarsest_level;
    if ( ML_CSolve_Check( &(ml_->csolve[i]) ) == 1 ) {
	sprintf(coarsest, "%s", ml_->csolve[i].label);
    }
    
    else {
      if (ml_->pre_smoother[i].ML_id != ML_EMPTY) 
	sprintf(coarsest, "%s", ml_->pre_smoother[i].label);
      if (ml_->post_smoother[i].ML_id != ML_EMPTY) 
	sprintf(coarsest, "%s/%s", coarsest,ml_->post_smoother[i].label);
    }
  }
  sprintf(Label_,"%d level AMG (%s, %s)", ml_->ML_num_actual_levels, finest, coarsest);
  
  return 0;
    
}

int Epetra_ML_Preconditioner::ApplyInverse(const Epetra_MultiVector& X,
					   Epetra_MultiVector& Y) const
{

  if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled
                              // in solver or if X = Y
  Y.PutScalar(0.0); // Always start with Y = 0

  // ML_iterate doesn't handle multivectors, so extract and iterate one at
  // a time on them.
  double **xvectors;
  double **yvectors;
  EPETRA_CHK_ERR(xtmp.ExtractView(&xvectors));
  EPETRA_CHK_ERR(Y.ExtractView(&yvectors));

  //note: solver_ is the ML handle

  for (int i=0; i < X.NumVectors(); i++) {
    switch(ml_->ML_scheme) {
    case(ML_MGFULLV):
      ML_Solve_MGFull(ml_, xvectors[i], yvectors[i]); 
      break;
    case(ML_SAAMG): //Marian Brezina's solver
      ML_Solve_AMGV(ml_, xvectors[i], yvectors[i]); 
      break;
    case(ML_ONE_LEVEL_DD):
      ML_DD_OneLevel(&(ml_->SingleLevel[ml_->ML_finest_level]),
		     yvectors[i], xvectors[i],
		     ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);
      break;
    case(ML_TWO_LEVEL_DD_ADD):
      ML_DD_Additive(&(ml_->SingleLevel[ml_->ML_finest_level]),
		     yvectors[i], xvectors[i],
		     ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);
      break;
    case(ML_TWO_LEVEL_DD_HYBRID):
      ML_DD_Hybrid(&(ml_->SingleLevel[ml_->ML_finest_level]),
		   yvectors[i], xvectors[i],
		   ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);    
      break;
    case(ML_TWO_LEVEL_DD_HYBRID_2):
      ML_DD_Hybrid_2(&(ml_->SingleLevel[ml_->ML_finest_level]),
		     yvectors[i], xvectors[i],
		     ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);    
      break;
    default:
      ML_Solve_MGV(ml_, xvectors[i], yvectors[i]); 
    }
  }
  
  return 0;
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
