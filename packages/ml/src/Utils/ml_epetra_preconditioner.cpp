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
  ML_Smoother *pre,  *post;
  ML_CSolve   *csolve;
  
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
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <cstring>
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"

#include "ml_epetra_preconditioner.h"
#include "ml_agg_ParMETIS.h"

#include "ml_anasazi.h"

#ifdef HAVE_ML_TRIUTILS
#include "Trilinos_Util_CommandLineParser.h"
#endif

using namespace Teuchos;
using namespace ML_Epetra;
  
// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::PrintLine() 
{
  cout << "--------------------------------------------------------------------------------" << endl;
}

// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::Destroy_ML_Preconditioner()
{

  char parameter[80];
  
  {
    int NumDestroy = OutputList_.get("number of destruction phases", 0);
    OutputList_.set("number of destruction phases", ++NumDestroy);
  }
  
  if( agg_ != 0 ) {
    ML_Aggregate_Destroy(&agg_); agg_ = 0;
  }
  if( ml_ != 0 ) {
    ML_Destroy(&ml_);
    ml_ = 0;
  }
  if( ml_nodes_ != 0 ) {
    ML_Destroy(&ml_nodes_);
    ml_nodes_ = 0;
  }
  if( ml_edges_ != 0 ) {
    ML_Destroy(&ml_edges_);
    ml_edges_ = 0;
  }

  if( Label_ ) delete Label_;
  
  if( LevelID_ != 0 ) delete [] LevelID_;
  
  if( ML_Get_PrintLevel() > 5 && Comm().MyPID() == 0 && NumApplications_ ) {

    double AvgTime = ApplicationTime_/NumApplications_;
    
    // stick data in OutputList

    int i;
    double t;
    
    sprintf(parameter,"%stime: total application", Prefix_);
    OutputList_.set(parameter, FirstApplicationTime_+ApplicationTime_);

    sprintf(parameter,"%stime: first application", Prefix_);
    OutputList_.set(parameter, FirstApplicationTime_);

    sprintf(parameter,"%stime: construction", Prefix_);
    OutputList_.set(parameter, ConstructionTime_);

    sprintf(parameter,"%snumber of applications", Prefix_);
    OutputList_.set(parameter,NumApplications_);
    
    // print on screen
    
    PrintLine();
    double TotalTime = FirstApplicationTime_ + ApplicationTime_;
    cout << PrintMsg_ << "Time for all applications = " << TotalTime << " (s)" << endl;
    cout << PrintMsg_ << "Time for first application = " << FirstApplicationTime_ << " (s)" << endl;
    cout << PrintMsg_ << "Each of " << NumApplications_
	 << " applications took  " << TotalTime/NumApplications_ << " (s)" << endl;
    PrintLine();
  }

  if( NullSpaceToFree_ != 0 ) delete [] NullSpaceToFree_;

  if( RowMatrixAllocated_ ) delete RowMatrixAllocated_;
  
  IsComputePreconditionerOK_ = false;

}

// ================================================ ====== ==== ==== == =

MultiLevelPreconditioner::MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
						   const bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  RowMatrixAllocated_(0)
{

  sprintf(Prefix_,"");
  
  ParameterList NewList;
  List_ = NewList;
  ML_Epetra::SetDefaults("DD",List_,Prefix_,SmootherOptions_,SmootherParams_);
    
  Initialize();

  // construct hierarchy
  if( ComputePrec == true ) ComputePreconditioner();
}
  
// ================================================ ====== ==== ==== == =

MultiLevelPreconditioner::MultiLevelPreconditioner( const Epetra_RowMatrix & RowMatrix,
						    const ParameterList & List, const bool ComputePrec,
						    const char Prefix[] ) :
  RowMatrix_(&RowMatrix),
  RowMatrixAllocated_(0)
{
  sprintf(Prefix_,"%s",Prefix);

  List_ = List;

  Initialize();

  // construct hierarchy
  if( ComputePrec == true ) ComputePreconditioner();
}

// ================================================ ====== ==== ==== == =

MultiLevelPreconditioner::MultiLevelPreconditioner( const Epetra_RowMatrix & EdgeMatrix,
						    const Epetra_RowMatrix & TMatrix,
						    const Epetra_RowMatrix & NodeMatrix,
						    const ParameterList & List,
						    const bool ComputePrec,
						    const char Prefix[] ) :
  RowMatrix_(&EdgeMatrix),
  RowMatrixAllocated_(0)
{
  sprintf(Prefix_,"%s",Prefix);

  List_ = List;

  Initialize();

  // set Maxwell here.
  // NOTE: RowMatrix_ and EdgeMatrix_ pointer to the same Epetra_RowMatrix
  SolvingMaxwell_ = true;
  NodeMatrix_ = & NodeMatrix;
  TMatrix_ = & TMatrix;
  EdgeMatrix_ = & EdgeMatrix;

  // construct hierarchy
  if( ComputePrec == true ) ComputePreconditioner();

}

// ================================================ ====== ==== ==== == =

MultiLevelPreconditioner::MultiLevelPreconditioner( ML_Operator * Operator,
						    const ParameterList & List, const bool ComputePrec,
						    const char Prefix[] )
{

  // need to wrap an Epetra_RowMatrix around Operator.
  // This is quite not the best approach for small matrices

  int MaxNumNonzeros;
  double CPUTime;
  
  ML_Operator2EpetraCrsMatrix(Operator,RowMatrixAllocated_,MaxNumNonzeros,
			      true,CPUTime);

  // this matrix must be freed by dtor. Keep trace of it in this pointer
  RowMatrix_ = RowMatrixAllocated_;

  // from now on as for the other constructors
  
  sprintf(Prefix_,"%s",Prefix);

  List_ = List;

  Initialize();

  // construct hierarchy
  if( ComputePrec == true ) ComputePreconditioner();
}

// ================================================ ====== ==== ==== == =

#ifdef HAVE_ML_TRIUTILS
MultiLevelPreconditioner::MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
						   Trilinos_Util_CommandLineParser & CLP,
						   const bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  RowMatrixAllocated_(0)
{
  Prefix_[0] = '\0';

  /* ********************************************************************** */
  /* Parse command line to get main options                                 */
  /* ********************************************************************** */

  string def = CLP.Get("-default","DD");
  SetDefaults(def,List_,"",SmootherOptions_,SmootherParams_);
  
  List_.set("max levels",CLP.Get("-num_levels",2));
  List_.set("increasing or decreasing",CLP.Get("-incr_or_decr","increasing"));

  List_.set("aggregation: type", CLP.Get("-aggr_scheme","Uncoupled"));
  List_.set("smoother: type", CLP.Get("-smoother_type","Gauss-Seidel"));
  List_.set("aggregation: nodes per aggregate", CLP.Get("-num_nodes_per_aggr",512));
  List_.set("aggregation: local aggregate", CLP.Get("-num_local_aggr",512));
  List_.set("smoother: pre or post", CLP.Get("-smoother_pre_or_post","both"));
  List_.set("smoother: damping factor", CLP.Get("-smoother_damping_factor",1.0));
  List_.set("coarse: type", CLP.Get("-coarse_type","Amesos_KLU"));
  List_.set("coarse: max processes", CLP.Get("-coarse_max_procs",4));
  List_.set("aggregation: damping factor", CLP.Get("-aggr_damping_factor",1.333));
  
  /*
    AZ_defaults(options,params);
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_ilu;
    List_.set("smoother: aztec options", options);
    List_.set("smoother: aztec params", params);
    List_.set("smoother: damping factor", 0.67);
  */

  List_.set("eigen-analysis: use symmetric algorithms", CLP.Has("-symmetricsssxx"));
  List_.set("eigen-analysis: tolerance", CLP.Get("-eigen_analysis_tol",1e-2));
  List_.set("compute null space", CLP.Has("-compute_null_space"));
  List_.set("null space dimension", CLP.Get("-null_space_dim",1));
  List_.set("add default null space", CLP.Has("-add_default_null_space"));
  List_.set("R and P smoothing: type", CLP.Get("-RP_smoothing","classic"));
  List_.set("R and P smoothing: damping", CLP.Get("-RP_damping","default"));

  /* ********************************************************************** */
  /* back to normal initialization                                          */
  /* ********************************************************************** */

  Initialize();

  // construct hierarchy
  if( ComputePrec == true ) ComputePreconditioner();
}
#endif

// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::Initialize()
{

  Comm_ = &(RowMatrix_->Comm());
  DomainMap_ = &(RowMatrix_->OperatorDomainMap());
  RangeMap_ = &(RowMatrix_->OperatorRangeMap());

  verbose_ = false;
  MaxLevels_ = 20;
  IsComputePreconditionerOK_ = false;

  NumPDEEqns_ = 1;
  
  NullSpaceToFree_ = 0;

  Label_ = 0;
  LevelID_ = new int[MaxLevels_];
  
  sprintf(ErrorMsg_,"ERROR (ML_Prec) : ");
  PrintMsg_ = "ML_Prec : ";
  
  AZ_defaults(SmootherOptions_,SmootherParams_);
  SmootherOptions_[AZ_precond] = AZ_dom_decomp;
  SmootherOptions_[AZ_subdomain_solve] = AZ_ilut;

  // Maxwell stuff is off by default
  SolvingMaxwell_ = false;
  NodeMatrix_ = 0;
  EdgeMatrix_ = 0;
  TMatrix_ = 0;
  ml_edges_ = 0;
  ml_nodes_ = 0;

  // timing
  NumApplications_ = 0;
  ApplicationTime_ = 0.0;
  FirstApplication_ = true;
  FirstApplicationTime_ = 0.0;
  NumConstructions_ = 0;
  ConstructionTime_ = 0.0;
  
  // some tracking here
  int NumInitializations = OutputList_.get("number of initialization phases", 0);
  OutputList_.set("number of initialization phases", ++NumInitializations);
}

// ================================================ ====== ==== ==== == =

int MultiLevelPreconditioner::ComputePreconditioner()
{

  Epetra_Time Time(Comm());
  
  {
    int NumCompute = OutputList_.get("number of construction phases", 0);
    OutputList_.set("number of construction phases", ++NumCompute);
  }
  
  char parameter[80];
  
  // get rid of what done before 
  
  if( IsComputePreconditionerOK_ == true ) {
    Destroy_ML_Preconditioner();
  }

  FirstApplication_ = true;

  if( Label_ ) delete Label_;
  
  Label_ = new char[80];
  
#ifdef HAVE_MPI
  const Epetra_MpiComm * MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());
  AZ_set_proc_config(ProcConfig_,MpiComm->Comm());
#else
  AZ_set_proc_config(ProcConfig_,AZ_NOT_MPI);
#endif

  sprintf(parameter,"%smax levels", Prefix_);
  NumLevels_ = List_.get(parameter,10);  

  sprintf(parameter,"%soutput", Prefix_);
  int OutputLevel = List_.get(parameter, 10);  
  ML_Set_PrintLevel(OutputLevel);

  verbose_ = ( 5 < ML_Get_PrintLevel() && ProcConfig_[AZ_node] == 0);

  if( verbose_ ) PrintLine();

  // user's defined output message
  sprintf(parameter,"%soutput prefix", Prefix_);
  PrintMsg_ = List_.get(parameter,PrintMsg_);

  // compute how to traverse levels (increasing of descreasing)
  // By default, use ML_INCREASING.
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  string IsIncreasing = List_.get(parameter,"increasing");

  int FinestLevel;
  
  if( IsIncreasing == "increasing" ) {
    FinestLevel = 0;
    for( int i=0 ; i<NumLevels_ ; ++i ) LevelID_[i] = FinestLevel+i;
  } else {
    FinestLevel = NumLevels_-1;
    for( int i=0 ; i<NumLevels_ ; ++i ) LevelID_[i] = FinestLevel-i;
  }
  
  // check no levels are negative
  for( int i=0 ; i<NumLevels_ ; ++i )
    if( LevelID_[i] <0 ) {
      cerr << ErrorMsg_ << "Level " << i << " has a negative ID" << endl;
      exit( EXIT_FAILURE );
    }  

  if( verbose_ ) {
    cout << PrintMsg_ << "Maximum number of levels = " << NumLevels_ << endl;
    if( IsIncreasing == "increasing" ) cout << PrintMsg_ << "Using increasing levels. ";
    else                               cout << PrintMsg_ << "Using decreasing levels. ";
    cout << "Finest level  = " << LevelID_[0];
    cout << ", coarsest level = " << LevelID_[NumLevels_-1] << endl;
  }


  int MaxCreationLevels = NumLevels_;
  if( IsIncreasing == "decreasing" )  MaxCreationLevels = FinestLevel+1;

  if( SolvingMaxwell_ == false ) {
    
    ML_Create(&ml_,MaxCreationLevels);

    int NumMyRows;
    
    NumMyRows = RowMatrix_->NumMyRows();
    int N_ghost = RowMatrix_->NumMyCols() - NumMyRows;
    
    if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix
    
    ML_Init_Amatrix(ml_,LevelID_[0],NumMyRows, NumMyRows, (void *) RowMatrix_);
    ML_Set_Amatrix_Getrow(ml_, LevelID_[0], Epetra_ML_getrow,
			  Epetra_ML_comm_wrapper, NumMyRows+N_ghost);
    
    ML_Set_Amatrix_Matvec(ml_, LevelID_[0], Epetra_ML_matvec);

  } else {

    ML_Create(&ml_edges_,MaxCreationLevels);
    ML_Create(&ml_nodes_,MaxCreationLevels);

    // DO SOMETHING HERE WITH MATRICES
    
  }
    
  ML_Aggregate_Create(&agg_);
  
  /* ********************************************************************** */
  /* set null space                                                         */
  /* ********************************************************************** */

  SetNullSpace();
  
  /* ********************************************************************** */
  /* pick up coarsening strategy. METIS and ParMETIS requires additional    */
  /* lines, as we have to set the number of aggregates                      */
  /* ********************************************************************** */

  SetAggregation();

  /* ********************************************************************** */
  /* minor settings                                                         */
  /* ********************************************************************** */
  
  int MaxCoarseSize = 50;
  sprintf(parameter,"%scoarse: max size", Prefix_);
  MaxCoarseSize = List_.get(parameter, MaxCoarseSize);
  ML_Aggregate_Set_MaxCoarseSize(agg_, MaxCoarseSize );

  double Threshold = 0.0;
  sprintf(parameter,"%saggregation: threshold", Prefix_);
  Threshold = List_.get(parameter, Threshold);
  ML_Aggregate_Set_Threshold(agg_,Threshold);
  
  int ReqAggrePerProc = 128;
  // compatibility with an older version
  sprintf(parameter,"%saggregation: req aggregates per process", Prefix_);
  if( List_.isParameter(parameter) ) 
    ReqAggrePerProc = List_.get(parameter, ReqAggrePerProc);
  else {
    sprintf(parameter,"%saggregation: next-level aggregates per process", Prefix_);
    ReqAggrePerProc = List_.get(parameter, ReqAggrePerProc);
  }
  ML_Aggregate_Set_ReqLocalCoarseSize( ml_, agg_, -1, ReqAggrePerProc);

  if( verbose_ ) {
    cout << PrintMsg_ << "Aggregation threshold = " << Threshold << endl;
    cout << PrintMsg_ << "Max coarse size = " << MaxCoarseSize << endl;
    cout << PrintMsg_ << "Requested next-level aggregates per process (ParMETIS) = " << ReqAggrePerProc << endl << endl;
  }

  /* ********************************************************************** */
  /* METIS and ParMETIS may suffer (that is, core dump) is the graph is     */
  /* highly non-symmetric. In this case, it is better to set this parameter */
  /* to ML_NO. It does affect METIS and ParMETIS only.                      */
  /* ********************************************************************** */

  bool UseDropping = true;
  sprintf(parameter,"%saggregation: use dropping", Prefix_);
  UseDropping = List_.get(parameter, UseDropping);
  if( UseDropping == true ) ML_Aggregate_Set_UseDropping( ML_YES );
  else                      ML_Aggregate_Set_UseDropping( ML_NO );  
 
  /* ********************************************************************** */
  /* Define scheme to determine damping parameter in prolongator and        */
  /* restriction smoother.                                                  */
  /* ********************************************************************** */

  SetSmoothingDamping();

  /************************************************************************/
  /* Build hierarchy using smoothed aggregation.                          */
  /* Then, retrive parameters for each level. Default values are given by */
  /* entries in parameter list without (level %d)                         */
  /*----------------------------------------------------------------------*/

  int Direction;
  if( IsIncreasing == "increasing" ) Direction = ML_INCREASING;
  else                               Direction = ML_DECREASING;

  if( SolvingMaxwell_ == false ) 
#ifdef MARZIO
    NumLevels_ = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml_, LevelID_[0], Direction, agg_);
#else
  NumLevels_ = ML_Gen_MGHierarchy_UsingAggregation(ml_, LevelID_[0], Direction, agg_);
#endif
  /* FIXME *

    put here creating for Maxwell
     NOTE: Tmatrix_ must be converted to an ML_Operator to be defined
     and TMatrixTranspose_ too.
  else
    // other parameters to be added to ParameterList ???
    NumLevels_ = ML_Gen_MGHierarchy_UsingReitzinger(ml_edges_, ml_nodes_,LevelID_[0],
						    Direction,agg_,TMatrix_,TMatrixTranspose_, 
						    &Tmat_array,&Tmat_trans_array, 
						    ML_YES, 1.5); 
  */
  if( verbose_ ) cout << PrintMsg_ << "Number of actual levels : " << NumLevels_ << endl;

  /* ********************************************************************** */
  /* Now cycling over all levels                                            */
  /* ********************************************************************** */

  if( SolvingMaxwell_ == false ) SetSmoothers();
  else                           SetSmoothersMaxwell();

  // MS ?? Isn't it already set somewhere ????
  sprintf(parameter,"%ssmoother: damping factor", Prefix_);
  double omega = List_.get(parameter,1.0);
  // what about parallel????  JJH

  /* ********************************************************************** */
  /* solution of the coarse problem                                         */
  /* ********************************************************************** */

  if( NumLevels_ > 1 ) SetCoarse();

  ownership_ = false;

  /* ********************************************************************** */
  /* Specific   preconditioners                                             */
  /* NOTE: the two-level DD preconditioners are kind of experimental!       */
  /* I suppose that the user knows what he/she is doing.... No real checks  */
  /* are performed. Note also that the coarse solver is somehow supposed to */
  /* be implemented as a post smoother (FIXME)                              */
  /* ********************************************************************** */

  ML_Gen_Solver(ml_, ML_MGV, LevelID_[0], LevelID_[NumLevels_-1]);

  CreateLabel();
  
  SetPreconditioner();
  if( verbose_ ) PrintLine();

  sprintf(parameter,"%sprint unused", Prefix_);
  if( List_.isParameter(parameter) ) {
    int ProcID = List_.get(parameter,-2);
    if( Comm().MyPID() == ProcID || ProcID == -1 ) PrintUnused();
  }

  /* ------------------- that's all folks --------------------------------- */

  IsComputePreconditionerOK_ = true;

  ConstructionTime_ += Time.ElapsedTime();
  
  return 0;
  
}

// ============================================================================

void MultiLevelPreconditioner::PrintUnused(int MyPID) 
{
  if( Comm().MyPID() == MyPID ) {
    PrintLine();
    cout << PrintMsg_ << "Unused parameters:" << endl;
    PrintUnused();
    PrintLine();
  }
}

// ============================================================================

void MultiLevelPreconditioner::PrintList(int MyPID) 
{
  if( Comm().MyPID() == MyPID ) {
    PrintLine();
    cout << List_;
    PrintLine();
  }
}

// ============================================================================

int MultiLevelPreconditioner::SetParameterList(const ParameterList & List) 
{
  if( IsComputePreconditionerOK_ == true ) DestroyPreconditioner();
  List_ = List;
  return 0;
  
}

// ============================================================================

int MultiLevelPreconditioner::CreateLabel()
{
  
  int i = ml_->ML_finest_level;
  char finest[80];
  char coarsest[80];
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
  sprintf(Label_,"%d level SA (%s, %s)", ml_->ML_num_actual_levels, finest, coarsest);
  
  return 0;
    
}

// ============================================================================

int MultiLevelPreconditioner::ApplyInverse(const Epetra_MultiVector& X,
					   Epetra_MultiVector& Y) const
{

  Epetra_Time Time(Comm());
  
  if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);
  if( !IsPreconditionerComputed() ) EPETRA_CHK_ERR(-10);

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

  MultiLevelPreconditioner * This = const_cast<MultiLevelPreconditioner *>(this);
  
  double t = Time.ElapsedTime();
  if( FirstApplication_ ) {
    This->FirstApplication_ = false;
    This->FirstApplicationTime_ += t;
  }
  
  This->ApplicationTime_ += t;
  
  ++(This->NumApplications_);

  return 0;
}

// ============================================================================

void MultiLevelPreconditioner::SetSmoothers() 
{

  char parameter[80];

  sprintf(parameter,"%ssmoother: sweeps", Prefix_);
  int num_smoother_steps = List_.get(parameter, 1);

  sprintf(parameter,"%ssmoother: damping factor", Prefix_);
  double omega = List_.get(parameter,1.0);

  sprintf(parameter,"%ssmoother: pre or post", Prefix_);
  int pre_or_post;
  string PreOrPostSmoother = List_.get(parameter,"post");

  sprintf(parameter,"%ssmoother: type", Prefix_);
  string Smoother = List_.get(parameter,"MLS");

  sprintf(parameter,"%ssmoother: aztec options", Prefix_);
  int * SmootherOptionsPtr = NULL;
  SmootherOptionsPtr = List_.get(parameter,SmootherOptionsPtr);

  sprintf(parameter,"%ssmoother: aztec params", Prefix_);
  double * SmootherParamsPtr = NULL;
  SmootherParamsPtr = List_.get(parameter,SmootherParamsPtr);

  sprintf(parameter,"%ssmoother: aztec as solver", Prefix_);
  bool AztecSmootherAsASolver = List_.get(parameter,false);
  int aztec_its;

  sprintf(parameter,"%ssmoother: MLS polynomial order", Prefix_);
  int MLSPolynomialOrder = List_.get(parameter,3);

  sprintf(parameter,"%ssmoother: MLS alpha",Prefix_);
  double MLSalpha = List_.get(parameter,30.0);;
  
  int SmootherLevels = (NumLevels_>1)?(NumLevels_-1):1;

  for( int level=0 ; level<SmootherLevels ; ++level ) {

    sprintf(parameter,"%ssmoother: sweeps (level %d)", Prefix_, LevelID_[level] );
    num_smoother_steps = List_.get(parameter,num_smoother_steps);

    sprintf(parameter,"%ssmoother: damping factor (level %d)", Prefix_, LevelID_[level] );
    omega = List_.get(parameter,omega);

    sprintf(parameter,"%ssmoother: pre or post (level %d)", Prefix_, LevelID_[level] );
    PreOrPostSmoother = List_.get(parameter, PreOrPostSmoother);
    
    if( PreOrPostSmoother == "post" ) pre_or_post = ML_POSTSMOOTHER;
    else if( PreOrPostSmoother == "pre" ) pre_or_post = ML_PRESMOOTHER;
    else if( PreOrPostSmoother == "both" ) pre_or_post = ML_BOTH;
    else 
      cerr << ErrorMsg_ << "smoother not recognized (" << PreOrPostSmoother << ")\n";
    
    sprintf(parameter,"%ssmoother: type (level %d)", Prefix_, LevelID_[level]);
    Smoother = List_.get(parameter,Smoother);

    if( Smoother == "Jacobi" ) {
      if( verbose_ ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : Jacobi (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_Jacobi(ml_, LevelID_[level], pre_or_post,
			     num_smoother_steps, omega);
    } else if( Smoother == "Gauss-Seidel" ) {
      if( verbose_ ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[level], pre_or_post,
				  num_smoother_steps, omega);
    } else if( Smoother == "symmetric Gauss-Seidel" ) {
      if( verbose_ ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : symmetric Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[level], pre_or_post,
				     num_smoother_steps, omega);
    } else if( Smoother == "block Gauss-Seidel" ) {
      if( verbose_ ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : block Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml_, LevelID_[level], pre_or_post,
				       num_smoother_steps, omega, NumPDEEqns_);
    } else if( Smoother == "MLS" ) {
      sprintf(parameter,"smoother: MLS polynomial order (level %d)", LevelID_[level]);
      if( verbose_ ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : MLS,"
			 << PreOrPostSmoother << endl;
      
      ML_Gen_Smoother_MLS(ml_, LevelID_[level], pre_or_post, 30.,
			  num_smoother_steps);
    } else if( Smoother == "aztec" ) {
      
      sprintf(parameter,"%ssmoother: aztec options (level %d)", Prefix_, LevelID_[level]);
      SmootherOptionsPtr = List_.get(parameter, SmootherOptionsPtr);
      sprintf(parameter,"%ssmoother: aztec params (level %d)", Prefix_, LevelID_[level]);
      SmootherParamsPtr = List_.get(parameter, SmootherParamsPtr);
      sprintf(parameter,"%ssmoother: aztec as solver (level %d)", Prefix_, LevelID_[level]);
      AztecSmootherAsASolver = List_.get(parameter,AztecSmootherAsASolver);
     
      if( AztecSmootherAsASolver == false ) aztec_its = AZ_ONLY_PRECONDITIONER;
      else                                  aztec_its = num_smoother_steps;
      
      if( SmootherOptionsPtr == NULL || SmootherParamsPtr == NULL ) {
      	SmootherOptionsPtr = SmootherOptions_;
      	SmootherParamsPtr = SmootherParams_;
      }

      if( verbose_ ) {
	cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : aztec";
	if( SmootherOptionsPtr[AZ_precond] == AZ_dom_decomp ) {
	  cout << " DD, overlap=" << SmootherOptionsPtr[AZ_overlap] << ", ";
	  if( SmootherOptionsPtr[AZ_reorder] == 1 ) cout << "reord, ";
	  else cout << "no reord, ";
	  switch( SmootherOptionsPtr[AZ_subdomain_solve] ) {
	  case AZ_lu: cout << " with LU"; break;
	  case AZ_ilu:
	    cout << "ILU(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_ilut:
	    cout << "ILUT(fill=" << SmootherParamsPtr[AZ_ilut_fill] << ",drop="
		 << SmootherParamsPtr[AZ_drop] << ")";
	    break;
	  case AZ_icc:
	    cout << "ICC(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_bilu:
	    cout << "BILU(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_rilu:
	    cout << "RILU(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ",omega="
		 << SmootherParamsPtr[AZ_omega] << ")";
	    break;
	  }
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_Jacobi ) {
	  cout << PrintMsg_ << " Jacobi preconditioner";
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_Neumann ) {
	  cout << PrintMsg_ << " Neumann preconditioner, order = " << SmootherOptionsPtr[AZ_poly_ord];
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_ls ) {
	  cout << PrintMsg_ << " LS preconditioner, order = " << SmootherOptionsPtr[AZ_poly_ord];
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_sym_GS ) {
	  cout << PrintMsg_ << " symmetric Gauss-Seidel preconditioner, sweeps = " << SmootherOptionsPtr[AZ_poly_ord];
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_none ) {
	  cout << PrintMsg_ << " with no preconditioning";
	}
	cout << ", "  << PreOrPostSmoother << endl;
      }
      
      ML_Gen_SmootherAztec(ml_, LevelID_[level], SmootherOptionsPtr, SmootherParamsPtr,
			   ProcConfig_, SmootherStatus_,
			   aztec_its, pre_or_post, NULL);
      
    } else if( Smoother == "ifpack" ) {
      if( verbose_ ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : Ifpack" << ","
			 << PreOrPostSmoother << endl;
      // get ifpack options from list ??? pass list ???
      ML_Gen_Smoother_Ifpack(ml_, LevelID_[level], pre_or_post, NULL, NULL);
    } else {
      if( ProcConfig_[AZ_node] == 0 )
	cerr << ErrorMsg_ << "Smoother not recognized!" << endl
	     << ErrorMsg_ << "(file " << __FILE__ << ",line " << __LINE__ << ")" << endl
	     << ErrorMsg_ << "Should be: " << endl
	     << ErrorMsg_ << "<Jacobi>/<Gauss-Seidel>/<Block Gauss-Seidel> / <MLS>" << endl
	     << ErrorMsg_ << "<aztec> " << endl;
      exit( EXIT_FAILURE );
    }
    
  } /* for */

  return;
}

// ============================================================================

void MultiLevelPreconditioner::SetSmoothersMaxwell()
{
  return;
}

// ============================================================================

void MultiLevelPreconditioner::SetCoarse() 
{

  char parameter[80];

  sprintf(parameter,"%scoarse: type", Prefix_);
  string CoarseSolution = List_.get(parameter, "Amesos_KLU");
  sprintf(parameter,"%scoarse: sweeps", Prefix_);
  int NumSmootherSteps = List_.get(parameter, 1);
  sprintf(parameter,"%scoarse: damping factor", Prefix_);
  double Omega = List_.get(parameter, 0.67);
    
  sprintf(parameter,"%scoarse: max processes", Prefix_);
  int MaxProcs = List_.get(parameter, -1);
  
  if( CoarseSolution == "Jacobi" ) 
    ML_Gen_Smoother_Jacobi(ml_, LevelID_[NumLevels_-1], ML_POSTSMOOTHER,
			   NumSmootherSteps, Omega);
  else if( CoarseSolution == "Gauss-Seidel" ) 
    ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[NumLevels_-1], ML_POSTSMOOTHER,
				NumSmootherSteps, Omega);
  else if( CoarseSolution == "SuperLU" ) 
    ML_Gen_CoarseSolverSuperLU( ml_, LevelID_[NumLevels_-1]);
  else if( CoarseSolution == "Amesos_KLU" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1], ML_AMESOS_KLU, MaxProcs);
  else if( CoarseSolution == "Amesos_UMFPACK" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1], ML_AMESOS_UMFPACK, MaxProcs);
  else if(  CoarseSolution == "Amesos_Superludist" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1], ML_AMESOS_SUPERLUDIST, MaxProcs);
  else if( CoarseSolution == "Amesos_MUMPS" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1], ML_AMESOS_MUMPS, MaxProcs);
  else
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1], ML_AMESOS_KLU, MaxProcs);
    
}

// ============================================================================

void MultiLevelPreconditioner::SetAggregation() 
{

  char parameter[80];
  
  int value = -777;
  sprintf(parameter,"%saggregation: type",Prefix_);
  string CoarsenScheme = List_.get(parameter,"Uncoupled");

  if ( CoarsenScheme == "Uncoupled-MIS" )
      ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_);
  else {
     for( int level=0 ; level<NumLevels_-1 ; ++level ) {  
   
       sprintf(parameter,"%saggregation: type (level %d)",Prefix_,LevelID_[level]);
       CoarsenScheme = List_.get(parameter,CoarsenScheme);
   
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
         if( Comm().MyPID() == 0 ) {
       cout << ErrorMsg_ << "specified options ("
            << CoarsenScheme << ") not valid. Should be:" << endl;
       cout << ErrorMsg_ << "<METIS> <ParMETIS> <MIS> <Uncoupled> <Coupled> <Hybrid>" << endl;
         }
         ML_Aggregate_Set_CoarsenSchemeLevel_METIS(LevelID_[level],NumLevels_,agg_);
       } 
   
       if( CoarsenScheme == "METIS" || CoarsenScheme == "ParMETIS" ) {
         
         bool isSet = false;
   
         // first look for parameters without any level specification
         
         sprintf(parameter,"%saggregation: global aggregates", Prefix_);
         if( List_.isParameter(parameter) ){
       value = -777; // simply means not set
       value = List_.get(parameter,value);
       if( value != -777 ) {
         ML_Aggregate_Set_GlobalNumber(ml_,agg_,LevelID_[level],value );
         isSet = true;
       }
         }
         
         sprintf(parameter,"%saggregation: local aggregates", Prefix_);
         if( List_.isParameter(parameter) ){
       value = -777;
       value = List_.get(parameter,value);
       if( value != -777 ) {
         ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value );
         isSet = true;
       }
         }
         
         sprintf(parameter,"%saggregation: nodes per aggregate", Prefix_);
         if( List_.isParameter(parameter) ){
       value = -777;
       value = List_.get(parameter,value);
       if( value != -777 ) {
         ML_Aggregate_Set_NodesPerAggr(ml_,agg_,LevelID_[level],value );
         isSet = true;
       }
         }
   
         // now for level-specific data
   
         sprintf(parameter,"%saggregation: global aggregates (level %d)", Prefix_, LevelID_[level]);
         if( List_.isParameter(parameter) ){
       value = -777; // simply means not set
       value = List_.get(parameter,value);
       if( value != -777 ) {
         ML_Aggregate_Set_GlobalNumber(ml_,agg_,LevelID_[level],value );
         isSet = true;
       }
         }
         
         sprintf(parameter,"%saggregation: local aggregates (level %d)", Prefix_, LevelID_[level]);
         if( List_.isParameter(parameter) ){
       value = -777;
       value = List_.get(parameter,value);
       if( value != -777 ) {
         ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value );
         isSet = true;
       }
         }
         
         sprintf(parameter,"%saggregation: nodes per aggregate (level %d)", Prefix_, LevelID_[level]);
         if( List_.isParameter(parameter) ){
       value = -777;
       value = List_.get(parameter,value);
       if( value != -777 ) {
         ML_Aggregate_Set_NodesPerAggr(ml_,agg_,LevelID_[level],value );
         isSet = true;
       }
         }
         
         if( isSet == false ) {
       // put default values
       sprintf(parameter,"%saggregation: local aggregates (level %d)", Prefix_, LevelID_[level]);
       value = List_.get(parameter,1);
       ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value);
         }
         
       } // if( CoarsenScheme == "METIS" || CoarsenScheme == "ParMETIS" )
       
     } /* for */
     } /* else */
  
}

// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::SetPreconditioner() 
{

  char parameter[80];
  
  sprintf(parameter,"%sprec type", Prefix_);
  string str = List_.get(parameter,"MGV");

  if( str == "1-level postsmoothing only" ) {
    
    sprintf(Label_, "1-level postsmoothing only");
    ml_->ML_scheme = ML_ONE_LEVEL_DD;

  } else if( str == "two-level additive DD" ) {

    sprintf(Label_, "two-level additive DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_ADD;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level additive DD' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level hybrid DD") {

    sprintf(Label_, "two-level hybrid DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level hybrid DD' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level hybrid DD (2)") {

    sprintf(Label_, "two-level hybrid DD (2)");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID_2;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level hybrid DD (2)' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "full MGV" ) {
    ml_->ML_scheme = ML_MGFULLV;

  } else if( str == "MGV" ) {
    // it is the default
    ml_->ML_scheme = ML_MGV;
  } else {

    if( Comm().MyPID() == 0 ) {
      cerr << ErrorMsg_ << "`prec type' has an incorrect value. It should be" << endl
	   << ErrorMsg_ << "<1-level postsmoothing only> / <two-level additive DD>" << endl
	   << ErrorMsg_ << "<two-level hybrid DD> / <two-level hybrid DD (2)>" << endl;
    }
    exit( EXIT_FAILURE );
    
  }
  
}

// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::SetNullSpace() 
{

  char parameter[80];

  const Epetra_VbrMatrix * VbrMatrix = dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
  if( VbrMatrix == 0 ) {
    sprintf(parameter,"%sPDE equations", Prefix_);
    NumPDEEqns_ = List_.get(parameter, 1);
  }
  else {
    int NumBlockRows = VbrMatrix->RowMap().NumGlobalElements();
    int NumRows = VbrMatrix->RowMap().NumGlobalPoints();
    if( NumRows % NumBlockRows ) {
      cerr << "# rows must be a multiple of # block rows ("
	   << NumRows << "," << NumBlockRows << ")" << endl;
      exit( EXIT_FAILURE );
    }
    NumPDEEqns_ = NumRows/NumBlockRows;
  }

  int NullSpaceDim = NumPDEEqns_;
  double * NullSpacePtr = NULL;
  
  sprintf(parameter,"%snull space dimension", Prefix_);
  NullSpaceDim = List_.get(parameter, NumPDEEqns_);
  sprintf(parameter,"%snull space vectors", Prefix_);
  NullSpacePtr = List_.get(parameter, NullSpacePtr);

  sprintf(parameter,"%scompute null space", Prefix_);
  bool ComputeNullSpace = List_.get(parameter, false);
  
  if( ComputeNullSpace == false ) {

    // sanity check for default null-space
    if( NullSpacePtr == NULL ) NullSpaceDim = NumPDEEqns_;
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NullSpaceDim,NullSpacePtr,
			       RowMatrix_->NumMyRows());
    
  } else {

#ifdef HAVE_ML_ANASAZI

    Epetra_Time Time(Comm());
    
    if( NullSpacePtr != NULL && Comm().MyPID() == 0 ) {
      cerr << ErrorMsg_ << "Null space vectors is not NULL!" << endl
	   << ErrorMsg_ << "Now reallocating memory and computing null space " << endl
	   << ErrorMsg_ << "using eigenvectors estimates, for " << NullSpaceDim << " vector(s)..." << endl;
    }

    sprintf(parameter,"%sadd default null space", Prefix_);
    bool UseDefaultVectors = List_.get(parameter, false);

    if( verbose_ ) {
      cout << PrintMsg_ << "Computing " << NullSpaceDim << " null space vector(s)" << endl;
      if( parameter ) cout << PrintMsg_ << "(plus " << NumPDEEqns_ << " constant vector(s))" << endl;
    }
    
    int TotalNullSpaceDim = NullSpaceDim;
    if( UseDefaultVectors ) TotalNullSpaceDim += NumPDEEqns_;

    // create a double vector hosting null space
    int LDA = NumMyRows();
    
    NullSpacePtr = new double[TotalNullSpaceDim*LDA];

    // and fill it with normal 0's and 1's
    if( UseDefaultVectors )
      for( int i=0 ; i<NumPDEEqns_ ; ++i )
	for( int j=0 ; j<LDA ; ++j )
	  if( j%NumPDEEqns_ == i ) NullSpacePtr[j+i*LDA] = 1.0;
	  else                     NullSpacePtr[j+i*LDA] = 0.0;

    double * start = NullSpacePtr;
    if( UseDefaultVectors ) start += NumPDEEqns_*LDA;
    
    Epetra_MultiVector EigenVectors(View,OperatorDomainMap(),start,LDA, NullSpaceDim);
    
    EigenVectors.Random();
    
    double * RealEigenvalues = new double[NullSpaceDim];
    double * ImagEigenvalues = new double[NullSpaceDim];
    
    // create List for Anasazi (kept separate from List_, I don't want to pollute it)
    // Also, I keep it local (not use EigenList_
    ParameterList AnasaziList(EigenList_);
    
    // new parameters specific for this function only (not set by the user)
    AnasaziList.set("matrix operation", "A");    
    AnasaziList.set("action", "SM");

    ML_Anasazi::Interface(RowMatrix_,EigenVectors,RealEigenvalues,
			  ImagEigenvalues, AnasaziList);
    
    NullSpaceToFree_ = NullSpacePtr; // this null space will be freed later
    
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,TotalNullSpaceDim,
			       NullSpacePtr,
			       NumMyRows());
   
    delete [] RealEigenvalues;
    delete [] ImagEigenvalues;
    
    if( verbose_ ) cout << PrintMsg_ << "Total time for eigen-analysis = " << Time.ElapsedTime() << " (s)\n";
    
#else
     cout << "ML_Anasazi ERROR: you must compile with --with-ml_anasazi "  << endl
       << "ML_Anasazi ERROR: for eigen-analysis." << endl;
     exit( EXIT_FAILURE );
#endif

  }
}

// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::SetEigenList() 
{

  char parameter[80];
  
  // eigen-analysis:
  sprintf(parameter,"%seigen-analysis: use symmetric algorithms", Prefix_);
  bool IsSymmetric = List_.get(parameter,false);
    
  if( IsSymmetric ) EigenList_.set("eigen-analysis: symmetric problem",true);
  else              EigenList_.set("eigen-analysis: symmetric problem",false);

  sprintf(parameter,"%seigen-analysis: tolerance", Prefix_);
  EigenList_.set("eigen-analysis: tolerance", List_.get(parameter, 1e-2));

  sprintf(parameter,"%seigen-analysis: use diagonal scaling", Prefix_);    
  EigenList_.set("eigen-analysis: use diagonal scaling", List_.get(parameter,false));
    
  sprintf(parameter,"%seigen-analysis: restart", Prefix_);
  int itemp = List_.get(parameter, 100);
  EigenList_.set("eigen-analysis: restart", itemp);

  sprintf(parameter,"%seigen-analysis: length", Prefix_);
  itemp =  List_.get(parameter, 20);
  EigenList_.set("eigen-analysis: length", itemp);

  sprintf(parameter,"%seigen-analysis: normalize eigenvectors", Prefix_);
  bool btemp =  List_.get(parameter, false);
  EigenList_.set("eigen-analysis: normalize eigenvectors",btemp);

  // field of values:

  sprintf(parameter,"%sfield-of-values: tolerance", Prefix_);
  EigenList_.set("field-of-values: tolerance", List_.get(parameter, 1e-2));

  sprintf(parameter,"%sfield-of-values: use diagonal scaling", Prefix_);    
  EigenList_.set("field-of-values: use diagonal scaling", List_.get(parameter,false));
    
  sprintf(parameter,"%sfield-of-values: restart", Prefix_);
  itemp = List_.get(parameter, 100);
  EigenList_.set("field-of-values: restart", itemp);

  sprintf(parameter,"%sfield-of-values: length", Prefix_);
  itemp =  List_.get(parameter, 20);
  EigenList_.set("field-of-values: ", itemp);

  sprintf(parameter,"%sfield-of-values: print current status", Prefix_);
  btemp =  List_.get(parameter, false);
  EigenList_.set("field-of-values: print current status", btemp);

  // general output
  
  sprintf(parameter,"%soutput", Prefix_);
  itemp =  List_.get(parameter, 10);
  EigenList_.set("output",itemp);
    
}

// ================================================ ====== ==== ==== == =

void MultiLevelPreconditioner::SetSmoothingDamping() 
{

#ifdef MARZIO

  char parameter[80];
  Epetra_Time Time(Comm());
  
  // almost everything here is experimental ;)
  
  sprintf(parameter,"%sR and P smoothing: type", Prefix_);
  string RandPSmoothing = List_.get(parameter, "symmetric");

  double DampingFactor = 1.333;
  if( SolvingMaxwell_ ) DampingFactor = 0.0;

  /* ********************************************************************** */
  /* For "classical" (non-Anasazi) approach to determine lambda_max only.   */
  /* ********************************************************************** */

  sprintf(parameter,"%seigen-analysis: use symmetric algorithms", Prefix_);
  bool IsSymmetric = List_.get(parameter,false);
  
  if( IsSymmetric ) ML_Aggregate_Set_SpectralNormScheme_Calc(agg_);
  else              ML_Aggregate_Set_SpectralNormScheme_Anorm(agg_);
  
  /* start looping over different options */
  
  if( RandPSmoothing == "classic" ) {

    /* ********************************************************************** */
    /* This is the standard approach followed by ML (R = P^T)                 */
    /* ********************************************************************** */

    sprintf(parameter,"%saggregation: damping factor", Prefix_);
    DampingFactor = List_.get(parameter, DampingFactor);
    ML_Aggregate_Set_DampingFactor( agg_, DampingFactor );
    
    agg_->Restriction_smoothagg_transpose = ML_FALSE;
    if( verbose_ )
      cout << PrintMsg_ << "R and P smoothing : Standard ML procedure" << endl;
    
  } else if( RandPSmoothing == "standard-use A in R" ) {

    /* ********************************************************************** */
    /* This is the ML way, but with A in R                                    */
    /* ********************************************************************** */

    sprintf(parameter,"%saggregation: damping factor", Prefix_);
    DampingFactor = List_.get(parameter, DampingFactor);
    ML_Aggregate_Set_DampingFactor( agg_, DampingFactor );

    if( verbose_ )
      cout << PrintMsg_ << "R and P smoothing : Standard ML procedure with A to smooth restriction"
	   << endl;
    agg_->Restriction_smoothagg_transpose = ML_TRUE;

  } else if( RandPSmoothing == "advanced" ) {

    /* ********************************************************************** */
    /* This is the new way, based on Anasazi to compute eigen-widgets         */
    /* ********************************************************************** */

    int MaxNumNonzeros;
    double CPUTime;
    
    if( verbose_ )
      cout << PrintMsg_ << "Use A to smooth restriction operator" << endl;
    agg_->Restriction_smoothagg_transpose = ML_TRUE;
    
    SetEigenList();
    
    struct ML_Field_Of_Values * field_of_values;

    // stick default values (undefined)
    field_of_values = (struct ML_Field_Of_Values *) malloc( sizeof(struct ML_Field_Of_Values) );
    field_of_values->eta     = 0.0;
    field_of_values->real_max= -1.0;
    field_of_values->imag_max= -1.0;
    field_of_values->poly_order = 0;
    // following values for choice:
    // -1 : undefined
    //  0 : do nothing, put eta = 0 (used in non-smoothed aggregation)
    //  1 : compute the box of the field of values (eta = imag_max/real_max)
    //  2 : compute ||lambda_max|| (eta = sqrt(imag_max^2 + real_max^2))
    
    field_of_values->choice     = -1;
    // and this is a pointer for the object's interal ParameterList
    // That's stilistically hugly, but I need this because ML is mainly C-coded
    field_of_values->EigenList = (void *) &EigenList_;

    // still to set up polynomial coeffiecients
    string DampingType =  List_.get("R and P smoothing: damping", "default");
  
    if( DampingType == "non_smoothed" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : non-smoothed aggregation" << endl;

      field_of_values->choice     =  0;
      field_of_values->poly_order =  0;
      // I don't really need them, smoothing will be set to zero
      field_of_values->R_coeff[0] =  0.0;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;
      
      field_of_values->P_coeff[0] =  0.0;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;
      
    } else if( DampingType == "almost_non_smoothed" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : almost non-smoothed aggregation" << endl;

      field_of_values->choice     =  0;
      field_of_values->poly_order =  0;
      // I don't really need them, smoothing will be set to zero
      field_of_values->R_coeff[0] =  0.000000001;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;
      
      field_of_values->P_coeff[0] =  0.000000001;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;
      
    } else if( DampingType == "default" || DampingType == "ray" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using default values" << endl;
      
      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
      
      field_of_values->R_coeff[0] =  1.107;
      field_of_values->R_coeff[1] =  0.285;
      field_of_values->R_coeff[2] =  0.718;
      
      field_of_values->P_coeff[0] =  1.878;
      field_of_values->P_coeff[1] = -2.515;
      field_of_values->P_coeff[2] =  0.942;
      
    } else if( DampingType == "marzio" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using marzio values" << endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
      
      field_of_values->R_coeff[0] =  1.138;
      field_of_values->R_coeff[1] =  1.162;
      field_of_values->R_coeff[2] = -2.384;
      
      field_of_values->P_coeff[0] =  2.143;
      field_of_values->P_coeff[1] = -2.179;
      field_of_values->P_coeff[2] =  0.101;

    } else if( DampingType == "random" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using random values" << endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
	    
      // initialize seed
      unsigned int s = (int) nearbyint(Time.ElapsedTime()*10000);
      srandom(s);
      
      // put random values
      field_of_values->R_coeff[0] =  (double)random()/RAND_MAX;
      field_of_values->R_coeff[1] =  (double)random()/RAND_MAX;
      field_of_values->R_coeff[2] =  (double)random()/RAND_MAX;
      
      field_of_values->P_coeff[0] =  (double)random()/RAND_MAX;
      field_of_values->P_coeff[1] =  (double)random()/RAND_MAX;
      field_of_values->P_coeff[2] =  (double)random()/RAND_MAX;

      if( verbose_ ) {
	cout << PrintMsg_ << "Random coefficients for R and P:" << endl
	     << PrintMsg_ << field_of_values->R_coeff[0] << "   "
	     << field_of_values->R_coeff[1] << "   "
	     << field_of_values->R_coeff[2] << endl
	     << PrintMsg_ << field_of_values->P_coeff[0] << "   "
	     << field_of_values->P_coeff[1] << "   "
	     << field_of_values->P_coeff[2] << "   (seed = "
	     << s << ")" << endl;
      }
      
    } else if( DampingType == "classic" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using 4/3" << endl;

      field_of_values->choice     =  2;
      field_of_values->poly_order =  2;
      
      field_of_values->R_coeff[0] =  1.333;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;
      
      field_of_values->P_coeff[0] =  1.333;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;

    } else {

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
      
      // get them from parameters' list.
      field_of_values->R_coeff[0] =  List_.get("R and P smoothing: c_0",  1.107);
      field_of_values->R_coeff[1] =  List_.get("R and P smoothing: c_1",  0.285);
      field_of_values->R_coeff[2] =  List_.get("R and P smoothing: c_2",  0.718);

      field_of_values->P_coeff[0] =  List_.get("R and P smoothing: g_0",  1.878);
      field_of_values->P_coeff[1] =  List_.get("R and P smoothing: g_1", -2.515);
      field_of_values->P_coeff[2] =  List_.get("R and P smoothing: g_2",  0.942);
      
    } 

    agg_->field_of_values = (void*) field_of_values;  

  } else {

    cerr << endl;
    cerr << ErrorMsg_ << "Parameter for `R and P smoothing : value for" << endl
	 << ErrorMsg_ << "`Standard ML procedure' not recognized (" << RandPSmoothing << ")" << endl
	 << ErrorMsg_ << "NO ACTION PERFORMED !!" << endl << endl;
    
  }

#else

  char parameter[80];
  
  sprintf(parameter,"%seigen-analysis: use symmetric algorithms", Prefix_);
  bool IsSymmetric = List_.get(parameter,false);
  
  if( IsSymmetric ) ML_Aggregate_Set_SpectralNormScheme_Calc(agg_);
  else              ML_Aggregate_Set_SpectralNormScheme_Anorm(agg_);

#endif       

}

// ============================================================================

int ML_Epetra::SetDefaults(string ProblemType, ParameterList & List, char * Prefix_,
			   int SmootherOptions[], double SmootherParams[] )
{
  
  int rv = 0;
  
  if( ProblemType == "SA" )
    return( ML_Epetra::SetDefaultsSA(List, Prefix_, SmootherOptions, SmootherParams ) );
  else if( ProblemType == "maxwell" )
    return( ML_Epetra::SetDefaultsMaxwell(List, Prefix_, SmootherOptions, SmootherParams ) );
  else if( ProblemType == "DD-ML" )
    return( ML_Epetra::SetDefaultsDD_3Levels(List, Prefix_, SmootherOptions, SmootherParams ) );
  else if( ProblemType == "DD" )
    return( ML_Epetra::SetDefaultsDD(List, Prefix_, SmootherOptions, SmootherParams ) );
  else {
    cerr << "ERROR: Wrong input parameter in `SetDefaults'. Should be: " << endl
	 << "ERROR: <SA> / <DD> / <DD-ML> / <maxwell>" << endl;
    rv = 1;
  }

  EPETRA_CHK_ERR(rv);

  return rv;
  
  
}

// ============================================================================

int ML_Epetra::SetDefaultsDD(ParameterList & List, char * Prefix,
			     int SmootherOptions[], double SmootherParams[]) 
{

  char parameter[80];

  sprintf(parameter,"%smax levels", Prefix);
  List.set(parameter,2);

  sprintf(parameter,"%soutput", Prefix);
  List.set(parameter,10);
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix);
  List.set(parameter,"increasing");

  sprintf(parameter,"%sPDE equations", Prefix);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type",Prefix);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: local aggregates",Prefix);
  List.set(parameter,1);
  
  sprintf(parameter,"%saggregation: damping factor",Prefix);
  List.set(parameter,0.01);

  sprintf(parameter,"%scoarse: max size",Prefix);
  List.set(parameter,128);

  sprintf(parameter,"%saggregation: threshold",Prefix);
  List.set(parameter,0.0);
  
  sprintf(parameter,"%ssmoother: sweeps",Prefix);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor",Prefix);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: pre or post",Prefix);
  List.set(parameter,"both");

  if( SmootherOptions != 0 && SmootherParams != 0 ) {
    
    sprintf(parameter,"%ssmoother: type",Prefix);
    List.set(parameter,"aztec");

    AZ_defaults(SmootherOptions,SmootherParams);
    SmootherOptions[AZ_precond] = AZ_dom_decomp;
    SmootherOptions[AZ_scaling] = AZ_none;
    SmootherOptions[AZ_subdomain_solve] = AZ_ilut;
    
    sprintf(parameter,"%ssmoother: aztec options",Prefix);
    List.set(parameter,SmootherOptions);
    
    sprintf(parameter,"%ssmoother: aztec params",Prefix);
    List.set(parameter,SmootherParams);
    
    sprintf(parameter,"%ssmoother: aztec as solver",Prefix);
    List.set(parameter,false);

  } else {

    sprintf(parameter,"%ssmoother: type",Prefix);
    List.set(parameter,"Gauss-Seidel");
    
  }
  
  sprintf(parameter,"%scoarse: type",Prefix);
  List.set(parameter,"Amesos_KLU");

  sprintf(parameter,"%sprec type",Prefix);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix);
  List.set(parameter,0);

  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_3Levels(ParameterList & List, char * Prefix_,
				     int SmootherOptions[], double SmootherParams[]) 
{

  char parameter[80];

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,3);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,10);
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"increasing");

  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type (level 0)",Prefix_);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: type (level 1)",Prefix_);
  List.set(parameter,"ParMETIS");

  sprintf(parameter,"%saggregation: nodes per aggregate (level 0)",Prefix_);
  List.set(parameter,512);

  sprintf(parameter,"%saggregation: nodes per aggregate (level 1)",Prefix_);
  List.set(parameter,32);
  
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,0.01);

  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,128);

  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);
  
  sprintf(parameter,"%ssmoother: sweeps (level 0)",Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor (level 0)",Prefix_);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: pre or post (level 0)",Prefix_);
  List.set(parameter,"both");

  if( SmootherOptions != 0 && SmootherParams != 0 ) {
    
    sprintf(parameter,"%ssmoother: type (level 0)",Prefix_);
    List.set(parameter,"aztec");
    
    AZ_defaults(SmootherOptions,SmootherParams);
    SmootherOptions[AZ_precond] = AZ_dom_decomp;
    SmootherOptions[AZ_subdomain_solve] = AZ_ilut;
    
    sprintf(parameter,"%ssmoother: aztec options (level 0)",Prefix_);
    List.set(parameter,SmootherOptions);
    
    sprintf(parameter,"%ssmoother: aztec params (level 0)",Prefix_);
    List.set(parameter,SmootherParams);
    
    sprintf(parameter,"%ssmoother: aztec as solver (level 0)",Prefix_);
    List.set(parameter,false);
  
  } else {

    sprintf(parameter,"%ssmoother: type",Prefix_);
    List.set(parameter,"Gauss-Seidel");
    
  }
  
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos_KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,0);
  
  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsMaxwell(ParameterList & List, char * Prefix_,
				  int SmootherOptions[], double SmootherParams[]) 
{

  // FIXME : here default values for Maxwell
  
  char parameter[80];
  int MaxLevels = 10;
  
  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,MaxLevels);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,10);
  
  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"decreasing");

  // aggregation: Uncoupled for first levels, then MIS
  sprintf(parameter,"%saggregation: type",Prefix_);
  List.set(parameter,"Hybrid");

  // optimal value for smoothed aggregation
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,1.3333);

  // relative small coarse size
  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,16);

  // don't forget any element
  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);

  // gauss-seidel for all levels
  sprintf(parameter,"%ssmoother: sweeps (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: type (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,"MLS");

  sprintf(parameter,"%ssmoother: MLS polynomial order", Prefix_);
  List.set(parameter,3);
  
  sprintf(parameter,"%ssmoother: pre or post (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,"both");
  
  // simplest solver on coarse problem
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"SuperLU");
//  Tim Davis' simple serial LU package.  It's part of Amesos
//  itself.
//  List.set(parameter,"Amesos_KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  // print unused parameters on proc 0
  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,0);

  return 0;
  
}

// ============================================================================

int ML_Epetra::SetDefaultsSA(ParameterList & List, char * Prefix_,
			     int SmootherOptions[], double SmootherParams[]) 
{

  char parameter[80];
  int MaxLevels = 16;
  
  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,MaxLevels);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,16);
  
  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"increasing");

  // aggregation: Uncoupled for first levels, then MIS
  sprintf(parameter,"%saggregation: type (level 0)",Prefix_);
  List.set(parameter,"Uncoupled");

  sprintf(parameter,"%saggregation: type (level 1)",Prefix_);
  List.set(parameter,"Uncoupled");

  sprintf(parameter,"%saggregation: type (level 2)",Prefix_);
  List.set(parameter,"MIS");
  
  // optimal value for smoothed aggregation
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,1.3333);

  // relative small coarse size
  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,16);

  // don't forget any element
  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);

  // guass-seidel for all levels
  sprintf(parameter,"%ssmoother: sweeps",Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor",Prefix_);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: type",Prefix_);
  List.set(parameter,"Gauss-Seidel");
  
  sprintf(parameter,"%ssmoother: pre or post",Prefix_);
  List.set(parameter,"both");
  
  // simplest solver on coarse problem
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos_KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  // print unused parameters on proc 0
  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,0);
  
  return 0;

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
