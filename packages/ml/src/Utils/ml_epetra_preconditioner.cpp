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
#include "Epetra_RowMatrix.h"
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

#include "ml_anasazi.h"

// ================================================ ====== ==== ==== == =

void Epetra_ML_Preconditioner::PrintLine() 
{
  cout << "--------------------------------------------------------------------------------" << endl;
}

// ================================================ ====== ==== ==== == =

void Epetra_ML_Preconditioner::Destroy_ML_Preconditioner()
{
  
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

  if( LevelID_ != 0 ) delete [] LevelID_;
  
  if( ML_Get_PrintLevel() > 5 && Comm().MyPID() == 0 && ML_Time_->NumApplications() ) {

    double AvgTime = (ML_Time_->Application()-ML_Time_->FirstApplication())/ML_Time_->NumApplications();
    
    // stick data in OutputList

    int i;
    double t;
    
    t = OutputList_.get("time for applications", 0.0) + ML_Time_->Application();
    OutputList_.set("application time", t);
    
    t = OutputList_.get("time for startup", 0.0) + ML_Time_->FirstApplication()-AvgTime;
    OutputList_.set("startup time", t);

    i = OutputList_.get("number of applications", 0);
    OutputList_.set("number of applications", (i+ML_Time_->NumApplications()));

    // print on screen
    
    PrintLine();
    cout << PrintMsg_ << "Total application time = " << ML_Time_->Application()
	 << " (s)" << endl;
    cout << PrintMsg_ << "Estimated startup time = " << ML_Time_->FirstApplication() - AvgTime << " (s)" << endl;
    cout << PrintMsg_ << "Each of " << ML_Time_->NumApplications()
	 << " applications took  " << AvgTime << " (s)" << endl;
    PrintLine();
  }

  if( NullSpaceToFree_ != 0 ) delete NullSpaceToFree_;
  
  delete ML_Time_;
  
  IsComputePreconditionerOK_ = false;

}

// ================================================ ====== ==== ==== == =

Epetra_ML_Preconditioner::Epetra_ML_Preconditioner(const Epetra_RowMatrix & RowMatrix,
						   char ProblemType[],
						   bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  Comm_(RowMatrix.Comm()),
  DomainMap_(RowMatrix.OperatorDomainMap()),
  RangeMap_(RowMatrix.OperatorRangeMap())
{

  sprintf(Prefix_,"");
  
  ParameterList NewList;
  List_ = NewList;
  SetDefaults(List_,ProblemType);
  
  Initialize(ComputePrec);
}

// ================================================ ====== ==== ==== == =

Epetra_ML_Preconditioner::Epetra_ML_Preconditioner(const Epetra_RowMatrix & RowMatrix,
						   char ProblemType[],
						   bool ComputePrec, char Prefix[] ) :
  RowMatrix_(&RowMatrix),
  Comm_(RowMatrix.Comm()),
  DomainMap_(RowMatrix.OperatorDomainMap()),
  RangeMap_(RowMatrix.OperatorRangeMap())
{

  sprintf(Prefix_,"%s",Prefix);
  
  ParameterList NewList;
  List_ = NewList;
  SetDefaults(List_,ProblemType);
  
  Initialize(ComputePrec);
}

// ================================================ ====== ==== ==== == =

Epetra_ML_Preconditioner::Epetra_ML_Preconditioner(const Epetra_RowMatrix & RowMatrix,
						   bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  Comm_(RowMatrix.Comm()),
  DomainMap_(RowMatrix.OperatorDomainMap()),
  RangeMap_(RowMatrix.OperatorRangeMap())
{

  sprintf(Prefix_,"");
  
  ParameterList NewList;
  List_ = NewList;
  SetDefaults(List_,"DD");
    
  Initialize(ComputePrec);
}

// ================================================ ====== ==== ==== == =

Epetra_ML_Preconditioner::Epetra_ML_Preconditioner( const Epetra_RowMatrix & RowMatrix,
						    ParameterList & List, bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  Comm_(RowMatrix.Comm()),
  DomainMap_(RowMatrix.OperatorDomainMap()),
  RangeMap_(RowMatrix.OperatorRangeMap())
{
  sprintf(Prefix_,"");

  List_ = List;
  
  Initialize(ComputePrec);
}

// ================================================ ====== ==== ==== == =

Epetra_ML_Preconditioner::Epetra_ML_Preconditioner( const Epetra_RowMatrix & RowMatrix,
						    ParameterList & List, bool ComputePrec,
						    char Prefix[] ) :
  RowMatrix_(&RowMatrix),
  Comm_(RowMatrix.Comm()),
  DomainMap_(RowMatrix.OperatorDomainMap()),
  RangeMap_(RowMatrix.OperatorRangeMap())
{
  sprintf(Prefix_,"%s",Prefix);

  List_ = List;

  Initialize(ComputePrec);
}

// ================================================ ====== ==== ==== == =

void Epetra_ML_Preconditioner::Initialize(bool ComputePrec)
{
  MaxLevels_ = 20;
  IsComputePreconditionerOK_ = false;
  ML_Time_ = 0;

  NullSpaceToFree_ = 0;
  
  LevelID_ = new int[MaxLevels_];
  
  sprintf(ErrorMsg_,"ERROR (ML_Prec) : ");
  PrintMsg_ = "ML_Prec : ";
  if( ComputePrec == true ) ComputePreconditioner();

  AZ_defaults(SmootherOptions_,SmootherParams_);
  SmootherOptions_[AZ_precond] = AZ_dom_decomp;
  SmootherOptions_[AZ_subdomain_solve] = AZ_ilut;

  int NumInitializations = OutputList_.get("number of initialization phases", 0);
  OutputList_.set("number of initialization phases", ++NumInitializations);
}

// ================================================ ====== ==== ==== == =

int Epetra_ML_Preconditioner::ComputePreconditioner()
{
  
  {
    int NumCompute = OutputList_.get("number of construction phases", 0);
    OutputList_.set("number of construction phases", ++NumCompute);
  }
  
  char parameter[80];
  
  // get rid of what done before 
  
  if( IsComputePreconditionerOK_ == true ) {
    Destroy_ML_Preconditioner();
  }

  // reset timing
  
  if( ML_Time_ ) delete ML_Time_;
  ML_Time_ = new Epetra_ML_Time();

  Label_ = new char[80];
  
#ifdef HAVE_MPI
  const Epetra_MpiComm * MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm_);
  AZ_set_proc_config(ProcConfig_,MpiComm->Comm());
#else
  AZ_set_proc_config(ProcConfig_,AZ_NOT_MPI);
#endif

  sprintf(parameter,"%smax levels", Prefix_);
  NumLevels_ = List_.get(parameter,10);  

  sprintf(parameter,"%soutput", Prefix_);
  int OutputLevel = List_.get(parameter, 10);  
  ML_Set_PrintLevel(OutputLevel);

  bool verbose = ( 5 < ML_Get_PrintLevel() && ProcConfig_[AZ_node] == 0);

  if( verbose ) PrintLine();

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

  if( verbose ) {
    cout << PrintMsg_ << "Maximum number of levels = " << NumLevels_ << endl;
    if( IsIncreasing == "increasing" ) cout << PrintMsg_ << "Using increasing levels. ";
    else                               cout << PrintMsg_ << "Using decreasing levels. ";
    cout << "Finest level  = " << LevelID_[0];
    cout << ", coarsest level = " << LevelID_[NumLevels_-1] << endl;
  }


  int MaxCreationLevels = NumLevels_;
  if( IsIncreasing == "decreasing" )  MaxCreationLevels = FinestLevel+1;

  ML_Create(&ml_,MaxCreationLevels);

  int NumMyRows, osize;

  NumMyRows = RowMatrix_->NumMyRows();
  int N_ghost = RowMatrix_->NumMyCols() - NumMyRows;

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  ML_Init_Amatrix(ml_,LevelID_[0],NumMyRows, NumMyRows, (void *) RowMatrix_);
  ML_Set_Amatrix_Getrow(ml_, LevelID_[0], Epetra_ML_getrow,
			Epetra_ML_comm_wrapper, NumMyRows+N_ghost);

  ML_Set_Amatrix_Matvec(ml_, LevelID_[0], Epetra_ML_matvec);

  ML_Aggregate_Create(&agg_);
  
  /* ********************************************************************** */
  /* set null space                                                         */
  /* ********************************************************************** */

  int NumPDEEqns;
  
  const Epetra_VbrMatrix * VbrMatrix = dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
  if( VbrMatrix == 0 ) {
    sprintf(parameter,"%sPDE equations", Prefix_);
    NumPDEEqns = List_.get(parameter, 1);
  }
  else {
    int NumBlockRows = VbrMatrix->RowMap().NumGlobalElements();
    int NumRows = VbrMatrix->RowMap().NumGlobalPoints();
    if( NumRows % NumBlockRows ) {
      cerr << "# rows must be a multiple of # block rows ("
	   << NumRows << "," << NumBlockRows << ")" << endl;
      exit( EXIT_FAILURE );
    }
    NumPDEEqns = NumRows/NumBlockRows;
  }

  int NullSpaceDim = NumPDEEqns;
  double * NullSpacePtr = NULL;
  
  sprintf(parameter,"%snull space dimension", Prefix_);
  NullSpaceDim = List_.get(parameter, NumPDEEqns);
  sprintf(parameter,"%snull space vectors", Prefix_);
  NullSpacePtr = List_.get(parameter, NullSpacePtr);

  sprintf(parameter,"%scompute null space", Prefix_);
  bool ComputeNullSpace = List_.get(parameter, false);
  
  if( ComputeNullSpace == true ) {
    
    if( NullSpacePtr != NULL && Comm_.MyPID() == 0 ) {
      cerr << ErrorMsg_ << "Null space vectors is not NULL!" << endl
	   << ErrorMsg_ << "Now reallocating memory and computing null space " << endl
	   << ErrorMsg_ << "using eigenvectors estimates, for " << NullSpaceDim << " vectors..." << endl;
    }

    sprintf(parameter,"%sanasazi", Prefix_);
    ParameterList & AnasaziList = List_.sublist(parameter);
    ML_Anasazi_Interface(RowMatrix_,NullSpaceDim,NullSpacePtr,AnasaziList);
    NullSpaceToFree_ = NullSpacePtr; // this null space will be freed later
    
  }
  
  ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns,NullSpaceDim,NullSpacePtr,NumMyRows);

  /* ********************************************************************** */
  /* pick up coarsening strategy. METIS and ParMETIS requires additional    */
  /* lines, as we have to set the number of aggregates                      */
  /* ********************************************************************** */

  int value = -777;
  sprintf(parameter,"%saggregation: type",Prefix_);
  string CoarsenScheme = List_.get(parameter,"Hybrid");
  
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
    else if ( CoarsenScheme == "Hybrid" )
      ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(level,NumLevels_,agg_);
    else if(  CoarsenScheme == "Coupled" ) 
      ML_Aggregate_Set_CoarsenSchemeLevel_Coupled(level,NumLevels_,agg_);
    else {
      if( Comm_.MyPID() == 0 ) {
	cout << ErrorMsg_ << "specified options ("
	     << CoarsenScheme << ") not valid. Should be:" << endl;
	cout << ErrorMsg_ << "<METIS> <ParMETIS> <MIS> <Uncoupled> <Coupled>" << endl;
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

  /* ********************************************************************** */
  /* minor settings                                                         */
  /* ********************************************************************** */

  double DampingFactor = 0.01;
  sprintf(parameter,"%saggregation: damping factor", Prefix_);
  DampingFactor = List_.get(parameter, DampingFactor);
  ML_Aggregate_Set_DampingFactor( agg_, DampingFactor ); 
  
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

  if( verbose ) {
    cout << PrintMsg_ << "Aggregation damping factor = " << DampingFactor << endl;
    cout << PrintMsg_ << "Aggregation threshold =" << Threshold << endl;
    cout << PrintMsg_ << "Max coarse size = " << MaxCoarseSize << endl;
    cout << PrintMsg_ << "Requested next-level aggregates per process (ParMETIS) = " << ReqAggrePerProc << endl << endl;
  }

  /************************************************************************/
  /* Build hierarchy using smoothed aggregation.                          */
  /* Then, retrive parameters for each level. Default values are given by */
  /* entries in parameter list without (level %d)                         */
  /*----------------------------------------------------------------------*/

  if( IsIncreasing == "increasing" ) 
    NumLevels_ = ML_Gen_MGHierarchy_UsingAggregation(ml_, LevelID_[0], ML_INCREASING, agg_);
  else
    NumLevels_ = ML_Gen_MGHierarchy_UsingAggregation(ml_, LevelID_[0], ML_DECREASING, agg_);
  
  if( verbose ) cout << PrintMsg_ << "Number of actual levels : " << NumLevels_ << endl;

  sprintf(parameter,"%ssmoother: sweeps", Prefix_);
  int num_smoother_steps = List_.get(parameter, 1);

  sprintf(parameter,"%ssmoother: damping factor", Prefix_);
  double omega = List_.get(parameter,1.0);
// what about parallel????  JJH

  sprintf(parameter,"%ssmoother: pre or post", Prefix_);
  int pre_or_post;
  string PreOrPostSmoother = List_.get(parameter,"post");

  sprintf(parameter,"%ssmoother: type", Prefix_);
  string Smoother = List_.get(parameter,"MLS");
// Gauss Seidel ...... JJH

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
  int MLSPolynomialOrder = List_.get(parameter,2);

  sprintf(parameter,"%ssmoother: MLS alpha",Prefix_);
  double MLSalpha = List_.get(parameter,30.0);;

  /* ********************************************************************** */
  /* Now cycling over all levels                                            */
  /* ********************************************************************** */

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
      if( verbose ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : Jacobi (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_Jacobi(ml_, LevelID_[level], pre_or_post,
			     num_smoother_steps, omega);
    } else if( Smoother == "Gauss-Seidel" ) {
      if( verbose ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[level], pre_or_post,
				  num_smoother_steps, omega);
    } else if( Smoother == "symmetric Gauss-Seidel" ) {
      if( verbose ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : symmetric Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[level], pre_or_post,
				     num_smoother_steps, omega);
    } else if( Smoother == "block Gauss-Seidel" ) {
      if( verbose ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : block Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml_, LevelID_[level], pre_or_post,
				       num_smoother_steps, omega, NumPDEEqns);
    } else if( Smoother == "MLS" ) {
      sprintf(parameter,"smoother: MLS polynomial order (level %d)", LevelID_[level]);
      if( verbose ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : MLS,"
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

      if( verbose ) {
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
      if( verbose ) cout << PrintMsg_ << "Smoother (level " << LevelID_[level] << ") : Ifpack" << ","
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

  /* ********************************************************************** */
  /* solution of the coarse problem                                         */
  /* ********************************************************************** */

  if( NumLevels_ > 1 ) {  

    sprintf(parameter,"%scoarse: type", Prefix_);
    string CoarseSolution = List_.get(parameter, "SuperLU");
//JJH SuperLU
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
  
  ML_Gen_Solver(ml_, ML_MGV, LevelID_[0], LevelID_[NumLevels_-1]);

  ownership_ = false;

  CreateLabel();
  
  /* ********************************************************************** */
  /* Specific   preconditioners                                             */
  /* NOTE: the two-level DD preconditioners are kind of experimental!       */
  /* I suppose that the user knows what he/she is doing.... No real checks  */
  /* are performed. Note also that the coarse solver is somehow supposed to */
  /* be implemented as a post smoother (FIXME)                              */
  /* ********************************************************************** */

  sprintf(parameter,"%sprec type", Prefix_);
  string str = List_.get(parameter,"MGV");

  if( str == "1-level postsmoothing only" ) {
    
    sprintf(Label_, "1-level postsmoothing only");
    ml_->ML_scheme = ML_ONE_LEVEL_DD;

  } else if( str == "two-level additive DD" ) {

    sprintf(Label_, "two-level additive DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_ADD;
    if( NumLevels_ != 2 ) {
      if( Comm_.MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level additive DD' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level hybrid DD") {

    sprintf(Label_, "two-level hybrid DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID;
    if( NumLevels_ != 2 ) {
      if( Comm_.MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level hybrid DD' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level hybrid DD (2)") {

    sprintf(Label_, "two-level hybrid DD (2)");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID_2;
    if( NumLevels_ != 2 ) {
      if( Comm_.MyPID() == 0 ) {
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

    if( Comm_.MyPID() == 0 ) {
      cerr << ErrorMsg_ << "`prec type' has an incorrect value. It should be" << endl
	   << ErrorMsg_ << "<1-level postsmoothing only> / <two-level additive DD>" << endl
	   << ErrorMsg_ << "<two-level hybrid DD> / <two-level hybrid DD (2)>" << endl;
    }
    exit( EXIT_FAILURE );
    
  }

  if( verbose ) PrintLine();

  sprintf(parameter,"%sprint unused", Prefix_);
  if( List_.isParameter(parameter) ) {
    int ProcID = List_.get(parameter,-2);
    if( Comm_.MyPID() == ProcID || ProcID == -1 ) PrintUnused();
  }

  /* ------------------- that's all folks --------------------------------- */

  IsComputePreconditionerOK_ = true;
  
  return 0;
  
}

// ============================================================================

void Epetra_ML_Preconditioner::PrintUnused(int MyPID) 
{
  if( Comm_.MyPID() == MyPID ) {
    PrintLine();
    cout << PrintMsg_ << "Unused parameters:" << endl;
    PrintUnused();
    PrintLine();
  }
}

// ============================================================================

void Epetra_ML_Preconditioner::PrintList(int MyPID) 
{
  if( Comm_.MyPID() == MyPID ) {
    PrintLine();
    cout << List_;
    PrintLine();
  }
}

// ============================================================================

int Epetra_ML_Preconditioner::SetParameterList(const ParameterList & List) 
{
  if( IsComputePreconditionerOK_ == true ) DestroyPreconditioner();
  List_ = List;
  return 0;
  
}

// ============================================================================

int Epetra_ML_Preconditioner::SetDefaults(ParameterList & List, const string ProblemType)
{

  int rv = 0;
  
  if( ProblemType == "SA" )              return( SetDefaultsSA(List) );
  else if( ProblemType == "maxwell" )     return( SetDefaultsMaxwell(List) );
  else if( ProblemType == "DD 3-levels" ) return( SetDefaultsDD_3Levels(List) );
  else if( ProblemType == "DD" )          return( SetDefaultsDD_3Levels(List) );
  else if( ProblemType == "empty" ) {
    // empty list
    ParameterList NewList;
    List_ = List;
  } else {
    cerr << ErrorMsg_ << "Wrong input parameter in `SetDefaults'. Should be: " << endl
	 << ErrorMsg_ << "<SA> / <DD> / <DD 3-levels> / <maxwell> / <emtpy>" << endl;
    rv = 1;
  }

  EPETRA_CHK_ERR(rv);

  return rv;
  
  
}

// ============================================================================

int Epetra_ML_Preconditioner::SetDefaultsDD(ParameterList & List) 
{

  char parameter[80];

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,10);
  
  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type (level 0)",Prefix_);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: local aggregates (level 0)",Prefix_);
  List.set(parameter,1);
  
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
  
  sprintf(parameter,"%ssmoother: type (level 0)",Prefix_);
  List.set(parameter,"aztec");

  AZ_defaults(SmootherOptions_,SmootherParams_);
  SmootherOptions_[AZ_precond] = AZ_dom_decomp;
  SmootherOptions_[AZ_scaling] = AZ_none;
  SmootherOptions_[AZ_subdomain_solve] = AZ_ilut;
  
  sprintf(parameter,"%ssmoother: aztec options (level 0)",Prefix_);
  List.set(parameter,SmootherOptions_);

  sprintf(parameter,"%ssmoother: aztec params (level 0)",Prefix_);
  List.set(parameter,SmootherParams_);

  sprintf(parameter,"%ssmoother: aztec as solver (level 0)",Prefix_);
  List.set(parameter,false);
  
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos_KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,0);

  return 0;

}

// ============================================================================

int Epetra_ML_Preconditioner::SetDefaultsDD_3Levels(ParameterList & List) 
{

  char parameter[80];

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,3);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,10);
  
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
  
  sprintf(parameter,"%ssmoother: type (level 0)",Prefix_);
  List.set(parameter,"aztec");

  AZ_defaults(SmootherOptions_,SmootherParams_);
  SmootherOptions_[AZ_precond] = AZ_dom_decomp;
  SmootherOptions_[AZ_subdomain_solve] = AZ_ilut;
  
  sprintf(parameter,"%ssmoother: aztec options (level 0)",Prefix_);
  List.set(parameter,SmootherOptions_);

  sprintf(parameter,"%ssmoother: aztec params (level 0)",Prefix_);
  List.set(parameter,SmootherParams_);

  sprintf(parameter,"%ssmoother: aztec as solver (level 0)",Prefix_);
  List.set(parameter,false);
  
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos_KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,0);
  
  return 0;

}

// ============================================================================

int Epetra_ML_Preconditioner::SetDefaultsMaxwell(ParameterList & List) 
{}

// ============================================================================

int Epetra_ML_Preconditioner::SetDefaultsSA(ParameterList & List) 
{

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
  sprintf(parameter,"%saggregation: type (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,"Uncoupled");

  sprintf(parameter,"%saggregation: type (level %d)",Prefix_,MaxLevels-2);
  List.set(parameter,"Uncoupled");

  sprintf(parameter,"%saggregation: type (level %d)",Prefix_,MaxLevels-3);
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
  sprintf(parameter,"%ssmoother: sweeps (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: type (level %d)",Prefix_,MaxLevels-1);
  List.set(parameter,"Gauss-Seidel");
  
  sprintf(parameter,"%ssmoother: pre or post (level %d)",Prefix_,MaxLevels-1);
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

// ============================================================================

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
  sprintf(Label_,"%d level SA (%s, %s)", ml_->ML_num_actual_levels, finest, coarsest);
  
  return 0;
    
}

// ============================================================================

int Epetra_ML_Preconditioner::ApplyInverse(const Epetra_MultiVector& X,
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
  
  bool flag = false;
  if( ML_Time_->NumApplications() == 0 ) flag = true;
  ML_Time_->Update( Time.ElapsedTime(), flag);
  
  return 0;
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
