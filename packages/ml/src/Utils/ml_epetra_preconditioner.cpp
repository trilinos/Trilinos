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

#if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_operator.h"

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

  for( int i=0 ; i<MaxLevels_ ; ++i )  AZ_defaults(SmootherOptions_[i],SmootherParams_[i]);
    
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
    
    sprintf(parameter,"num nodes per aggregate, level %d",level );
    value = List_.getParameter(parameter,value);
    
    if( value != -777 )  ML_Aggregate_Set_NodesPerAggr(ml_,agg_,level,value );

  } /* for */

  /* ********************************************************************** */
  /* minor settings                                                         */
  /* ********************************************************************** */

  int MaxCoarseSize = List_.getParameter("max coarse size", 128);
  ML_Aggregate_Set_MaxCoarseSize(agg_,MaxCoarseSize );

  double Threshold = List_.getParameter("aggregation threshold", 0.0);
  ML_Aggregate_Set_Threshold(agg_,Threshold );

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
  cout<< CoarseSolution << endl;
  
  int MaxProcs = List_.getParameter("max processes for coarse solution", -1);
  
  if( CoarseSolution == "Amesos_KLU" )
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
  // fix it with more info ???

  CreateLabel();
  
}

int Epetra_ML_Preconditioner::CreateLabel()
{
  
  int i = ml_->ML_finest_level;
  char finest[80];
  char coarsest[80];
  char label[80];
  finest[0] = '\0';
  coarsest[0] = '\0';
  label[0] = '\0';
  
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
  sprintf(label,"%d level AMG (%s, %s)", ml_->ML_num_actual_levels, finest, coarsest);
  
  Label_ = label;
  
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
      ML_Solve_MGFull(ml_,
		      xvectors[i],  //rhs
		      yvectors[i]); //solution
      break;
    case(ML_SAAMG): //Marian Brezina's solver
      ML_Solve_AMGV(ml_,
		    xvectors[i],  //rhs
		    yvectors[i]); //solution
      break;
    default:
      ML_Solve_MGV(ml_,
		   xvectors[i],  //rhs
		   yvectors[i]); //solution
    }
  }
  
  return 0;
}

#else

  /*noop for certain compilers*/
int ML_EPETRA_EMPTY;
int ML_TEUCHOS_EMTPY;

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/

