/*!
 *  \file ml_MultiLevelPreconditioner.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */

#include "ml_common.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_memory.h"
#include "ml_DD_prec.h"
#include <iomanip>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//#include <cstring>
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"

#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_ParMETIS.h"

#include "ml_anasazi.h"

#ifdef HAVE_ML_TRIUTILS
#include "Trilinos_Util_CommandLineParser.h"
#endif

using namespace Teuchos;

// ================================================ ====== ==== ==== == =
void ML_Epetra::MultiLevelPreconditioner::PrintMem(char *fmt, int min, int avg, int max)
{

  if( Comm().MyPID() == 0 ) printf(fmt,min, avg, max);
  puts(" (Mb)");

  return;
}

// ================================================ ====== ==== ==== == =
int ML_Epetra::MultiLevelPreconditioner::DestroyPreconditioner()
{

  {
    int NumDestroy = OutputList_.get("number of destruction phases", 0);
    OutputList_.set("number of destruction phases", ++NumDestroy);
  }
  
  // may need to clean up after visualization and statistics
  ML_Aggregate_Stats_CleanUp_Info(ml_, agg_);
  ML_Aggregate_VizAndStats_Clean( agg_, NumLevels_);

  // destroy main objects
  if (agg_ != 0) { ML_Aggregate_Destroy(&agg_); agg_ = 0; }
  if (ml_ != 0) { ML_Destroy(&ml_); ml_ = 0; }
  if (ml_nodes_ != 0) { ML_Destroy(&ml_nodes_); ml_nodes_ = 0; }
  if (ml_edges_ != 0) { ML_Destroy(&ml_edges_); ml_edges_ = 0; }

  if (TMatrixML_ != 0) {
    ML_Operator_Destroy(&TMatrixML_);
    TMatrixML_ = 0;
  }

  if (TMatrixTransposeML_ != 0) {
    ML_Operator_Destroy(&TMatrixTransposeML_);
    TMatrixTransposeML_ = 0;
  }

  if (Tmat_array != 0) {
    ML_MGHierarchy_ReitzingerDestroy(LevelID_[0]-1, &Tmat_array,
				     &Tmat_trans_array);
    Tmat_array = 0;
    Tmat_trans_array = 0;
  }

  if (nodal_args_ != 0) {
    ML_Smoother_Arglist_Delete(&nodal_args_);
    nodal_args_ = 0;
  }
  
  if (edge_args_  != 0) {
    ML_Smoother_Arglist_Delete(&edge_args_);
    edge_args_ = 0;
  }

  if (Label_) { 
    delete [] Label_; 
    Label_ = 0; 
  }
  
  if (LevelID_) { 
    delete [] LevelID_; 
    LevelID_ = 0; 
  }
  
  if (Scaling_) {
    delete Scaling_;
    Scaling_ = 0;
  }

  if (InvScaling_) {
    delete InvScaling_;
    InvScaling_ = 0;
  }
  // stick data in OutputList

  OutputList_.set(Prefix_ + "time: total",
		  FirstApplicationTime_+ApplicationTime_);

  OutputList_.set(Prefix_ + "time: first application",
		  FirstApplicationTime_);

  OutputList_.set(Prefix_ + "time: construction",
		  ConstructionTime_);

  OutputList_.set(Prefix_ + "number of applications",
		  NumApplications_);

  int min[ML_MEM_SIZE], max[ML_MEM_SIZE], avg[ML_MEM_SIZE];
  for( int i=0 ; i<ML_MEM_SIZE ; ++i ) avg[i] = 0;

  if( AnalyzeMemory_ ) {
    memory_[ML_MEM_TOT2] = memory_[ML_MEM_TOT1] + memory_[ML_MEM_PREC_FIRST];
    memory_[ML_MEM_TOT2_USED] = memory_[ML_MEM_TOT1_USED] + memory_[ML_MEM_PREC_FIRST_USED];
    Comm().MinAll(memory_,min,ML_MEM_SIZE);
    Comm().MaxAll(memory_,max,ML_MEM_SIZE);
    Comm().SumAll(memory_,avg,ML_MEM_SIZE);
  }
  
  for( int i=0 ; i<ML_MEM_SIZE ; ++i ) avg[i] /= Comm().NumProc();

  if( verbose_ && NumApplications_ ) {

    // print on screen
    
    ML_print_line("-",78);
    double TotalTime = FirstApplicationTime_ + ApplicationTime_;
    cout << PrintMsg_ << "   ML time information" << endl << endl
	 << "   1- Construction time             = " << ConstructionTime_ << " (s)" << endl;
    cout << PrintMsg_ << "   2- Time for all applications     = " << TotalTime << " (s)" << endl;
    cout << PrintMsg_ << "   3- Time for first application(s) = " << FirstApplicationTime_ << " (s)" << endl;
    cout << PrintMsg_ << "   4- Each of " << NumApplications_
	 << " applications took " << TotalTime/NumApplications_ << " (s)" << endl;
    
  }

  if( verbose_ && AnalyzeMemory_ ) {
    
    // print memory usage

    cout << endl;
    cout << "   ML memory information:                 min     avg     max" << endl
	 << "   (warning: data may be incorrect for small runs, or" << endl
	 << "    if more processes share the same physical memory)" << endl
	 << endl
	 << "   1- max allocatable contiguous chunk, using malloc()" << endl;
    PrintMem("      before the ML construction     = %5d   %5d   %5d",
	     min[ML_MEM_INITIAL],avg[ML_MEM_INITIAL],max[ML_MEM_INITIAL]);
    PrintMem("      after the ML construction      = %5d   %5d   %5d",
	     min[ML_MEM_FINAL],avg[ML_MEM_FINAL],max[ML_MEM_FINAL]);
    cout << endl
	 << "   2- estimated ML memory usage, using malloc()" << endl;
    PrintMem("      for the hierarchy              = %5d   %5d   %5d",
	     min[ML_MEM_HIERARCHY],avg[ML_MEM_HIERARCHY],max[ML_MEM_HIERARCHY]);
    PrintMem("      for the smoother(s)            = %5d   %5d   %5d",
	     min[ML_MEM_SMOOTHER],avg[ML_MEM_SMOOTHER],max[ML_MEM_SMOOTHER]);
    PrintMem("      for the coarse solver          = %5d   %5d   %5d",
	     min[ML_MEM_COARSE],avg[ML_MEM_COARSE],max[ML_MEM_COARSE]);
    PrintMem("      preconditioning                = %5d   %5d   %5d",
	     min[ML_MEM_PREC_FIRST],avg[ML_MEM_PREC_FIRST],max[ML_MEM_PREC_FIRST]);
    PrintMem("      total (w/o other prec data)    = %5d   %5d   %5d",
	     min[ML_MEM_TOT1],avg[ML_MEM_TOT1],max[ML_MEM_TOT1]);
    PrintMem("      total (w/  other prec data)    = %5d   %5d   %5d",
	     min[ML_MEM_TOT2],avg[ML_MEM_TOT2],max[ML_MEM_TOT2]);
    cout << endl;
#ifdef ML_MALLINFO    
    cout << "   3- estimated ML memory usage, using mallinfo()" << endl;
    PrintMem("      for the hierarchy              = %5d   %5d   %5d",
	     min[ML_MEM_HIERARCHY_USED],avg[ML_MEM_HIERARCHY_USED],max[ML_MEM_HIERARCHY_USED]);
    PrintMem("      for the smoother(s)            = %5d   %5d   %5d",
	     min[ML_MEM_SMOOTHER_USED],avg[ML_MEM_SMOOTHER_USED],max[ML_MEM_SMOOTHER_USED]);
    PrintMem("      for the coarse solver          = %5d   %5d   %5d",
	     min[ML_MEM_COARSE_USED],avg[ML_MEM_COARSE_USED],max[ML_MEM_COARSE_USED]);
    PrintMem("      preconditioning                = %5d   %5d   %5d",
	     min[ML_MEM_PREC_FIRST_USED],avg[ML_MEM_PREC_FIRST_USED],max[ML_MEM_PREC_FIRST_USED]);
    PrintMem("      total (w/o other prec data)    = %5d   %5d   %5d",
	     min[ML_MEM_TOT1_USED],avg[ML_MEM_TOT1_USED],max[ML_MEM_TOT1_USED]);
    PrintMem("      total (w/  other prec data)    = %5d   %5d   %5d",
	     min[ML_MEM_TOT2_USED],avg[ML_MEM_TOT2_USED],max[ML_MEM_TOT2_USED]);
#endif
  }

  if (verbose_ && NumApplications_) 
    ML_print_line("-",78);
    
  if (NullSpaceToFree_ != 0) { 
    delete [] NullSpaceToFree_;  
    NullSpaceToFree_ = 0; 
  }
  
  if (RowMatrixAllocated_) { 
    delete RowMatrixAllocated_; 
    RowMatrixAllocated_ = 0; 
  }

  // filtering stuff

  if (flt_R_) { 
    delete flt_R_; 
    flt_R_ = 0; 
  }
  
  if (flt_NullSpace_) { 
    delete [] flt_NullSpace_; 
    flt_NullSpace_ = 0; 
  }
  
  if (flt_MatrixData_) { 
    delete flt_MatrixData_; 
    flt_MatrixData_ = 0; 
  }
  
  if (flt_ml_) { 
    ML_Destroy(&flt_ml_); 
    flt_ml_ = 0; 
  }
  
  if (flt_agg_) { 
    ML_Aggregate_Destroy(&flt_agg_); 
    flt_agg_ = 0; 
  }
  
  // CheckPreconditioner stuff
  if( SchurDecomposition_ ) { 
    delete SchurDecomposition_; 
    SchurDecomposition_ = 0; 
  }

  IsComputePreconditionerOK_ = false;

#ifdef ML_MEM_CHECK
  // print out allocated memory. It should be zero.
  if( Comm().MyPID() == 0 ) 
    cout << PrintMsg_ << "Calling ML_Print_it()..." << endl
         << PrintMsg_ << "no memory should be allocated by the ML preconditioner" 
	 << " at this point." << endl;
  ML_print_it();
#endif

  return 0;
}

// ================================================ ====== ==== ==== == =

ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
			 const bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  RowMatrixAllocated_(0)
{

  Prefix_ = "";
  
  ParameterList NewList;
  List_ = NewList;
  ML_Epetra::SetDefaults("SA",List_,(int *)0, (double *)0, Prefix_);
    
  ML_CHK_ERRV(Initialize());

  // construct hierarchy
  if (ComputePrec == true) 
    ML_CHK_ERRV(ComputePreconditioner());
}
  
// ================================================ ====== ==== ==== == =

ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner( const Epetra_RowMatrix & RowMatrix,
			 const ParameterList & List, const bool ComputePrec,
			 const char Prefix[] ) :
  RowMatrix_(&RowMatrix),
  RowMatrixAllocated_(0)
{
  Prefix_ = Prefix;

  List_ = List;

  ML_CHK_ERRV(Initialize());

  // construct hierarchy
  if (ComputePrec == true) 
    ML_CHK_ERRV(ComputePreconditioner());
}

// ================================================ ====== ==== ==== == =

/*! The constructor for the Maxwell equations requires three Epetra_RowMatrix's. 
 * Two conditions are required on their maps:
 * - TMatrix.OperatorDomainMap() == NodeMatrix.OperatorRangeMap()
 * - TMatrix.OperatorRangeMap()  == EdgeMatrix.OperatorDomainMap()
 */
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(const Epetra_RowMatrix & EdgeMatrix,
			 const Epetra_RowMatrix & TMatrix,
			 const Epetra_RowMatrix & NodeMatrix,
			 const ParameterList & List,
			 const bool ComputePrec,
			 const char Prefix[] ) :
  RowMatrix_(&EdgeMatrix),
  RowMatrixAllocated_(0)
{

  // some sanity checks
  if (! TMatrix.OperatorDomainMap().SameAs(NodeMatrix.OperatorRangeMap()) ) {
    cerr << ErrorMsg_ << "discrete grad DomainMap != node RangeMap..." << endl;
    ML_CHK_ERRV(-1); // error on discrete grad
  }

  if (! TMatrix.OperatorRangeMap().SameAs(EdgeMatrix.OperatorDomainMap()) ) {
    cerr << ErrorMsg_ << "discrete grad RangeMap != edge DomainMap..." <<endl;
    ML_CHK_ERRV(-2); // error on discrete grad
  }

  Prefix_ = Prefix;

  List_ = List;

  ML_CHK_ERRV(Initialize());

  // set Maxwell here.
  // NOTE: RowMatrix_ and EdgeMatrix_ pointer to the same Epetra_RowMatrix
  SolvingMaxwell_ = true;
  NodeMatrix_ = & NodeMatrix;
  TMatrix_ = & TMatrix;
  EdgeMatrix_ = & EdgeMatrix;

  // construct hierarchy
  if (ComputePrec == true) 
    ML_CHK_ERRV(ComputePreconditioner());

}

// ================================================ ====== ==== ==== == =
// FIXME: should I be deleted??
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(ML_Operator * Operator,
			 const ParameterList & List, const bool ComputePrec,
			 const char Prefix[] )
{

  // need to wrap an Epetra_RowMatrix around Operator.
  // This is quite not the best approach for small matrices

  int MaxNumNonzeros;
  double CPUTime;
  
  // FIXME: substitute with a better version
  ML_Operator2EpetraCrsMatrix(Operator,RowMatrixAllocated_,MaxNumNonzeros,
			      true,CPUTime);

  // this matrix must be freed by dtor. Keep trace of it in this pointer
  RowMatrix_ = RowMatrixAllocated_;

  // from now on as for the other constructors
  
  Prefix_ = Prefix;

  List_ = List;

  ML_CHK_ERRV(Initialize());

  // construct hierarchy
  if (ComputePrec == true) 
    ML_CHK_ERRV(ComputePreconditioner());
}

// ================================================ ====== ==== ==== == =
/*! - set to 0 all allocatable pointers
 *  - put default values in Aztec vectors.
 *  - zero-out timing
 */
int ML_Epetra::MultiLevelPreconditioner::Initialize()
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
  LevelID_ = 0;

  ml_ = 0;
  agg_ = 0;
  
  sprintf(ErrorMsg_,"*ML*ERR* : ");
  PrintMsg_ = "";
  
  AZ_defaults(SmootherOptions_,SmootherParams_);
  SmootherOptions_[AZ_precond] = AZ_dom_decomp;
  SmootherOptions_[AZ_subdomain_solve] = AZ_ilut;
  SmootherOptions_[AZ_overlap] = 0;

  // Maxwell stuff is off by default
  SolvingMaxwell_ = false;
  NodeMatrix_ = 0;
  EdgeMatrix_ = 0;
  TMatrix_ = 0;
  ml_edges_ = 0;
  ml_nodes_ = 0;
  TMatrixML_ = 0;
  TMatrixTransposeML_ = 0;
  Tmat_array = 0;
  Tmat_trans_array = 0;
  nodal_args_ = 0;
  edge_args_ = 0;
  
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

  // memory
  AnalyzeMemory_ = false;

  for (int i = 0 ; i < ML_MEM_SIZE ; ++i) 
    memory_[i] = 0;

  // filtering vectors
  flt_R_ = 0;
  flt_NullSpace_ = 0;
  flt_MatrixData_ = 0;
  flt_ml_ = 0;
  flt_agg_ = 0;

  // CheckPreconditioner stuff
  SchurDecomposition_ = 0;
  RateOfConvergence_ = -1.0;

  // scaling stuff
  Scaling_ = 0;
  InvScaling_ = 0;

  return 0;
}


// ================================================ ====== ==== ==== == =
// FIXME: test me of delete me??
int ML_Epetra::MultiLevelPreconditioner::ComputeFilteringPreconditioner()
{

  if (IsComputePreconditionerOK_ == true ){
    DestroyPreconditioner();
  }

  // 1.- disable filtering/GGB in ComputePreconditioner()
  List_.set(Prefix_ + "filtering: enable", false);

  ComputePreconditioner();

  // 2.- now enable, and call the function to compute the "bad-modes"
  //     (that is, the bad boys. Go Detroit !!!! Beat LA !!!)

  List_.set(Prefix_ + "filtering: enable", true);
  
  List_.set(Prefix_ + "filtering: type", "let ML be my master");

  if( NullSpaceToFree_ ) delete [] NullSpaceToFree_; NullSpaceToFree_ = 0;
  int NullSpaceDim = SetFiltering();
  assert( NullSpaceDim > 0 );
  
  List_.set(Prefix_ + "null space: type", "pre-computed");

  List_.set(Prefix_ + "null space: dimension", NullSpaceDim);
  
  List_.set(Prefix_ + "null space: vectors", NullSpaceToFree_);

  NullSpaceToFree_ = 0;
  
  DestroyPreconditioner();

  // 4.- recompute preconditioner with new options
  List_.set(Prefix_ + "filtering: enable", false);

  ComputePreconditioner();

  return 0;

}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::
ComputePreconditioner(const bool CheckPreconditioner)
{

  BreakForDebugger();

  // ============================================================== //
  // check whether the old filtering is still ok for the new matrix //
  // ============================================================== //

  if (CheckPreconditioner == true && RateOfConvergence_ != -1.0) {

    // If the previous preconditioner was computed with option
    // "adaptive: enable" == true, we know the rate of convergence
    // with the previous matrix (and this preconditioner. Now, 
    // we recompute this ratio, and compare it with the previous one.
    // This requires an AztecOO object to be defined (ML must have been
    // configured with --enable-aztecoo)

    if (CheckPreconditionerKrylov() == false) {
      ML_CHK_ERR(DestroyPreconditioner());
    }
    else                                       
      return 0;
    
  } else if( CheckPreconditioner == true && SchurDecomposition_ ) {
 
    // If the previous preconditioner was computed with option
    // "filtering: enable" == true, we have some eigen-information to use.
    // In this case, we
    // check whether the old preconditioner is still ok for the new matrix
    
    if (CheckPreconditionerFiltering() == false) {
      ML_CHK_ERR(DestroyPreconditioner());
    }
    else
      return 0;
    
  } else if( CheckPreconditioner == false && IsComputePreconditionerOK_ == true ) {
  
    // get rid of what done before 
    ML_CHK_ERR(DestroyPreconditioner());
    
  } // nothing else if left

  // ======================== //
  // build the preconditioner //
  // ======================== //

  Epetra_Time Time(Comm());
  {
    int NumCompute = OutputList_.get("number of construction phases", 0);
    OutputList_.set("number of construction phases", ++NumCompute);
  }
  
#ifdef HAVE_MPI
  const Epetra_MpiComm * MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());
  AZ_set_proc_config(ProcConfig_,MpiComm->Comm());
#else
  AZ_set_proc_config(ProcConfig_,AZ_NOT_MPI);
#endif

  // user's defined output message
  PrintMsg_ = List_.get(Prefix_ + "output prefix",PrintMsg_);
  
  NumLevels_ = List_.get(Prefix_ + "max levels",10);  

  int OutputLevel = List_.get(Prefix_ + "output", 10);  
  ML_Set_PrintLevel(OutputLevel);

  verbose_ = (5 < ML_Get_PrintLevel() && ProcConfig_[AZ_node] == 0);

  if( verbose_ ) 
    ML_print_line("-",78);
  
  FirstApplication_ = true;

  int call1 = 0, call2 = 0, call1_used = 0, call2_used = 0;

  AnalyzeMemory_ = List_.get(Prefix_ + "analyze memory", false);  

  if( AnalyzeMemory_ ) {
    memory_[ML_MEM_INITIAL] = ML_MaxAllocatableSize();
    call1 = memory_[ML_MEM_INITIAL];
    call1_used = ML_MaxMemorySize();
    if( verbose_ ) cout << "Memory : max allocatable block = " << call1 << " Mbytes" << endl;
  }

  if( Label_ ) delete [] Label_;
  
  Label_ = new char[80];
  
  // compute how to traverse levels (increasing of descreasing)
  // By default, use ML_INCREASING.
  
  string IsIncreasing = List_.get(Prefix_ + "increasing or decreasing",
				  "increasing");

  if( SolvingMaxwell_ == true ) IsIncreasing = "decreasing";
  
  int FinestLevel;

  LevelID_ = new int[NumLevels_];
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
    cout << PrintMsg_ << "*** " << endl;
    cout << PrintMsg_ << "*** ML_Epetra::MultiLevelPreconditioner" << endl;
    cout << PrintMsg_ << "***" << endl;
    cout << PrintMsg_ << "Matrix has " << RowMatrix_->NumGlobalRows()
	 << " rows, distributed over " << Comm().NumProc() << " process(es)" << endl;
    {
      const Epetra_CrsMatrix * dummy = dynamic_cast<const Epetra_CrsMatrix *>(RowMatrix_);
      if( dummy ) cout << "The linear system matrix is an Epetra_CrsMatrix" << endl;
    }	
    {
      const Epetra_VbrMatrix * dummy = dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
      if( dummy ) cout << "The linear system matrix is an Epetra_VbrMatrix" << endl;
    }	
    
    if( List_.isParameter("default values") ){
      cout << PrintMsg_ << "Default values for `" << List_.get("default values", "DD") << "'" << endl;
    }
    cout << PrintMsg_ << "Maximum number of levels = " << NumLevels_ << endl;
    if( IsIncreasing == "increasing" ) cout << PrintMsg_ << "Using increasing levels. ";
    else                               cout << PrintMsg_ << "Using decreasing levels. ";
    cout << "Finest level  = " << LevelID_[0];
    cout << ", coarsest level = " << LevelID_[NumLevels_-1] << endl;
  }


  int MaxCreationLevels = NumLevels_;
  if( IsIncreasing == "decreasing" )  MaxCreationLevels = FinestLevel+1;

  if( SolvingMaxwell_ == false ) {

    /* ********************************************************************** */
    /* create hierarchy for classical equations (all but Maxwell)             */
    /* ********************************************************************** */

    ML_Create(&ml_,MaxCreationLevels);
    ml_->output_level = OutputLevel;

    int NumMyRows;
    
    NumMyRows = RowMatrix_->NumMyRows();
    int N_ghost = RowMatrix_->NumMyCols() - NumMyRows;
    
    if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix
    
    ML_Init_Amatrix(ml_,LevelID_[0],NumMyRows, NumMyRows, (void *) RowMatrix_);
    ML_Set_Amatrix_Getrow(ml_, LevelID_[0], Epetra_ML_getrow,
			  Epetra_ML_comm_wrapper, NumMyRows+N_ghost);
    
    ML_Set_Amatrix_Matvec(ml_, LevelID_[0], Epetra_ML_matvec);

  } else {

    /* ********************************************************************** */
    /* create hierarchy for Maxwell. Needs to define ml_edges_ and ml_nodes_  */
    /* The only way to activate SolvingMaxwell_ == true is through the        */
    /* constructor for Maxwell. I suppose that the matrices are called not    */
    /* Matrix_, but moreover NodeMatrix_ and EdgeMatrix_.                     */
    /* ********************************************************************** */

    if( verbose_ ) cout << PrintMsg_ << "Solving Maxwell Equations..." << endl;

    int NumMyRows, N_ghost;
    
    // create hierarchy for edges

    ML_Create(&ml_edges_,MaxCreationLevels);

    NumMyRows = EdgeMatrix_->NumMyRows();
    N_ghost   = EdgeMatrix_->NumMyCols() - NumMyRows;
    
    if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix
    
    ML_Init_Amatrix(ml_edges_,LevelID_[0],NumMyRows,
		    NumMyRows, (void *) EdgeMatrix_);
    ML_Set_Amatrix_Getrow(ml_edges_, LevelID_[0], Epetra_ML_getrow,
			  Epetra_ML_comm_wrapper, NumMyRows+N_ghost);

    ML_Set_Amatrix_Matvec(ml_edges_, LevelID_[0], Epetra_ML_matvec);

    // create hierarchy for nodes
    
    ML_Create(&ml_nodes_,MaxCreationLevels);

    NumMyRows = NodeMatrix_->NumMyRows();
    N_ghost   = NodeMatrix_->NumMyCols() - NumMyRows;
    
    if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix
    
    ML_Init_Amatrix(ml_nodes_,LevelID_[0],NumMyRows, NumMyRows,
		    (void *) NodeMatrix_);
    ML_Set_Amatrix_Getrow(ml_nodes_, LevelID_[0], Epetra_ML_getrow,
			  Epetra_ML_comm_wrapper, NumMyRows+N_ghost);
    
    ML_Set_Amatrix_Matvec(ml_nodes_, LevelID_[0], Epetra_ML_matvec);
    
  }
  
  ML_Aggregate_Create(&agg_);
  
  /* **********************************************************************
   * visualize aggregate shape and other statistics. 
   * ********************************************************************** */
  
  bool viz = List_.get(Prefix_ + "viz: enable",false);
  if( viz == true )
    ML_Aggregate_VizAndStats_Setup(agg_,NumLevels_);

  /* ********************************************************************** */
  /* pick up coarsening strategy. METIS and ParMETIS requires additional    */
  /* lines, as we have to set the number of aggregates                      */
  /* ********************************************************************** */

  SetAggregation();

  /* ********************************************************************** */
  /* minor settings                                                         */
  /* ********************************************************************** */

  double Threshold = 0.0;
  Threshold = List_.get(Prefix_ + "aggregation: threshold", Threshold);
  ML_Aggregate_Set_Threshold(agg_,Threshold);
    
  int MaxCoarseSize = 50;
  MaxCoarseSize = List_.get(Prefix_ + "coarse: max size", MaxCoarseSize);
  ML_Aggregate_Set_MaxCoarseSize(agg_, MaxCoarseSize );

  int ReqAggrePerProc = 128;
  // FIXME: delete me???
  // compatibility with an older version
  if( List_.isParameter(Prefix_ + "aggregation: req aggregates per process") ) 
    ReqAggrePerProc = List_.get(Prefix_ + "aggregation: req aggregates per proces", ReqAggrePerProc);
  else {
    ReqAggrePerProc = List_.get(Prefix_ + "aggregation: next-level aggregates per process", ReqAggrePerProc);
  }

  if( SolvingMaxwell_ == false ) { 
    ML_Aggregate_Set_ReqLocalCoarseSize( ml_, agg_, -1, ReqAggrePerProc);
  } else {
    // Jonathan, is it right ???
    ML_Aggregate_Set_ReqLocalCoarseSize( ml_edges_, agg_, -1, ReqAggrePerProc);
    ML_Aggregate_Set_ReqLocalCoarseSize( ml_nodes_, agg_, -1, ReqAggrePerProc);
  }
  
  if( verbose_ ) {
    cout << PrintMsg_ << "Aggregation threshold = " << Threshold << endl;
    cout << PrintMsg_ << "Max coarse size = " << MaxCoarseSize << endl;
    
  }

  /* ********************************************************************** */
  /* METIS and ParMETIS may suffer (that is, core dump) is the graph is     */
  /* highly non-symmetric. In this case, it is better to set this parameter */
  /* to ML_NO. It does affect METIS and ParMETIS only.                      */
  /* ********************************************************************** */

  bool UseDropping = true;
  UseDropping = List_.get(Prefix_ + "aggregation: use dropping", UseDropping);
  if( UseDropping == true ) ML_Aggregate_Set_UseDropping( ML_YES );
  else                      ML_Aggregate_Set_UseDropping( ML_NO );  
 
  /* ********************************************************************** */
  /* alternatively, one can use dropping on a symmetrized matrix            */
  /* (although, this can be quite expensive for the finest levels)          */
  /* ********************************************************************** */
     
  if( SolvingMaxwell_ == false ) {
    bool UseSymmetrize = false;
    UseSymmetrize = List_.get(Prefix_ + "aggregation: symmetrize",
			      UseSymmetrize);
    if( UseSymmetrize == true ) ML_Set_Symmetrize(ml_, ML_YES );
    else                        ML_Set_Symmetrize(ml_, ML_NO );  
  }

  /* ********************************************************************** */
  /* Define scheme to determine damping parameter in prolongator and        */
  /* restriction smoother. Only for non-Maxwell.                            */
  /* ********************************************************************** */

  if( SolvingMaxwell_ == false ) {
    ML_CHK_ERR(SetSmoothingDamping());
  }
  else 
    ML_Aggregate_Set_DampingFactor( agg_, 0.0);

  /* ********************************************************************** */
  /* set scaling                                                            */
  /* ********************************************************************** */

  if (SolvingMaxwell_ == false) {
    ML_CHK_ERR(SetScaling());
  }

  /* ********************************************************************** */
  /* set null space                                                         */
  /* ********************************************************************** */

  if( SolvingMaxwell_ == false ) {
    ML_CHK_ERR(SetNullSpace());
  }
  else {
    ML_CHK_ERR(SetNullSpaceMaxwell());
  }
  
  /************************************************************************/
  /* Build hierarchy using smoothed aggregation.                          */
  /* Then, retrive parameters for each level. Default values are given by */
  /* entries in parameter list without (level %d)                         */
  /*----------------------------------------------------------------------*/

  int Direction;
  if (IsIncreasing == "increasing")
    Direction = ML_INCREASING;
  else
    Direction = ML_DECREASING;

  if (SolvingMaxwell_ == false) {

    bool CreateFakeProblem = 
      List_.get(Prefix_ + "aggregation: use auxiliary matrix", false);

    // I use FE matrix because I can set off-process elements. This may
    // help to create symmetric graphs (undirected graphs) from 
    // non-symmetric matrices (after dropping). 

    Epetra_FECrsMatrix * FakeMatrix = 0;

    if( CreateFakeProblem == true ) {
      ML_CHK_ERR(CreateAuxiliaryMatrix(FakeMatrix));
    }
    
    NumLevels_ = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml_, LevelID_[0], Direction, agg_);

    if( CreateFakeProblem == true ) {

      delete FakeMatrix;
      
      ml_->Amat[LevelID_[0]].data = (void *)RowMatrix_;

      // generate new hierarchy with "good" matrix
      if( verbose_ ) {
	cout << endl;
	cout << PrintMsg_ 
	     << "*** Now re-building the ML hierarchy with the actual matrix..." 
	     << endl;
	cout << endl;
      }

      //ML_Gen_MGHierarchy_UsingSmoothedAggr_ReuseExistingAgg(ml_, agg_);
      ML_Gen_MultiLevelHierarchy_UsingSmoothedAggr_ReuseExistingAgg(ml_, agg_);

    } // nothing special to be done if CreateFakeProblem is false
      
  } else {

    if( TMatrixML_ == 0 ) {

      // convert TMatrix to ML_Operator

      TMatrixML_ = ML_Operator_Create(ml_edges_->comm);

      Epetra2MLMatrix(const_cast<Epetra_RowMatrix*>(TMatrix_),TMatrixML_);

    }

    TMatrixTransposeML_ = ML_Operator_Create(ml_edges_->comm);
    ML_Operator_Transpose_byrow(TMatrixML_,TMatrixTransposeML_);
    
    NumLevels_ = ML_Gen_MGHierarchy_UsingReitzinger(ml_edges_, &ml_nodes_,LevelID_[0],
						    Direction,agg_,TMatrixML_,TMatrixTransposeML_, 
						    &Tmat_array,&Tmat_trans_array, 
						    ML_YES, 1.5); 


  }

  {
    int NL2 = OutputList_.get("max number of levels", 0);
    OutputList_.set("max number of levels", NL2+NumLevels_);
  }
  
  if( AnalyzeMemory_ ) {
    call2 = ML_MaxAllocatableSize();
    if( verbose_ ) cout << "Memory : max allocatable block = " << call2 << " Mbytes" << endl;    
    memory_[ML_MEM_HIERARCHY] = call1-call2;
    call1 = call2;

    call2_used = ML_MaxMemorySize();
    memory_[ML_MEM_HIERARCHY_USED] = call2_used-call1_used;
    call1_used = call2_used;
  }
  
  if( verbose_ ) cout << PrintMsg_ << "Number of actual levels : " << NumLevels_ << endl;

  /* ********************************************************************** */
  /* Now cycling over all levels                                            */
  /* ********************************************************************** */

  if( SolvingMaxwell_ == false ) {
    ML_CHK_ERR(SetSmoothers());
  }
  else {
    ML_CHK_ERR(SetSmoothersMaxwell());
  }

  if( AnalyzeMemory_ ) {
    call2 = ML_MaxAllocatableSize();
    if( verbose_ ) cout << "Memory : max allocatable block = " << call2 << " Mbytes" << endl;        
    memory_[ML_MEM_SMOOTHER] = call1 - call2;
    call1 = call2;
    call2_used = ML_MaxMemorySize();
    memory_[ML_MEM_SMOOTHER_USED] = call2_used - call1_used;
    call1_used = call2_used;
  }
  
  /* ********************************************************************** */
  /* solution of the coarse problem                                         */
  /* ********************************************************************** */

  if( NumLevels_ > 1 ) {
    ML_CHK_ERR(SetCoarse());
  }

  ownership_ = false;

  /* ********************************************************************** */
  /* Specific   preconditioners                                             */
  /* NOTE: the two-level DD preconditioners are kind of experimental!       */
  /* I suppose that the user knows what he/she is doing.... No real checks  */
  /* are performed. Note also that the coarse solver is somehow supposed to */
  /* be implemented as a post smoother (FIXME)                              */
  /* ********************************************************************** */

  if( AnalyzeMemory_ ) {
    call2 = ML_MaxAllocatableSize();
    if( verbose_ ) cout << "Memory : max allocatable block = " << call2 << " Mbytes" << endl;        
    memory_[ML_MEM_COARSE] = call1 - call2;
    call1 = call2;
    call2_used = ML_MaxMemorySize();
    memory_[ML_MEM_COARSE_USED] = call2_used - call1_used;
    call1_used = call2_used;
  }
  
  if( SolvingMaxwell_ == false ) 
    ML_Gen_Solver(ml_, ML_MGV, LevelID_[0], LevelID_[NumLevels_-1]);
  else {
    ML_Gen_Solver(ml_edges_, ML_MGV, LevelID_[0], LevelID_[NumLevels_-1]);
  }

  /* ********************************************************************** */
  /* Use of filtering functions, here called `filtering: enable'            */
  /* If this option is true, the code detects the non-converging modes of   */
  /* I - ML^{-1}A (where ML is the preconditioner we have just built), and  */
  /* creates a new V cycle to be added to the preconditioner. This part is  */
  /* equivalent to the GGB files of Haim Waisman (files ml_struct.c and     */
  /* ml_ggb.c, in the Main subdirectory).                                   */
  /* ********************************************************************** */

  ML_CHK_ERR(SetFiltering());
  
  /* ********************************************************************** */
  /* One may decide to print out the entire hierarchy (to be analyzed in    */
  /* MATLAB, for instance).                                                 */
  /* ********************************************************************** */

  bool PrintHierarchy = List_.get(Prefix_ + "print hierarchy", false);
  
  if( Comm().NumProc() > 1 && PrintHierarchy == true ) {
    if( Comm().MyPID() == 0 ) {
      cerr << endl;
      cerr << ErrorMsg_ << "Option `print hierarchy' == `true' is available" << endl
	<< ErrorMsg_ << "only for serial runs." << endl;
      cerr << endl;
    }
  }
  
  if( PrintHierarchy == true && Comm().NumProc() == 1 ) {
    if( Comm().MyPID() == 0 ) {
      cout << endl;
      cout << PrintMsg_ << "You are printing the entire hierarchy," << endl
	   << PrintMsg_ << "from finest level (" << LevelID_[0] 
	                << ") to coarsest (" << LevelID_[NumLevels_-1] << ")." << endl
	   << PrintMsg_ << "MATLAB can be used to load the matrices, using spconvert()" << endl;
      cout << endl;
    }

    // Amat (one for each level)
    for( int i=0 ; i<NumLevels_ ; ++i ) {
      char name[80];
      sprintf(name,"Amat_%d", LevelID_[i]);
      ML_Operator_Print(&(ml_->Amat[LevelID_[i]]), name);
    }
    
    // Pmat (one for each level, except last)
    for( int i=1 ; i<NumLevels_ ; ++i ) {
      char name[80];
      sprintf(name,"Pmat_%d", LevelID_[i]);
      ML_Operator_Print(&(ml_->Pmat[LevelID_[i]]), name);
    }

    // Rmat (one for each level, except first)
    for( int i=0 ; i<NumLevels_-1 ; ++i ) {
      char name[80];
      sprintf(name,"Rmat_%d", LevelID_[i]);
      ML_Operator_Print(&(ml_->Rmat[LevelID_[i]]), name);
    }
  
  }

  /* ********************************************************************** */
  /* Scaling was used on the matrix, need to scale back                     */
  /* ********************************************************************** */

  if (Scaling_) {
    // need to do the worst to get the best?
    // is C++ ethically correct?
    Epetra_RowMatrix* RM = const_cast<Epetra_RowMatrix*>(RowMatrix_);
    RM->RightScale(*InvScaling_);
  }

  /* ********************************************************************** */
  /* Other minor settings                                                   */
  /* ********************************************************************** */
  
  CreateLabel();
  
  if( SolvingMaxwell_ == false ) 
    SetPreconditioner();

  IsComputePreconditionerOK_ = true;
  
  // ====================================================================== //
  // Compute the rate of convergence (for adaptive preconditioners)         //
  // ====================================================================== //

  if( List_.get(Prefix_ + "adaptive: enable", false) == true ) 
    CheckPreconditionerKrylov();
  
  if( AnalyzeMemory_ ) {
    memory_[ML_MEM_FINAL] = ML_MaxAllocatableSize();
    memory_[ML_MEM_TOT1] = memory_[ML_MEM_INITIAL] - memory_[ML_MEM_FINAL];
    memory_[ML_MEM_FINAL_USED] = ML_MaxMemorySize();
    memory_[ML_MEM_TOT1_USED] = memory_[ML_MEM_FINAL_USED] - memory_[ML_MEM_INITIAL_USED];
  }
  
  // print unused parameters
  if( List_.isParameter(Prefix_ + "print unused") ) {
    int ProcID = List_.get(Prefix_ + "print unused",-2);
    if( Comm().MyPID() == ProcID || ProcID == -1 ) PrintUnused();
  }

  // ===================================================================== //
  // compute the coordinates for each level (that is, the center of        //
  // gravity, as no mesh is really available for coarser levels)           //
  // ===================================================================== //
 
  if (viz == true) {
    double * x_coord = List_.get(Prefix_ + "viz: x-coordinates", (double *)0);
    double * y_coord = List_.get(Prefix_ + "viz: y-coordinates", (double *)0);
    double * z_coord = List_.get(Prefix_ + "viz: z-coordinates", (double *)0);
    ML_Aggregate_Stats_ComputeCoordinates(ml_, agg_,
					  x_coord, y_coord, z_coord);
  }

  /* ------------------- that's all folks --------------------------------- */

  if( verbose_ )
    ML_print_line("-",78);

  ConstructionTime_ += Time.ElapsedTime();
  
  return 0;
  
}

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::PrintUnused(const int MyPID) const
{
  if( Comm().MyPID() == MyPID ) {
    ML_print_line("-",78);
    cout << PrintMsg_ << "Unused parameters:" << endl;
    PrintUnused();
    ML_print_line("-",78);
  }
}

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::PrintList(int MyPID) 
{
  if( Comm().MyPID() == MyPID ) {
    ML_print_line("-",78);
    cout << List_;
    ML_print_line("-",78);
  }
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::
SetParameterList(const ParameterList & List) 
{
  if( IsComputePreconditionerOK_ == true ) DestroyPreconditioner();
  List_ = List;
  return 0;
  
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::CreateLabel()
{  

  char finest[80];
  char coarsest[80];
  finest[0] = '\0';
  coarsest[0] = '\0';
  char * label;

  ML * ml_ptr;
  if( SolvingMaxwell_ == false ) ml_ptr = ml_;
  else                           ml_ptr = ml_edges_;
 
  int i = ml_ptr->ML_finest_level;
     
  if (ml_ptr->pre_smoother[i].smoother->func_ptr != NULL) {
    label = ml_ptr->pre_smoother[i].label;
    if( strncmp(label,"PreS_",4) == 0 ) sprintf(finest, "~");
    else                                sprintf(finest, "%s", label);
  } else                                sprintf(finest, "~"); 

  if (ml_ptr->post_smoother[i].smoother->func_ptr != NULL) {
    label = ml_ptr->pre_smoother[i].label;
    if( strncmp(label,"PostS_", 5) == 0 ) sprintf(finest, "%s/~", finest);
    else                                  sprintf(finest, "%s/%s", finest, label);
  } else                                  sprintf(finest, "%s/~", finest);  

  if (i != ml_ptr->ML_coarsest_level) {
    i = ml_ptr->ML_coarsest_level;
    if ( ML_CSolve_Check( &(ml_ptr->csolve[i]) ) == 1 ) {
	sprintf(coarsest, "%s", ml_ptr->csolve[i].label);
    }
    
    else {
      if (ml_ptr->pre_smoother[i].smoother->func_ptr != NULL) {
	label = ml_ptr->pre_smoother[i].label;
	if( strncmp(label,"PreS_",4) == 0 ) sprintf(coarsest, "~");
	else                                sprintf(coarsest, "%s", label);
      } else                                sprintf(coarsest, "~"); 
      if (ml_ptr->post_smoother[i].smoother->func_ptr != NULL) {
	label = ml_ptr->post_smoother[i].label;
	if( strncmp(label,"PostS_", 5) == 0 ) sprintf(coarsest, "%s/~", coarsest); 
	else                                  sprintf(coarsest, "%s/%s",coarsest, label);
      } else                                  sprintf(coarsest, "%s/~", coarsest); 
    }
  }

  if( SolvingMaxwell_ == false ) 
    sprintf(Label_,"ML (L=%d, %s, %s)", 
	    ml_->ML_num_actual_levels, finest, coarsest);
  else
    sprintf(Label_,"ML (Maxwell, L=%d, %s, %s)", 
	    ml_ptr->ML_num_actual_levels, finest, coarsest);
  
  return 0;
    
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::
ApplyInverse(const Epetra_MultiVector& X,
	     Epetra_MultiVector& Y) const
{

  // FIXME: What the hell am I doing here?
  if (Scaling_) {
    // need to do the worst to get the best?
    // is C++ ethically correct?
    Epetra_RowMatrix* RM = const_cast<Epetra_RowMatrix*>(RowMatrix_);
    RM->RightScale(*Scaling_);
  }

  int before = 0, after = 0, before_used = 0, after_used = 0;
  if( AnalyzeMemory_ ) {
    before = ML_MaxAllocatableSize();
    before_used = ML_MaxMemorySize();
  }
    
  Epetra_Time Time(Comm());
  
  if (!X.Map().SameAs(OperatorDomainMap())) 
    ML_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) 
    ML_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) 
    ML_CHK_ERR(-3);
  if( !IsPreconditionerComputed() ) 
    ML_CHK_ERR(-10);

  // FIXME: allocate me before, just once???
  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled
                              // in solver or if X = Y

  if (Scaling_ && false)
    xtmp.Multiply(1.0,xtmp,*Scaling_,0.0);

  Y.PutScalar(0.0); // Always start with Y = 0

  // ML_iterate doesn't handle multivectors, so extract and iterate one at
  // a time on them.
  double** xvectors;
  double** yvectors;
  ML_CHK_ERR(xtmp.ExtractView(&xvectors));
  ML_CHK_ERR(Y.ExtractView(&yvectors));

  ML * ml_ptr;
  
  if( SolvingMaxwell_ == false ) ml_ptr = ml_;
  else                           ml_ptr = ml_edges_;

  for (int i=0; i < X.NumVectors(); i++) {
    switch(ml_ptr->ML_scheme) {
    case(ML_MGFULLV):
      ML_Solve_MGFull(ml_ptr, xvectors[i], yvectors[i]); 
      break;
    case(ML_SAAMG): //Marian Brezina's solver
      ML_Solve_AMGV(ml_ptr, xvectors[i], yvectors[i]); 
      break;
    case(ML_ONE_LEVEL_DD):
      ML_DD_OneLevel(&(ml_ptr->SingleLevel[ml_ptr->ML_finest_level]),
		     yvectors[i], xvectors[i],
		     ML_ZERO, ml_ptr->comm, ML_NO_RES_NORM, ml_ptr);
      break;
    case(ML_TWO_LEVEL_DD_ADD):
      ML_DD_Additive(&(ml_ptr->SingleLevel[ml_ptr->ML_finest_level]),
		     yvectors[i], xvectors[i],
		     ML_ZERO, ml_ptr->comm, ML_NO_RES_NORM, ml_ptr);
      break;
    case(ML_TWO_LEVEL_DD_HYBRID):
      ML_DD_Hybrid(&(ml_ptr->SingleLevel[ml_ptr->ML_finest_level]),
		   yvectors[i], xvectors[i],
		   ML_ZERO, ml_ptr->comm, ML_NO_RES_NORM, ml_ptr);    
      break;
    case(ML_TWO_LEVEL_DD_HYBRID_2):
      ML_DD_Hybrid_2(&(ml_ptr->SingleLevel[ml_ptr->ML_finest_level]),
		     yvectors[i], xvectors[i],
		     ML_ZERO, ml_ptr->comm, ML_NO_RES_NORM, ml_ptr);    
      break;
    default:
      ML_Solve_MGV(ml_ptr, xvectors[i], yvectors[i]); 
    }
    
    if( flt_ml_ ) {
      ML_Cycle_MG(&(flt_ml_->SingleLevel[flt_ml_->ML_finest_level]),
		  yvectors[i], xvectors[i],
		  ML_NONZERO, flt_ml_->comm, ML_NO_RES_NORM, flt_ml_);
    }
  }

  if (Scaling_)
    Y.Multiply(1.0,*Scaling_,Y,0.0);

  if (Scaling_) {
    // need to do the worst to get the best?
    // is C++ ethically correct?
    Epetra_RowMatrix* RM = const_cast<Epetra_RowMatrix*>(RowMatrix_);
    RM->RightScale(*InvScaling_);
  }

  /* ********************************************************************** */
  /* filtering stuff if required                                            */
  /* ********************************************************************** */

  if (flt_R_) {
    
    // apply matvec
    if( Y.NumVectors() != 1 ) {
      if( verbose_ ) 
	cerr << ErrorMsg_ << "My dear user, filtering currently works only with one vector," << endl
	     << ErrorMsg_ << "I am very sorry for this, now I give up..." << endl;
      exit( EXIT_FAILURE );
    }

    xtmp.PutScalar(0.0);
    RowMatrix_->Multiply(false,Y,xtmp);
    
    xtmp.Update(1.0,X,-1.0);

    //    Epetra_MultiVector xtmp2(xtmp);
    /*
    xtmp2.PutScalar(0.0);
    RowMatrix_->Multiply(false,xtmp,xtmp2);
    */
    int size = flt_A_.N();

    for( int i=0 ; i<size ; ++i ) {
      double val = 0.0;
      (*flt_R_)(i)->Dot(*(xtmp(0)),&val);
      flt_rhs_(i) = val;
    }
    
    flt_lhs_.Multiply('N','N',1.0, flt_A_, flt_rhs_, 0.0);
    /*
      flt_solver_.SetVectors(flt_lhs_, flt_rhs_);
      flt_solver_.SolveToRefinedSolution(true);
      flt_solver_.Solve();
    */
    
    xtmp.PutScalar(0.0);
    for( int i=0 ; i<size ; ++i ) {
      double val = flt_lhs_(i);
      xtmp(0)->Update(val,*(*flt_R_)(i),1.0); 
    }
    
    Y.Update(1.0,xtmp,1.0);
    
  }
  
  /* ********************************************************************** */
  /* track timing                                                           */
  /* ********************************************************************** */

  MultiLevelPreconditioner * This = const_cast<MultiLevelPreconditioner *>(this);

  if( AnalyzeMemory_ ) {
    after = ML_MaxAllocatableSize();
    after_used = ML_MaxMemorySize();
  }
  
  double t = Time.ElapsedTime();
  if( FirstApplication_ ) {
    This->FirstApplication_ = false;
    This->FirstApplicationTime_ += t;
    This->memory_[ML_MEM_PREC_FIRST] = before - after;
    This->memory_[ML_MEM_PREC_FIRST_USED] = after_used - before_used;
  } else {
    This->memory_[ML_MEM_PREC_OTHER] = before - after;
    This->memory_[ML_MEM_PREC_OTHER_USED] = after_used - before_used;
  }
  
  This->ApplicationTime_ += t;
  
  ++(This->NumApplications_);

  return 0;
}

// ============================================================================
/*! Values for \c "coarse: type"
 * - \c Jacobi
 * - \c Gauss-Seidel
 * - \c SuperLU (deprecated)
 * - \c Amesos-KLU
 * - \c Amesos-UMFPACK
 * - \c Amesos-Superludist
 * - \c Amesos-MUMPS
 * - \c Amesos-ScALAPACK (under development in Amesos)
 * - \c do-nothing
 */
int ML_Epetra::MultiLevelPreconditioner::SetCoarse() 
{

  string CoarseSolution = List_.get(Prefix_ + "coarse: type", "Amesos-KLU");
  int NumSmootherSteps = List_.get(Prefix_ + "coarse: sweeps", 1);
  double Omega = List_.get(Prefix_ + "coarse: damping factor", 0.67);
    
  int MaxProcs = List_.get(Prefix_ + "coarse: max processes", -1);

  ML * ml_ptr;
  
  if( SolvingMaxwell_ == false ) ml_ptr = ml_;
  else                           ml_ptr = ml_edges_;
    
  if( CoarseSolution == "Jacobi" ) 
    ML_Gen_Smoother_Jacobi(ml_ptr, LevelID_[NumLevels_-1], ML_POSTSMOOTHER,
			   NumSmootherSteps, Omega);
  else if( CoarseSolution == "Gauss-Seidel" ) 
    ML_Gen_Smoother_GaussSeidel(ml_ptr, LevelID_[NumLevels_-1], ML_POSTSMOOTHER,
				NumSmootherSteps, Omega);
  else if( CoarseSolution == "SuperLU" ) 
    ML_Gen_CoarseSolverSuperLU( ml_ptr, LevelID_[NumLevels_-1]);
  else if( CoarseSolution == "Amesos-KLU" ) {
    ML_Gen_Smoother_Amesos(ml_ptr, LevelID_[NumLevels_-1], ML_AMESOS_KLU, MaxProcs);
  } else if( CoarseSolution == "Amesos-UMFPACK" )
    ML_Gen_Smoother_Amesos(ml_ptr, LevelID_[NumLevels_-1], ML_AMESOS_UMFPACK, MaxProcs);
  else if(  CoarseSolution == "Amesos-Superludist" )
    ML_Gen_Smoother_Amesos(ml_ptr, LevelID_[NumLevels_-1], ML_AMESOS_SUPERLUDIST, MaxProcs);
  else if( CoarseSolution == "Amesos-MUMPS" )
    ML_Gen_Smoother_Amesos(ml_ptr, LevelID_[NumLevels_-1], ML_AMESOS_MUMPS, MaxProcs);
  else if( CoarseSolution == "Amesos-ScaLAPACK" )
    ML_Gen_Smoother_Amesos(ml_ptr, LevelID_[NumLevels_-1], ML_AMESOS_SCALAPACK, MaxProcs);
  else if( CoarseSolution == "do-nothing" ) {
    // do nothing, ML will not use any coarse solver 
  } else {
    ML_Gen_Smoother_Amesos(ml_ptr, LevelID_[NumLevels_-1], ML_AMESOS_KLU, MaxProcs);
    ML_EXIT(-1); // Amesos is only the default...
  }
    
  return 0;
}

// ============================================================================
/*! Values for \c "aggregation: type"
 * - \c Uncoupled-MIS
 * - \c METIS
 * - \c ParMETIS
 * - \c Uncoupled
 * - \c Coupled (deprecated)
 * - \c MIS
 */
int ML_Epetra::MultiLevelPreconditioner::SetAggregation() 
{

  char parameter[80];
  
  int value = -777;
  string CoarsenScheme = List_.get(Prefix_ + "aggregation: type","Uncoupled");

  if ( CoarsenScheme == "Uncoupled-MIS" )
      ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_);
  else {
     for( int level=0 ; level<NumLevels_-1 ; ++level ) {  
   
       sprintf(parameter,"%saggregation: type (level %d)",
	       Prefix_.c_str(),LevelID_[level]);
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
	   cout << ErrorMsg_ << "<METIS> / <ParMETIS> / <MIS> / <Uncoupled> / <Coupled> / <Uncoupled-MIS>" << endl;
	 }
	 exit( EXIT_FAILURE );
       } 
   
       if( CoarsenScheme == "METIS" || CoarsenScheme == "ParMETIS" ) {
         
         bool isSet = false;

	 // first look for parameters without any level specification

	 sprintf(parameter,"%saggregation: global aggregates", 
		 Prefix_.c_str());
	 if( List_.isParameter(parameter) ){
	   value = -777; // simply means not set
	   value = List_.get(parameter,value);
	   if( value != -777 ) {
	     ML_Aggregate_Set_GlobalNumber(ml_,agg_,LevelID_[level],value );
	     isSet = true;
	   }
	 }

	 sprintf(parameter,"%saggregation: local aggregates", 
		 Prefix_.c_str());
	 if( List_.isParameter(parameter) ){
	   value = -777;
	   value = List_.get(parameter,value);
	   if( value != -777 ) {
	     ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value );
	     isSet = true;
       }
         }

	 sprintf(parameter,"%saggregation: nodes per aggregate", 
		 Prefix_.c_str());
	 if( List_.isParameter(parameter) ){
	   value = -777;
	   value = List_.get(parameter,value);
	   if( value != -777 ) {
	     ML_Aggregate_Set_NodesPerAggr(ml_,agg_,LevelID_[level],value );
	     isSet = true;
	   }
	 }

	 // now for level-specific data

	 sprintf(parameter,"%saggregation: global aggregates (level %d)", 
		 Prefix_.c_str(), LevelID_[level]);
	 if( List_.isParameter(parameter) ){
	   value = -777; // simply means not set
	   value = List_.get(parameter,value);
	   if( value != -777 ) {
	     ML_Aggregate_Set_GlobalNumber(ml_,agg_,LevelID_[level],value );
	     isSet = true;
	   }
	 }

	 sprintf(parameter,"%saggregation: local aggregates (level %d)", 
		 Prefix_.c_str(), LevelID_[level]);
	 if( List_.isParameter(parameter) ){
	   value = -777;
	   value = List_.get(parameter,value);
	   if( value != -777 ) {
	     ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value );
	     isSet = true;
	   }
	 }

	 sprintf(parameter,"%saggregation: nodes per aggregate (level %d)", 
		 Prefix_.c_str(), LevelID_[level]);
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
	   sprintf(parameter,"%saggregation: local aggregates (level %d)", 
		   Prefix_.c_str(), LevelID_[level]);
	   value = List_.get(parameter,1);
	   ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value);
	 }

       } // if( CoarsenScheme == "METIS" || CoarsenScheme == "ParMETIS" )

     } /* for */
  } /* else */

  return 0;
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetPreconditioner() 
{

  string str = List_.get(Prefix_ + "prec type","MGV");

  if( str == "one-level-postsmoothing" ) {
    
    sprintf(Label_, "1-level postsmoothing only");
    ml_->ML_scheme = ML_ONE_LEVEL_DD;

  } else if( str == "two-level-additive" ) {

    sprintf(Label_, "two-level additive DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_ADD;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level additive DD' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level-hybrid") {

    sprintf(Label_, "two-level hybrid DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level hybrid DD' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "two-level-hybrid2") {

    sprintf(Label_, "two-level hybrid DD (2)");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID_2;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
	cerr << ErrorMsg_ << "You asked for `two-level hybrid DD (2)' but you don't have" << endl
	     << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << endl;
      }
    }

  } else if( str == "full-MGV" ) {
    ml_->ML_scheme = ML_MGFULLV;

  } else if( str == "MGV" ) {
    // it is the default
    ml_->ML_scheme = ML_MGV;
  } else {

    if( Comm().MyPID() == 0 ) {
      cerr << ErrorMsg_ << "`prec type' has an incorrect value. It should be" << endl
	   << ErrorMsg_ << "<one-level-postsmoothing> / <two-level-additive>" << endl
	   << ErrorMsg_ << "<two-level-hybrid> / <two-level-hybrid2>" << endl;
    }
    exit( EXIT_FAILURE );
    
  }
  
  return 0;
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetSmoothingDamping() 
{

  Epetra_Time Time(Comm());
  
  /* ********************************************************************** */
  /* Strategies to determine the field-of-values.                           */
  /* almost everything here is experimental ;)                              */
  /* ********************************************************************** */

  string RandPSmoothing = List_.get(Prefix_ + "R and P smoothing: type", 
				    "classic");

  /* start looping over different options */

  if( RandPSmoothing == "classic" ) {

    /* ********************************************************************** */
    /* For "classical" approach to determine lambda_max only.                 */
    /* ********************************************************************** */

    ML_CHK_ERR(SetSmoothingDampingClassic());
    
  } else if( RandPSmoothing == "advanced" ) {

    /* ********************************************************************** */
    /* This is the new way, based on Anasazi to compute eigen-widgets         */
    /* ********************************************************************** */

    //    if( verbose_ )
    //      cout << PrintMsg_ << "Use A to smooth restriction operator" << endl;
    agg_->Restriction_smoothagg_transpose = ML_TRUE;
    
    SetEigenList();
    
    struct ML_Field_Of_Values * field_of_values;

    // stick default values (undefined)
    field_of_values = (struct ML_Field_Of_Values *) ML_allocate( sizeof(struct ML_Field_Of_Values) );
    field_of_values->eta     = 0.0;
    field_of_values->real_max= -1.0;
    field_of_values->imag_max= -1.0;
    field_of_values->poly_order = 0;

    if( List_.get(Prefix_ + "aggregation: compute field of values",true) )
      field_of_values->compute_field_of_values = ML_YES;
    else
      field_of_values->compute_field_of_values = ML_NO;
    
    if( List_.get(Prefix_ + "aggreation: compute field of values for non-scaled",false) )
      field_of_values->compute_field_of_values_non_scaled = ML_YES;
    else
      field_of_values->compute_field_of_values_non_scaled = ML_NO;

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
  
    if( DampingType == "non-smoothed" ) {

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

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

      agg_->Restriction_smoothagg_transpose = ML_FALSE;      
      
    } else if( DampingType == "almost-non-smoothed" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : almost non-smoothed aggregation" << endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  0;
      // I don't really need them, smoothing will be set to zero
      field_of_values->R_coeff[0] =  0.000000001;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;
      
      field_of_values->P_coeff[0] =  0.000000001;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "fov-1" ) {

      // those are the coefficients proposed by Ray
      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `fov-1' values" << endl;
      
      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
      
      field_of_values->R_coeff[0] =  1.107;
      field_of_values->R_coeff[1] =  0.285;
      field_of_values->R_coeff[2] =  0.718;
      
      field_of_values->P_coeff[0] =  1.878;
      field_of_values->P_coeff[1] = -2.515;
      field_of_values->P_coeff[2] =  0.942;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "fov-2" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `fov-2' values" << endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
      
      field_of_values->R_coeff[0] =  1.138;
      field_of_values->R_coeff[1] =  1.162;
      field_of_values->R_coeff[2] = -2.384;
      
      field_of_values->P_coeff[0] =  2.143;
      field_of_values->P_coeff[1] = -2.179;
      field_of_values->P_coeff[2] =  0.101;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "fov-5" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `fov-5' values" << endl;
   
      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;

      field_of_values->R_coeff[0] = 1.631;
      field_of_values->R_coeff[1] = -3.9015;
      field_of_values->R_coeff[2] = 2.5957;
      
      field_of_values->P_coeff[0] = 1.00145;
      field_of_values->P_coeff[1] = -1.4252;
      field_of_values->P_coeff[2] = 0.6627;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "fov-10" ) {

      // those are Marzio's best values ;^)
      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `fov-10' values" << endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;


      field_of_values->R_coeff[0] = 1.768909e+00;
      field_of_values->R_coeff[1] = -4.132227e+00;
      field_of_values->R_coeff[2] = 2.669318e+00;
      field_of_values->P_coeff[0] = 1.619455e+00;
      field_of_values->P_coeff[1] = -2.347773e+00;
      field_of_values->P_coeff[2] = 8.652273e-01;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "random" ) {

      // only to play with
      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `random' values" << endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
	    
      // initialize seed
      unsigned int s = (int)(Time.ElapsedTime()*10000);
      srand(s);
      
      // put random values
      field_of_values->R_coeff[0] =  (double)rand()/RAND_MAX;
      field_of_values->R_coeff[1] =  (double)rand()/RAND_MAX;
      field_of_values->R_coeff[2] =  (double)rand()/RAND_MAX;
      
      field_of_values->P_coeff[0] =  (double)rand()/RAND_MAX;
      field_of_values->P_coeff[1] =  (double)rand()/RAND_MAX;
      field_of_values->P_coeff[2] =  (double)rand()/RAND_MAX;

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

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "classic" ) {

      // This must be as with classic ML approach.
      // It can be used to estimate the field of value, tough.
      
      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `classic'" << endl;

      // First set damping as usual
      SetSmoothingDampingClassic();
      
      agg_->Restriction_smoothagg_transpose = ML_FALSE;      
      //      ml_->symmetrize_matrix == ML_TRUE;

    }  else if( DampingType == "classic-use-A" ) {

      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `classic-use-A'" << endl;

      field_of_values->choice     =  2;
      field_of_values->poly_order =  2;
      
      field_of_values->R_coeff[0] =  1.333;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;
      
      field_of_values->P_coeff[0] =  1.333;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else if( DampingType == "user-defined") {
      
      if( verbose_ )
	cout << PrintMsg_ << "R and P smoothing : Using `user-defined'" << endl;

      // user may specify each coefficient. Default values as for fov-1
      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;
      
      // get them from parameters' list.
      field_of_values->R_coeff[0] =  List_.get("R and P smoothing: c_0",  1.107);
      field_of_values->R_coeff[1] =  List_.get("R and P smoothing: c_1",  0.285);
      field_of_values->R_coeff[2] =  List_.get("R and P smoothing: c_2",  0.718);

      field_of_values->P_coeff[0] =  List_.get("R and P smoothing: g_0",  1.878);
      field_of_values->P_coeff[1] =  List_.get("R and P smoothing: g_1", -2.515);
      field_of_values->P_coeff[2] =  List_.get("R and P smoothing: g_2",  0.942);

      ML_Aggregate_Set_DampingFactor(agg_,0.0);
      
    } else {

      if( Comm().MyPID() == 0 ) {
	cerr << endl;
	cerr << ErrorMsg_ << "Parameter for `R and P smoothing : damping' not recognized" << endl
	     << ErrorMsg_ << "It is: `" << DampingType << "'. It should be one of:" << endl
	     << ErrorMsg_ << "<fov-1> / <fov-2> / <fov-5> / <fov-10> / <user-defined>" << endl;
      }

      ML_EXIT(-10); // wrong input parameter
    }

    agg_->field_of_values = (void*) field_of_values;  

  } else {

    if( Comm().MyPID() == 0 ) {
      cerr << endl;
      cerr << ErrorMsg_ << "Parameter for `R and P smoothing : type' not recognized" << endl
	   << ErrorMsg_ << "It is: `" << RandPSmoothing << "'. It should be one of:" << endl
	   << ErrorMsg_ << "<classic> / <advanced>" << endl;
    }
    
    ML_EXIT(-11); // wrong input parameter
    
  }

  return 0;
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::SetSmoothingDampingClassic()
{
  
  double DampingFactor = 1.333;
  if( SolvingMaxwell_ ) DampingFactor = 0.0;

  DampingFactor = List_.get(Prefix_ + "aggregation: damping factor", 
			    DampingFactor);
  ML_Aggregate_Set_DampingFactor( agg_, DampingFactor );
  
  if( verbose_ ) {
    cout << PrintMsg_ << "R and P smoothing : P = (I-\\omega A) P_t, R = P^T" << endl;
    cout << PrintMsg_ << "R and P smoothing : \\omega = " << DampingFactor << "/lambda_max" <<endl;
  }
    
  string str = List_.get(Prefix_ + "eigen-analysis: type","Anorm");
  
  if( verbose_ ) cout << PrintMsg_ << "Using `" << str << "' scheme for eigen-computations" << endl;
  
  if( str == "cg" )                ML_Aggregate_Set_SpectralNormScheme_Calc(agg_);
  else if( str == "Anorm" )        ML_Aggregate_Set_SpectralNormScheme_Anorm(agg_);
  else if( str == "Anasazi" )      ML_Aggregate_Set_SpectralNormScheme_Anasazi(agg_);
  else if( str == "power-method" ) ML_Aggregate_Set_SpectralNormScheme_PowerMethod(agg_);
  else {
    if( Comm().MyPID() == 0 ) {
      cerr << ErrorMsg_ << "parameter `eigen-analysis: type' has an incorrect value"
	   << "(" << str << ")" << endl;
      cerr << ErrorMsg_ << "It should be: " << endl
	   << ErrorMsg_ << "<cg> / <Anorm> / <Anasazi> / <power-method>" << endl;
    }
    ML_EXIT(-10); // wrong input parameter
  }
    
  return 0;
}
 
// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::
PrintStencil2D(const int nx, const int ny, 
	       int NodeID,
	       const int EquationID)
{

  if (nx <= 0)
    ML_CHK_ERR(-1); // need nodes along the X-axis
 
  if (ny <= 0)
    ML_CHK_ERR(-2); // need nodes along Y-axis


  if (RowMatrix_ == 0) 
    ML_CHK_ERR(-3); // matrix still not set

  // automatically compute NodeID, somewhere in the middle of the grid
  if (NodeID == -1) {
    if (ny == 1) 
      NodeID = (int)(nx/2);
    else 
      NodeID = (int)(ny*(nx/2) + nx/2);
  }
  
  // need to convert from NodeID (BlockRowID) to PointRowID
  int GID = NodeID * NumPDEEqns_;
  
  int LID = RowMatrix_->RowMatrixRowMap().LID(GID);

  // only processor having this node will go on
  if (LID == -1)
    return(0);

  int MaxPerRow = RowMatrix_->MaxNumEntries();
  int NumEntriesRow;   // local entries on each row
  vector<double> Values;  Values.resize(MaxPerRow);
  vector<int>    Indices; Indices.resize(MaxPerRow);
  
  int ierr = RowMatrix_->ExtractMyRowCopy(LID, MaxPerRow, NumEntriesRow,
					  &Values[0], &Indices[0]);

  if (ierr) 
    ML_CHK_ERR(-4);

  // cycle over nonzero elements, look for elements in positions that we
  // can understand
  
  Epetra_IntSerialDenseMatrix StencilInd(3,3);
  Epetra_SerialDenseMatrix StencilVal(3,3);
  for( int i=0 ; i<3 ; ++i ) 
    for( int j=0 ; j<3 ; ++j ) {
      StencilVal(i,j) = 0.0;
    }
  
  // look for the following positions
  StencilInd(0,0) = RowMatrix_->RowMatrixColMap().LID(NodeID-1-nx);
  StencilInd(1,0) = RowMatrix_->RowMatrixColMap().LID(NodeID-nx);
  StencilInd(2,0) = RowMatrix_->RowMatrixColMap().LID(NodeID+1-nx);
  StencilInd(0,1) = RowMatrix_->RowMatrixColMap().LID(NodeID-1);
  StencilInd(1,1) = RowMatrix_->RowMatrixColMap().LID(NodeID);
  StencilInd(2,1) = RowMatrix_->RowMatrixColMap().LID(NodeID+1);
  StencilInd(0,2) = RowMatrix_->RowMatrixColMap().LID(NodeID-1+nx);
  StencilInd(1,2) = RowMatrix_->RowMatrixColMap().LID(NodeID+nx);
  StencilInd(2,2) = RowMatrix_->RowMatrixColMap().LID(NodeID+1+nx);

  for( int i=0 ; i<NumEntriesRow ; ++i ) {
    // get only the required equation
    if( Indices[i]%NumPDEEqns_ ) continue;
    // convert into block row
    int LocalColID = Indices[i]/NumPDEEqns_;
    // look for known positions
    for( int ix=0 ; ix<3 ; ++ix ) {
      for( int iy=0 ; iy<3 ; ++iy ) {
	if( StencilInd(ix,iy) == LocalColID ) {
	  StencilVal(ix,iy) = Values[i];
	}
      }
    }
  }
  
  cout << "2D computational stencil for equation " << EquationID << " at node " << NodeID
       << " (grid is " << nx << " x " << ny << ")" << endl;
  cout << endl;
  for( int iy=0 ; iy<3 ; ++iy ) {
    cout << "\t";
    for( int ix=0 ; ix<3 ; ++ix ) {
      cout << " " << std::setw(15) << StencilVal(ix,iy);
    }
    cout << endl;
  }
  cout << endl;
  
  return 0;
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::BreakForDebugger()
{
  // print out some junk for debugging (copied from code in
  // Utils/ml_utils.c, suggested by Jonathan)
  // LAM/MPI has some difficulties related to environmental variables.
  // The problem is that LAM/MPI uses ssh to log in, and on some
  // machine "export ML_BREAK_FOR_DEBUGGER" does not work. So, we
  // allow two options:
  // 1.) export ML_BREAK_FOR_DEBUGGER=1
  // 2.) create a file in the executable directory, called ML_debug_now

  char * str = (char *) getenv("ML_BREAK_FOR_DEBUGGER");
  int i = 0, j = 0;
  char buf[80];
  char go = ' ';
  char hostname[80];
  if (str != NULL) i++;

  FILE * ML_capture_flag;
  ML_capture_flag = fopen("ML_debug_now","r");
  if(ML_capture_flag) {
    i++;
    fclose(ML_capture_flag);
  }

  Comm().SumAll(&i, &j, 1);

  if (j != 0)
  {
    if (Comm().MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
    for (i = 0; i <Comm().NumProc() ; i++) {
      if (i == Comm().MyPID() ) {
#ifdef COUGAR
	sprintf(buf, "Host: %s   PID: %d", "janus", getpid());
#else
	gethostname(hostname, sizeof(hostname));
	sprintf(buf, "Host: %s\tComm().MyPID(): %d\tPID: %d", 
		hostname, Comm().MyPID(), getpid());
#endif
	printf("%s\n",buf);
	fflush(stdout);
	sleep(1);
      }
    }
    if(Comm().MyPID() == 0) {
      printf("\n");
      printf("** Pausing because environment variable ML_BREAK_FOR_DEBUGGER has been set,\n");
      puts("** or file ML_debug_now has been created");
      printf("**\n");
      printf("** You may now attach debugger to the processes listed above.\n");
      printf( "**\n");
      printf( "** Enter a character to continue > "); fflush(stdout);
      scanf("%c",&go);
    }
  }

  return 0;
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::CreateAuxiliaryMatrix(Epetra_FECrsMatrix * & FakeMatrix)
{

  int NumMyRows = RowMatrix_->NumMyRows();

  // create the auxiliary matrix as VBR matrix. This should help
  // to save memory with respect to the creation of "pure" VBR
  // matrices.

  // FIXME
  Epetra_Map Map(-1,NumMyRows,0,Comm());
  FakeMatrix = new Epetra_FECrsMatrix(Copy,Map,
				      2*RowMatrix_->MaxNumEntries());
  if (FakeMatrix == 0) {
    ML_CHK_ERR(-1); // something went wrong
  }

  int NumDimensions = 0;

  double* x_coord = List_.get(Prefix_ + "aggregation: x-coordinates",
			      (double *)0);
  if( x_coord != 0 ) ++NumDimensions;

  // at least x-coordinates myst be not null
  if( NumDimensions == 0 ) {
    cerr << ErrorMsg_ << "Option `aggregation: use auxiliary matrix' == true" << endl
      << ErrorMsg_ << "requires x-, y-, or z-coordinates." << endl
      << ErrorMsg_ << "You must specify them using options" << endl
      << ErrorMsg_ << "`aggregation: x-coordinates' (and equivalently for" << endl
      << ErrorMsg_ << "y- and z-." << endl;
    ML_CHK_ERR(-2); // wrong parameters
  }

  double * y_coord = List_.get(Prefix_ + "aggregation: y-coordinates",
			       (double *)0);
  if( y_coord != 0 ) ++NumDimensions;

      double * z_coord = List_.get(Prefix_ + "aggregation: z-coordinates",
				   (double *)0);
      if( z_coord != 0 ) ++NumDimensions;

      // small check to avoid strange behavior
      if( z_coord != 0 && y_coord == 0 ) {
        cerr << ErrorMsg_ << "Something wrong: `aggregation: y-coordinates'" << endl
             << ErrorMsg_ << "is null, while `aggregation: z-coordinates' is null" << endl;
	ML_CHK_ERR(-3); // something went wrong
      }
  
  double theta = List_.get(Prefix_ + "aggregation: theta",0.0);

  bool SymmetricPattern = List_.get(Prefix_ + "aggregation: use symmetric pattern",false);

  // usual crap to clutter the output
  if( verbose_ ) {
    cout << endl;
    cout << PrintMsg_ << "*** Using auxiliary matrix to create the aggregates" << endl
      << PrintMsg_ << "*** Number of dimensions = " << NumDimensions << endl
      << PrintMsg_ << "*** theta = " << theta;
    if( SymmetricPattern ) cout << ", using symmetric pattern" << endl;
    else                   cout << ", using original pattern" << endl;
    cout << endl;
  }

  // create the auxiliary matrix

  int allocated = 128;
  int* colInd = new int[allocated];
  double* colVal = new double[allocated];
  int NnzRow;

  double coord_i[3], coord_j[3];
  for( int i = 0; i<3 ; ++i ) {
    coord_i[i] = 0.0; coord_j[i] = 0.0;
  }

  // get global column number
  int isize  = RowMatrix_->OperatorDomainMap().NumMyElements();
  int Nghost = RowMatrix_->RowMatrixColMap().NumMyElements() - isize;
  if (Nghost < 0) Nghost = 0; 
  assert(isize == NumMyRows);

  vector<double> global_isize;        global_isize.resize(isize+Nghost+1);
  vector<int>    global_isize_as_int; global_isize_as_int.resize(isize+Nghost+1);

  int isize_offset;
  Comm().ScanSum(&isize,&isize_offset,1); 
  isize_offset -= isize;

  for (int i = 0 ; i < isize ; ++i) {
    global_isize[i] = (double) (isize_offset + i);
  }

  for (int i = 0 ; i < Nghost; ++i) 
    global_isize[i+isize] = -1;

  ML_exchange_bdry(&global_isize[0],(&(ml_->Amat[LevelID_[0]]))->getrow->pre_comm,
                  (&(ml_->Amat[LevelID_[0]]))->invec_leng,ml_->comm,ML_OVERWRITE,NULL);

  for ( int j = 0; j < isize+Nghost; ++j) {
    global_isize_as_int[j] = (int) global_isize[j];
  }
  
  // =================== //
  // cycle over all rows //
  // =================== //

  for (int i = 0; i < NumMyRows ; i += NumPDEEqns_) {

    int GlobalRow = global_isize_as_int[i];

    assert(GlobalRow != -1);

    if( i%NumPDEEqns_ == 0 ) { // do it just once for each block row
      switch( NumDimensions ) {
      case 3:
	coord_i[2] = z_coord[i/NumPDEEqns_];
      case 2:
	coord_i[1] = y_coord[i/NumPDEEqns_];
      case 1:
	coord_i[0] = x_coord[i/NumPDEEqns_];

      }

      int ierr = ML_Operator_Getrow(&(ml_->Amat[LevelID_[0]]),1,&i,
				    allocated,colInd,colVal,&NnzRow);

      if (ierr == 0) {
	do {
	  delete [] colInd;
	  delete [] colVal;
	  allocated *= 2;
	  colInd = new int[allocated];
	  colVal = new double[allocated];
	  ierr = ML_Operator_Getrow(&(ml_->Amat[LevelID_[0]]),1,&i,
				    allocated,colInd,colVal,&NnzRow);
	} while (ierr == 0);
      }

      // NOTE: for VBR matrices, the "real" value that will be used in
      // the subsequent part of the code is only the one for the first
      // equations. For each block, I replace values with the sum of
      // the abs of each block entry.

      for (int j = 0 ; j < NnzRow ; j += NumPDEEqns_) {
	colVal[j] = abs(colVal[j]);
	for (int k = 1 ; k < NumPDEEqns_ ; ++k) {
	  colVal[j] += abs(colVal[j+k]);
	}
      }

      // work only on the first equations. Theta will blend the
      // coordinate part with the sub of abs of row elements.

      int GlobalCol;
      double total = 0.0;

      for (int j = 0 ; j < NnzRow ; j += NumPDEEqns_) {

	if (colInd[j]%NumPDEEqns_ == 0) { 

	  // insert diagonal later
	  if( colInd[j] != i ) {

	    // get coorinates of this node
	    switch( NumDimensions ) {
	    case 3:
	      coord_j[2] = z_coord[colInd[j]/NumPDEEqns_];
	    case 2:
	      coord_j[1] = y_coord[colInd[j]/NumPDEEqns_];
	    case 1:
	      coord_j[0] = x_coord[colInd[j]/NumPDEEqns_];
	    }

	    // d2 is the square of the distance between node `i' and
	    // node `j'
	    double d2 = (coord_i[0] - coord_j[0]) * (coord_i[0] - coord_j[0]) +
	      (coord_i[1] - coord_j[1]) * (coord_i[1] - coord_j[1]) +		     
	      (coord_i[2] - coord_j[2]) * (coord_i[2] - coord_j[2]);

	    if (d2 == 0.0) {
	      cerr << endl;
	      cerr << ErrorMsg_ << "distance between node " << i/NumPDEEqns_ << " and node " 
		<< colInd[j]/NumPDEEqns_ << endl
		<< ErrorMsg_ << "is zero. Coordinates of these nodes are" << endl
		<< ErrorMsg_ << "x_i = " << coord_i[0] << ", x_j = " << coord_j[0] << endl  
		<< ErrorMsg_ << "y_i = " << coord_i[1] << ", y_j = " << coord_j[1] << endl  
		<< ErrorMsg_ << "z_i = " << coord_i[2] << ", z_j = " << coord_j[2] << endl  
		<< ErrorMsg_ << "Now proceeding with distance = 1.0" << endl;
	      cerr << endl;
	      d2 = 1.0;
	    }

	    // blend d2 with the actual values of the matrix
	    // FIXME: am I useful?
	    double val = -(1.0-theta)*(1.0/d2) + theta*(colVal[j]);

	    GlobalCol = global_isize_as_int[colInd[j]];
	    if (GlobalCol == -1)
	      ML_CHK_ERR(-9);

	    // insert this value on all rows
	    for (int k = 0 ; k < NumPDEEqns_ ; ++k) {
	      int row = GlobalRow+k;
	      int col = GlobalCol+k;
	      if (row >= NumGlobalRows())
		ML_CHK_ERR(-9);
	      if (col >= NumGlobalCols())
		ML_CHK_ERR(-9);
/*
	      assert(FakeMatrix->InsertGlobalValues(row,1,&val,&col) == 0);
*/
	      if( FakeMatrix->SumIntoGlobalValues(1,&row,1,&col,&val) != 0 ) {
		assert(FakeMatrix->InsertGlobalValues(1,&row,1,&col,&val) == 0);
	      }

	    }

	    total -= val;

	    // put (j,i) element as well, only for in-process stuff.
	    // I have some problems with off-processor elements
	    // FIXME: leave the FECrsMatrix???
	    if (SymmetricPattern == true && colInd[j] < NumMyRows ) {

	      for( int k=0 ; k<NumPDEEqns_ ; ++k ) {
		int row = GlobalCol+k;
		int col = GlobalRow+k;
		if (row >= NumGlobalRows())
		  ML_CHK_ERR(-9);
		if (col >= NumGlobalCols())
		  ML_CHK_ERR(-9);
/*
		assert(FakeMatrix->InsertGlobalValues(row,1,&val,&col) == 0);
*/
		if( FakeMatrix->SumIntoGlobalValues(1,&row,1,&col,&val) != 0 ) { 
		  assert(FakeMatrix->InsertGlobalValues(1,&row,1,&col,&val) == 0);
		}

	      }
	      total -= val;
	    }
	  } 
	}
      }

      // create lines with zero-row sum
      for (int k = 0 ; k < NumPDEEqns_ ; ++k) {
	int row = GlobalRow + k;
	if (row >= NumGlobalRows())
	  ML_CHK_ERR(-9); // something went wrong
	/*
	assert(FakeMatrix->InsertGlobalValues(row,1,&total,&row) == 0);
	*/
	if (FakeMatrix->SumIntoGlobalValues(1,&row,1,&row,&total) != 0) {
	  if (FakeMatrix->InsertGlobalValues(1,&row,1,&row,&total) != 0)
	    ML_CHK_ERR(-9); // something went wrong
	}

      }
    }
  }

  if (FakeMatrix->FillComplete())
    ML_CHK_ERR(-5); // something went wrong

  delete [] colInd;
  delete [] colVal;

  // stick pointer in Amat for level 0 (finest level)
  ml_->Amat[LevelID_[0]].data = (void *)FakeMatrix;

  // tell ML to keep the tentative prolongator
  ML_Aggregate_Set_Reuse(agg_);

  return 0;
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
