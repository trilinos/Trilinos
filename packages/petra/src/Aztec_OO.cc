#include <stdlib.h>
#include "Aztec_OO.h"


//=============================================================================
Aztec_OO::Aztec_OO(Petra_RDP_RowMatrix * A, 
                   Petra_RDP_MultiVector * X,
                   Petra_RDP_MultiVector * B) 
  : A_(A), X_(X), B_(B), Prec_(0), Scaling_(0), azVarsAllocated_(false), condest_(-1)
{
  AllocAzArrays();
  SetAztecDefaults();

  if (SetAztecVariables() != 0) {
     cerr << __FILE__ << ", line "<<__LINE__<<
         ": error in SetAztecVariables. aborting." << endl;
     abort();
  }
}

//=============================================================================
Aztec_OO::Aztec_OO(const Petra_RDP_LinearProblem& prob) 
 : A_(prob.GetOperator()), X_(prob.GetLHS()), B_(prob.GetRHS()), Prec_(0), Scaling_(0), 
   azVarsAllocated_(false), condest_(-1)
{
  AllocAzArrays();
  SetAztecDefaults();

  if (SetAztecVariables() != 0) {
     cerr << __FILE__ << ", line "<<__LINE__<<
         ": error in SetAztecVariables. aborting." << endl;
     abort();
  }
  SetProblemOptions(prob.GetPDL(), prob.IsOperatorSymmetric());
}

//=============================================================================
Aztec_OO::Aztec_OO() 
 : A_(NULL), X_(NULL), B_(NULL), Prec_(0), Scaling_(0), azVarsAllocated_(false), condest_(-1)
{
  AllocAzArrays();
  SetAztecDefaults();
}

//=============================================================================
Aztec_OO::Aztec_OO(const Aztec_OO& azOO) 
  : A_(azOO.A_), X_(azOO.X_), B_(azOO.B_), Prec_(azOO.Prec_), Scaling_(azOO.Scaling_), 
    azVarsAllocated_(false), condest_(azOO.condest_)
{
  AllocAzArrays();
  SetAztecDefaults();

  if (SetAztecVariables() != 0) {
     cerr << __FILE__ << ", line "<<__LINE__<<
         ": error in SetAztecVariables. aborting." << endl;
     abort();
  }
}

//=============================================================================
Aztec_OO::~Aztec_OO(void)  
{
  DeleteMemory();
  DeleteAzArrays();
}

//=============================================================================
void Aztec_OO::DeleteMemory()
{
  if (!azVarsAllocated_) return;

  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_); 
    Prec_ = 0;
  }

  if (Amat_ != 0) {
    AZ_matrix_destroy(&Amat_);
    Amat_ = 0;
  }

  delete [] update_;
}

//=============================================================================
int Aztec_OO::SetAztecDefaults() {

 AZ_defaults(options_, params_);
 useAdaptiveDefaults_ = true;
 NumTrials_ = 0;
 maxFill_ = 0;
 maxKspace_ = 0;
 athresholds_ = 0;
 rthresholds_ = 0;
 condestThreshold_ = 0.0;
 return(0);

}

//=============================================================================
int Aztec_OO::SetAztecVariables() {
  DeleteMemory();

  if (X_ == NULL || B_ == NULL || A_ == NULL) PETRA_CHK_ERR(-1);

  N_update_ = B_->MyLength();
  N_local_ = N_update_;
  update_ = new int[N_update_];

  B_->Map().MyGlobalElements(update_);
  X_->ExtractView(&x_, &x_LDA_);
  B_->ExtractView(&b_, &b_LDA_);

  if (update_ == NULL) PETRA_CHK_ERR(-1);

#ifdef AZTEC_MPI
  AZ_set_proc_config(proc_config_, A_->Comm().Comm());
#else
  AZ_set_proc_config(proc_config_, 0);
#endif


  /* initialize Aztec matrix to be solved */

   Amat_ = AZ_matrix_create(N_local_);
   AZ_set_MATFREE(Amat_, (void *) A_, Petra_Aztec_matvec);
   AZ_set_MATFREE_matrix_norm(Amat_, A_->NormInf()); /* Aztec needs upper bound for  */
                                                     /* matrix norm                  */

   int N_ghost = A_->NumMyCols() - A_->NumMyRows();
   AZ_set_MATFREE_getrow(Amat_, A_, Petra_Aztec_getrow,
			 Petra_Aztec_comm_wrapper,N_ghost,proc_config_);

     

   Prec_ = NULL;    /* When preconditioning structure is NULL, AZ_iterate()  */
                   /* applies Aztec's preconditioners to application matrix */
                   /* (i.e. user does not supply a preconditioning routine  */
                   /* or an additional matrix for preconditioning.          */

  azVarsAllocated_ = true;

  return(0);
}

//=============================================================================
int Aztec_OO::AllocAzArrays()
{
  proc_config_ = new int[AZ_PROC_SIZE];

  options_ = new int[AZ_OPTIONS_SIZE];
  params_ = new double[AZ_PARAMS_SIZE];
  status_ = new double[AZ_STATUS_SIZE];

  if (options_ == NULL || params_ == NULL || status_ == NULL) PETRA_CHK_ERR(-1);

  return(0);
}

//=============================================================================
void Aztec_OO::DeleteAzArrays()
{
  delete [] proc_config_;
  delete [] options_;
  delete [] params_;
  delete [] status_;
  if (athresholds_!=0) delete [] athresholds_;
  if (rthresholds_!=0) delete [] rthresholds_;
}

//=============================================================================
int Aztec_OO::SetProblemOptions(ProblemDifficultyLevel PDL,
				bool ProblemSymmetric)
{
  if (ProblemSymmetric)
    {
      SetAztecOption(AZ_solver, AZ_cg);
      switch (PDL)
      {
      case easy:
        SetAztecOption(AZ_poly_ord, 1);
        SetAztecOption(AZ_precond, AZ_Jacobi);
        break;

      case moderate:

        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_icc);
        break;

      case hard:
      case unsure:
      default:
        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_icc);
        SetAztecParam(AZ_omega, 1.2);
      }
    }
  else
    {
      switch (PDL)
      {
      case easy:
        SetAztecOption(AZ_poly_ord, 1);
        SetAztecOption(AZ_precond, AZ_Jacobi);
        SetAztecOption(AZ_solver, AZ_bicgstab);
        break;

      case moderate:

        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_ilu);
        SetAztecOption(AZ_solver, AZ_gmres);
        break;

      case hard:
      case unsure:
      default:
        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_ilut);
        SetAztecOption(AZ_overlap, 1);
        SetAztecParam(AZ_ilut_fill, 3.0);
        SetAztecParam(AZ_drop, 0.01);
        SetAztecOption(AZ_kspace, 1000);
      }
    }

  return(0);
}

//=============================================================================
int Aztec_OO::SetPreconditioner(void  (*prec_function)
					      (double *, int *, int *, double *,
					       struct AZ_MATRIX_STRUCT  *,
					       struct AZ_PREC_STRUCT *),
					      void *p_data)
{
  Prec_ = AZ_precond_create(Amat_,prec_function, p_data);
  options_[AZ_precond] = AZ_user_precond;

  return(0);
}
//
//=============================================================================
int Aztec_OO::ConstructPreconditioner(double & condest) {

  if (A_==0) return(-1); // No matrix yet

  int precond_flag = options_[AZ_precond];

  if (precond_flag) {

  // Create default Aztec preconditioner if one not defined
  if (Prec_==0) Prec_ = AZ_precond_create(Amat_, AZ_precondition, NULL);

  AZ_mk_context(options_, params_, Amat_->data_org, Prec_, proc_config_);


    int NN = A_->NumMyCols();
    double * condvec = new double[NN];
    for (int i = 0 ; i < N_local_ ; i++ ) condvec[i] = 1.0;
    Prec_->prec_function(condvec,options_,proc_config_,params_,Amat_,Prec_);
    condest_ = 0.0;
    for (int i=0; i<N_local_; i++)
      if (fabs(condvec[i]) > condest_)
	condest_ = fabs(condvec[i]);
    delete [] condvec;
    options_[AZ_pre_calc] = AZ_reuse;
    double tmp_condest = condest_;
    A_->Comm().MaxAll(&tmp_condest, &condest_, 1); // Get the worst of all condition estimates

    condest = condest_;
  }
  return(0);
}
//
//=============================================================================
int Aztec_OO::DestroyPreconditioner() {

  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_);
    Prec_ = 0;
    options_[AZ_pre_calc] = AZ_calc;
  }
  return(0);
}
//=============================================================================
int Aztec_OO::SetMatrixName(int label)

{
  Amat_->data_org[AZ_name] = label;

  return(0);
}

//=============================================================================
int Aztec_OO::recursiveIterate(int MaxIters, double Tolerance)
{
  SetAztecOption(AZ_max_iter, MaxIters);
  SetAztecParam(AZ_tol, Tolerance);

  int prec_allocated = 0;
  if (Prec_ == 0) {
    if (options_[AZ_precond] == AZ_user_precond) {
      if (proc_config_[AZ_node] == 0) {
	printf("AZ_iterate: Can not use NULL for precond argument when\n");
	printf("            options[AZ_precond] == AZ_user_precond.\n");
      }
      exit(1);
    }
    Prec_ = AZ_precond_create(Amat_, AZ_precondition, NULL);
    prec_allocated = 1;
  }

  options_[AZ_recursion_level]++;
  AZ_oldsolve(x_, b_, options_, params_, status_, proc_config_,
	      Amat_, Prec_, Scaling_);
  options_[AZ_recursion_level]--;
  if (prec_allocated == 1) {
    AZ_precond_destroy(&Prec_);
    Prec_ = NULL;
    prec_allocated = 0;
  }
          


  // Determine end status

  int ierr = 0;
  if (status_[AZ_why]==AZ_normal) ierr = 0;
  else if (status_[AZ_why]==AZ_param) ierr = -1;
  else if (status_[AZ_why]==AZ_breakdown) ierr = -2;
  else if (status_[AZ_why]==AZ_loss) ierr = -3;
  else if (status_[AZ_why]==AZ_ill_cond) ierr = -4;
  else if (status_[AZ_why]==AZ_maxits) return(1);
  else cerr << "Internal Aztec_OO Error" << endl;
  abort();

  PETRA_CHK_ERR(ierr);
}
//=============================================================================
int Aztec_OO::Iterate(int MaxIters, double Tolerance)
{
  if (X_ == 0 || B_ == 0 || A_ == 0) PETRA_CHK_ERR(-1);

  SetAztecOption(AZ_max_iter, MaxIters);
  SetAztecParam(AZ_tol, Tolerance);

  
  int prec_allocated = 0;
  if (Prec_ == 0) {
    if (options_[AZ_precond] == AZ_user_precond) {
      if (proc_config_[AZ_node] == 0) {
	cerr << "Aztec_OO::Iterate: Can not use NULL for precond argument when\n";
	cerr << "                   options[AZ_precond] == AZ_user_precond.\n";
      }
      PETRA_CHK_ERR(-2);
    }
    Prec_ = AZ_precond_create(Amat_, AZ_precondition, NULL);
    prec_allocated = 1;
  }

  AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
	     Amat_, Prec_, Scaling_);

  if (prec_allocated == 1) {
    AZ_precond_destroy(&Prec_);
    Prec_ = NULL;
    prec_allocated = 0;
  }
  

  // Determine end status

  int ierr = 0;
  if (status_[AZ_why]==AZ_normal) ierr = 0;
  else if (status_[AZ_why]==AZ_param) ierr = -1;
  else if (status_[AZ_why]==AZ_breakdown) ierr = -2;
  else if (status_[AZ_why]==AZ_loss) ierr = -3;
  else if (status_[AZ_why]==AZ_ill_cond) ierr = -4;
  else if (status_[AZ_why]==AZ_maxits) return(1);
  else {
    cerr << "Internal Aztec_OO Error" << endl;
    abort();
  }

  PETRA_CHK_ERR(ierr);
}

//=============================================================================
int Aztec_OO::Iterate(Petra_RDP_RowMatrix * A,
                      Petra_RDP_MultiVector * X,
                      Petra_RDP_MultiVector * B,
                      int MaxIters, double Tolerance)
{
  A_ = A;
  X_ = X;
  B_ = B;
  SetAztecVariables();

  PETRA_CHK_ERR( Iterate(MaxIters, Tolerance) );
}
//=============================================================================
int Aztec_OO::SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds,
				double condestThreshold, double maxFill, int maxKspace) {

  if (athresholds_!=0) delete [] athresholds_;
  if (rthresholds_!=0) delete [] rthresholds_;

  NumTrials_ = NumTrials;
  maxFill_ = maxFill;
  maxKspace_ = maxKspace;
  athresholds_ = new double[NumTrials];
  rthresholds_ = new double[NumTrials];
  for (int i=0; i<NumTrials; i++) {
    athresholds_[i] = athresholds[i];
    rthresholds_[i] = rthresholds[i];
  }
  if (condestThreshold>0) condestThreshold_ = condestThreshold;
  useAdaptiveDefaults_ = false;
  return(0);
}
//=============================================================================
int Aztec_OO::AdaptiveIterate(int MaxIters, double Tolerance) {

  // Check if adaptive strategy is appropriate (only works for domain decomp and if subdomain
  // solve is not AZ_lu.  If not, call standard Iterate() method

  if (options_[AZ_precond] != AZ_dom_decomp) PETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
  if (options_[AZ_subdomain_solve] == AZ_lu) PETRA_CHK_ERR(Iterate(MaxIters,Tolerance));


  SetAztecOption(AZ_max_iter, MaxIters);
  SetAztecParam(AZ_tol, Tolerance);

  // Make sure we are using IFPACK BILU
  if (options_[AZ_subdomain_solve] == AZ_bilu) options_[AZ_subdomain_solve] = AZ_bilu_ifp;

  // Construct adaptive strategy if necessary
  if (useAdaptiveDefaults_) {

    if (options_[AZ_subdomain_solve] == AZ_bilu_ifp) {
      int NumTrials = 3;
      double athresholds[] = {0.0, 1.0E-14, 1.0E-3};
      double rthresholds[] = {0.0, 1.0E-14, 1.0E-3};
      double condestThreshold = 1.0E16;
      double maxFill = 4.0;
      int maxKspace = 4*options_[AZ_kspace];
      SetAdaptiveParams(NumTrials, athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
    }
    else {
      int NumTrials = 7;
      double athresholds[] = {0.0, 1.0E-12, 1.0E-12, 1.0E-5, 1.0E-5, 1.0E-2, 1.0E-2};
      double rthresholds[] = {1.0, 1.0,     1.01,    1.0,    1.01,   1.01,   1.1   };
      double condestThreshold = 1.0E16;
      double maxFill = 4.0;
      int maxKspace = 4*options_[AZ_kspace];
      SetAdaptiveParams(NumTrials, athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
    }
  }

  // If no trials were defined, just call regular Iterate method

  if (NumTrials_<1) PETRA_CHK_ERR(Iterate(MaxIters,Tolerance));


  bool FirstCallToSolver = true;
  
  // ********************************************************************************
  //  Phase:  Tweak fill level and drop tolerances
  // ********************************************************************************
  
  
  double fill;
  if (options_[AZ_subdomain_solve] == AZ_ilut) fill = params_[AZ_ilut_fill];
  else fill = (double) options_[AZ_graph_fill];

  double curMaxFill = maxfn(fill, maxFill_);
  int curMaxKspace = maxfn(options_[AZ_kspace], maxKspace_);
  if (options_[AZ_solver]!=AZ_gmres) curMaxKspace = options_[AZ_kspace]; // GMRES is only solver sensitive to kspace
  int kspace =options_[AZ_kspace];
 
  while ((status_[AZ_why]!=AZ_normal || FirstCallToSolver) && kspace<=curMaxKspace) {
    SetAztecOption(AZ_kspace,kspace);

    if (options_[AZ_subdomain_solve] == AZ_ilut) params_[AZ_ilut_fill] = fill;
    else options_[AZ_graph_fill] = (int) fill;

    // ********************************************************************************
    //  Phase:  Find a preconditioner whose condest is below the condest threshold
    // ********************************************************************************

    // Start with first trial
    int curTrial = 0;
    // Initial guess for conditioner number
    condest_ = condestThreshold_; 

    // While not converged or doing at least one trial AND still have trials to execute
    while ((status_[AZ_why]!=AZ_normal || FirstCallToSolver) && curTrial<NumTrials_) {

      // Current condest threshold number (force one while iteration)
      double curCondestThreshold = condest_;
    
      // While condest is too large
  
      while ( condest_>=curCondestThreshold) {
	if (Prec_!=0) DestroyPreconditioner(); // Get rid of any existing preconditioner
      
	SetAztecParam(AZ_athresh, athresholds_[curTrial]); // Set threshold values
	SetAztecParam(AZ_rthresh, rthresholds_[curTrial]);
      
	// Preconstruct preconditioner and get a condition number estimate
	ConstructPreconditioner(condest_);
	curTrial++;
      }

  
      // ********************************************************************************
      // Phase:  Solve using preconditioner from first phase
      // ********************************************************************************
    
   	Petra_RDP_MultiVector Xold(*X_);
	double oldResid = status_[AZ_r];
       
      AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
		 Amat_, Prec_, Scaling_);
      
      if (!FirstCallToSolver) {
	if (!(oldResid>status_[AZ_r])) *X_ = Xold;
      }
      FirstCallToSolver = false;

      if (status_[AZ_why]==AZ_maxits) { // This means we need more robust preconditioner
	if (fill<curMaxFill) {
	  fill = minfn(2*fill, curMaxFill); // double fill
	  curTrial = NumTrials_; // force exit of current trial loop
	}
	else if (options_[AZ_subdomain_solve] == AZ_ilut && params_[AZ_drop]>0) {
	  params_[AZ_drop] = 0.0; // if nonzero drop used with ILUT, try one more time with drop = 0
	  curTrial = NumTrials_; // force exit of current trial loop
	}
	else {
	  kspace *= 2; // double kspace and try again
	  curTrial = 0;
	}
      }
    }

  }
  DestroyPreconditioner(); // Delete preconditioner
  // Determine end status

  int ierr = 0;
  if (status_[AZ_why]==AZ_normal) ierr = 0;
  else if (status_[AZ_why]==AZ_param) ierr = -1;
  else if (status_[AZ_why]==AZ_breakdown) ierr = -2;
  else if (status_[AZ_why]==AZ_loss) ierr = -3;
  else if (status_[AZ_why]==AZ_ill_cond) ierr = -4;
  else if (status_[AZ_why]==AZ_maxits) return(1);
  else {
    cerr << "Internal Aztec_OO Error" << endl;
    abort();
  }

  PETRA_CHK_ERR(ierr);
}



//=============================================================================
void Petra_Aztec_matvec(double *x, double *y, AZ_MATRIX *Amat, int proc_config[]) {

  Petra_RDP_RowMatrix * A = (Petra_RDP_RowMatrix *) AZ_get_matvec_data(Amat);
  Petra_RDP_Vector X(View, A->BlockRowMap(), x);
  Petra_RDP_Vector Y(View, A->BlockRowMap(), y);

  //cout << X << endl;
  //cout << Y << endl;

  A->Multiply(false, X, Y);
  //cout << X << endl;
  //cout << Y << endl;

}

//=============================================================================
int Petra_Aztec_getrow(int columns[], double values[],
     int row_lengths[], AZ_MATRIX *Amat, int N_requested_rows,
		       int requested_rows[], int allocated_space){
/*
 * Supply local matrix (without ghost node columns) for rows given by
 * requested_rows[0 ... N_requested_rows-1].  Return this information in
 * 'row_lengths, columns, values'.  If there is not enough space to complete
 * this operation, return 0. Otherwise, return 1.
 *
 * Parameters
 * ==========
 * Amat             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be
 *                  returned.
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the
 *                  row indices of the rows for which nonzero values are
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the
 *                  row 'requested_rows[i]'
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and
 *                  'values' for storing nonzeros. If more space is needed,
 *                  return 0.
 */
  Petra_RDP_RowMatrix * A = (Petra_RDP_RowMatrix *) AZ_get_matvec_data(Amat);
  int nz_ptr = 0;
  int NumRows = A->NumMyRows();
  int NumEntries;

  double *Values = values;
  int * Indices = columns;
  int Length = allocated_space;

  for (int i = 0; i < N_requested_rows; i++) {
    int LocalRow = requested_rows[i];
    // Do copy, return if copy failed
    if (A->ExtractMyRowCopy(LocalRow, Length, NumEntries, Values, Indices)!=0) return(0);
    // update row lengths and pointers 
    row_lengths[i] = NumEntries;
    Values += NumEntries;
    Indices+= NumEntries;
    Length -= NumEntries;
  }
  /*
  for (int i = 0; i < N_requested_rows; i++) {
     int LocalRow = requested_rows[i];
     A->ExtractMyRowView (LocalRow, NumEntries, Values, Indices);
     row_lengths[i] = NumEntries;
     if (nz_ptr+NumEntries>allocated_space) return(0);
     for (int j=0; j<NumEntries; j++) {
       columns[nz_ptr] = Indices[j];
       values[nz_ptr++] = Values[j];
    }
   }
  */
  return(1);
}
//=============================================================================
int Petra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) {
/*
 * Update vec's ghost node via communication. Note: the length of vec is
 * given by N_local + N_ghost where Amat was created via
 *                 AZ_matrix_create(N_local);
 * and a 'getrow' function was supplied via
 *                 AZ_set_MATFREE_getrow(Amat, , , , N_ghost, );
 *
 * Parameters
 * ==========
 * vec              On input, vec contains data. On output, ghost values
 *                  are updated.
 *
 * Amat             On input, points to user's data containing matrix values.
 *                  and communication information.
 */


  Petra_RDP_RowMatrix * A = (Petra_RDP_RowMatrix *) AZ_get_matvec_data(Amat);


  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode

  //Petra_RDP_Vector & X = *new Petra_RDP_Vector(View, A->BlockImportMap(), vec);
  Petra_RDP_Vector X(View, A->BlockImportMap(), vec);

  assert(X.Import(X, *(A->Importer()),Insert)==0);

  return(1);
}

#ifdef AZTEC_OO_WITH_ML
//=============================================================================
int Aztec_OO::PetraMatrix2MLMatrix(ML *ml_handle, int level, 
						Petra_RDP_RowMatrix * A)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = osize;
  int N_ghost = A->NumMyCols() - A->NumMyRows();

  
  ML_Init_Amatrix(ml_handle, level,isize, osize, (void *) A);
  ML_Set_Amatrix_Getrow(ml_handle, level, Petra_ML_getrow,
			Petra_ML_comm_wrapper,	isize+N_ghost);

  ML_Set_Amatrix_Matvec(ml_handle,  level, Petra_ML_matvec);

  return 1;
}
int Petra_ML_matvec(void *data, int in, double *p, int out, double *ap) {

  Petra_RDP_RowMatrix * A = (Petra_RDP_RowMatrix *) data;
  Petra_RDP_Vector X(View, A->BlockRowMap(), p);
  Petra_RDP_Vector Y(View, A->BlockRowMap(), ap);
  
  A->Multiply(false, X, Y);
  

  return 1;
}

//=============================================================================
int Petra_ML_getrow(void *data, int N_requested_rows, int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[]) {
  /*
 * Supply local matrix (without ghost node columns) for rows given by
 * requested_rows[0 ... N_requested_rows-1].  Return this information in
 * 'row_lengths, columns, values'.  If there is not enough space to complete
 * this operation, return 0. Otherwise, return 1.
 *
 * Parameters
 * ==========
 * Amat             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be
 *                  returned.
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the
 *                  row indices of the rows for which nonzero values are
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the
 *                  row 'requested_rows[i]'
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and
 *                  'values' for storing nonzeros. If more space is needed,
 *                  return 0.
 */
  Petra_RDP_RowMatrix * A = (Petra_RDP_RowMatrix *) data;
  int nz_ptr = 0;
  int NumRows = A->NumMyRows();
  int NumEntries;
  double * Values;
  int * Indices;
  for (int i = 0; i < N_requested_rows; i++) {
    int LocalRow = requested_rows[i];
    A->ExtractMyRowView (LocalRow, NumEntries, Values, Indices);
    row_lengths[i] = NumEntries;
    if (nz_ptr+NumEntries>allocated_space) return(0);
    for (int j=0; j<NumEntries; j++) {
      columns[nz_ptr] = Indices[j];
      values[nz_ptr++] = Values[j];
    }
  }

  return(1);
}
//=============================================================================
int Petra_ML_comm_wrapper(double vec[], void *data) {
  /*
   * Update vec's ghost node via communication. Note: the length of vec is
   * given by N_local + N_ghost where Amat was created via
   *                 AZ_matrix_create(N_local);
   * and a 'getrow' function was supplied via
   *                 AZ_set_MATFREE_getrow(Amat, , , , N_ghost, );
 *
 * Parameters
 * ==========
 * vec              On input, vec contains data. On output, ghost values
 *                  are updated.
 *
 * Amat             On input, points to user's data containing matrix values.
 *                  and communication information.
 */


  Petra_RDP_RowMatrix * A = (Petra_RDP_RowMatrix *) data;


  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode

  Petra_RDP_Vector X(View, A->BlockImportMap(), vec);

  assert(X.Import(X, *(A->Importer()),Insert)==0);

  return(1);
}
#endif // AZTEC_OO_WITH_ML
