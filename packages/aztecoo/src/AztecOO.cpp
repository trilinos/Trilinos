
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "AztecOO.h"
#ifdef AZTEC_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"



//=============================================================================
AztecOO::AztecOO(Epetra_Operator * A, 
                   Epetra_MultiVector * X,
                   Epetra_MultiVector * B) {
  AllocAzArrays();
  SetAztecDefaults();

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while

  Epetra_RowMatrix * UserMatrix = dynamic_cast<Epetra_RowMatrix *>(A); // Try to cast operator to matrix 
  if (UserMatrix!=0) 
    SetUserMatrix(UserMatrix);
  else
    SetUserOperator(A);

  SetLHS(X);
  SetRHS(B);
  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO(Epetra_RowMatrix * A, 
		 Epetra_MultiVector * X,
		 Epetra_MultiVector * B) {
  AllocAzArrays();
  SetAztecDefaults();

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  SetUserMatrix(A);
  
  SetLHS(X);
  SetRHS(B);
  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO(const Epetra_LinearProblem& prob) {
  AllocAzArrays();
  SetAztecDefaults();

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  // Try to cast operator to matrix 
  Epetra_RowMatrix * UserMatrix = dynamic_cast<Epetra_RowMatrix *>(prob.GetOperator());
  if (UserMatrix!=0) 
    SetUserMatrix(UserMatrix);
  else
    SetUserOperator(prob.GetOperator());
  
  SetLHS(prob.GetLHS());
  SetRHS(prob.GetRHS());

  SetProblemOptions(prob.GetPDL(), prob.IsOperatorSymmetric());
  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO() {
  AllocAzArrays();
  SetAztecDefaults();
}

//=============================================================================
AztecOO::AztecOO(const AztecOO& source) {
  AllocAzArrays();
  SetAztecDefaults();

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  SetUserMatrix(source.GetUserMatrix());
  SetUserOperator(source.GetUserOperator());
  SetPrecMatrix(source.GetPrecMatrix());   // Assume user want to base preconditioner on this matrix
  SetPrecOperator(source.GetPrecOperator());   // Assume user want to base preconditioner on this matrix
  SetLHS(source.GetLHS());
  SetRHS(source.GetRHS());
  SetAllAztecOptions(source.options_);
  SetAllAztecParams(source.params_);
  inConstructor_ = false;
}

//=============================================================================
AztecOO::~AztecOO(void) {
  DeleteMemory();
  DeleteAzArrays();
}

//=============================================================================
int AztecOO::AllocAzArrays()
{
  proc_config_ = new int[AZ_PROC_SIZE];

  options_ = new int[AZ_OPTIONS_SIZE];
  params_ = new double[AZ_PARAMS_SIZE];
  status_ = new double[AZ_STATUS_SIZE];

  if (proc_config_ ==0 ||options_ == 0 || params_ == 0 || status_ == 0) EPETRA_CHK_ERR(-1);

  return(0);
}
//=============================================================================
void AztecOO::DeleteMemory() {
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
int AztecOO::SetAztecDefaults() {

 AZ_defaults(options_, params_);
 UserOperator_ = 0;
 UserMatrix_ = 0;
 PrecOperator_ = 0;
 PrecMatrix_ = 0;
 X_ = 0;
 B_ = 0;
 
 N_update_ = 0;
 N_local_ = 0;
 update_ = 0;
 x_LDA_ = 0;
 x_ = 0;
 b_LDA_ = 0;
 b_ = 0;
 Amat_ = 0;
 Pmat_ = 0;
 Prec_ = 0;
 Scaling_ = 0;
 
 condest_ = (-1.0); 
 useAdaptiveDefaults_ = true;
 NumTrials_ = 0;
 maxFill_ = 0;
 maxKspace_ = 0;
 athresholds_ = 0;
 rthresholds_ = 0;
 condestThreshold_ = 0.0;
 procConfigSet_ = false;
 return(0);

}

//=============================================================================
int AztecOO::SetUserOperator(Epetra_Operator * UserOperator) {

  if (UserOperator == 0 && inConstructor_ == true) return(0);
  if (UserOperator == 0) EPETRA_CHK_ERR(-1);

  if (Amat_ != 0) {
    AZ_matrix_destroy(&Amat_);
    Amat_ = 0;
    delete [] update_; update_ = 0;
  }
  UserOperator_ = UserOperator;

  N_update_ = UserOperator_->RangeMap().NumMyElements();
  N_local_ = N_update_;
  update_ = new int[N_update_];
  if (update_ == 0) EPETRA_CHK_ERR(-1);

  UserOperator_->RangeMap().MyGlobalElements(update_);

  if (update_ == 0) EPETRA_CHK_ERR(-1);


   Amat_ = AZ_matrix_create(N_local_);
   AZ_set_MATFREE(Amat_, (void *) UserOperator_, Epetra_Aztec_operatorvec);

   // Aztec needs upper bound for matrix norm if doing polynomial preconditioning
   if (UserOperator_->HasNormInf())
     AZ_set_MATFREE_matrix_norm(Amat_, UserOperator_->NormInf()); 
  return(0);
}

//=============================================================================
int AztecOO::SetUserMatrix(Epetra_RowMatrix * UserMatrix) {

  if (UserMatrix == 0 && inConstructor_ == true) return(0);
  if (UserMatrix == 0) EPETRA_CHK_ERR(-1);

  UserMatrix_ = UserMatrix;

  SetProcConfig(UserMatrix_->Comm());
  EPETRA_CHK_ERR(SetUserOperator(UserMatrix));
  AZ_set_MATFREE(Amat_, (void *) UserMatrix_, Epetra_Aztec_matvec);
  int N_ghost = UserMatrix_->NumMyCols() - UserMatrix_->NumMyRows();
  AZ_set_MATFREE_getrow(Amat_, (void *) UserMatrix_, Epetra_Aztec_getrow,
			Epetra_Aztec_comm_wrapper,N_ghost,proc_config_);

  // If preconditioner not defined, set up to possibly use native Aztec precons
  if (Prec_==0) EPETRA_CHK_ERR(SetPrecMatrix(UserMatrix));
  return(0);
}



//=============================================================================
int AztecOO::SetPrecMatrix(Epetra_RowMatrix * PrecMatrix) {

  if (PrecMatrix == 0 && inConstructor_ == true) return(0);
  if (PrecMatrix == 0) EPETRA_CHK_ERR(-1);
  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_); 
    Prec_ = 0;
  }
  if (Pmat_ != 0) {
    AZ_matrix_destroy(&Pmat_);
    Pmat_ = 0;
  }

  PrecMatrix_ = PrecMatrix;
  SetProcConfig(PrecMatrix_->Comm());
  Pmat_ = AZ_matrix_create(N_local_);
  AZ_set_MATFREE(Pmat_, (void *) PrecMatrix_, Epetra_Aztec_matvec);
  
  // Aztec needs upper bound for matrix norm if doing polynomial preconditioning
  AZ_set_MATFREE_matrix_norm(Pmat_, PrecMatrix_->NormInf()); 
  int N_ghost = PrecMatrix_->NumMyCols() - PrecMatrix_->NumMyRows();
  AZ_set_MATFREE_getrow(Pmat_, (void *) PrecMatrix_, Epetra_Aztec_getrow,
			Epetra_Aztec_comm_wrapper,N_ghost,proc_config_);
     

   Prec_ = 0;      /* When preconditioning structure is NULL, AZ_iterate()  */
                   /* applies Aztec's preconditioners to the Pmat_ matrix   */
                   /* (i.e. user does not supply a preconditioning routine  */
                   /* or an additional matrix for preconditioning.          */


  return(0);
}


//=============================================================================
int AztecOO::SetPrecOperator(Epetra_Operator * PrecOperator) {

  if (PrecOperator == 0 && inConstructor_ == true) return(0);
  if (PrecOperator == 0) EPETRA_CHK_ERR(-1);

  // Get rid of any other preconditioner
  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_); 
    Prec_ = 0;
  }
  if (Pmat_ != 0) {
    AZ_matrix_destroy(&Pmat_);
    Pmat_ = 0;
  }

  PrecOperator_ = PrecOperator;

  if (Amat_==0) EPETRA_CHK_ERR(-2); // UserOperator must be defined first

   Prec_ = AZ_precond_create(Amat_, Epetra_Aztec_precond, (void *) PrecOperator_);

   options_[AZ_precond] = AZ_user_precond;
   char * label = PrecOperator_->Label();
   if (label==0)
     AZ_set_precond_print_string(Prec_,"User-defined preconditioner");
   else
     AZ_set_precond_print_string(Prec_,label);
   return(0);
}
//=============================================================================
int AztecOO::SetLHS(Epetra_MultiVector * X) {

  if (X == 0 && inConstructor_ == true) return(0);
  if (X == 0) EPETRA_CHK_ERR(-1);
  X_ = X;
  X_->ExtractView(&x_, &x_LDA_);
  return(0);
}
//=============================================================================
int AztecOO::SetRHS(Epetra_MultiVector * B) {

  if (B == 0 && inConstructor_ == true) return(0);
  if (B == 0) EPETRA_CHK_ERR(-1);
  B_ = B;
  B_->ExtractView(&b_, &b_LDA_);

  return(0);
}
//=============================================================================
int AztecOO::SetProcConfig(const Epetra_Comm & Comm) {

  if (!procConfigSet_) {
#ifdef AZTEC_MPI
    const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
    AZ_set_proc_config(proc_config_, comm1.Comm());
#else
    AZ_set_proc_config(proc_config_, 0);
#endif
    procConfigSet_ = true;
  }
  return(0);
}
//=============================================================================
void AztecOO::DeleteAzArrays() {
  if (proc_config_!=0) {delete [] proc_config_; proc_config_ = 0;}
  if (options_!=0) {delete [] options_; options_ = 0;}
  if (params_!=0) {delete [] params_; params_ = 0;}
  if (status_!=0) {delete [] status_; status_ = 0;}
  if (athresholds_!=0) {delete [] athresholds_; athresholds_ = 0;}
  if (rthresholds_!=0) {delete [] rthresholds_; rthresholds_ = 0;}
}

//=============================================================================
int AztecOO::SetProblemOptions(ProblemDifficultyLevel PDL,
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
int AztecOO::SetPreconditioner(void  (*prec_function)
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
int AztecOO::ConstructPreconditioner(double & condest) {

  if (PrecMatrix_==0) return(-1); // No matrix yet

  int precond_flag = options_[AZ_precond];

  if (precond_flag) {

  // Create default Aztec preconditioner if one not defined
  if (Prec_==0) Prec_ = AZ_precond_create(Amat_, AZ_precondition, NULL);

  AZ_mk_context(options_, params_, Amat_->data_org, Prec_, proc_config_);


    int NN = PrecMatrix_->NumMyCols();
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
    // if any processor has a tmp_condest==0.0, then it has a singular preconditioner, check for that first
    PrecMatrix_->Comm().MinAll(&tmp_condest, &condest_, 1); // Get the min of all condition estimates
    if (condest_!=0.0)
      PrecMatrix_->Comm().MaxAll(&tmp_condest, &condest_, 1); // Get the worst of all condition estimates

    condest = condest_;
  }
  return(0);
}
//
//=============================================================================
int AztecOO::DestroyPreconditioner() {

  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_);
    Prec_ = 0;
    options_[AZ_pre_calc] = AZ_calc;
  }
  if (Pmat_ != 0) {
    AZ_matrix_destroy(&Pmat_);
    Pmat_ = 0;
  }
  return(0);
}
//=============================================================================
int AztecOO::SetMatrixName(int label)

{
  Amat_->data_org[AZ_name] = label;

  return(0);
}

//=============================================================================
int AztecOO::recursiveIterate(int MaxIters, double Tolerance)
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
  else throw B_->ReportError("Internal AztecOO Error", -5);

  EPETRA_CHK_ERR(ierr);
  return(0);
}
//=============================================================================
int AztecOO::Iterate(int MaxIters, double Tolerance)
{
  if (X_ == 0 || B_ == 0 || UserOperator_ == 0) EPETRA_CHK_ERR(-1);

  SetAztecOption(AZ_max_iter, MaxIters);
  SetAztecParam(AZ_tol, Tolerance);

  
  int prec_allocated = 0;
  if (Prec_ == 0) {
    if (options_[AZ_precond] == AZ_user_precond) {
      if (proc_config_[AZ_node] == 0) {
	cerr << "AztecOO::Iterate: Can not use NULL for precond argument when\n";
	cerr << "                   options[AZ_precond] == AZ_user_precond.\n";
      }
      EPETRA_CHK_ERR(-2);
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
  else throw B_->ReportError("Internal AztecOO Error", -5);

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int AztecOO::Iterate(Epetra_RowMatrix * A,
                      Epetra_MultiVector * X,
                      Epetra_MultiVector * B,
                      int MaxIters, double Tolerance)
{
  SetUserMatrix(A);
  
  SetLHS(X);
  SetRHS(B);

  EPETRA_CHK_ERR( Iterate(MaxIters, Tolerance) );
  return(0);
}
//=============================================================================
int AztecOO::SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds,
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
int AztecOO::AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double Tolerance) {

  // Check if adaptive strategy is appropriate (only works for domain decomp and if subdomain
  // solve is not AZ_lu.  If not, call standard Iterate() method

  if (options_[AZ_precond] != AZ_dom_decomp) {
    EPETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
    return(0);
  }
  if (options_[AZ_subdomain_solve] == AZ_lu) {
    EPETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
    return(0);
  }

  int NumSolveAttempts = 0;
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

  if (NumTrials_<1) {
    EPETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
    return(0);
  }


  bool FirstCallToSolver = true;
  double oldResid;

  // Keep copy of best old solution
  Epetra_MultiVector Xold(*X_);
  
  // ********************************************************************************
  //  Phase:  Tweak fill level and drop tolerances
  // ********************************************************************************
  
  
  double fill;
  if (options_[AZ_subdomain_solve] == AZ_ilut) fill = params_[AZ_ilut_fill];
  else fill = (double) options_[AZ_graph_fill];

  double curMaxFill = EPETRA_MAX(fill, maxFill_);
  int curMaxKspace = EPETRA_MAX(options_[AZ_kspace], maxKspace_);

  // GMRES is only solver sensitive to kspace
  if (options_[AZ_solver]!=AZ_gmres) curMaxKspace = options_[AZ_kspace]; 
  int kspace =options_[AZ_kspace];

  int maxLossTries = 6;
  int numLossTries = 0;
 

  // *****************************************************************************
  // Master "while" loop 
  // *****************************************************************************

  while ((FirstCallToSolver || status_[AZ_why]!=AZ_normal) && NumSolveAttempts<MaxSolveAttempts) {

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
    
    // While condest is too large or iterative solver says the preconditioned subspace is too
    // ill-conditioned, try to form a better-conditioned preconditioner
    
    while ( curTrial<NumTrials_ && (condest_>=condestThreshold_ || condest_==0.0)) {
      if (Prec_!=0) DestroyPreconditioner(); // Get rid of any existing preconditioner
      
      SetAztecParam(AZ_athresh, athresholds_[curTrial]); // Set threshold values
      SetAztecParam(AZ_rthresh, rthresholds_[curTrial]);
      
      // Preconstruct preconditioner and get a condition number estimate
      ConstructPreconditioner(condest_);
      curTrial++;
    }

    // ********************************************************************************
    // Try solver now, but test if solver says  preconditioned subspace is too
    // ill-conditioned, if so try to form a better-conditioned preconditioner
    // ********************************************************************************

    if (!FirstCallToSolver) {
      // keep old residual and solution to see if progress is made (will keep old solution if not)
      oldResid = status_[AZ_r];
      Xold = *X_;
      // Adjust the tolerance to account for any progress made in the solution
      // The test for "oldResid>0.0" insures that NaN will be ignored
      if (options_[AZ_conv]==AZ_r0 && oldResid>0.0) {
	Tolerance = Tolerance/oldResid;
	SetAztecParam(AZ_tol, Tolerance);
      }
    }
    AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
	       Amat_, Prec_, Scaling_);

    if (!FirstCallToSolver) {
      if (!(oldResid>status_[AZ_r])) *X_ = Xold;  // keep old solution
    }

    NumSolveAttempts++;
    FirstCallToSolver = false;
    
    while ( curTrial<NumTrials_ && 
	    (status_[AZ_why]==AZ_ill_cond || status_[AZ_why]==AZ_breakdown)) {
      
      if (Prec_!=0) DestroyPreconditioner(); // Get rid of any existing preconditioner
      
      SetAztecParam(AZ_athresh, athresholds_[curTrial]); // Set next threshold values
      SetAztecParam(AZ_rthresh, rthresholds_[curTrial]);
      
      // Preconstruct preconditioner and get a condition number estimate
      ConstructPreconditioner(condest_);
      curTrial++;

      // keep old residual and solution to see if progress is made (will keep old solution if not)
      oldResid = status_[AZ_r];
      Xold = *X_;
      // If a relative residual is being used as a stopping criterion,
      // adjust the tolerance to account for any progress made in the solution
      // The test for "status_[AZ_scaled_r]>0.0" insures that NaN will be ignored
      if (options_[AZ_conv]==AZ_r0 && status_[AZ_scaled_r]>0.0) {
	Tolerance = Tolerance/status_[AZ_scaled_r];
	SetAztecParam(AZ_tol, Tolerance);
      }

      AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
		 Amat_, Prec_, Scaling_);
      NumSolveAttempts++;
      
      if (!(oldResid>status_[AZ_r])) *X_ = Xold;  // keep old solution
    }

    if (status_[AZ_why]==AZ_loss) { // This means we should try to solve one more time only
      if (numLossTries<maxLossTries) {
	numLossTries++;
      }
    }
    else if (status_[AZ_why]==AZ_maxits) { // This means we need more robust preconditioner
      if (fill<curMaxFill) {
	fill = EPETRA_MIN(2*fill, curMaxFill); // double fill
	curTrial = 0; // force new pass through preconditioner setup
      }
      else if (options_[AZ_subdomain_solve] == AZ_ilut && params_[AZ_drop]>0) {
	params_[AZ_drop] = 0.0; // if nonzero drop used with ILUT, try one more time with drop = 0
	curTrial = 0; // force new pass through preconditioner setup
      }
      else if (kspace<curMaxKspace && kspace<MaxIters) {
	kspace = EPETRA_MIN(curMaxKspace,2*kspace); // double kspace and try again
      }
    }
    else break;
  }
  if (Prec_!=0) DestroyPreconditioner(); // Delete preconditioner
  // Determine end status

  int ierr = 0;
  if (status_[AZ_why]==AZ_normal) ierr = 0;
  else if (status_[AZ_why]==AZ_param) ierr = -1;
  else if (status_[AZ_why]==AZ_breakdown) ierr = -2;
  else if (status_[AZ_why]==AZ_loss) ierr = -3;
  else if (status_[AZ_why]==AZ_ill_cond) ierr = -4;
  else if (status_[AZ_why]==AZ_maxits) return(1);
  else throw B_->ReportError("Internal AztecOO Error", -5);

  EPETRA_CHK_ERR(ierr);
  return(0);
}



//=============================================================================
void Epetra_Aztec_matvec(double x[], double y[], AZ_MATRIX *Amat, int proc_config[]) {

  Epetra_RowMatrix * A = (Epetra_RowMatrix *) AZ_get_matvec_data(Amat);
  Epetra_Vector X(View, A->DomainMap(), x);
  Epetra_Vector Y(View, A->RangeMap(), y);

  //cout << X << endl;
  //cout << Y << endl;

  A->Apply(X, Y);
  //cout << X << endl;
  //cout << Y << endl;

}

//=============================================================================
void Epetra_Aztec_operatorvec(double x[], double y[], AZ_MATRIX *Amat, int proc_config[]) {

  Epetra_Operator * A = (Epetra_Operator *) AZ_get_matvec_data(Amat);
  Epetra_Vector X(View, A->DomainMap(), x);
  Epetra_Vector Y(View, A->RangeMap(), y);

  //cout << X << endl;
  //cout << Y << endl;

  A->Apply(X, Y);
  //cout << X << endl;
  //cout << Y << endl;

}

//=============================================================================
void Epetra_Aztec_precond(double x[], int input_options[],
			  int proc_config[], double input_params[], AZ_MATRIX *Amat,
			  AZ_PRECOND *prec) {

  Epetra_Operator * A = (Epetra_Operator *) AZ_get_precond_data(prec);
  Epetra_Vector X(View, A->DomainMap(), x);
  Epetra_Vector Y(View, A->RangeMap(), x);

  //cout << X << endl;
  //cout << Y << endl;

  A->ApplyInverse(X, Y);
  //cout << X << endl;
  //cout << Y << endl;

}

//=============================================================================
int Epetra_Aztec_getrow(int columns[], double values[],
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
  Epetra_RowMatrix * A = (Epetra_RowMatrix *) AZ_get_matvec_data(Amat);
  int nz_ptr = 0;
  int NumRows = A->NumMyRows();
  int NumEntries;

  double *Values = values;
  int * Indices = columns;
  int Length = allocated_space;

  for (int i = 0; i < N_requested_rows; i++) {
    int LocalRow = requested_rows[i];
    // Do copy, return if copy failed
    A->NumMyRowEntries(LocalRow, NumEntries);
    if (NumEntries>Length) return(0);  // This check should eliminate error reports from next line
    if (A->ExtractMyRowCopy(LocalRow, Length, NumEntries, Values, Indices)!=0) return(0);
    // update row lengths and pointers 
    row_lengths[i] = NumEntries;
    Values += NumEntries;
    Indices+= NumEntries;
    Length -= NumEntries;
  }
  return(1);
}
//=============================================================================
int Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) {
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


  Epetra_RowMatrix * A = (Epetra_RowMatrix *) AZ_get_matvec_data(Amat);


  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode

  Epetra_Vector X_target(View, A->Importer()->TargetMap(), vec);
  Epetra_Vector X_source(View, A->Importer()->SourceMap(), vec);

  assert(X_target.Import(X_source, *(A->Importer()),Insert)==0);

  return(1);
}

#ifdef AZTEC_OO_WITH_ML
//=============================================================================
int AztecOO::PetraMatrix2MLMatrix(ML *ml_handle, int level, 
						Epetra_RowMatrix * A)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = osize;
  int N_ghost = A->NumMyCols() - A->NumMyRows();

  
  ML_Init_Amatrix(ml_handle, level,isize, osize, (void *) A);
  ML_Set_Amatrix_Getrow(ml_handle, level, Epetra_ML_getrow,
			Epetra_ML_comm_wrapper,	isize+N_ghost);

  ML_Set_Amatrix_Matvec(ml_handle,  level, Epetra_ML_matvec);

  return 1;
}
int Epetra_ML_matvec(void *data, int in, double *p, int out, double *ap) {

  Epetra_RowMatrix * A = (Epetra_RowMatrix *) data;
  Epetra_Vector X(View, A->DomainMap(), p);
  Epetra_Vector Y(View, A->RangeMap(), ap);
  
  A->Multiply(false, X, Y);
  

  return 1;
}

//=============================================================================
int Epetra_ML_getrow(void *data, int N_requested_rows, int requested_rows[], 
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
  Epetra_RowMatrix * A = (Epetra_RowMatrix *) data;
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
int Epetra_ML_comm_wrapper(double vec[], void *data) {
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


  Epetra_RowMatrix * A = (Epetra_RowMatrix *) data;


  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode

  Epetra_Vector X(View, A->ImportMap(), vec);

  assert(X.Import(X, *(A->Importer()),Insert)==0);

  return(1);
}
#endif // AZTEC_OO_WITH_ML
