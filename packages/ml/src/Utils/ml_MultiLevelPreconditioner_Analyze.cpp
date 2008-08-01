/*!
 *  \file ml_MultiLevelPreconditioner_Analyze.cpp
 *
 *  \brief Visualization utilities for MultiLevelPreconditioner class
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update to Doxygen: 09-Aug-04
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_operator.h"
#include "ml_op_utils.h"
#include "ml_utils.h"
#include "ml_RowMatrix.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include <iomanip>
#ifdef HAVE_ML_IFPACK
#include "Ifpack_Utils.h"
#endif

using namespace std;

// ============================================================================
void ML_Epetra::MultiLevelPreconditioner::VectorNorms(double* vector, 
						      int size, 
						      double* L_inf, 
						      double* L_2)
{
 
  double* Linf = new double[NumPDEEqns_];
  double* L2   = new double[NumPDEEqns_];

  for (int i = 0 ; i < NumPDEEqns_ ; ++i) {
    Linf[i] = 0.0;
    L2[i] = 0.0;
  }

  for (int i = 0 ; i < size ; ++i) {
    // Linf norm 
    if (fabs(vector[i]) > Linf[i % NumPDEEqns_]) 
      Linf[i % NumPDEEqns_] = fabs(vector[i]);
    // L2 norm
    L2[i % NumPDEEqns_] += vector[i] * vector[i];
  }

  Comm().SumAll(Linf,L_inf,NumPDEEqns_);
  Comm().SumAll(L2,L_2,NumPDEEqns_);

  for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
    L_inf[eq] = sqrt(L_inf[eq]);
    L_2[eq] = sqrt(L_2[eq]);
  }

  delete [] Linf;
  delete [] L2;

  return;
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
AnalyzeHierarchy(const bool AnalyzeMatrices, 
                 const int PreCycles, const int PostCycles,
                 const int MLCycles)
{

  // sanity checks

  if (RowMatrix_ == 0) 
    ML_CHK_ERR(-1); // Matrix not yet set
  if (ml_ == 0) 
    ML_CHK_ERR(-4); // at present does not work with Maxwell (easy fix?)

  // execution begins

  double time;

  time = GetClock();

  if (Comm().MyPID() == 0) {
    cout << endl;
    ML_print_line("-",80);
  }

  if (AnalyzeMatrices) {
    for (int i = 0 ; i < NumLevels_ ; ++i) {
      char name[80];
      sprintf(name,"A matrix level %d", LevelID_[i]);
#ifdef HAVE_ML_IFPACK
      ML_Epetra::RowMatrix wrapper(&(ml_->Amat[LevelID_[i]]), &(Comm()));
      if (verbose_) cout << endl;
      wrapper.SetLabel(name);
      Ifpack_Analyze(wrapper,true,ml_->Amat[LevelID_[i]].num_PDEs);
      if (verbose_) cout << endl;
#else
      cerr << ErrorMsg_ << "AnalyzeHierarchy(true) requires IFPACK" << endl;
      cerr << ErrorMsg_ << "Please configure ML with the (at least) the following:" << endl;
      cerr << ErrorMsg_ << "--enable-epetra" << endl;
      cerr << ErrorMsg_ << "--enable-teuchos" << endl;
      cerr << ErrorMsg_ << "--enable-ifpack" << endl;
      ML_CHK_ERR(-1);
#endif
    }
  }

  if (Comm().MyPID() == 0) {
    cout << endl;
    cout << "Solving Ae = 0, with a random initial guess" << endl;
    cout << "- number of pre-smoother cycle(s)  = " << PreCycles << endl;
    cout << "- number of post-smoother cycle(s) = " << PostCycles << endl;
    cout << "- number of ML cycle(s)            = " << MLCycles << endl;
    cout << "- all reported data are scaled with respect to the corresponding" << endl
         << "  value before the application of the solver" << endl;
    cout << "  (0 == perfect solution, 1 == no effect)" << endl;
    cout << endl;
    cout.width(40); cout.setf(ios::left); 
    cout << "Solver";
    cout.width(10); cout.setf(ios::left); 
    cout << "  Linf";
    cout.width(10); cout.setf(ios::left); 
    cout << "   L2";
    cout << endl;
    cout.width(40); cout.setf(ios::left); 
    cout << "------";
    cout.width(10); cout.setf(ios::left); 
    cout << "  ----";
    cout.width(10); cout.setf(ios::left); 
    cout << "   --";
    cout << endl;
    cout << endl;
  }

  if (PreCycles > 0 || PostCycles > 0)
    AnalyzeSmoothers(PreCycles, PostCycles);

  AnalyzeCoarse();

  if (MLCycles > 0)
    AnalyzeCycle(MLCycles);

  if (Comm().MyPID() == 0) {
    cout << endl;
    cout << "*** Total time for analysis = " 
         << GetClock() - time << " (s)" << endl;
    ML_print_line("-",80);
    cout << endl;
  }

  return(0);

}

// ============================================================================
// date: Aug-17
int ML_Epetra::
MultiLevelPreconditioner::AnalyzeSmoothers(const int NumPreCycles,
                                           const int NumPostCycles)
{

  // sanity checks

  if (IsPreconditionerComputed() == false) 
    ML_CHK_ERR(-1); // need preconditioner to do this job

  if( ml_ == 0 ) {
    ML_CHK_ERR(-2); // Does not work with Maxwell
  }

  // =============================================================== //
  // Cycle over all levels.                                          //
  // =============================================================== //

  for (int ilevel = 0 ; ilevel<NumLevels_ -1 ; ++ilevel) {

    // ============ //
    // pre-smoother //
    // ============ //

    ML_Smoother* ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).pre_smoother);

    vector<double> before_Linf(NumPDEEqns_);
    vector<double> before_L2(NumPDEEqns_);
    vector<double> after_Linf(NumPDEEqns_);
    vector<double> after_L2(NumPDEEqns_);

    int n = ml_->Amat[LevelID_[ilevel]].outvec_leng;
    vector<double> tmp_rhs(n);
    vector<double> tmp_sol(n);

    if (ptr != NULL) {

      RandomAndZero(&tmp_sol[0], &tmp_rhs[0],
                    ml_->Amat[LevelID_[ilevel]].outvec_leng);
      VectorNorms(&tmp_sol[0], n, &before_Linf[0], &before_L2[0]);

      // increase the number of applications of the smoother
      // and run the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumPreCycles;
      ML_Smoother_Apply(ptr, n, &tmp_sol[0], n, &tmp_rhs[0], ML_NONZERO);
      ptr->ntimes = old_ntimes;

      VectorNorms(&tmp_sol[0], n, &after_Linf[0], &after_L2[0]);

      if (Comm().MyPID() == 0) {
	for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
	  cout << "Presmoother  (level " << LevelID_[ilevel] 
	       << ", eq " << eq << ")\t\t";
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_Linf[eq] / before_Linf[eq];
	  cout << ' ';
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_L2[eq] / before_L2[eq] << endl;
	}
        cout << endl;
      }
    }

    // ============= //
    // post-smoother //
    // ============= //

    ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).post_smoother);
    if (ptr != NULL) {

      // random solution and 0 rhs
      RandomAndZero(&tmp_sol[0], &tmp_rhs[0], n);

      VectorNorms(&tmp_sol[0], n, &before_Linf[0], &before_L2[0]);

      // increase the number of applications of the smoother
      // and run the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumPostCycles;
      ML_Smoother_Apply(ptr, n, &tmp_sol[0], n, &tmp_rhs[0], ML_ZERO);
      ptr->ntimes = old_ntimes;

      VectorNorms(&tmp_sol[0], n, &after_Linf[0], &after_L2[0]);

      if (Comm().MyPID() == 0) {
	for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
	  cout << "Postsmoother (level " << LevelID_[ilevel] 
	       << ", eq " << eq << ")\t\t";
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_Linf[eq] / before_Linf[eq];
	  cout << ' ';
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_L2[eq] / before_L2[eq] << endl;;
	}
        cout << endl;
      }
    } 
  } // for( ilevel )

  if (Comm().MyPID() == 0) cout << endl;

  return(0);

}

// ============================================================================
// NOTE: requires Amesos
int ML_Epetra::
MultiLevelPreconditioner::AnalyzeCoarse()
{

  // sanity checks

  if (IsPreconditionerComputed() == false) 
    ML_CHK_ERR(-1); // need preconditioner to do this job

  if( ml_ == 0 ) {
    ML_CHK_ERR(-2); // Does not work with Maxwell
  }

  // execution begins

  vector<double> before_Linf(NumPDEEqns_);
  vector<double> before_L2(NumPDEEqns_);
  vector<double> after_Linf(NumPDEEqns_);
  vector<double> after_L2(NumPDEEqns_);

  int level = ml_->ML_coarsest_level;

  int n = ml_->Amat[level].outvec_leng;

  vector<double> tmp_rhs(n);
  vector<double> tmp_sol(n);

  ML_Smoother* ptr;

  ptr = ((ml_->SingleLevel[level]).post_smoother);

  if (ptr != NULL) {

    RandomAndZero(&tmp_sol[0], &tmp_rhs[0], ml_->Amat[level].outvec_leng);
    VectorNorms(&tmp_sol[0], n, &before_Linf[0], &before_L2[0]);

    ML_Smoother_Apply(ptr, n, &tmp_sol[0], n, &tmp_rhs[0], ML_NONZERO);

    VectorNorms(&tmp_sol[0], n, &after_Linf[0], &after_L2[0]);

    if (Comm().MyPID() == 0) {
      for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
        cout << "Coarse Solver (level " << level
          << ", eq " << eq << ")\t\t";
        cout.width(10); cout.setf(ios::left); 
        cout << after_Linf[eq] / before_Linf[eq];
        cout << ' ';
        cout.width(10); cout.setf(ios::left); 
        cout << after_L2[eq] / before_L2[eq] << endl;
      }
    }
  }

  if (Comm().MyPID() == 0) cout << endl;

  return(0);

}

// ============================================================================
// run ML cycle on a random vector
// date: Aug-17
int ML_Epetra::MultiLevelPreconditioner::AnalyzeCycle(const int NumCycles)
{

  // sanity checks

  if (IsPreconditionerComputed() == false) 
    ML_CHK_ERR(-1); // need preconditioner to do this job

  if( ml_ == 0 ) {
    ML_CHK_ERR(-2); // Does not work with Maxwell (yet)
  }

  // execution begins

  double* before_Linf = new double[NumPDEEqns_];
  double* before_L2   = new double[NumPDEEqns_];
  double* after_Linf  = new double[NumPDEEqns_];
  double* after_L2    = new double[NumPDEEqns_];

  assert(NumMyRows() == ml_->Amat[LevelID_[0]].outvec_leng);
  int Nghost = RowMatrix_->RowMatrixColMap().NumMyElements() - NumMyRows();
  if (Nghost < 0) Nghost = 0; 

  double * tmp_rhs = new double[NumMyRows()]; 
  double * tmp_sol = new double[NumMyRows() + Nghost]; 

  // random solution and zero rhs
  RandomAndZero(tmp_sol, tmp_rhs,NumMyRows());

  VectorNorms(tmp_sol, NumMyRows(), before_Linf, before_L2);

  // run the cycle
  for (int i=0 ; i < NumCycles ; ++i) 
    ML_Cycle_MG(&(ml_->SingleLevel[ml_->ML_finest_level]),
		tmp_sol, tmp_rhs,
		ML_NONZERO, ml_->comm, ML_NO_RES_NORM, ml_);

  VectorNorms(tmp_sol, NumMyRows(), after_Linf, after_L2);

  if (Comm().MyPID() == 0) {
    for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
      cout << "complete ML cycle (eq" << eq << ")\t\t\t";
      cout.width(10); cout.setf(ios::left); 
      cout << after_Linf[eq] / before_Linf[eq];
      cout << ' ';
      cout.width(10); cout.setf(ios::left); 
      cout << after_L2[eq] / before_L2[eq] << endl;
    }
  }

  delete [] before_Linf;
  delete [] after_Linf;
  delete [] before_L2;
  delete [] after_L2;

  delete [] tmp_sol;
  delete [] tmp_rhs;

  return(0);
}

// ============================================================================
#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
