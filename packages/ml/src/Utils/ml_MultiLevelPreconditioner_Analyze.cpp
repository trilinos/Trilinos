/*!
 *  \file ml_MultiLevelPreconditioner_Analyze.cpp
 *
 *  \brief Visualization utilities for MultiLevelPreconditioner class
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 09-Aug-04
 *
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_operator.h"
#include "ml_op_utils.h"
#include "ml_utils.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_epetra_preconditioner.h"
#include <iomanip>

// ============================================================================
// vector must hold space for external nodes too
void ML_Epetra::MultiLevelPreconditioner::SmoothnessFactor(ML_Operator* Op,
							   double* vector,
							   double* S_F)
{
  double* SF = new double[NumPDEEqns_];
  for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq)
    SF[eq] = 0.0;
 
  int allocated = 10;

  int* colInd = new int[allocated];
  double* colVal = new double[allocated];
  int NumNonzeros;

  // update external nodes
  ML_exchange_bdry(vector,Op->getrow->pre_comm,
		   Op->invec_leng,Op->comm,ML_OVERWRITE,NULL);

  // cycle over all rows
  for (int i = 0 ; i < Op->invec_leng ; ++i) {

    // get nonzero elements
    int ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NumNonzeros);

    if( ierr == 0 ) {
      do {
	delete [] colInd;
	delete [] colVal;
	
	allocated *= 2;
	colInd = new int[allocated];
	colVal = new double[allocated];

	ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NumNonzeros);
      } while( ierr == 0 );
    }

    // compute the smoothness factor for this row
    for (int j = 0 ; j < NumNonzeros ; ++j) {
      if (colInd[j] % NumPDEEqns_ == i % NumPDEEqns_) {
	double val = (vector[i] - vector[colInd[j]]) / vector[i];
	SF[i % NumPDEEqns_] += val * val;
      }
    }
    
  } // for each row of Op

  delete [] colInd;
  delete [] colVal;
   
  Comm().SumAll(SF,S_F,NumPDEEqns_);

  for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq)
    S_F[eq] = sqrt(S_F[eq]);

  return;
}

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
    if (abs(vector[i]) > Linf[i % NumPDEEqns_]) 
      Linf[i % NumPDEEqns_] = abs(vector[i]);
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
int ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixCheap()
{

  // sanity checks

  if (RowMatrix_ == 0) 
    return(-1); // Matrix not yet set
  if (ml_ == 0) 
    return(-4); // at present does not work with Maxwell (easy fix?)

  // execution begins

  double time;

  time = GetClock();

  if( Comm().MyPID() == 0 ) {
    cout << endl;
    ML_print_line("-",78);
  }

  // Analyze the Amat for each level.
  // I do not use on Epetra matrices because I would like to analyze all
  // the levels, and not only the finest one.
  for (int i = 0 ; i < NumLevels_ ; ++i) {
    char name[80];
    sprintf(name,"A matrix level %d", LevelID_[i]);
    ML_Operator_Analyze(&(ml_->Amat[LevelID_[i]]), name);
  }

  if (Comm().MyPID() == 0) {
    cout << endl;
    cout << "*** Total time for analysis = " 
         << GetClock() - time << " (s)" << endl;
    ML_print_line("-",78);
    cout << endl;
  }

  return(0);

}

// ============================================================================
#include "ml_anasazi.h"
#include "float.h"
#include <fstream>
int ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixEigenvaluesSparse(char* MatVec,
									bool IsSymmetric)
{

  // sanity checks

  if (RowMatrix_ == 0) 
    return(-1); // Matrix not yet set
  if ( (MatVec[0] != 'A') && (IsPreconditionerComputed() != true) )
    return(-2); // need preconditioner for this job
  if (ml_ == 0) 
    return(-4); // at present does not work with Maxwell (easy fix?)

  // execution begins

  double time;

  time = GetClock();

  if( Comm().MyPID() == 0 ) {
    cout << endl;
    ML_print_line("-",78);
    cout << "*** ***************************************************** ***" << endl;
    cout << "*** Analysis of the spectral properties using Anasazi ***" << endl;
    cout << "*** ***************************************************** ***" << endl;
    cout << "*** Operator = " << MatVec << endl;
    cout << endl;
  }

  int def = EPETRA_MIN(NumMyRows()/2,128);
  char filename[80];
  char MATLABname[80];

  if( *MatVec == 'A' ) {
    sprintf(filename,"seig_A.m");
    sprintf(MATLABname,"seig_A");
  }
  else {
    sprintf(filename,"seig_PA.m");
    sprintf(MATLABname,"seig_PA");
  }

  int NumEigenvalues = List_.get("analysis: eigenvalues", def);
  int length = List_.get("analysis: Anasazi length", def);
  double tol = List_.get("analysis: Anasazi tolerance", 1e-2);
  int restarts = List_.get("analysis: Anasazi restart", 2);
  int output = List_.get("analysis: Anasazi output",0);
  bool PrintStatus = List_.get("analysis: print Anasazi status", false);
  

  // 1.- set parameters for Anasazi
  Teuchos::ParameterList AnasaziList;
  // MatVec should be either "A" or "ML^{-1}A"
  AnasaziList.set("eigen-analysis: matrix operation", MatVec);
  AnasaziList.set("eigen-analysis: use diagonal scaling", false);
  AnasaziList.set("eigen-analysis: symmetric problem", false);
  AnasaziList.set("eigen-analysis: length", length);
  AnasaziList.set("eigen-analysis: block-size",1);
  AnasaziList.set("eigen-analysis: tolerance", tol);
  AnasaziList.set("eigen-analysis: restart", restarts);
  AnasaziList.set("eigen-analysis: output", output);
  AnasaziList.get("eigen-analysis: print current status",PrintStatus);

  // data to hold real and imag for eigenvalues and eigenvectors
  double * RealEigenvalues = new double[NumEigenvalues];
  double * ImagEigenvalues = new double[NumEigenvalues];

  // this is the starting value -- random
  Epetra_MultiVector EigenVectors(OperatorDomainMap(),NumEigenvalues);
  EigenVectors.Random();

  int NumRealEigenvectors, NumImagEigenvectors;

  AnasaziList.set("eigen-analysis: action", "LM");

#ifdef HAVE_ML_ANASAZI
  // 2.- call Anasazi and store the results in eigenvectors      
  if (Comm().MyPID() == 0) {
    cout << "\t*** Computing LM eigenvalues of " << MatVec << endl;
    cout << "\t*** using Anasazi. NumEigenvalues = " << NumEigenvalues 
         << ", tolerance = " << tol << endl;
  }
  ML_Anasazi::Interface(RowMatrix_,EigenVectors,RealEigenvalues,
			ImagEigenvalues, AnasaziList, 0, 0,
			&NumRealEigenvectors, &NumImagEigenvectors, ml_);
#else
  if (Comm().MyPID() == 0) {
    cerr << ErrorMsg_ << "ML has been configure without the Anasazi interface" << endl
      << ErrorMsg_ << "You must add the option --enable-anasazi to use" << endl
      << ErrorMsg_ << "filtering and Anasazi" << endl;
  }
  exit(EXIT_FAILURE);
#endif

  // 3.- some output, to print out how many vectors we are actually using

  double max = 0;

  if (verbose_) {

    cout << "\t*** Done. Results are in file `" << filename << "'." << endl;
    std::ofstream seig_A(filename);

    seig_A << "%largest magnitude eigenvalues" << endl;
    seig_A << "%computed eigenvalues = " << NumEigenvalues << endl;
    seig_A << "%tolerance = " << tol << endl;
    seig_A << MATLABname << " = [" << endl;

    for (int i=0 ; i<NumEigenvalues ; ++i) {
      seig_A <<	RealEigenvalues[i] << " + i * (" << ImagEigenvalues[i] << ')' << endl;
      double eig =sqrt(pow(RealEigenvalues[i],2.0) + pow(ImagEigenvalues[i],2.0));
      if( eig>max ) max = eig;
    }
    seig_A.close();

    cout << PrintMsg_ << "\tmax |lambda_{" << MatVec << "}|   = " << max << endl;
  }

  // ========================================= //
  // Now as before, but for smallest magnitude //
  // ========================================= //

  EigenVectors.Random();

  AnasaziList.set("eigen-analysis: action", "SM");

#ifdef HAVE_ML_ANASAZI
  if (Comm().MyPID() == 0) {
    cout << "\t*** Computing SM eigenvalues of " << MatVec << endl;
    cout << "\t*** using Anasazi. NumEigenvalues = " << NumEigenvalues 
         << ", tolerance = " << tol << endl;
  }

  ML_Anasazi::Interface(RowMatrix_,EigenVectors,RealEigenvalues,
			ImagEigenvalues, AnasaziList, 0, 0,
			&NumRealEigenvectors, &NumImagEigenvectors, ml_);
#else
  if (Comm().MyPID() == 0) {
    cerr << ErrorMsg_ << "ML has been configure without the Anasazi interface" << endl
      << ErrorMsg_ << "You must add the option --enable-anasazi to use" << endl
      << ErrorMsg_ << "filtering and Anasazi" << endl;
  }
  exit(EXIT_FAILURE);
#endif

  // 3.- some output, to print out how many vectors we are actually using

  double min = DBL_MAX;

  if (verbose_) {

    cout << "\t*** Done. Results are in file `" << filename << "'." << endl;
    std::ofstream seig_A(filename, std::ios::app);

    seig_A << "%smallest magnitude eigenvalues" << endl;
    seig_A << "%computed eigenvalues = " << NumEigenvalues << endl;
    seig_A << "%tolerance = " << tol << endl;

    for (int i=0 ; i<NumEigenvalues ; ++i) {
      seig_A <<	RealEigenvalues[i] << " + i * (" << ImagEigenvalues[i] << ')' << endl;
      double eig =sqrt(pow(RealEigenvalues[i],2.0) + pow(ImagEigenvalues[i],2.0));
      if (eig < min) min = eig;
    }
    seig_A << "];" << endl;
    seig_A.close();

    cout << PrintMsg_ << "\tmin |lambda_{" << MatVec <<"}|   = " << min << endl;

    cout << endl;
    cout << PrintMsg_ << "\tSpectral condition number of " << MatVec
         << " = " << max/min << endl;
    cout << endl;
  }

  if (Comm().MyPID() == 0) {
    cout << endl;
    cout << "*** Total time for analysis = " 
         << GetClock() - time << " (s)" << endl;
    ML_print_line("-",78);
    cout << endl;
  }

  return(0);
}

// ============================================================================
// This function works only for serial runs.
// I convert the matrix into dense format, than I call
// LAPACK to do the job.
// This has advantages and disadvantages:
// A) eigenvalues are "exact", and I get easily real
//    and imaginary part
// A) I can compute the eigenvalues of ML^{-1} A
// D) only serial runs
// D) only small problems
// ============================================================================
#include "ml_lapack.h"
int ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixEigenvaluesDense(char* MatVec,
									bool IsSymmetric)
{

  // sanity checks

  if (RowMatrix_ == 0) 
    return(-1); // Matrix not yet set
  if ( (MatVec[0] != 'A') && (IsPreconditionerComputed() != true) )
    return(-2); // need preconditioner for this job
  if (Comm().NumProc() != 1)
    return(-3); // only serial runs

  // execution begins
 
  double time;

  time = GetClock();

  if( Comm().MyPID() == 0 ) {
    cout << endl;
    ML_print_line("-",78);
    cout << "*** ************************************************ ***" << endl;
    cout << "*** Analysis of the spectral properties using LAPACK ***" << endl;
    cout << "*** ************************************************ ***" << endl;
    cout << "*** Operator = " << MatVec << endl;
    cout << endl;
  }

  // FIXME: I don't consider at all IsSymmetric (might save some time
  // doing so)
  int n = NumMyRows();
  ML_Operator* Amat = &(ml_->Amat[LevelID_[0]]);
  char jobvl, jobvr;
  int lda,  ldvl, ldvr, lwork, info;
  double* a = 0;
  double* vl = 0;
  double* vr = 0;
  double* work = 0;
  int *ipiv;
  double *Er, *Ei;
  double min = DBL_MAX, max = DBL_MIN;

  jobvl = 'N'; /* V/N to calculate/not calculate left eigenvectors
                  of matrix H.*/

  jobvr = 'N'; /* As above, but for right eigenvectors. */

  lda = n; /* The leading dimension of matrix a. */

#ifdef NOT_DEF
  if( n * n * 8 > 536870912 ) {
    cerr << ErrorMsg_ << "LAPACK analysis of finest matrix would" << endl
         << ErrorMsg_ << "require " << n * n * 8 /1024 << " Kbytes. This seems too" << endl
	 << ErrorMsg_ << "much to me. Now I return; maybe you can change the" << endl
	 << ErrorMsg_ << "source code (file " << __FILE__ << ", line "
	 << __LINE__ << ")" << endl;
    return(-1);
  }
#endif
    
  info = 0;

  ipiv = new int[n];
  ldvl = n;
  ldvr = n;
  work = new double[4*n];
  lwork = 4*n;

  Er = new double[n];
  Ei = new double[n];

  for( int i=0 ; i<n ; i++ ) {
    Er[i] = 0.0, Ei[i] = 0.0;
  }

  a = new double[NumMyRows()*NumMyRows()];
  if( a == 0 ) {
    cerr << ErrorMsg_ << "Not enough memory to allocate"
         << 8*NumMyRows()*NumMyRows() << "bytes. Now" << endl
	 << ErrorMsg_ << "skipping the analysis of the eigenvalues." << endl;
    return(-10);
  }
  
  int allocated = 1;
  int * colInd = new int[allocated];
  double * colVal = new double[allocated];
  int ierr;
  int ncnt;

  if (MatVec[0] == 'A') {

    // set to zero the elements
    for (int i=0; i<NumMyRows(); ++i)
      for( int j=0 ; j<NumMyRows() ; ++j ) 
	a[i+j*NumMyRows()] = 0.0;

    // now insert the nonzero elements, row by row

    for (int i = 0; i < n ; i++) {

      ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);

      if( ierr == 0 ) {
	do {
	  delete [] colInd;
	  delete [] colVal;
	  allocated *= 2;
	  colInd = new int[allocated];
	  colVal = new double[allocated];
	  ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);
	} while( ierr == 0 );
      }
      for (int j = 0 ; j < ncnt ; ++j) {
	a[i+n*colInd[j]] = colVal[j];
      }
    }

    /* ================================================ */
    /* largest and smallest eigenvalue (in module) of A */
    /* ------------------------------------------------ */

    cout << "\t*** Computing eigenvalues of finest-level matrix" << endl;
    cout << "\t*** using LAPACK. This may take some time..." << endl;

    DGEEV_F77(CHAR_MACRO(jobvl), CHAR_MACRO(jobvr), &n, a, &n, Er, Ei, vl,
	      &ldvl, vr, &ldvr, work, &lwork, &info);

    cout << "\t*** results are on file `eig_A.m'." << endl;

    std::ofstream eig_A("eig_A.m");
    eig_A << "eig_A = [" << endl;
    min = DBL_MAX, max = DBL_MIN;

    for (int i=0 ; i<n ; i++) {
      double eig = sqrt(Er[i] * Er[i] + Ei[i] * Ei[i]);
      if( eig < min ) min = eig;
      if( eig > max ) max = eig;
      eig_A << Er[i] << " + i * " << Ei[i] << endl;
    }
    eig_A << "];";
    eig_A.close();

    cout << endl << "\tmin |lambda_i(A)|         = " << min << endl;
    cout << "\tmax |lambda_i(A)|         = " << max << endl;
    if( min != 0.0 ) 
      cout << "\tspectral condition number = " << max/min << endl;

  }
  else {

    // ================================================= //
    // Now the same job on ML^{-1}A                      //
    // 1.- extract one column of A                       //
    // 2.- apply the ML preconditioner to it             //
    // 3.- substitute the old column with the new one    //
    // 4.- compute the eigenvalues of ML^{-1}A           //
    // ================================================= //

    // LAPACK overwrites the matrix, need to getrow() again
    for (int i = 0; i < n ; i++) {

      ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);

      if (ierr == 0) {
	do {
	  delete [] colInd;
	  delete [] colVal;
	  allocated *= 2;
	  colInd = new int[allocated];
	  colVal = new double[allocated];
	  ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);
	} while( ierr == 0 );
      }
      for (int j = 0 ; j < ncnt ; ++j) {
	a[i+n*colInd[j]] = colVal[j];
      }
    }

    Epetra_Vector x(RowMatrixRowMap());
    Epetra_Vector y(RowMatrixRowMap());


    for (int j = 0 ; j < n ; ++j) {
      for (int i = 0 ; i < n ; ++i) {
	x[i] = a[i+j*n];
      }
      // y is put to zero in ApplyInverse()
      // apply the preconditioner of choice
      ApplyInverse(x,y);
      for (int i = 0 ; i < n ; ++i) {
	a[i+j*n] = y[i];
      }
    }

    cout << "\t*** Computing eigenvalues of ML^{-1}A" << endl;
    cout << "\t*** using LAPACK. This may take some time..." << endl;

    DGEEV_F77(CHAR_MACRO(jobvl), CHAR_MACRO(jobvr), &n, a, &n, Er, Ei, vl,
	      &ldvl, vr, &ldvr, work, &lwork, &info);

    cout << "\t*** results are on file `eig_PA.m'." << endl;

    std::ofstream eig_PA("eig_PA.m");
    eig_PA << "eig_PA = [" << endl;
    min = DBL_MAX, max = DBL_MIN;

    for (int i=0 ; i<n ; i++) {
      double eig = sqrt(Er[i] * Er[i] + Ei[i] * Ei[i]);
      if( eig < min ) min = eig;
      if( eig > max ) max = eig;
      eig_PA << Er[i] << " + i * " << Ei[i] << endl;
    }
    eig_PA << "];";
    eig_PA.close();

    cout << endl << "\tmin |lambda_i(ML^{-1}A)|  = " << min << endl;
    cout << "\tmax |lambda_i(ML^{-1}A)|  = " << max << endl;
    if( min != 0.0 ) 
      cout << "\tspectral condition number = " << max/min << endl;

  }
  
  // free memory

  delete [] colInd;
  delete [] colVal;
  delete [] a;
  delete [] ipiv;
  delete [] work;
  delete [] Er;
  delete [] Ei;

  cout << endl << "** Total time = " << GetClock() - time << "(s)" << endl;
  ML_print_line("-",78);
  cout << endl;

  return(0);
}

// ============================================================================
// date: Aug-17
int ML_Epetra::MultiLevelPreconditioner::AnalyzeSmoothersSparse(int NumPreCycles,
								int NumPostCycles)
{

  // sanity checks

  if (IsPreconditionerComputed() == false) 
    return(-1); // need preconditioner to do this job

  if( ml_ == 0 ) {
    return(-2); // Does not work with Maxwell
  }

  // execution begins

  double* before_Linf = new double[NumPDEEqns_];
  double* before_L2   = new double[NumPDEEqns_];
  double* after_Linf  = new double[NumPDEEqns_];
  double* after_L2    = new double[NumPDEEqns_];
  double* before_SF   = new double[NumPDEEqns_];
  double* after_SF    = new double[NumPDEEqns_];

  double * tmp_rhs = new double[NumMyRows()]; 
  double * tmp_sol = new double[NumMyRows()]; 
  
  ML_Smoother* ptr;

  if (Comm().MyPID() == 0) {
    cout << endl;
    cout << "Solving Ae = 0, with a random initial guess" << endl;
    cout << "- number of pre-smoother cycle(s) = " << NumPreCycles << endl;
    cout << "- number of post-smoother cycle(s) = " << NumPostCycles << endl;
    cout << "- all reported data are scaled with their values" << endl
         << "  before the application of the solver" << endl;
    cout << "  (0 == perfect solution, 1 == no effect)" << endl;
    cout << "- SF is the smoothness factor" << endl;
    cout << endl;
    cout.width(40); cout.setf(ios::left); 
    cout << "Solver";
    cout.width(10); cout.setf(ios::left); 
    cout << "Linf";
    cout.width(10); cout.setf(ios::left); 
    cout << " L2";
    cout.width(10); cout.setf(ios::left); 
    cout << "  SF" << endl;
  }

  // =============================================================== //
  // Cycle over all levels.                                          //
  // =============================================================== //

  for (int ilevel=0 ; ilevel<NumLevels_ - 1 ; ++ilevel) {

    // ============ //
    // pre-smoother //
    // ============ //

    ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).pre_smoother);

    if (ptr != NULL) {

      RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[LevelID_[ilevel]].outvec_leng);
      VectorNorms(tmp_sol, NumMyRows(), before_Linf, before_L2);
      SmoothnessFactor(&(ml_->Amat[LevelID_[ilevel]]), tmp_sol, before_SF);

      // increase the number of applications of the smoother
      // and run the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumPreCycles;
      ML_Smoother_Apply(ptr, 
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_sol,
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_rhs, ML_NONZERO);
      ptr->ntimes = old_ntimes;

      VectorNorms(tmp_sol, NumMyRows(), after_Linf, after_L2);
      SmoothnessFactor(&(ml_->Amat[LevelID_[ilevel]]), tmp_sol, after_SF);

      if (Comm().MyPID() == 0) {
	for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
	  cout << "Presmoother  (level " << LevelID_[ilevel] 
	       << ", eq " << eq << ")\t\t";
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_Linf[eq] / before_Linf[eq];
	  cout << ' ';
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_L2[eq] / before_L2[eq];
	  cout << ' ';
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_SF[eq] / before_SF[eq] << endl;
	}
      }
    }

    // ============= //
    // post-smoother //
    // ============= //

    ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).post_smoother);
    if (ptr != NULL) {

      // random solution and 0 rhs
      RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[LevelID_[ilevel]].outvec_leng);

      VectorNorms(tmp_sol, NumMyRows(), before_Linf, before_L2);
      SmoothnessFactor(&(ml_->Amat[LevelID_[ilevel]]), tmp_sol, before_SF);

      // increase the number of applications of the smoother
      // and run the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumPostCycles;
      ML_Smoother_Apply(ptr, 
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_sol,
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_rhs, ML_ZERO);
      ptr->ntimes = old_ntimes;

      VectorNorms(tmp_sol, NumMyRows(), after_Linf, after_L2);
      SmoothnessFactor(&(ml_->Amat[LevelID_[ilevel]]), tmp_sol, after_SF);

      if (Comm().MyPID() == 0) {
	for (int eq = 0 ; eq < NumPDEEqns_ ; ++eq) {
	  cout << "Postsmoother (level " << LevelID_[ilevel] 
	       << ", eq " << eq << ")\t\t";
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_Linf[eq] / before_Linf[eq];
	  cout << ' ';
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_L2[eq] / before_L2[eq];
	  cout << ' ';
	  cout.width(10); cout.setf(ios::left); 
	  cout << after_SF[eq] / before_SF[eq] << endl;
	}
      }
    } 
  } // for( ilevel )

  if (Comm().MyPID() == 0) 
    cout << endl;

  delete [] before_Linf;
  delete [] after_Linf;
  delete [] before_L2;
  delete [] after_L2;
  delete [] before_SF;
  delete [] after_SF;

  delete [] tmp_sol;
  delete [] tmp_rhs;

  return(0);

}

// ============================================================================
// run ML cycle on a random vector
// date: Aug-17
int ML_Epetra::MultiLevelPreconditioner::AnalyzeCycle(int NumCycles)
{

  // sanity checks

  if (IsPreconditionerComputed() == false) 
    return(-1); // need preconditioner to do this job

  if( ml_ == 0 ) {
    return(-2); // Does not work with Maxwell (yet)
  }

  // execuction begins

  double* before_Linf = new double[NumPDEEqns_];
  double* before_L2   = new double[NumPDEEqns_];
  double* after_Linf  = new double[NumPDEEqns_];
  double* after_L2    = new double[NumPDEEqns_];
  double* before_SF   = new double[NumPDEEqns_];
  double* after_SF    = new double[NumPDEEqns_];

  double * tmp_rhs = new double[NumMyRows()]; 
  double * tmp_sol = new double[NumMyRows()]; 

  // random solution and zero rhs
  RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[LevelID_[0]].outvec_leng);

  VectorNorms(tmp_sol, NumMyRows(), before_Linf, before_L2);
  SmoothnessFactor(&(ml_->Amat[LevelID_[0]]), tmp_sol, before_SF);

  // run the cycle
  for (int i=0 ; i < NumCycles ; ++i) 
    ML_Cycle_MG(&(ml_->SingleLevel[ml_->ML_finest_level]),
		tmp_sol, tmp_rhs,
		ML_NONZERO, ml_->comm, ML_NO_RES_NORM, ml_);

  VectorNorms(tmp_sol, NumMyRows(), after_Linf, after_L2);
  SmoothnessFactor(&(ml_->Amat[LevelID_[0]]), tmp_sol, after_SF);

  if (Comm().MyPID() == 0) {
    cout << endl;
    cout << PrintMsg_ << "Solving Ae = 0, with a random initial guess" << endl;
    cout << PrintMsg_ << "using " << NumCycles << " ML cycle(s)." << endl;
    for (int i = 0 ; i < NumPDEEqns_ ; ++i) {
      cout << "- (eq " << i << ") scaled Linf norm after application(s) = " 
	<< after_Linf[i] / before_Linf[i] << endl;
      cout << "- (eq " << i << ") scaled L2 norm after application(s)   = " 
	<< after_L2[i] / before_L2[i] << endl;
      cout << "- (eq " << i << ") scaled smoothness factor              = "
	<< after_SF[i] / before_SF[i] << endl;
      cout << endl;
    }
  }

  delete [] before_Linf;
  delete [] after_Linf;
  delete [] before_L2;
  delete [] after_L2;
  delete [] before_SF;
  delete [] after_SF;

  delete [] tmp_sol;
  delete [] tmp_rhs;

  return(0);
}

// ============================================================================
// date: Aug-17
int ML_Epetra::MultiLevelPreconditioner::AnalyzeSmoothersDense(int NumPreCycles,
							       int NumPostCycles,
							       int MaxSize)
{

  // sanity checks

  if (IsPreconditionerComputed() == false) 
    return(-1); // need preconditioner to do this job

  if( ml_ == 0 ) {
    return(-2); // Does not work with Maxwell (yet)
  }

  if (Comm().NumProc() != 1) 
    return(-3); // only serial runs

  // execution begins

  double before_L2;
  double after_L2;

  double* tmp_rhs = new double[NumMyRows()]; 
  double* tmp_sol = new double[NumMyRows()]; 
  
  ML_Smoother* ptr;

  // =============================================================== //
  // Cycle over all levels.                                          //
  // =============================================================== //

  for (int ilevel=0 ; ilevel<NumLevels_ - 1 ; ++ilevel) {

    ML_Operator* Op;
    int n;
    double* EigenVectors = 0;
    double* RealEigenValues = 0;
    double* ImagEigenValues = 0;
    
    Op = &(ml_->Amat[LevelID_[ilevel]]); 
    n = Op->outvec_leng;

    // you'd better skip big matrices
    if (n > MaxSize) 
      continue;
    
    // only square matrices
    if (n != Op->invec_leng)
      return(-4);

    EigenVectors = new double[n * n];
    RealEigenValues = new double[n];
    ImagEigenValues = new double[n];

    // compute the eigenvalues and eigenvectors for this operator

    int ierr = ML_Operator_Eigensolver_Dense(Op,RealEigenValues,
					     ImagEigenValues,EigenVectors);

    // define the following variable to plot on file
    // the eigenvectors of the operator
#ifdef NOT_DEF
    double * x_coord = List_.get("viz: x-coordinates", (double *)0);
    double * y_coord = List_.get("viz: y-coordinates", (double *)0);
    
    if ((x_coord != 0) && (y_coord != 0) ) {

      assert(x_coord != 0);
      assert(y_coord != 0);

      char filename[80];
      sprintf(filename,"eigenvalues.xyz");
      std::ofstream evec(filename);
      for (int k = 0 ; k < n ; ++k) {
	evec << RealEigenValues[k] << "   ";
	evec << ImagEigenValues[k] << "   ";
	evec << sqrt(pow(RealEigenValues[k],2) + pow(ImagEigenValues[k],2));
	evec << endl;
      }
      evec.close();

      for (int kk = 0 ; kk < n ; ++kk) {
	char filename[80];
	sprintf(filename,"eigenvectors%d.xyz",kk);
	std::ofstream evec(filename);
	for (int k = 0 ; k < n ; ++k) {
	  evec << x_coord[k] << "   ";
	  evec << y_coord[k] << "   ";
	  evec << EigenVectors[k + kk * n];
	  evec << endl;
	} // filename
	evec.close();
      }
    }
#endif 

    assert(ierr == 0);

    char filename[80];
    sprintf(filename,"presmoother-reduction-level%d.xyz",LevelID_[ilevel]);
    std::ofstream pre(filename);

    if (Comm().MyPID() == 0) {
      cout << "Effect of presmoother is reported on file `" << filename
	   << "'" << endl;

    }

    sprintf(filename,"postsmoother-reduction-level%d.xyz",LevelID_[ilevel]);
    std::ofstream post(filename);

    if (Comm().MyPID() == 0) {
      cout << "Effect of postsmoother is reported on file `" << filename
	   << "'" << endl;

    }

    // =========================== //
    // cycle over all eigenvectors //
    // =========================== //

    for (int k = 0 ; k < n ; ++k) {

      // ============ //
      // pre-smoother //
      // ============ //

      ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).pre_smoother);

      if (ptr != NULL) {

	before_L2 = 0.0;
	for (int j = 0 ; j < n ; ++j) {
	  tmp_sol[j] = EigenVectors[j + k * n];
	  tmp_rhs[j] = 0.0;
	  before_L2 += tmp_sol[j] * tmp_sol[j];
	}

	before_L2 = sqrt(before_L2);

	// increase the number of applications of the smoother
	// and run the smoother
	int old_ntimes = ptr->ntimes;
	ptr->ntimes = NumPreCycles;
	ML_Smoother_Apply(ptr, n, tmp_sol, n, tmp_rhs, ML_NONZERO);
	ptr->ntimes = old_ntimes;

	after_L2 = 0.0;
	for (int j = 0 ; j < n ; ++j) {
	  after_L2 += tmp_sol[j] * tmp_sol[j];
	}

	after_L2 = sqrt(after_L2);

	pre << RealEigenValues[k] << "   ";
	pre << ImagEigenValues[k] << "   ";
	pre << after_L2 / before_L2 << endl;

      } // ptr != NULL

      // ============ //
      // post-smoother //
      // ============ //

      ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).post_smoother);

      if (ptr != NULL) {

	before_L2 = 0.0;
	for (int j = 0 ; j < n ; ++j) {
	  tmp_sol[j] = EigenVectors[j + k * n];
	  tmp_rhs[j] = 0.0;
	  before_L2 += tmp_sol[j] * tmp_sol[j];
	}

	before_L2 = sqrt(before_L2);

	// increase the number of applications of the smoother
	// and run the smoother
	int old_ntimes = ptr->ntimes;
	ptr->ntimes = NumPostCycles;
	ML_Smoother_Apply(ptr, n, tmp_sol, n, tmp_rhs, ML_NONZERO);
	ptr->ntimes = old_ntimes;

	after_L2 = 0.0;
	for (int j = 0 ; j < n ; ++j) {
	  after_L2 += tmp_sol[j] * tmp_sol[j];
	}

	after_L2 = sqrt(after_L2);

	post << RealEigenValues[k] << "   ";
	post << ImagEigenValues[k] << "   ";
	post << after_L2 / before_L2 << endl;
      } // ptr != NULL
    } // for (k)

    pre.close();
    post.close();

    delete [] EigenVectors;
    delete [] RealEigenValues;
    delete [] ImagEigenValues;

  } // for (ilevel)

  delete [] tmp_sol;
  delete [] tmp_rhs;

  return(0);

}
// ============================================================================
#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
