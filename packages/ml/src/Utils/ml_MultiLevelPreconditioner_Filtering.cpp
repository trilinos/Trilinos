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

/* utilities for filtering (or GGB).
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 19-Jan-05
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "ml_amesos_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_agg_ParMETIS.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_anasazi.h"

using namespace Teuchos;
using namespace ML_Epetra;

#ifdef HAVE_ML_AZTECOO
#include "AztecOO.h"
#endif

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::SetFiltering() 
{

  Epetra_Time Time(Comm());

  if (List_.get("filtering: enable",false) == false) 
    return(0);

  int restarts = List_.get("eigen-analysis: restart", 50);
  int NumEigenvalues = List_.get("filtering: eigenvalues to compute", 5);
  int length = List_.get("eigen-analysis: length", NumEigenvalues);
  int BlockSize = List_.get("eigen-analysis: block-size", 1);
  double tol = List_.get("eigen-analysis: tolerance", 1e-5);

  if (length <= NumEigenvalues) 
    length = NumEigenvalues+1;

  if (verbose_) {
    cout << endl;
    cout << PrintMsg_ << "\tFiltering the preconditioner: computing low-convergent modes..." << endl;
    cout << PrintMsg_ << "\t- number of eigenvectors to compute = " << NumEigenvalues << endl;
    cout << PrintMsg_ << "\t- tolerance = " << tol << endl;
    cout << PrintMsg_ << "\t- block size = " << BlockSize << endl;
    cout << PrintMsg_ << "\t- length     = " << length << endl;
    cout << PrintMsg_ << "\t  (Note that the resulting precondition is non-symmetric)" << endl;
  }    

  // set parameters for Anasazi
  Teuchos::ParameterList AnasaziList;
  AnasaziList.set("eigen-analysis: matrix operation", "I-ML^{-1}A");
  AnasaziList.set("eigen-analysis: use diagonal scaling", false);
  AnasaziList.set("eigen-analysis: symmetric problem", false);
  AnasaziList.set("eigen-analysis: length", length);
  AnasaziList.set("eigen-analysis: block-size", BlockSize);
  AnasaziList.set("eigen-analysis: tolerance", tol);
  AnasaziList.set("eigen-analysis: action", "LM");
  AnasaziList.set("eigen-analysis: restart", restarts);
  AnasaziList.set("eigen-analysis: output", 0);

  // data to hold real and imag for eigenvalues and eigenvectors
  vector<double> RealEigenvalues(NumEigenvalues);
  vector<double> ImagEigenvalues(NumEigenvalues);

  vector<double> RealEigenvectors(NumEigenvalues * NumMyRows());
  vector<double> ImagEigenvectors(NumEigenvalues * NumMyRows());

  // this is the starting value -- random
  Epetra_MultiVector EigenVectors(OperatorDomainMap(),NumEigenvalues);
  EigenVectors.Random();

  int NumRealEigenvectors = 0, NumImagEigenvectors = 0;

  // FIXME: this is not a genius stroke, I should simply get
  // the Schur decomposition from Anasazi
#ifdef HAVE_ML_ANASAxI
  // call Anasazi and store the results in eigenvectors      
  ML_Anasazi::Interface(RowMatrix_,EigenVectors,&RealEigenvalues[0],
                        &ImagEigenvalues[0], AnasaziList, &RealEigenvectors[0],
                        &ImagEigenvectors[0],
                        &NumRealEigenvectors, &NumImagEigenvectors, ml_);
#else
  /*
  if (Comm().MyPID() == 0) {
    cerr << ErrorMsg_ << "ML has been configure without the Anasazi interface" << endl
      << ErrorMsg_ << "You must add the option --enable-anasazi to use" << endl
      << ErrorMsg_ << "filtering and Anasazi" << endl;
  }
  */
  if (Comm().MyPID() == 0) {
    cerr << ErrorMsg_ << "The Anasazi support is no longer maintained in ML." << endl;
    cerr << ErrorMsg_ << "Please check with the developers for more details." << endl;
  }

  ML_EXIT(EXIT_FAILURE);
#endif

  if (NumRealEigenvectors + NumImagEigenvectors == 0) {
    cerr << ErrorMsg_ << "Anasazi has computed no nonzero eigenvalues" << endl
      << ErrorMsg_ << "This sometimes happen because your fine grid matrix" << endl
      << ErrorMsg_ << "is too small. In this case, try to change the Anasazi input" << endl
      << ErrorMsg_ << "parameters, or drop the filtering correction." << endl;
    ML_EXIT(EXIT_FAILURE);
  }

  // some output, to print out how many vectors we are actually using
  if (verbose_) {
    cout << PrintMsg_ << "\t- Computed largest magnitude eigenvalues of I - ML^{-1}A:" << endl;
    for( int i=0 ; i<NumEigenvalues ; ++i ) {
      cout << PrintMsg_ << "\t  z = " << std::setw(10) << RealEigenvalues[i]
        << " + i(" << std::setw(10) << ImagEigenvalues[i] << " ),  |z| = "
        << sqrt(RealEigenvalues[i]*RealEigenvalues[i] + ImagEigenvalues[i]*ImagEigenvalues[i]) << endl;
    }
    cout << PrintMsg_ << "\t- Using " << NumRealEigenvectors << " real and "
      << NumImagEigenvectors << " imaginary eigenvector(s)" << endl;
  }

  int size = NumRealEigenvectors+NumImagEigenvectors;

  assert (size < 2 * NumEigenvalues + 1);

  // this is equivalent to the "fattening" of Haim. I build a new ML
  // hierarchy, using the null space just computed. I have at least one
  // aggregate per subdomain, using METIS.

  // copy the null space in a double vector, as now I have real
  // and imag null space in two different vectors (that may be only
  // partially populated)

  flt_NullSpace_.resize((NumRealEigenvectors + NumImagEigenvectors) * NumMyRows());

  int count = 0;
  for (int i = 0 ; i < NumRealEigenvectors ; ++i)
    for (int j = 0 ; j < NumMyRows() ; ++j) 
      flt_NullSpace_[count++] = RealEigenvectors[j + i * NumMyRows()];

  for (int i = 0 ; i < NumImagEigenvectors ; ++i)
    for (int j = 0 ; j < NumMyRows() ; ++j) 
      flt_NullSpace_[count++] = ImagEigenvectors[j + i * NumMyRows()];

  int NumAggr = List_.get("filtering: local aggregates",1);

  // create a new ML hierarchy, whose null space has been computed
  // with the iteration matrix modes

  if (verbose_) 
    cout << endl << PrintMsg_
         << "Building ML hierarchy for filtering" << endl << endl;

  ML_Create(&flt_ml_,2); // now only 2-level methods
  ML_Operator_halfClone_Init(&(flt_ml_->Amat[1]),
                             &(ml_->Amat[ml_->ML_finest_level]));

  ML_Aggregate_Create(&flt_agg_);
  ML_Aggregate_Set_CoarsenScheme_METIS(flt_agg_);
  ML_Aggregate_Set_LocalNumber(flt_ml_,flt_agg_,-1,NumAggr);
  ML_Aggregate_Set_DampingFactor(flt_agg_,0.0);
  ML_Aggregate_Set_Threshold(flt_agg_, 0.0);
  ML_Aggregate_Set_MaxCoarseSize(flt_agg_, 1);
  ML_Aggregate_Set_NullSpace(flt_agg_, NumPDEEqns_,
                             NumRealEigenvectors + NumImagEigenvectors,
                             &flt_NullSpace_[0], NumMyRows());
  int CL = ML_Gen_MultiLevelHierarchy_UsingAggregation(flt_ml_,1, 
                                                       ML_DECREASING, 
                                                       flt_agg_);

  assert (CL == 2);

  ML_Gen_Smoother_Amesos(flt_ml_, 0, ML_AMESOS_KLU, -1, 0.0);
  ML_Gen_Solver(flt_ml_, ML_MGV, 1, 0);

  if (verbose_) 
    cout << endl;

  if (verbose_) 
    cout << PrintMsg_ << "\t- Total Time for filtering setup = " 
         << Time.ElapsedTime() << " (s)" << endl;

  return(0);

}

// ============================================================================
// FIXME: into another file
bool ML_Epetra::MultiLevelPreconditioner::CheckPreconditionerKrylov() 
{
#ifdef HAVE_ML_AZTECOO

  Epetra_Time Time(Comm());

  if (verbose_) 
    cout << PrintMsg_ << endl << "\tComputing the rate of convergence..." << endl;

  int MaxIters = List_.get("reuse: max iters",(int)5);
  double Ratio = List_.get("reuse: ratio",(double)0.5);
  int Output = List_.get("reuse: output",-1);
  
  Epetra_Vector LHS(Map());
  Epetra_Vector RHS(Map());

  LHS.PutScalar(0.0);
  RHS.Random();

  Epetra_LinearProblem Problem(const_cast<Epetra_RowMatrix*>(RowMatrix_),&LHS, &RHS);
  
  AztecOO solver(Problem);

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_kspace, MaxIters);
  solver.SetAztecOption(AZ_conv, AZ_r0);
  if (Output == -1) solver.SetAztecOption(AZ_output, AZ_none);
  else              solver.SetAztecOption(AZ_none, Output);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(this);

  solver.Iterate(MaxIters, 1e-15);
  
  double status[AZ_STATUS_SIZE];
  solver.GetAllAztecStatus(status);

  double NewRateOfConvergence = status[AZ_scaled_r];
#ifdef LATER
  if( (status[AZ_why] != AZ_normal ||
       status[AZ_why] != AZ_maxits ) && verbose_ ) {
    cerr << endl;
    cerr << ErrorMsg_ << "\tConvergence tests did not converge normally:" << endl;
         << ErrorMsg_ << "\tstatus[AZ_why] = " << status[AZ_why] << endl;
    cerr << endl;
  }
#endif 

  // print out the computed rate of convergence.
  if (RateOfConvergence_ == -1.0) {

    // This is the first time we are computing this
    RateOfConvergence_ = NewRateOfConvergence;
    if( verbose_ ) {
      cout << PrintMsg_ << "\tRate of convergence : current = " 
	   << RateOfConvergence_ << endl;
      cout << PrintMsg_ << "\tTime to check convergence rate = " 
	   << Time.ElapsedTime() << " (s)" << endl;
    }
    return(true);

  } else {
    // this is NOT the first time this function is called; we have to compare the
    // current rate of convergence with the previous one
    if (verbose_) {
      cout << PrintMsg_ << "\tRate of convergence : previous = " << RateOfConvergence_ << endl;
      cout << PrintMsg_ << "\tRate of convergence : current  = " << NewRateOfConvergence << endl;
    }
    // if the current rate of convergence is not too different from the previous
    // one, the already computed preconditioner is good enough, and we return true
    // (that is, keep what we have). Otherwise, return false (that is, recompute
    // all the stuff).
    
    bool rv;
    if (Ratio * NewRateOfConvergence >= RateOfConvergence_) rv = false;
    else                                                    rv = true;

    if (rv == true && verbose_ ) 
      cout << PrintMsg_ << endl << "\tTest passed: keep old preconditioner" << endl;
    if (rv == false && verbose_) 
      cout << PrintMsg_ << endl << "\tTest failed: now recompute the preconditioner" << endl;
    cout << PrintMsg_ << "\tTime to check convergence rate = " 
      << Time.ElapsedTime() << " (s)" << endl;

    RateOfConvergence_ = NewRateOfConvergence;
    return(rv);

  }
  
#else
  cerr << ErrorMsg_ << "reuse preconditioner requires ML to be configured with" << endl
       << ErrorMsg_ << "--enable-aztecoo." << endl;
  ML_EXIT(EXIT_FAILURE);
  return(false);
#endif
  
}

#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS*/
