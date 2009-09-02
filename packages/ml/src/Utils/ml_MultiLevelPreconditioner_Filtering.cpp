/* utilities for filtering (or GGB).
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 19-Jan-05
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

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
    std::cout << std::endl;
    std::cout << PrintMsg_ << "\tFiltering the preconditioner: computing low-convergent modes..." << std::endl;
    std::cout << PrintMsg_ << "\t- number of eigenvectors to compute = " << NumEigenvalues << std::endl;
    std::cout << PrintMsg_ << "\t- tolerance = " << tol << std::endl;
    std::cout << PrintMsg_ << "\t- block size = " << BlockSize << std::endl;
    std::cout << PrintMsg_ << "\t- length     = " << length << std::endl;
    std::cout << PrintMsg_ << "\t  (Note that the resulting precondition is non-symmetric)" << std::endl;
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
  std::vector<double> RealEigenvalues(NumEigenvalues);
  std::vector<double> ImagEigenvalues(NumEigenvalues);

  std::vector<double> RealEigenvectors(NumEigenvalues * NumMyRows());
  std::vector<double> ImagEigenvectors(NumEigenvalues * NumMyRows());

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
    std::cerr << ErrorMsg_ << "ML has been configure without the Anasazi interface" << std::endl
      << ErrorMsg_ << "You must add the option --enable-anasazi to use" << std::endl
      << ErrorMsg_ << "filtering and Anasazi" << std::endl;
  }
  */
  if (Comm().MyPID() == 0) {
    std::cerr << ErrorMsg_ << "The Anasazi support is no longer maintained in ML." << std::endl;
    std::cerr << ErrorMsg_ << "Please check with the developers for more details." << std::endl;
  }

  ML_EXIT(EXIT_FAILURE);
#endif

  if (NumRealEigenvectors + NumImagEigenvectors == 0) {
    std::cerr << ErrorMsg_ << "Anasazi has computed no nonzero eigenvalues" << std::endl
      << ErrorMsg_ << "This sometimes happen because your fine grid matrix" << std::endl
      << ErrorMsg_ << "is too small. In this case, try to change the Anasazi input" << std::endl
      << ErrorMsg_ << "parameters, or drop the filtering correction." << std::endl;
    ML_EXIT(EXIT_FAILURE);
  }

  // some output, to print out how many vectors we are actually using
  if (verbose_) {
    std::cout << PrintMsg_ << "\t- Computed largest magnitude eigenvalues of I - ML^{-1}A:" << std::endl;
    for( int i=0 ; i<NumEigenvalues ; ++i ) {
      std::cout << PrintMsg_ << "\t  z = " << std::setw(10) << RealEigenvalues[i]
        << " + i(" << std::setw(10) << ImagEigenvalues[i] << " ),  |z| = "
        << std::sqrt(RealEigenvalues[i]*RealEigenvalues[i] + ImagEigenvalues[i]*ImagEigenvalues[i]) << std::endl;
    }
    std::cout << PrintMsg_ << "\t- Using " << NumRealEigenvectors << " real and "
      << NumImagEigenvectors << " imaginary eigenvector(s)" << std::endl;
  }

  int size = NumRealEigenvectors+NumImagEigenvectors;

  assert (size < 2 * NumEigenvalues + 1);

  // this is equivalent to the "fattening" of Haim. I build a new ML
  // hierarchy, using the null space just computed. I have at least one
  // aggregate per subdomain, using METIS.

  // copy the null space in a double std::vector, as now I have real
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
    std::cout << std::endl << PrintMsg_
         << "Building ML hierarchy for filtering" << std::endl << std::endl;

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
    std::cout << std::endl;

  if (verbose_) 
    std::cout << PrintMsg_ << "\t- Total Time for filtering setup = " 
         << Time.ElapsedTime() << " (s)" << std::endl;

  return(0);

}

// ============================================================================
// FIXME: into another file
bool ML_Epetra::MultiLevelPreconditioner::CheckPreconditionerKrylov() 
{
#ifdef HAVE_ML_AZTECOO

  Epetra_Time Time(Comm());

  if (verbose_) 
    std::cout << PrintMsg_ << std::endl << "\tComputing the rate of convergence..." << std::endl;

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
    std::cerr << std::endl;
    std::cerr << ErrorMsg_ << "\tConvergence tests did not converge normally:" << std::endl;
         << ErrorMsg_ << "\tstatus[AZ_why] = " << status[AZ_why] << std::endl;
    std::cerr << std::endl;
  }
#endif 

  // print out the computed rate of convergence.
  if (RateOfConvergence_ == -1.0) {

    // This is the first time we are computing this
    RateOfConvergence_ = NewRateOfConvergence;
    if( verbose_ ) {
      std::cout << PrintMsg_ << "\tRate of convergence : current = " 
	   << RateOfConvergence_ << std::endl;
      std::cout << PrintMsg_ << "\tTime to check convergence rate = " 
	   << Time.ElapsedTime() << " (s)" << std::endl;
    }
    return(true);

  } else {
    // this is NOT the first time this function is called; we have to compare the
    // current rate of convergence with the previous one
    if (verbose_) {
      std::cout << PrintMsg_ << "\tRate of convergence : previous = " << RateOfConvergence_ << std::endl;
      std::cout << PrintMsg_ << "\tRate of convergence : current  = " << NewRateOfConvergence << std::endl;
    }
    // if the current rate of convergence is not too different from the previous
    // one, the already computed preconditioner is good enough, and we return true
    // (that is, keep what we have). Otherwise, return false (that is, recompute
    // all the stuff).
    
    bool rv;
    if (Ratio * NewRateOfConvergence >= RateOfConvergence_) rv = false;
    else                                                    rv = true;

    if (rv == true && verbose_ ) 
      std::cout << PrintMsg_ << std::endl << "\tTest passed: keep old preconditioner" << std::endl;
    if (rv == false && verbose_) 
      std::cout << PrintMsg_ << std::endl << "\tTest failed: now recompute the preconditioner" << std::endl;
    std::cout << PrintMsg_ << "\tTime to check convergence rate = " 
      << Time.ElapsedTime() << " (s)" << std::endl;

    RateOfConvergence_ = NewRateOfConvergence;
    return(rv);

  }
  
#else
  std::cerr << ErrorMsg_ << "reuse preconditioner requires ML to be configured with" << std::endl
       << ErrorMsg_ << "--enable-aztecoo." << std::endl;
  ML_EXIT(EXIT_FAILURE);
  return(false);
#endif
  
}

#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS*/
