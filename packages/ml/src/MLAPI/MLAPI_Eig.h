#ifndef MLAPI_EIG_H
#define MLAPI_EIG_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include "MLAPI_Error.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"

namespace MLAPI {

//! Computes the maximum eigenvalue of \c Op using the A-norm of the operator.
double MaxEigAnorm(const Operator& Op, const bool DiagonalScaling = false) 
{
  return(ML_Operator_MaxNorm(Op.GetML_Operator(), DiagonalScaling));
}

//! Computes the maximum eigenvalue of \c Op using the CG method.
double MaxEigCG(const Operator& Op, const bool DiagonalScaling = false) 
{

  ML_Krylov *kdata;
  double MaxEigen;

  kdata = ML_Krylov_Create(GetML_Comm());
  if (DiagonalScaling == false)
    kdata->ML_dont_scale_by_diag = ML_TRUE;
  else
    kdata->ML_dont_scale_by_diag = ML_FALSE;
  ML_Krylov_Set_PrintFreq(kdata, 0);
  ML_Krylov_Set_ComputeEigenvalues( kdata );
  ML_Krylov_Set_Amatrix(kdata, Op.GetML_Operator());
  ML_Krylov_Solve(kdata, Op.GetML_Operator()->outvec_leng, NULL, NULL);
  MaxEigen = ML_Krylov_Get_MaxEigenvalue(kdata);
  ML_Krylov_Destroy(&kdata);

  return(MaxEigen);
}

//! Computes the maximum eigenvalue of \c Op using the power method.
double MaxEigPowerMethod(const Operator& Op, const bool DiagonalScaling = false) 
{

  ML_Krylov *kdata;
  double MaxEigen;

    kdata = ML_Krylov_Create(GetML_Comm());
    if (DiagonalScaling == false)
      kdata->ML_dont_scale_by_diag = ML_TRUE;
    else
      kdata->ML_dont_scale_by_diag = ML_FALSE;
    ML_Krylov_Set_PrintFreq(kdata, 0);
    ML_Krylov_Set_ComputeNonSymEigenvalues(kdata);
    ML_Krylov_Set_Amatrix(kdata, Op.GetML_Operator());
    ML_Krylov_Solve(kdata, Op.GetML_Operator()->outvec_leng, NULL, NULL);
    MaxEigen = ML_Krylov_Get_MaxEigenvalue(kdata);
  ML_Krylov_Destroy(&kdata);

  return(MaxEigen);
}

//! Computes the maximum eigenvalue of \c Op using Anasazi
double MaxEigAnasazi(const Operator& Op, const bool DiagonalScaling = false) 
{

  double MaxEigen = 0.0;

  bool DiagScal;
  if (DiagonalScaling)
    DiagScal = ML_TRUE;
  else
    DiagScal = ML_FALSE;

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAZI) && defined(HAVE_ML_TEUCHOS)
  ML_Anasazi_Get_SpectralNorm_Anasazi(Op.GetML_Operator(), 0, 10, 1e-5,
                                      ML_FALSE, DiagScal, &MaxEigen);
#else
  ML_THROW("Configure w/ --enable-epetra --enable-anasazi --enable-teuchos",
           -1);
#endif
  return(MaxEigen);
}

//! Computes eigenvalues and eigenvectors using LAPACK (w/ one process only).
void Eig(const Operator& Op, MultiVector& ER, 
         MultiVector& EI, MultiVector& V)
{
  if (GetNumProcs() != 1)
    ML_THROW("Function Eig() works only w/ one process, use Eigs() instead.", -1);

  int ierr;
  if (Op.GetDomainSpace() != Op.GetRangeSpace())
    ML_THROW("Matrix is not square", -1);

  V.Reshape(Op.GetDomainSpace(), Op.GetDomainSpace().GetNumGlobalElements());
  ER.Reshape(Op.GetDomainSpace());
  EI.Reshape(Op.GetDomainSpace());

  ierr = ML_Operator_Eigensolver_Dense(Op.GetML_Operator(), ER.GetValues(), 
                                       EI.GetValues(), V.GetValues());
  if (ierr)
    ML_THROW("Error occurred", -1);
  
}

#include "ml_anasazi.h"
#include "float.h"
#include <fstream>
// FIXME: Add List
void Eigs(const Operator& A, int NumEigenvalues, 
          MultiVector& ER, MultiVector& EI)
{

  if (A.GetDomainSpace() != A.GetRangeSpace())
    ML_THROW("Input Operator is not square", -1);

  double time;

  time = GetClock();

  int length = NumEigenvalues;
  double tol = 1e-3;
  int restarts = 1;
  int output = 10;
  bool PrintStatus = true;

  // 1.- set parameters for Anasazi
  Teuchos::ParameterList AnasaziList;
  // MatVec should be either "A" or "ML^{-1}A"
  AnasaziList.set("eigen-analysis: matrix operation", "A");
  AnasaziList.set("eigen-analysis: use diagonal scaling", false);
  AnasaziList.set("eigen-analysis: symmetric problem", false);
  AnasaziList.set("eigen-analysis: length", length);
  AnasaziList.set("eigen-analysis: block-size",1);
  AnasaziList.set("eigen-analysis: tolerance", tol);
  AnasaziList.set("eigen-analysis: restart", restarts);
  AnasaziList.set("eigen-analysis: output", output);
  AnasaziList.get("eigen-analysis: print current status",PrintStatus);

  // data to hold real and imag for eigenvalues and eigenvectors
  Space ESpace(-1, NumEigenvalues);
  ER.Reshape(ESpace);
  EI.Reshape(ESpace);

  // this is the starting value -- random
  Epetra_MultiVector EigenVectors(A.GetRowMatrix()->OperatorDomainMap(),
                                  NumEigenvalues);
  EigenVectors.Random();
#ifdef HAVE_ML_ANASAZI
  int NumRealEigenvectors, NumImagEigenvectors;
#endif

  AnasaziList.set("eigen-analysis: action", "LM");

#ifdef HAVE_ML_ANASAZI
  ML_Anasazi::Interface(A.GetRowMatrix(),EigenVectors,ER.GetValues(),
			EI.GetValues(), AnasaziList, 0, 0,
			&NumRealEigenvectors, &NumImagEigenvectors, 0);
#else
  ML_THROW("Configure ML with --enable-anasazi to use Eigs()", -1);
#endif

  cout << EI;
  cout << ER;
  exit(0);

  return;
}

} // namespace MLAPI

#endif // HAVE_ML_MLAPI

#endif // MLAPI_EIG_H
