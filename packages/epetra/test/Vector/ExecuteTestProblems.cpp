#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "BuildTestProblems.h"
#include "Epetra_Comm.h"
  int MatrixTests(const Epetra_BlockMap & Map, const Epetra_LocalMap & LocalMap, 
		      bool verbose)
  {
    int NumVectors = 1;
    const Epetra_Comm & Comm = Map.Comm();
    int ierr = 0, i;
    int IndexBase = 0;
    double *residual = new double[NumVectors];

    /* get ID of this processor */

    int MyPID   = Comm.MyPID();


    // Test GEMM first.  7 cases:
    
    //                                       Num
    //     OPERATIONS                        case  Notes
    // 1) C(local) = A^X(local) * B^X(local)  4   (X=Trans or Not, No Comm needed) 
    // 2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)
    // 3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no Comm needed)

    // ==================================================================
    // Case 1 through 4 (A, B, C all local) Strided and non-strided cases
    // ==================================================================

    // Construct Vectors

  {
    Epetra_Vector A(LocalMap);
    Epetra_Vector B(LocalMap);
    Epetra_LocalMap  Map2d(NumVectors, IndexBase, Comm);
    Epetra_Vector C(Map2d);
    Epetra_Vector C_GEMM(Map2d);

    double **App, **Bpp, **Cpp;
    
    Epetra_Vector *Ap, *Bp, *Cp;

    // For testing non-strided mode, create Vectors that are scattered throughout memory

    App = new double *[NumVectors];
    Bpp = new double *[NumVectors];
    Cpp = new double *[NumVectors];
    for (i=0; i<NumVectors; i++) App[i] = new double[A.MyLength()+i];
    for (i=0; i<NumVectors; i++) Bpp[i] = new double[B.MyLength()+i];
    for (i=0; i<NumVectors; i++) Cpp[i] = new double[C.MyLength()+i];
    
    Epetra_Vector A1(View, LocalMap, App[0]);
    Epetra_Vector B1(View, LocalMap, Bpp[0]);
    Epetra_Vector C1(View, Map2d, Cpp[0]);

    for (int strided = 0; strided<2; strided++){
    int ierr;
    // Loop through all trans cases using a variety of values for alpha and beta
    for (i=0; i<4; i++){
	ierr = 0;
	char transa = 'N'; if (i>1) transa = 'T';
	char transb = 'N'; if (i%2!=0) transb = 'T';
	double alpha = (double) i+1;
	double beta  = (double) (i/2);
	EPETRA_TEST_ERR(C.Random(),ierr);  // Fill C with random numbers
	ierr += BuildMatrixTests(C,transa, transb, alpha, A, B, beta, C_GEMM );
	if (strided)
	  {
	    Ap = &A; Bp = &B; Cp = &C;
	  }
	else
	  {
	    A.ExtractCopy(App[0]); Ap = &A1;
	    B.ExtractCopy(Bpp[0]); Bp = &B1;
	    C.ExtractCopy(Cpp[0]); Cp = &C1;
	  }
	  
	ierr += Cp->Multiply(transa, transb, alpha, *Ap, *Bp, beta);
	ierr += Cp->Update(-1.0, C_GEMM, 1.0);
	ierr += Cp->Norm2(residual);

	if (verbose && ierr==0)
	  {
	    cout << "XXXXX Replicated Local Vector GEMM tests";
	    if (strided)
	    cout << " (Strided Multivectors)" << endl;
	    else
	    cout << " (Non-Strided Multivectors)" << endl;
	    cout << "  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
		 <<", transb = " << transb;
	  }
	if (ierr==0 && BadResidual(verbose,residual)) return(-1);
      }

      }
    for (i=0; i<NumVectors; i++)
      {
	delete [] App[i];
	delete [] Bpp[i];
	delete [] Cpp[i];
      }
    delete [] App;
    delete [] Bpp;
    delete [] Cpp;
  }
      
    // ====================================
    // Case 5  (A, B distributed C  local)
    // ====================================

    // Construct Vectors
  {
    Epetra_Vector A(Map);
    Epetra_Vector B(Map);
    Epetra_LocalMap Map2d(NumVectors, IndexBase, Comm);
    Epetra_Vector C(Map2d);
    Epetra_Vector C_GEMM(Map2d);

    char transa = 'T';
    char transb = 'N';
    double alpha = 2.0;
    double beta  = 1.0;
    EPETRA_TEST_ERR(C.Random(),ierr);  // Fill C with random numbers
    ierr += BuildMatrixTests(C, transa, transb, alpha, A, B, beta, C_GEMM );
    ierr += C.Multiply(transa, transb, alpha, A, B, beta);
    ierr += C.Update(-1.0, C_GEMM, 1.0);
    ierr += C.Norm2(residual);

    if (verbose && ierr==0)
      {
	cout << "XXXXX Generalized 2D dot product via GEMM call     " << endl;
	cout << "  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
	     <<", transb = " << transb;
      }
    if (BadResidual(verbose,residual)) return(-1);
    
    
  }      
    // ====================================
    // Case 6-7  (A, C distributed, B local)
    // ====================================

    // Construct Vectors
  {
    Epetra_Vector A(Map);
    Epetra_LocalMap Map2d(NumVectors, IndexBase, Comm);
    Epetra_Vector B(Map2d);
    Epetra_Vector C(Map);
    Epetra_Vector C_GEMM(Map);

    for (i=0; i<2; i++)
      {
	char transa = 'N';
	char transb = 'N'; if (i>0) transb = 'T';
	double alpha = 2.0;
	double beta  = 1.1;
	EPETRA_TEST_ERR(C.Random(),ierr);  // Fill C with random numbers
	ierr += BuildMatrixTests(C,transa, transb, alpha, A, B, beta, C_GEMM );
	ierr += C.Multiply(transa, transb, alpha, A, B, beta);
	ierr += C.Update(-1.0, C_GEMM, 1.0);
	ierr += C.Norm2(residual);
	
	if (verbose)
	  {
	    cout << "XXXXX Generalized 2D vector update via GEMM call     " << endl;
	    cout << "  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
		 <<", transb = " << transb;
	  }
	if (BadResidual(verbose,residual)) return(-1);
      }

    delete [] residual;
    
    return(ierr);
  }
  }

int VectorTests(const Epetra_BlockMap & Map, bool verbose)
{
  int NumVectors = 1;
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0, i;
  double *residual = new double[NumVectors];
  
  Epetra_BLAS BLAS;
  /* get number of processors and the name of this processor */
  
  // int NumProc = Comm.getNumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct Vectors
  
  Epetra_Vector A(Map);
  Epetra_Vector sqrtA(Map);
  Epetra_Vector B(Map);
  Epetra_Vector C(Map);
  Epetra_Vector C_alphaA(Map);
  Epetra_Vector C_alphaAplusB(Map);
  Epetra_Vector C_plusB(Map);
  Epetra_Vector Weights(Map);
  
  // Construct double vectors
  double *dotvec_AB   = new double[NumVectors];
  double *norm1_A     = new double[NumVectors];
  double *norm2_sqrtA = new double[NumVectors];
  double *norminf_A = new double[NumVectors];
  double *normw_A = new double[NumVectors];
  double *minval_A = new double[NumVectors];
  double *maxval_A = new double[NumVectors];
  double *meanval_A = new double[NumVectors];
  
  // Generate data 

  
  EPETRA_TEST_ERR(C.Random(),ierr); // Fill C with random numbers.
  double alpha = 2.0;
  BuildVectorTests (C,alpha, A, sqrtA, B, C_alphaA, C_alphaAplusB,
			     C_plusB, dotvec_AB, norm1_A, norm2_sqrtA, norminf_A, 
			     normw_A, Weights, minval_A, maxval_A, meanval_A);

  int err = 0;
  if (verbose) cout << "XXXXX Testing alpha * A     ";
  // Test alpha*A
  Epetra_Vector alphaA(A); // Copy of A
  EPETRA_TEST_ERR(alphaA.Scale(alpha),err);
  EPETRA_TEST_ERR(alphaA.Update(-1.0, C_alphaA, 1.0),err);
  EPETRA_TEST_ERR(alphaA.Norm2(residual),err);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }

  err = 0;
  if (verbose) cout << "XXXXX Testing C = alpha * A + B      ";
  // Test alpha*A + B
  Epetra_Vector alphaAplusB(A); // Copy of A
  EPETRA_TEST_ERR(alphaAplusB.Update(1.0, B, alpha, A, 0.0),err);
  EPETRA_TEST_ERR(alphaAplusB.Update(-1.0, C_alphaAplusB, 1.0),err);
  EPETRA_TEST_ERR(alphaAplusB.Norm2(residual),err);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }

  err = 0;
  if (verbose) cout << "XXXXX Testing C += B      ";
  // Test + B
  Epetra_Vector plusB(C); // Copy of C
  EPETRA_TEST_ERR(plusB.Update(1.0, B, 1.0),err);
  EPETRA_TEST_ERR(plusB.Update(-1.0, C_plusB, 1.0),err);
  EPETRA_TEST_ERR(plusB.Norm2(residual),err);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }

  err = 0;
  if (verbose) cout << "XXXXX Testing A.dotProd(B)     ";
  // Test A.dotvec(B)
  double *dotvec = residual;
  EPETRA_TEST_ERR(A.Dot(B,dotvec),err);
  BLAS.AXPY(NumVectors,-1.0,dotvec_AB,dotvec);
  
  if (err) ierr += err;
  else {
  EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
  
  err = 0;
  if (verbose) cout << "XXXXX Testing norm1_A      ";
  // Test A.norm1()
  double *norm1 = residual;
  EPETRA_TEST_ERR(A.Norm1(norm1),err);
  BLAS.AXPY(NumVectors,-1.0,norm1_A,norm1);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
	  
  err = 0;
  if (verbose) cout << "XXXXX Testing norm2_sqrtA     ";
  // Test sqrtA.norm2()
  double *norm2 = residual;
  EPETRA_TEST_ERR(sqrtA.Norm2(norm2),err);
  BLAS.AXPY(NumVectors,-1.0,norm2_sqrtA,norm2);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
	
  err = 0;
  if (verbose) cout << "XXXXX Testing norminf_A     ";
  // Test A.norminf()
  double *norminf = residual;
  EPETRA_TEST_ERR(A.NormInf(norminf),ierr);
  BLAS.AXPY(NumVectors,-1.0,norminf_A,norminf);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
	
  err = 0;
  if (verbose) cout << "XXXXX Testing normw_A     ";
  // Test A.NormWeighted()
  double *normw = residual;
  EPETRA_TEST_ERR(A.NormWeighted(Weights, normw),err);
  BLAS.AXPY(NumVectors,-1.0,normw_A,normw);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
	
  err = 0;
  if (verbose) cout << "XXXXX Testing minval_A     ";
  // Test A.MinValue()
  double *minval = residual;
  EPETRA_TEST_ERR(A.MinValue(minval),err);
  BLAS.AXPY(NumVectors,-1.0,minval_A,minval);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
	
  err = 0;
  if (verbose) cout << "XXXXX Testing maxval_A     ";
  // Test A.MaxValue()
  double *maxval = residual;
  EPETRA_TEST_ERR(A.MaxValue(maxval),err);
  BLAS.AXPY(NumVectors,-1.0,maxval_A,maxval);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }

  err = 0;
  if (verbose) cout << "XXXXX Testing meanval_A     ";
  // Test A.MeanValue()
  double *meanval = residual;
  EPETRA_TEST_ERR(A.MeanValue(meanval),err);
  BLAS.AXPY(NumVectors,-1.0,meanval_A,meanval);
  
  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }
	
  err = 0;
  if (verbose) cout << "XXXXX Testing abs_A     ";
  // Test A.Abs()
  Epetra_Vector Abs_A = A;
  EPETRA_TEST_ERR(Abs_A.Abs(A),err);
  EPETRA_TEST_ERR(Abs_A.Update(1.0, A, -1.0),err); // Abs_A = A - Abs_A (should be zero since A > 0)
  EPETRA_TEST_ERR(Abs_A.Norm2(residual),err);

  if (err) ierr += err;
  else {
    EPETRA_TEST_ERR(BadResidual(verbose,residual),ierr);
  }	
  
  // Delete everything
  
  delete [] dotvec_AB;
  delete [] norm1_A;
  delete [] norm2_sqrtA;
  delete [] norminf_A;
  delete [] normw_A;
  delete [] minval_A;
  delete [] maxval_A;
  delete [] meanval_A;
  delete [] residual;
  
  return(ierr);
}

int BadResidual(bool verbose, double * Residual)
{
  double threshold = 5.0E-6;
  int ierr = 0;
    if (Residual[0]>threshold) {
      ierr = 1;// Output will be more useful after returning from method
      if (verbose) cout << endl << "     Residual = " << Residual[0];
    }
  if (verbose)
    if (ierr==0) cout << "\t Checked OK" << endl;
  
  return(ierr);
}
