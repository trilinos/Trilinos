#include "Petra_Petra.h"
#include "Petra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "BuildTestProblems.h"
  int RDP_MatrixTests(const Petra_BlockMap & Map, const Petra_LocalMap & LocalMap, int NumVectors,
		      bool verbose)
  {
    const Petra_Comm & Comm = Map.Comm();
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

    // Construct MultiVectors

  {
    Petra_RDP_MultiVector& A          = *new Petra_RDP_MultiVector(LocalMap, NumVectors);
    Petra_RDP_MultiVector& B          = *new Petra_RDP_MultiVector(LocalMap, NumVectors);
    Petra_LocalMap & Map2d = *new Petra_LocalMap(NumVectors, IndexBase, Comm);
    Petra_RDP_MultiVector& C          = *new Petra_RDP_MultiVector(Map2d, NumVectors);
    Petra_RDP_MultiVector& C_GEMM     = *new Petra_RDP_MultiVector(Map2d, NumVectors);

    double **App, **Bpp, **Cpp;
    
    Petra_RDP_MultiVector *Ap, *Bp, *Cp;

    // For testing non-strided mode, create MultiVectors that are scattered throughout memory

    App = new double *[NumVectors];
    Bpp = new double *[NumVectors];
    Cpp = new double *[NumVectors];
    for (i=0; i<NumVectors; i++) App[i] = new double[A.MyLength()+i];
    for (i=0; i<NumVectors; i++) Bpp[i] = new double[B.MyLength()+i];
    for (i=0; i<NumVectors; i++) Cpp[i] = new double[C.MyLength()+i];
    
    Petra_RDP_MultiVector& A1 = *new Petra_RDP_MultiVector(View, LocalMap, App, NumVectors);
    Petra_RDP_MultiVector& B1 = *new Petra_RDP_MultiVector(View, LocalMap, Bpp, NumVectors);
    Petra_RDP_MultiVector& C1 = *new Petra_RDP_MultiVector(View, Map2d, Cpp, NumVectors);

    for (int strided = 0; strided<2; strided++){
    int ierr;
    // Loop through all trans cases using a variety of values for alpha and beta
    for (i=0; i<4; i++){
	ierr = 0;
	char transa = 'N'; if (i>1) transa = 'T';
	char transb = 'N'; if (i%2!=0) transb = 'T';
	double alpha = (double) i+1;
	double beta  = (double) (i/2);
	ierr += C.Random();  // Fill C with random numbers
	ierr += BuildRDP_MatrixTests(C,transa, transb, alpha, A, B, beta, C_GEMM );
	if (strided)
	  {
	    Ap = &A; Bp = &B; Cp = &C;
	  }
	else
	  {
	    A.ExtractCopy(App); Ap = &A1;
	    B.ExtractCopy(Bpp); Bp = &B1;
	    C.ExtractCopy(Cpp); Cp = &C1;
	  }
	  
	ierr += Cp->Multiply(transa, transb, alpha, *Ap, *Bp, beta);
	ierr += Cp->Update(-1.0, C_GEMM, 1.0);
	ierr += Cp->Norm2(residual);

	if (verbose && ierr==0)
	  {
	    cout << "\n\nXXXXX Replicated Local MultiVector GEMM tests XXXXX\n" << endl;
	    if (strided)
	    cout << "\n\nXXXXX              Strided Multivectors       XXXXX\n" << endl;
	    else
	    cout << "\n\nXXXXX          Non-Strided Multivectors       XXXXX\n" << endl;
	    cout << "\n  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
		 <<", transb = " << transb << endl;
	    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
	  }
	if (ierr==0 && BadResidual(residual, NumVectors)) return(-1);
      }

      }
    for (i=0; i<NumVectors; i++)
      {
	delete [] App[i];
	delete [] Bpp[i];
	delete [] Cpp[i];
      }
    delete &A;
    delete &B;
    delete &C;
    delete &C_GEMM;
    delete &Map2d;
    delete [] App;
    delete [] Bpp;
    delete [] Cpp;
    
    delete &A1;
    delete &B1;
    delete &C1;
  }
      
    // ====================================
    // Case 5  (A, B distributed C  local)
    // ====================================

    // Construct MultiVectors
  {
    Petra_RDP_MultiVector& A          = *new Petra_RDP_MultiVector(Map, NumVectors);
    Petra_RDP_MultiVector& B          = *new Petra_RDP_MultiVector(Map, NumVectors);
    Petra_LocalMap & Map2d = *new Petra_LocalMap(NumVectors, IndexBase, Comm);
    Petra_RDP_MultiVector& C          = *new Petra_RDP_MultiVector(Map2d, NumVectors);
    Petra_RDP_MultiVector& C_GEMM     = *new Petra_RDP_MultiVector(Map2d, NumVectors);

    char transa = 'T';
    char transb = 'N';
    double alpha = 2.0;
    double beta  = 1.0;
    ierr += C.Random();  // Fill C with random numbers
    ierr += BuildRDP_MatrixTests(C, transa, transb, alpha, A, B, beta, C_GEMM );
    ierr += C.Multiply(transa, transb, alpha, A, B, beta);
    ierr += C.Update(-1.0, C_GEMM, 1.0);
    ierr += C.Norm2(residual);

    if (verbose && ierr==0)
      {
	cout << "\n\nXXXXX   Generalized 2D dot product via GEMM call  XXXXX\n" << endl;
	cout << "\n  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
	     <<", transb = " << transb << endl;
	for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
      }
    if (BadResidual(residual, NumVectors)) return(-1);
    
    
    delete &A;
    delete &B;
    delete &C;
    delete &C_GEMM;
    delete &Map2d;
  }      
    // ====================================
    // Case 6-7  (A, C distributed, B local)
    // ====================================

    // Construct MultiVectors
  {
    Petra_RDP_MultiVector& A          = *new Petra_RDP_MultiVector(Map, NumVectors);
    Petra_LocalMap & Map2d = *new Petra_LocalMap(NumVectors, IndexBase, Comm);
    Petra_RDP_MultiVector& B          = *new Petra_RDP_MultiVector(Map2d, NumVectors);
    Petra_RDP_MultiVector& C          = *new Petra_RDP_MultiVector(Map, NumVectors);
    Petra_RDP_MultiVector& C_GEMM     = *new Petra_RDP_MultiVector(Map, NumVectors);

    for (i=0; i<2; i++)
      {
	char transa = 'N';
	char transb = 'N'; if (i>0) transb = 'T';
	double alpha = 2.0;
	double beta  = 1.1;
	ierr += C.Random();  // Fill C with random numbers
	ierr += BuildRDP_MatrixTests(C,transa, transb, alpha, A, B, beta, C_GEMM );
	ierr += C.Multiply(transa, transb, alpha, A, B, beta);
	ierr += C.Update(-1.0, C_GEMM, 1.0);
	ierr += C.Norm2(residual);
	
	if (verbose)
	  {
	    cout << "\n\nXXXXX   Generalized 2D vector update via GEMM call  XXXXX\n" << endl;
	    cout << "\n  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
		 <<", transb = " << transb << endl;
	    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
	  }
	if (BadResidual(residual, NumVectors)) return(-1);
      }

    delete &A;
    delete &B;
    delete &C;
    delete &C_GEMM; 
    delete &Map2d;
    delete [] residual;
    
    return(ierr);
  }
  }

int RDP_MultiVectorTests(const Petra_BlockMap & Map, int NumVectors, bool verbose)
{
  const Petra_Comm & Comm = Map.Comm();
  int ierr = 0, i;
  double *residual = new double[NumVectors];
  
  Petra_BLAS BLAS;
  /* get number of processors and the name of this processor */
  
  // int NumProc = Comm.getNumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct MultiVectors
  
  Petra_RDP_MultiVector& A             = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& sqrtA         = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& B             = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& C             = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& C_alphaA      = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& C_alphaAplusB = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& C_plusB       = *new Petra_RDP_MultiVector(Map, NumVectors);
  Petra_RDP_MultiVector& Weights       = *new Petra_RDP_MultiVector(Map, NumVectors);
  
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

  
  C.Random(); // Fill C with random numbers.
  double alpha = 2.0;
  BuildRDP_MultiVectorTests (C,alpha, A, sqrtA, B, C_alphaA, C_alphaAplusB,
			     C_plusB, dotvec_AB, norm1_A, norm2_sqrtA, norminf_A, 
			     normw_A, Weights, minval_A, maxval_A, meanval_A);
  
  if (verbose) cout << "\n\nXXXXX   Testing alpha * A  XXXXX\n";
  // Test alpha*A
  Petra_RDP_MultiVector& alphaA = *new Petra_RDP_MultiVector(A); // Copy of A
  ierr += alphaA.Scale(alpha);
  ierr += alphaA.Update(-1.0, C_alphaA, 1.0);
  ierr += alphaA.Norm2(residual);
  
  if (ierr!=0 && verbose) 
    cout << "Error in alpha * A MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
	{
	  cout << "\n  alpha = " << alpha << endl;
	  for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
	}
  delete &alphaA;
  
  if (verbose) cout << "\n\nXXXXX   Testing C = alpha * A + B   XXXXX\n";
  // Test alpha*A + B
  Petra_RDP_MultiVector& alphaAplusB = *new Petra_RDP_MultiVector(A); // Copy of A
  ierr += alphaAplusB.Update(1.0, B, alpha, A, 0.0);
  ierr += alphaAplusB.Update(-1.0, C_alphaAplusB, 1.0);
  ierr += alphaAplusB.Norm2(residual);
  
  if (ierr!=0 && verbose) 
    cout << "Error in alpha * A + B MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    {
      cout << "\n  alpha = " << alpha << endl;
      for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
    }
  if (BadResidual(residual, NumVectors)) return(-1);
  delete &alphaAplusB;
  
  if (verbose) cout << "\n\nXXXXX   Testing C += B   XXXXX\n\n";
  // Test + B
  Petra_RDP_MultiVector& plusB = *new Petra_RDP_MultiVector(C); // Copy of C
  ierr += plusB.Update(1.0, B, 1.0);
  ierr += plusB.Update(-1.0, C_plusB, 1.0);
  ierr += plusB.Norm2(residual);
  
  if (ierr!=0 && verbose) 
    cout << "Error in + B MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
  
  delete &plusB;
  
  if (verbose) cout << "\n\nXXXXX  Testing A.dotProd(B)  XXXXX\n\n";
  // Test A.dotvec(B)
  double *dotvec = residual;
  ierr += A.Dot(B,dotvec);
  BLAS.AXPY(NumVectors,-1.0,dotvec_AB,dotvec);
  
  if (ierr!=0 && verbose) 
    cout << "Error dotvec MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
  
  
  if (verbose) cout << "\n\nXXXXX   Testing norm1_A   XXXXX\n\n";
  // Test A.norm1()
  double *norm1 = residual;
  ierr += A.Norm1(norm1);
  BLAS.AXPY(NumVectors,-1.0,norm1_A,norm1);
  
  if (ierr!=0 && verbose)
    cout << "Error in norm1 MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	  
  
  if (verbose) cout << "\n\nXXXXX Testing norm2_sqrtA  XXXXX\n\n";
  // Test sqrtA.norm2()
  double *norm2 = residual;
  ierr += sqrtA.Norm2(norm2);
  BLAS.AXPY(NumVectors,-1.0,norm2_sqrtA,norm2);
  
  if (ierr!=0 && verbose)
    cout << "Error in norm2 MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	
  
  if (verbose) cout << "\n\nXXXXX Testing norminf_A  XXXXX\n\n";
  // Test A.norminf()
  double *norminf = residual;
  ierr += A.NormInf(norminf);
  BLAS.AXPY(NumVectors,-1.0,norminf_A,norminf);
  
  if (ierr!=0 && verbose)
    cout << "Error in NormInf MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	
  
  if (verbose) cout << "\n\nXXXXX Testing normw_A  XXXXX\n\n";
  // Test A.NormWeighted()
  double *normw = residual;
  ierr += A.NormWeighted(Weights, normw);
  BLAS.AXPY(NumVectors,-1.0,normw_A,normw);
  
  if (ierr!=0 && verbose)
    cout << "Error in NormWeighted MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	
  
  if (verbose) cout << "\n\nXXXXX Testing minval_A  XXXXX\n\n";
  // Test A.MinValue()
  double *minval = residual;
  ierr += A.MinValue(minval);
  BLAS.AXPY(NumVectors,-1.0,minval_A,minval);
  
  if (ierr!=0 && verbose)
    cout << "Error in MinValue MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	
  
  if (verbose) cout << "\n\nXXXXX Testing maxval_A  XXXXX\n\n";
  // Test A.MaxValue()
  double *maxval = residual;
  ierr += A.MaxValue(maxval);
  BLAS.AXPY(NumVectors,-1.0,maxval_A,maxval);
  
  if (ierr!=0 && verbose)
    cout << "Error in MaxValue MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	  
  if (verbose) cout << "\n\nXXXXX Testing meanval_A  XXXXX\n\n";
  // Test A.MeanValue()
  double *meanval = residual;
  ierr += A.MeanValue(meanval);
  BLAS.AXPY(NumVectors,-1.0,meanval_A,meanval);
  
  if (ierr!=0 && verbose)
    cout << "Error in MeanValue MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	
  
  if (verbose) cout << "\n\nXXXXX Testing abs_A  XXXXX\n\n";
  // Test A.Abs()
  Petra_RDP_MultiVector& Abs_A = A;
  ierr += Abs_A.Abs(A);
  ierr += Abs_A.Update(1.0, A, -1.0); // Abs_A = A - Abs_A (should be zero since A > 0)
  ierr += Abs_A.Norm2(residual);
  if (ierr!=0 && verbose)
    cout << "Error in Absolute value MultiVector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual, NumVectors)) return(-1);
	
  
  // Delete everything
  
  delete &A;
  delete &sqrtA;
  delete &B;
  delete &C;
  delete &C_alphaA;
  delete &C_alphaAplusB;
  delete &C_plusB;
  delete &Weights;
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

int BadResidual(double * Residual, int NumVectors)
{
  double threshold = 5.0E-6;
  for (int i=0; i<NumVectors; i++) if (Residual[i]>threshold) return(1);
  return(0);
}
