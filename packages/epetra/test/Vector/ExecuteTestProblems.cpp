#include "Epetra_Object.h"
#include "ExecuteTestProblems.h"
#include "BuildTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_BLAS.h"
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
    Epetra_Vector& A          = *new Epetra_Vector(LocalMap);
    Epetra_Vector& B          = *new Epetra_Vector(LocalMap);
    Epetra_LocalMap & Map2d = *new Epetra_LocalMap(NumVectors, IndexBase, Comm);
    Epetra_Vector& C          = *new Epetra_Vector(Map2d);
    Epetra_Vector& C_GEMM     = *new Epetra_Vector(Map2d);

    double *App, *Bpp, *Cpp;
    
    Epetra_Vector *Ap, *Bp, *Cp;

    // For testing non-strided mode, create Vectors that are scattered throughout memory

    App = new double [NumVectors];
    Bpp = new double [NumVectors];
    Cpp = new double [NumVectors];
    Epetra_Vector& A1 = *new Epetra_Vector(View, LocalMap, App);
    Epetra_Vector& B1 = *new Epetra_Vector(View, LocalMap, Bpp);
    Epetra_Vector& C1 = *new Epetra_Vector(View, Map2d, Cpp);

    int ierr;
    // Loop through all trans cases using a variety of values for alpha and beta
    for (i=0; i<4; i++){
	ierr = 0;
	char transa = 'N'; if (i/2) transa = 'T';
	char transb = 'N'; if (i%2) transb = 'N';
	double alpha = (double) i+1;
	double beta  = (double) (i/2);
	ierr += C.Random();  // Fill C with random numbers
	ierr += BuildMatrixTests(C,transa, transb, alpha, A, B, beta, C_GEMM );
	Ap = &A; Bp = &B; Cp = &C;
	ierr += Cp->Multiply(transa, transb, alpha, *Ap, *Bp, beta);
	ierr += Cp->Update(-1.0, C_GEMM, 1.0);
	ierr += Cp->Norm2(residual);

	if (verbose && ierr==0)
	  {
	    cout << "\n\nXXXXX Replicated Local Vector GEMM tests XXXXX\n" << endl;
	    cout << "\n  alpha = " << alpha << ",  beta = " << beta <<", transa = "<<transa
		 <<", transb = " << transb << endl;
	    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
	  }
	if (ierr==0 && BadResidual(residual)) return(-1);
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

    // Construct Vectors
  {
    Epetra_Vector& A          = *new Epetra_Vector(Map);
    Epetra_Vector& B          = *new Epetra_Vector(Map);
    Epetra_LocalMap & Map2d = *new Epetra_LocalMap(NumVectors, IndexBase, Comm);
    Epetra_Vector& C          = *new Epetra_Vector(Map2d);
    Epetra_Vector& C_GEMM     = *new Epetra_Vector(Map2d);

    char transa = 'T';
    char transb = 'N';
    double alpha = 2.0;
    double beta  = 1.0;
    ierr += C.Random();  // Fill C with random numbers
    ierr += BuildMatrixTests(C, transa, transb, alpha, A, B, beta, C_GEMM );
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
    if (BadResidual(residual)) return(-1);
    
    
    delete &A;
    delete &B;
    delete &C;
    delete &C_GEMM;
    delete &Map2d;
  }      
    // ====================================
    // Case 6-7  (A, C distributed, B local)
    // ====================================

    // Construct Vectors
  {
    Epetra_Vector& A          = *new Epetra_Vector(Map);
    Epetra_LocalMap & Map2d = *new Epetra_LocalMap(NumVectors, IndexBase, Comm);
    Epetra_Vector& B          = *new Epetra_Vector(Map2d);
    Epetra_Vector& C          = *new Epetra_Vector(Map);
    Epetra_Vector& C_GEMM     = *new Epetra_Vector(Map);

    for (i=0; i<2; i++)
      {
	char transa = 'N';
	char transb = 'N'; if (i) transb = 'T';
	double alpha = 2.0;
	double beta  = 1.1;
	ierr += C.Random();  // Fill C with random numbers
	ierr += BuildMatrixTests(C,transa, transb, alpha, A, B, beta, C_GEMM );
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
	if (BadResidual(residual)) return(-1);
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

int VectorTests(const Epetra_BlockMap & Map, bool verbose)
{
  int NumVectors = 1;
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0, i;
  double *residual = new double[NumVectors];
  Epetra_BLAS Blas;
  
  /* get number of processors and the name of this processor */
  
  // int NumProc = Comm.getNumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct Vectors
  
  Epetra_Vector& A             = *new Epetra_Vector(Map);
  Epetra_Vector& sqrtA         = *new Epetra_Vector(Map);
  Epetra_Vector& B             = *new Epetra_Vector(Map);
  Epetra_Vector& C             = *new Epetra_Vector(Map);
  Epetra_Vector& C_alphaA      = *new Epetra_Vector(Map);
  Epetra_Vector& C_alphaAplusB = *new Epetra_Vector(Map);
  Epetra_Vector& C_plusB       = *new Epetra_Vector(Map);
  Epetra_Vector& Weights       = *new Epetra_Vector(Map);
  
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
  BuildMultiVectorTests (C,alpha, A, sqrtA, B, C_alphaA, C_alphaAplusB,
			     C_plusB, dotvec_AB, norm1_A, norm2_sqrtA, norminf_A, 
			     normw_A, Weights, minval_A, maxval_A, meanval_A);
  
  if (verbose) cout << "\n\nXXXXX   Testing alpha * A  XXXXX\n";
  // Test alpha*A
  Epetra_Vector& alphaA = *new Epetra_Vector(A); // Copy of A
  ierr += alphaA.Scale(alpha);
  ierr += alphaA.Update(-1.0, C_alphaA, 1.0);
  ierr += alphaA.Norm2(residual);
  
  if (ierr!=0 && verbose) 
    cout << "Error in alpha * A Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
	{
	  cout << "\n  alpha = " << alpha << endl;
	  for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
	}
  delete &alphaA;
  
  if (verbose) cout << "\n\nXXXXX   Testing C = alpha * A + B   XXXXX\n";
  // Test alpha*A + B
  Epetra_Vector& alphaAplusB = *new Epetra_Vector(A); // Copy of A
  ierr += alphaAplusB.Update(1.0, B, alpha, A, 0.0);
  ierr += alphaAplusB.Update(-1.0, C_alphaAplusB, 1.0);
  ierr += alphaAplusB.Norm2(residual);
  
  if (ierr!=0 && verbose) 
    cout << "Error in alpha * A + B Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    {
      cout << "\n  alpha = " << alpha << endl;
      for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
    }
  if (BadResidual(residual)) return(-1);
  delete &alphaAplusB;
  
  if (verbose) cout << "\n\nXXXXX   Testing C += B   XXXXX\n\n";
  // Test + B
  Epetra_Vector& plusB = *new Epetra_Vector(C); // Copy of C
  ierr += plusB.Update(1.0, B, 1.0);
  ierr += plusB.Update(-1.0, C_plusB, 1.0);
  ierr += plusB.Norm2(residual);
  
  if (ierr!=0 && verbose) 
    cout << "Error in + B Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
  
  delete &plusB;
  
  if (verbose) cout << "\n\nXXXXX  Testing A.dotProd(B)  XXXXX\n\n";
  // Test A.dotvec(B)
  double *dotvec = new double[NumVectors];
  ierr += A.Dot(B,dotvec);
  Blas.AXPY(NumVectors,-1.0,dotvec_AB,dotvec);
  
  if (ierr!=0 && verbose) 
    cout << "Error dotvec Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
  
  delete [] dotvec;
  
  if (verbose) cout << "\n\nXXXXX   Testing norm1_A   XXXXX\n\n";
  // Test A.norm1()
  double *norm1 = new double[NumVectors];
  ierr += A.Norm1(norm1);
  Blas.AXPY(NumVectors,-1.0,norm1_A,norm1);
  
  if (ierr!=0 && verbose)
    cout << "Error in norm1 Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] norm1;
  
  
  if (verbose) cout << "\n\nXXXXX Testing norm2_sqrtA  XXXXX\n\n";
  // Test sqrtA.norm2()
  double *norm2 = new double[NumVectors];
  ierr += sqrtA.Norm2(norm2);
  Blas.AXPY(NumVectors,-1.0,norm2_sqrtA,norm2);
  
  if (ierr!=0 && verbose)
    cout << "Error in norm2 Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] norm2;
  
  if (verbose) cout << "\n\nXXXXX Testing norminf_A  XXXXX\n\n";
  // Test A.norminf()
  double *norminf = new double[NumVectors];
  ierr += A.NormInf(norminf);
  Blas.AXPY(NumVectors,-1.0,norminf_A,norminf);
  
  if (ierr!=0 && verbose)
    cout << "Error in NormInf Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] norminf;
  
  if (verbose) cout << "\n\nXXXXX Testing normw_A  XXXXX\n\n";
  // Test A.NormWeighted()
  double *normw = new double[NumVectors];
  ierr += A.NormWeighted(Weights, normw);
  Blas.AXPY(NumVectors,-1.0,normw_A,normw);
  
  if (ierr!=0 && verbose)
    cout << "Error in NormWeighted Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] normw;
  
  if (verbose) cout << "\n\nXXXXX Testing minval_A  XXXXX\n\n";
  // Test A.MinValue()
  double *minval = new double[NumVectors];
  ierr += A.MinValue(minval);
  Blas.AXPY(NumVectors,-1.0,minval_A,minval);
  
  if (ierr!=0 && verbose)
    cout << "Error in MinValue Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] minval;
  
  if (verbose) cout << "\n\nXXXXX Testing maxval_A  XXXXX\n\n";
  // Test A.MaxValue()
  double *maxval = new double[NumVectors];
  ierr += A.MaxValue(maxval);
  Blas.AXPY(NumVectors,-1.0,maxval_A,maxval);
  
  if (ierr!=0 && verbose)
    cout << "Error in MaxValue Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] maxval;
  
  if (verbose) cout << "\n\nXXXXX Testing meanval_A  XXXXX\n\n";
  // Test A.MeanValue()
  double *meanval = new double[NumVectors];
  ierr += A.MeanValue(meanval);
  Blas.AXPY(NumVectors,-1.0,meanval_A,meanval);
  
  if (ierr!=0 && verbose)
    cout << "Error in MeanValue Vector testing\n";
  if (ierr) return(-2);
  if (verbose)
    for (int j=0; j< NumVectors; j++) cout << "     Residual[" << j <<"] = " << residual[j] << endl;
  if (BadResidual(residual)) return(-1);
	
  delete [] meanval;
  
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

int BadResidual(double * Residual)
{
  int NumVectors = 1;
  double threshold = 1.0E-7;
  for (int i=0; i<NumVectors; i++) if (Residual[i]>threshold) return(1);
  return(0);
}
