// tpetra/test/DenseMatrix/cc_main.cc
// 16-May-2002 - Changed names to use Tpetra instead of TPetra
// 17-May-2002 - Switched from Petra_Comm, Petra_Time, and Petra_Map to Epetra's versions
// 20-May-2002 - Changed formatting for readability, no real changes

// Used to easily change scalar type
#define SCALARTYPE float

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Time.h" 
#include "Tpetra_DenseMatrix.h"

// Local prototypes
template<class scalarType> bool check(int M, int N, Tpetra::DenseMatrix<scalarType> A, Tpetra::DenseMatrix<scalarType> B, Tpetra::DenseMatrix<scalarType>& C);

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI 

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

  bool verbose1 = verbose;
  verbose = (MyPID==0);  // Only print most results on PE 0

  if (verbose) cout << endl 
		    << "#################################################"  << endl
		    << "Testing ScalarType: " << Tpetra::ScalarTraits<SCALARTYPE>::name() << endl
		    << "#################################################" << endl;

  if (verbose && Tpetra::ScalarTraits<SCALARTYPE>::haveMachineParameters()) cout << endl << endl
		    << "Some properties of this ScalarType:" << endl
		    << "zero            = " << Tpetra::ScalarTraits<SCALARTYPE>::zero() << endl
		    << "one             = " << Tpetra::ScalarTraits<SCALARTYPE>::one() << endl
		    << "eps             = " << Tpetra::ScalarTraits<SCALARTYPE>::eps() << endl
		    << "sfmin           = " << Tpetra::ScalarTraits<SCALARTYPE>::sfmin() << endl
		    << "base            = " << Tpetra::ScalarTraits<SCALARTYPE>::base() << endl
		    << "prec            = " << Tpetra::ScalarTraits<SCALARTYPE>::prec() << endl
		    << "t               = " << Tpetra::ScalarTraits<SCALARTYPE>::t() << endl
		    << "rnd             = " << Tpetra::ScalarTraits<SCALARTYPE>::rnd() << endl
		    << "emin            = " << Tpetra::ScalarTraits<SCALARTYPE>::emin() << endl
		    << "rmin            = " << Tpetra::ScalarTraits<SCALARTYPE>::rmin() << endl
		    << "emax            = " << Tpetra::ScalarTraits<SCALARTYPE>::emax() << endl
		    << "rmax            = " << Tpetra::ScalarTraits<SCALARTYPE>::rmax() << endl << endl;

  Tpetra::DenseMatrix<SCALARTYPE> A;
  Tpetra::DenseMatrix<SCALARTYPE> B;
  Tpetra::DenseMatrix<SCALARTYPE> C1;
  Tpetra::DenseMatrix<SCALARTYPE> C2;
  

  int M = 600;
  int N = 700;
  SCALARTYPE zero = Tpetra::ScalarTraits<SCALARTYPE>::zero();
  SCALARTYPE one = Tpetra::ScalarTraits<SCALARTYPE>::one(); 
  A.shape(M,N);
  B.shape(N,M);
  C1.shape(M,M);
  C2.shape(M,M);

  for (int i=0; i<M; i++)
    for (int j=0; j<N; j++)
      {
	SCALARTYPE curValue = ((Tpetra::ScalarTraits<SCALARTYPE>::magnitudeType) ((i+1)*(j+1)))* Tpetra::ScalarTraits<SCALARTYPE>::one();
	A(i,j) = (SCALARTYPE) curValue;
	B(j,i) = (SCALARTYPE) A(i,j);
      }
  bool smallProblem = (M*N<200);
  if (verbose && smallProblem) cout << "\nContents of A:\n" << A << endl;
  if (verbose && smallProblem) cout << "\nContents of B:\n" << B << endl;

  Epetra_Flops flop_counter;
  C1.SetFlopCounter(flop_counter);
  Epetra_Time timer(Comm);
  
  double startFlops = C1.Flops();
  double startTime = timer.ElapsedTime();
  C1.multiply('N', 'N', one, A, B, zero);
  double time = timer.ElapsedTime() - startTime;
  double flops = C1.Flops() - startFlops;
  double MFLOPS = flops/time/1000000.0;

  if (verbose) cout << "Statistics for multiply method C(MxM) = A(MxN) * B(NxM) " << "M = " << M
		    << " N = " << N << " is:" << endl 
		    << " Total FLOPS  = " << flops << endl 
		    << " Total Time   = " << time  << endl 
		    << " Total MFLOPS = " << MFLOPS << endl;

  if (verbose && smallProblem) cout << "\nContents of C1:\n" << C1 << endl;

  if (!check(M, N, A, B, C2)) return 1; // return 1 if tests fail

      
  if (verbose && smallProblem) cout << "\nContents of C2:\n" << C2 << endl;

  Tpetra::ScalarTraits<SCALARTYPE>::magnitudeType c1Norm = C1.oneNorm();
  Tpetra::ScalarTraits<SCALARTYPE>::magnitudeType c2Norm = C2.oneNorm();
  //double c1Norm = 0;
  //double c2Norm = 0;

  if (verbose) {
    if (c1Norm!=c2Norm)
      cout << "Norms of C1 and C2 differ:" << endl;
    cout  << "Norm of C1 = " << c1Norm << endl
	  << "Norm of C2 = " << c2Norm << endl;
  if (c1Norm==c2Norm) cout << "Tpetra::DenseMatrix check OK" << endl;
  }
	  
  return 0; // All done
}
 
template<class scalarType> bool check(int M, int N, Tpetra::DenseMatrix<scalarType> A, Tpetra::DenseMatrix<scalarType> B, Tpetra::DenseMatrix<scalarType>& C)
{

  // Confirm correct dimensions

  assert(A.numRows()==M);
  assert(A.numCols()==N);
  assert(B.numRows()==N);
  assert(B.numCols()==M);
  assert(C.numRows()==M);
  assert(C.numCols()==M);

  scalarType zero = Tpetra::ScalarTraits<scalarType>::zero();
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      C(i,j) = zero; // Initialize C to zero
    
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
    {
      for (int k = 0; k<N; k++)
        C(i,j) += A(i,k)*B(k,j);
    }

  return(true);
}
