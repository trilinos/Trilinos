//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


#define EPETRA_HAVE_JADMATRIX
#define EPETRA_VERY_SHORT_PERFTEST
#define EPETRA_HAVE_STATICPROFILE
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"
#ifdef EPETRA_HAVE_JADMATRIX
#include "Epetra_JadMatrix.h"
#endif

// prototypes

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact, bool StaticProfile, bool MakeLocalOnly);

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, int nrhs,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact, bool StaticProfile, bool MakeLocalOnly);
 
void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			int nsizes, int * sizes,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact, bool StaticProfile, bool MakeLocalOnly);

void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, 
			int nsizes, int * sizes, int nrhs,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact, bool StaticProfile, bool MakeLocalOnly);

void GenerateMyGlobalElements(int numNodesX, int numNodesY, int numProcsX, int numProcs,
			      int myPID, long long * & myGlobalElements);

void runMatrixTests(Epetra_CrsMatrix * A,  Epetra_MultiVector * b, Epetra_MultiVector * bt,
		    Epetra_MultiVector * xexact, bool StaticProfile, bool verbose, bool summary);
#ifdef EPETRA_HAVE_JADMATRIX
void runJadMatrixTests(Epetra_JadMatrix * A,  Epetra_MultiVector * b, Epetra_MultiVector * bt,
		    Epetra_MultiVector * xexact, bool StaticProfile, bool verbose, bool summary);
#endif
void runLUMatrixTests(Epetra_CrsMatrix * L,  Epetra_MultiVector * bL, Epetra_MultiVector * btL, Epetra_MultiVector * xexactL, 
		      Epetra_CrsMatrix * U,  Epetra_MultiVector * bU, Epetra_MultiVector * btU, Epetra_MultiVector * xexactU, 
		      bool StaticProfile, bool verbose, bool summary);
int main(int argc, char *argv[])
{
  int ierr = 0;
  double elapsed_time;
  double total_flops;
  double MFLOPs;
    

#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif

  bool verbose = false;
  bool summary = false;

  // Check if we should print verbose results to standard out
  if (argc>6) if (argv[6][0]=='-' && argv[6][1]=='v') verbose = true;

  // Check if we should print verbose results to standard out
  if (argc>6) if (argv[6][0]=='-' && argv[6][1]=='s') summary = true;

  if(argc < 6) {
    cerr << "Usage: " << argv[0]
         << " NumNodesX NumNodesY NumProcX NumProcY NumPoints [-v|-s]" << endl
         << "where:" << endl
         << "NumNodesX         - Number of mesh nodes in X direction per processor" << endl
         << "NumNodesY         - Number of mesh nodes in Y direction per processor" << endl
         << "NumProcX          - Number of processors to use in X direction" << endl
         << "NumProcY          - Number of processors to use in Y direction" << endl
         << "NumPoints         - Number of points to use in stencil (5, 9 or 25 only)" << endl
         << "-v|-s             - (Optional) Run in verbose mode if -v present or summary mode if -s present" << endl
         << " NOTES: NumProcX*NumProcY must equal the number of processors used to run the problem." << endl << endl
	 << " Serial example:" << endl
         << argv[0] << " 16 12 1 1 25 -v" << endl
	 << " Run this program in verbose mode on 1 processor using a 16 X 12 grid with a 25 point stencil."<< endl <<endl
	 << " MPI example:" << endl
         << "mpirun -np 32 " << argv[0] << " 10 12 4 8 9 -v" << endl
	 << " Run this program in verbose mode on 32 processors putting a 10 X 12 subgrid on each processor using 4 processors "<< endl
	 << " in the X direction and 8 in the Y direction.  Total grid size is 40 points in X and 96 in Y with a 9 point stencil."<< endl
         << endl;
    return(1);

  }
    //char tmp;
    //if (comm.MyPID()==0) cout << "Press any key to continue..."<< endl;
    //if (comm.MyPID()==0) cin >> tmp;
    //comm.Barrier();

  comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  if (verbose && comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;
  if (summary && comm.MyPID()==0) {
    if (comm.NumProc()==1)
      cout << Epetra_Version() << endl << endl;
    else
      cout << endl << endl; // Print two blank line to keep output columns lined up
  }

  if (verbose) cout << comm <<endl;


  // Redefine verbose to only print on PE 0

  if (verbose && comm.MyPID()!=0) verbose = false;
  if (summary && comm.MyPID()!=0) summary = false;

  int numNodesX = atoi(argv[1]);
  int numNodesY = atoi(argv[2]);
  int numProcsX = atoi(argv[3]);
  int numProcsY = atoi(argv[4]);
  int numPoints = atoi(argv[5]);

  if (verbose || (summary && comm.NumProc()==1)) {
    cout << " Number of local nodes in X direction  = " << numNodesX << endl
	 << " Number of local nodes in Y direction  = " << numNodesY << endl
	 << " Number of global nodes in X direction = " << numNodesX*numProcsX << endl
	 << " Number of global nodes in Y direction = " << numNodesY*numProcsY << endl
	 << " Number of local nonzero entries       = " << numNodesX*numNodesY*numPoints << endl
	 << " Number of global nonzero entries      = " << numNodesX*numNodesY*numPoints*numProcsX*numProcsY << endl
	 << " Number of Processors in X direction   = " << numProcsX << endl
	 << " Number of Processors in Y direction   = " << numProcsY << endl
	 << " Number of Points in stencil           = " << numPoints << endl << endl;
  }
  // Print blank line to keep output columns lined up
  if (summary && comm.NumProc()>1)
    cout << endl << endl << endl << endl << endl << endl << endl << endl<< endl << endl;

  if (numProcsX*numProcsY!=comm.NumProc()) {
    cerr << "Number of processors = " << comm.NumProc() << endl
	 << " is not the product of " << numProcsX << " and " << numProcsY << endl << endl;
    return(1);
  }

  if (numPoints!=5 && numPoints!=9 && numPoints!=25) {
    cerr << "Number of points specified = " << numPoints << endl
	 << " is not 5, 9, 25" << endl << endl;
    return(1);
  }

  if (numNodesX*numNodesY<=0) {
    cerr << "Product of number of nodes is <= zero" << endl << endl;
    return(1);
  }

  Epetra_IntSerialDenseVector Xoff, XLoff, XUoff;
  Epetra_IntSerialDenseVector Yoff, YLoff, YUoff;
  if (numPoints==5) {

     // Generate a 5-point 2D Finite Difference matrix
    Xoff.Size(5);
    Yoff.Size(5);
    Xoff[0] = -1; Xoff[1] = 1; Xoff[2] = 0; Xoff[3] = 0;  Xoff[4] = 0; 
    Yoff[0] = 0;  Yoff[1] = 0; Yoff[2] = 0; Yoff[3] = -1; Yoff[4] = 1; 

     // Generate a 2-point 2D Lower triangular Finite Difference matrix
    XLoff.Size(2);
    YLoff.Size(2);
    XLoff[0] = -1; XLoff[1] =  0; 
    YLoff[0] =  0; YLoff[1] = -1;

     // Generate a 3-point 2D upper triangular Finite Difference matrix
    XUoff.Size(3);
    YUoff.Size(3);
    XUoff[0] =  0; XUoff[1] =  1; XUoff[2] = 0; 
    YUoff[0] =  0; YUoff[1] =  0; YUoff[2] = 1;
  }
  else if (numPoints==9) {
    // Generate a 9-point 2D Finite Difference matrix
    Xoff.Size(9);
    Yoff.Size(9);
    Xoff[0] = -1;  Xoff[1] =  0; Xoff[2] =  1; 
    Yoff[0] = -1;  Yoff[1] = -1; Yoff[2] = -1; 
    Xoff[3] = -1;  Xoff[4] =  0; Xoff[5] =  1; 
    Yoff[3] =  0;  Yoff[4] =  0; Yoff[5] =  0; 
    Xoff[6] = -1;  Xoff[7] =  0; Xoff[8] =  1; 
    Yoff[6] =  1;  Yoff[7] =  1; Yoff[8] =  1; 

    // Generate a 5-point lower triangular 2D Finite Difference matrix
    XLoff.Size(5);
    YLoff.Size(5);
    XLoff[0] = -1;  XLoff[1] =  0; Xoff[2] =  1; 
    YLoff[0] = -1;  YLoff[1] = -1; Yoff[2] = -1; 
    XLoff[3] = -1;  XLoff[4] =  0; 
    YLoff[3] =  0;  YLoff[4] =  0;

    // Generate a 4-point upper triangular 2D Finite Difference matrix
    XUoff.Size(4);
    YUoff.Size(4);
    XUoff[0] =  1; 
    YUoff[0] =  0; 
    XUoff[1] = -1;  XUoff[2] =  0; XUoff[3] =  1; 
    YUoff[1] =  1;  YUoff[2] =  1; YUoff[3] =  1; 

  }
  else {
    // Generate a 25-point 2D Finite Difference matrix
    Xoff.Size(25);
    Yoff.Size(25);
    int xi = 0, yi = 0;
    int xo = -2, yo = -2;
    Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
    Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
    xo = -2, yo++;
    Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
    Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
    xo = -2, yo++;
    Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
    Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
    xo = -2, yo++;
    Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
    Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
    xo = -2, yo++;
    Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
    Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 

    // Generate a 13-point lower triangular 2D Finite Difference matrix
    XLoff.Size(13);
    YLoff.Size(13);
    xi = 0, yi = 0;
    xo = -2, yo = -2;
    XLoff[xi++] = xo++;  XLoff[xi++] = xo++; XLoff[xi++] = xo++; XLoff[xi++] = xo++; XLoff[xi++] = xo++;
    YLoff[yi++] = yo  ;  YLoff[yi++] = yo  ; YLoff[yi++] = yo  ; YLoff[yi++] = yo  ; YLoff[yi++] = yo  ; 
    xo = -2, yo++;
    XLoff[xi++] = xo++;  XLoff[xi++] = xo++; XLoff[xi++] = xo++; XLoff[xi++] = xo++; XLoff[xi++] = xo++;
    YLoff[yi++] = yo  ;  YLoff[yi++] = yo  ; YLoff[yi++] = yo  ; YLoff[yi++] = yo  ; YLoff[yi++] = yo  ; 
    xo = -2, yo++;
    XLoff[xi++] = xo++;  XLoff[xi++] = xo++; XLoff[xi++] = xo++;
    YLoff[yi++] = yo  ;  YLoff[yi++] = yo  ; YLoff[yi++] = yo  ;

    // Generate a 13-point upper triangular 2D Finite Difference matrix
    XUoff.Size(13);
    YUoff.Size(13);
    xi = 0, yi = 0;
    xo = 0, yo = 0;
    XUoff[xi++] = xo++;  XUoff[xi++] = xo++; XUoff[xi++] = xo++;
    YUoff[yi++] = yo  ;  YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; 
    xo = -2, yo++;
    XUoff[xi++] = xo++;  XUoff[xi++] = xo++; XUoff[xi++] = xo++; XUoff[xi++] = xo++; XUoff[xi++] = xo++;
    YUoff[yi++] = yo  ;  YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; 
    xo = -2, yo++;
    XUoff[xi++] = xo++;  XUoff[xi++] = xo++; XUoff[xi++] = xo++; XUoff[xi++] = xo++; XUoff[xi++] = xo++;
    YUoff[yi++] = yo  ;  YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; YUoff[yi++] = yo  ; 

  }

  Epetra_Map * map;
  Epetra_Map * mapL;
  Epetra_Map * mapU;
  Epetra_CrsMatrix * A;
  Epetra_CrsMatrix * L;
  Epetra_CrsMatrix * U;
  Epetra_MultiVector * b;
  Epetra_MultiVector * bt;
  Epetra_MultiVector * xexact;
  Epetra_MultiVector * bL;
  Epetra_MultiVector * btL;
  Epetra_MultiVector * xexactL;
  Epetra_MultiVector * bU;
  Epetra_MultiVector * btU;
  Epetra_MultiVector * xexactU;
  Epetra_SerialDenseVector resvec(0);

  //Timings
  Epetra_Flops flopcounter;
  Epetra_Time timer(comm);

#ifdef EPETRA_VERY_SHORT_PERFTEST
  int jstop = 1;
#elif EPETRA_SHORT_PERFTEST
  int jstop = 1;
#else
  int jstop = 2;
#endif
  for (int j=0; j<jstop; j++) {
    for (int k=1; k<17; k++) {
#ifdef EPETRA_VERY_SHORT_PERFTEST
      if (k<3 || (k%4==0 && k<9)) {
#elif EPETRA_SHORT_PERFTEST
      if (k<6 || k%4==0) {
#else
      if (k<7 || k%2==0) {
#endif
      int nrhs=k;
      if (verbose) cout << "\n*************** Results for " << nrhs << " RHS with ";

      bool StaticProfile = (j!=0);
      if (verbose) {
        if (StaticProfile) cout << " static profile\n";
        else cout << " dynamic profile\n";
      }
      GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints,
			 Xoff.Values(), Yoff.Values(), nrhs, comm, verbose, summary,
			 map, A, b, bt, xexact, StaticProfile, false);

      
#ifdef EPETRA_HAVE_JADMATRIX
      
      timer.ResetStartTime();
      Epetra_JadMatrix JA(*A);
      elapsed_time = timer.ElapsedTime();
      if (verbose) cout << "Time to create Jagged diagonal matrix = " << elapsed_time << endl;

      //cout << "A = " << *A << endl;
      //cout << "JA = " << JA << endl;

      runJadMatrixTests(&JA, b, bt, xexact, StaticProfile, verbose, summary);

#endif
      runMatrixTests(A, b, bt, xexact, StaticProfile, verbose, summary);

      delete A;
      delete b;
      delete bt; 
      delete xexact;

      GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, XLoff.Length(),
			 XLoff.Values(), YLoff.Values(), nrhs, comm, verbose, summary,
			 mapL, L, bL, btL, xexactL, StaticProfile, true);
      

      GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, XUoff.Length(),
			 XUoff.Values(), YUoff.Values(), nrhs, comm, verbose, summary,
			 mapU, U, bU, btU, xexactU, StaticProfile, true);
      

      runLUMatrixTests(L, bL, btL, xexactL, U, bU, btU, xexactU, StaticProfile, verbose, summary);

      delete L;
      delete bL;
      delete btL; 
      delete xexactL;
      delete mapL;

      delete U;
      delete bU;
      delete btU; 
      delete xexactU;
      delete mapU;

      Epetra_MultiVector q(*map, nrhs);
      Epetra_MultiVector z(q);
      Epetra_MultiVector r(q);
      
      delete map;
      q.SetFlopCounter(flopcounter);
      z.SetFlopCounter(q);
      r.SetFlopCounter(q);

      resvec.Resize(nrhs);
      
    
      flopcounter.ResetFlops();
      timer.ResetStartTime();

      //10 norms
      for( int i = 0; i < 10; ++i )
	q.Norm2( resvec.Values() );

      elapsed_time = timer.ElapsedTime();
      total_flops = q.Flops();
      MFLOPs = total_flops/elapsed_time/1000000.0;
      if (verbose) cout << "\nTotal MFLOPs for 10 Norm2's= " << MFLOPs << endl;
      
      if (summary) {
	if (comm.NumProc()==1) cout << "Norm2" << '\t';
	cout << MFLOPs << endl;
      }
      
      flopcounter.ResetFlops();
      timer.ResetStartTime();
      
      //10 dot's
      for( int i = 0; i < 10; ++i )
	q.Dot(z, resvec.Values());
      
      elapsed_time = timer.ElapsedTime();
      total_flops = q.Flops();
      MFLOPs = total_flops/elapsed_time/1000000.0;
      if (verbose) cout << "Total MFLOPs for 10 Dot's  = " << MFLOPs << endl;
      
      if (summary) {
	if (comm.NumProc()==1) cout << "DotProd" << '\t';
	cout << MFLOPs << endl;
      }
      
      flopcounter.ResetFlops();
      timer.ResetStartTime();
      
      //10 dot's
      for( int i = 0; i < 10; ++i )
	q.Update(1.0, z, 1.0, r, 0.0);
      
      elapsed_time = timer.ElapsedTime();
      total_flops = q.Flops();
      MFLOPs = total_flops/elapsed_time/1000000.0;
      if (verbose) cout << "Total MFLOPs for 10 Updates= " << MFLOPs << endl;
      
      if (summary) {
	if (comm.NumProc()==1) cout << "Update" << '\t';
	cout << MFLOPs << endl;
      }
    }
    }
  }
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return ierr ;
}

// Constructs a 2D PDE finite difference matrix using the list of x and y offsets.
// 
// nx      (In) - number of grid points in x direction
// ny      (In) - number of grid points in y direction
//   The total number of equations will be nx*ny ordered such that the x direction changes
//   most rapidly: 
//      First equation is at point (0,0)
//      Second at                  (1,0)
//       ...
//      nx equation at             (nx-1,0)
//      nx+1st equation at         (0,1)

// numPoints (In) - number of points in finite difference stencil
// xoff    (In) - stencil offsets in x direction (of length numPoints)
// yoff    (In) - stencil offsets in y direction (of length numPoints)
//   A standard 5-point finite difference stencil would be described as:
//     numPoints = 5
//     xoff = [-1, 1, 0,  0, 0]
//     yoff = [ 0, 0, 0, -1, 1]

// nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed

// comm    (In) - an Epetra_Comm object describing the parallel machine (numProcs and my proc ID)
// map    (Out) - Epetra_Map describing distribution of matrix and vectors/multivectors
// A      (Out) - Epetra_CrsMatrix constructed for nx by ny grid using prescribed stencil
//                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
//                diagonal will be slightly diag dominant.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// bt     (Out) - Generated RHS.  Values satisfy b = A'*xexact
// xexact (Out) - Generated exact solution to Ax = b and b' = A'xexact

// Note: Caller of this function is responsible for deleting all output objects.

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact, bool StaticProfile, bool MakeLocalOnly) {

  Epetra_MultiVector * b1, * bt1, * xexact1;
	
  GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints, 
		     xoff, yoff, 1, comm, verbose, summary, 
		     map, A, b1, bt1, xexact1, StaticProfile, MakeLocalOnly);

  b = dynamic_cast<Epetra_Vector *>(b1);
  bt = dynamic_cast<Epetra_Vector *>(bt1);
  xexact = dynamic_cast<Epetra_Vector *>(xexact1);

  return;
}

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, int nrhs,
			const Epetra_Comm  &comm, bool verbose, bool summary,
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact, bool StaticProfile, bool MakeLocalOnly) {
  
  Epetra_Time timer(comm);
  // Determine my global IDs
  long long * myGlobalElements;
  GenerateMyGlobalElements(numNodesX, numNodesY, numProcsX, numProcsY, comm.MyPID(), myGlobalElements);

  int numMyEquations = numNodesX*numNodesY;
  
  map = new Epetra_Map((long long)-1, numMyEquations, myGlobalElements, 0, comm); // Create map with 2D block partitioning.
  delete [] myGlobalElements;

  long long numGlobalEquations = map->NumGlobalElements64();

  int profile = 0; if (StaticProfile) profile = numPoints;

#ifdef EPETRA_HAVE_STATICPROFILE

  if (MakeLocalOnly) 
    A = new Epetra_CrsMatrix(Copy, *map, *map, profile, StaticProfile); // Construct matrix with rowmap=colmap
  else
    A = new Epetra_CrsMatrix(Copy, *map, profile, StaticProfile); // Construct matrix

#else

  if (MakeLocalOnly) 
    A = new Epetra_CrsMatrix(Copy, *map, *map, profile); // Construct matrix with rowmap=colmap
  else
    A = new Epetra_CrsMatrix(Copy, *map, profile); // Construct matrix

#endif

  long long * indices = new long long[numPoints];
  double * values = new double[numPoints];

  double dnumPoints = (double) numPoints;
  int nx = numNodesX*numProcsX;

  for (int i=0; i<numMyEquations; i++) {

    long long rowID = map->GID64(i);
    int numIndices = 0;

    for (int j=0; j<numPoints; j++) {
      long long colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numGlobalEquations) {
	indices[numIndices] = colID;
	double value = - ((double) rand())/ ((double) RAND_MAX);
	if (colID==rowID)
	  values[numIndices++] = dnumPoints - value; // Make diagonal dominant
	else
	  values[numIndices++] = value;
      }
    }
    //cout << "Building row " << rowID << endl;
    A->InsertGlobalValues(rowID, numIndices, values, indices);
  }

  delete [] indices;
  delete [] values;
  double insertTime = timer.ElapsedTime();
  timer.ResetStartTime();
  A->FillComplete(false);
  double fillCompleteTime = timer.ElapsedTime();

  if (verbose)
    cout << "Time to insert matrix values = " << insertTime << endl
	 << "Time to complete fill        = " << fillCompleteTime << endl;
  if (summary) {
    if (comm.NumProc()==1) cout << "InsertTime" << '\t';
    cout << insertTime << endl;
    if (comm.NumProc()==1) cout << "FillCompleteTime" << '\t';
    cout << fillCompleteTime << endl;
  }

  if (nrhs<=1) {  
    b = new Epetra_Vector(*map);
    bt = new Epetra_Vector(*map);
    xexact = new Epetra_Vector(*map);
  }
  else {
    b = new Epetra_MultiVector(*map, nrhs);
    bt = new Epetra_MultiVector(*map, nrhs);
    xexact = new Epetra_MultiVector(*map, nrhs);
  }

  xexact->Random(); // Fill xexact with random values

  A->Multiply(false, *xexact, *b);
  A->Multiply(true, *xexact, *bt);

  return;
}


// Constructs a 2D PDE finite difference matrix using the list of x and y offsets.
// 
// nx      (In) - number of grid points in x direction
// ny      (In) - number of grid points in y direction
//   The total number of equations will be nx*ny ordered such that the x direction changes
//   most rapidly: 
//      First equation is at point (0,0)
//      Second at                  (1,0)
//       ...
//      nx equation at             (nx-1,0)
//      nx+1st equation at         (0,1)

// numPoints (In) - number of points in finite difference stencil
// xoff    (In) - stencil offsets in x direction (of length numPoints)
// yoff    (In) - stencil offsets in y direction (of length numPoints)
//   A standard 5-point finite difference stencil would be described as:
//     numPoints = 5
//     xoff = [-1, 1, 0,  0, 0]
//     yoff = [ 0, 0, 0, -1, 1]

// nsizes  (In) - Length of element size list used to create variable block map and matrix
// sizes   (In) - integer list of element sizes of length nsizes
//    The map associated with this VbrMatrix will be created by cycling through the sizes list.
//    For example, if nsize = 3 and sizes = [ 2, 4, 3], the block map will have elementsizes
//    of 2, 4, 3, 2, 4, 3, ...

// nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed

// comm    (In) - an Epetra_Comm object describing the parallel machine (numProcs and my proc ID)
// map    (Out) - Epetra_Map describing distribution of matrix and vectors/multivectors
// A      (Out) - Epetra_VbrMatrix constructed for nx by ny grid using prescribed stencil
//                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
//                diagonal will be slightly diag dominant.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// bt     (Out) - Generated RHS.  Values satisfy b = A'*xexact
// xexact (Out) - Generated exact solution to Ax = b and b' = A'xexact

// Note: Caller of this function is responsible for deleting all output objects.

void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			int nsizes, int * sizes,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact, bool StaticProfile, bool MakeLocalOnly) {
	
  Epetra_MultiVector * b1, * bt1, * xexact1;
	
  GenerateVbrProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints,
		     xoff, yoff, nsizes, sizes,
		     1, comm, verbose, summary, map, A, b1, bt1, xexact1, StaticProfile, MakeLocalOnly);

  b = dynamic_cast<Epetra_Vector *>(b1);
  bt = dynamic_cast<Epetra_Vector *>(bt1);
  xexact = dynamic_cast<Epetra_Vector *>(xexact1);

  return;
}

void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, 
			int nsizes, int * sizes, int nrhs,
			const Epetra_Comm  &comm, bool verbose, bool summary, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact, bool StaticProfile, bool MakeLocalOnly) {

  int i;

  // Determine my global IDs
  long long * myGlobalElements;
  GenerateMyGlobalElements(numNodesX, numNodesY, numProcsX, numProcsY, comm.MyPID(), myGlobalElements);

  int numMyElements = numNodesX*numNodesY;
  
  Epetra_Map ptMap((long long)-1, numMyElements, myGlobalElements, 0, comm); // Create map with 2D block partitioning.
  delete [] myGlobalElements;

  Epetra_IntVector elementSizes(ptMap); // This vector will have the list of element sizes
  for (i=0; i<numMyElements; i++) 
    elementSizes[i] = sizes[ptMap.GID64(i)%nsizes]; // cycle through sizes array

  map = new Epetra_BlockMap((long long)-1, numMyElements, ptMap.MyGlobalElements64(), elementSizes.Values(),
			    ptMap.IndexBase64(), ptMap.Comm());

  int profile = 0; if (StaticProfile) profile = numPoints;
  
// FIXME: Won't compile until Epetra_VbrMatrix is modified.
#if 0
  int j;
  long long numGlobalEquations = ptMap.NumGlobalElements64();

  if (MakeLocalOnly) 
    A = new Epetra_VbrMatrix(Copy, *map, *map, profile); // Construct matrix rowmap=colmap
  else
    A = new Epetra_VbrMatrix(Copy, *map, profile); // Construct matrix

  long long * indices = new long long[numPoints];

  // This section of code creates a vector of random values that will be used to create
  // light-weight dense matrices to pass into the VbrMatrix construction process.

  int maxElementSize = 0;
  for (i=0; i< nsizes; i++) maxElementSize = EPETRA_MAX(maxElementSize, sizes[i]);

  Epetra_LocalMap lmap((long long)maxElementSize*maxElementSize, ptMap.IndexBase(), ptMap.Comm());
  Epetra_Vector randvec(lmap);
  randvec.Random();
  randvec.Scale(-1.0); // Make value negative
  int nx = numNodesX*numProcsX;


  for (i=0; i<numMyElements; i++) {
    long long rowID = map->GID64(i);
    int numIndices = 0;
    int rowDim = sizes[rowID%nsizes];
    for (j=0; j<numPoints; j++) {
      long long colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numGlobalEquations)
	indices[numIndices++] = colID;
    }
			
    A->BeginInsertGlobalValues(rowID, numIndices, indices);
		
    for (j=0; j < numIndices; j++) {
      int colDim = sizes[indices[j]%nsizes];
      A->SubmitBlockEntry(&(randvec[0]), rowDim, rowDim, colDim);
    }
    A->EndSubmitEntries();
  }

  delete [] indices;

  A->FillComplete();

  // Compute the InvRowSums of the matrix rows
  Epetra_Vector invRowSums(A->RowMap());
  Epetra_Vector rowSums(A->RowMap());
  A->InvRowSums(invRowSums);
  rowSums.Reciprocal(invRowSums);

  // Jam the row sum values into the diagonal of the Vbr matrix (to make it diag dominant)
  int numBlockDiagonalEntries;
  int * rowColDims;
  int * diagoffsets = map->FirstPointInElementList();
  A->BeginExtractBlockDiagonalView(numBlockDiagonalEntries, rowColDims);
  for (i=0; i< numBlockDiagonalEntries; i++) {
    double * diagVals;
    int diagLDA;
    A->ExtractBlockDiagonalEntryView(diagVals, diagLDA);
    int rowDim = map->ElementSize(i);
    for (j=0; j<rowDim; j++) diagVals[j+j*diagLDA] = rowSums[diagoffsets[i]+j];
  }

  if (nrhs<=1) {  
    b = new Epetra_Vector(*map);
    bt = new Epetra_Vector(*map);
    xexact = new Epetra_Vector(*map);
  }
  else {
    b = new Epetra_MultiVector(*map, nrhs);
    bt = new Epetra_MultiVector(*map, nrhs);
    xexact = new Epetra_MultiVector(*map, nrhs);
  }

  xexact->Random(); // Fill xexact with random values

  A->Multiply(false, *xexact, *b);
  A->Multiply(true, *xexact, *bt);

#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

  return;
}

void GenerateMyGlobalElements(int numNodesX, int numNodesY, int numProcsX, int numProcs,
			      int myPID, long long * & myGlobalElements) {

  myGlobalElements = new long long[numNodesX*numNodesY];
  int myProcX = myPID%numProcsX;
  int myProcY = myPID/numProcsX;
  int curGID = myProcY*(numProcsX*numNodesX)*numNodesY+myProcX*numNodesX;
  for (int j=0; j<numNodesY; j++) {
    for (int i=0; i<numNodesX; i++) {
      myGlobalElements[j*numNodesX+i] = curGID+i;
    }
    curGID+=numNodesX*numProcsX;
  }
  //for (int i=0; i<numNodesX*numNodesY; i++) cout << "MYPID " << myPID <<" GID "<< myGlobalElements[i] << endl;
  
  return;
}

void runMatrixTests(Epetra_CrsMatrix * A,  Epetra_MultiVector * b, Epetra_MultiVector * bt,
		    Epetra_MultiVector * xexact, bool StaticProfile, bool verbose, bool summary) {

  Epetra_MultiVector z(*b);
  Epetra_MultiVector r(*b);
  Epetra_SerialDenseVector resvec(b->NumVectors());

  //Timings
  Epetra_Flops flopcounter;
  A->SetFlopCounter(flopcounter);
  Epetra_Time timer(A->Comm());
  std::string statdyn =        "dynamic";
  if (StaticProfile) statdyn = "static ";

  for (int j=0; j<4; j++) { // j = 0/2 is notrans, j = 1/3 is trans
    
    bool TransA = (j==1 || j==3);
    std::string contig = "without";
    if (j>1) contig =    "with   ";
    
#ifdef EPETRA_SHORT_PERFTEST
    int kstart = 1;
#else
    int kstart = 0;
#endif
    for (int k=kstart; k<2; k++) { // Loop over old multiply vs. new multiply
      
      std::string oldnew = "old";
      if (k>0) oldnew =    "new";

      if (j==2) A->OptimizeStorage();

      flopcounter.ResetFlops();
      timer.ResetStartTime();

      if (k==0) {
	//10 matvecs
#ifndef EPETRA_SHORT_PERFTEST
	for( int i = 0; i < 10; ++i )
	  A->Multiply1(TransA, *xexact, z); // Compute z = A*xexact or z = A'*xexact using old Multiply method
#endif
      }
      else {
	//10 matvecs
	for( int i = 0; i < 10; ++i )
	  A->Multiply(TransA, *xexact, z); // Compute z = A*xexact or z = A'*xexact
      }
      
      double elapsed_time = timer.ElapsedTime();
      double total_flops = A->Flops();

      // Compute residual
      if (TransA)
	r.Update(-1.0, z, 1.0, *bt, 0.0); // r = bt - z
      else
	r.Update(-1.0, z, 1.0, *b, 0.0); // r = b - z

      r.Norm2(resvec.Values());
      
      if (verbose) cout << "ResNorm = " << resvec.NormInf() << ": ";
      double MFLOPs = total_flops/elapsed_time/1000000.0;
      if (verbose) cout << "Total MFLOPs for 10 " << oldnew << " MatVec's with " << statdyn << " Profile (Trans = " << TransA
			<< ")  and " << contig << " optimized storage = " << MFLOPs << " (" << elapsed_time << " s)" <<endl;
      if (summary) {
	if (A->Comm().NumProc()==1) {
	  if (TransA) cout << "TransMv" << statdyn<< "Prof" << contig << "OptStor" << '\t';
	  else cout << "NoTransMv" << statdyn << "Prof" << contig << "OptStor" << '\t';
	}
	cout << MFLOPs << endl;
      }
    }
  }
  return;
}
#ifdef EPETRA_HAVE_JADMATRIX
void runJadMatrixTests(Epetra_JadMatrix * A,  Epetra_MultiVector * b, Epetra_MultiVector * bt,
		    Epetra_MultiVector * xexact, bool StaticProfile, bool verbose, bool summary) {

  Epetra_MultiVector z(*b);
  Epetra_MultiVector r(*b);
  Epetra_SerialDenseVector resvec(b->NumVectors());

  //Timings
  Epetra_Flops flopcounter;
  A->SetFlopCounter(flopcounter);
  Epetra_Time timer(A->Comm());

  for (int j=0; j<2; j++) { // j = 0 is notrans, j = 1 is trans
    
    bool TransA = (j==1);
    A->SetUseTranspose(TransA);
    flopcounter.ResetFlops();
    timer.ResetStartTime();

    //10 matvecs
    for( int i = 0; i < 10; ++i )
      A->Apply(*xexact, z); // Compute z = A*xexact or z = A'*xexact
    
    double elapsed_time = timer.ElapsedTime();
    double total_flops = A->Flops();
    
    // Compute residual
    if (TransA)
      r.Update(-1.0, z, 1.0, *bt, 0.0); // r = bt - z
    else
      r.Update(-1.0, z, 1.0, *b, 0.0); // r = b - z
    
    r.Norm2(resvec.Values());
    
    if (verbose) cout << "ResNorm = " << resvec.NormInf() << ": ";
    double MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Total MFLOPs for 10 " << " Jagged Diagonal MatVec's with (Trans = " << TransA
		      << ") " << MFLOPs << " (" << elapsed_time << " s)" <<endl;
    if (summary) {
      if (A->Comm().NumProc()==1) {
	if (TransA) cout << "TransMv" << '\t';
	else cout << "NoTransMv" << '\t';
      }
      cout << MFLOPs << endl;
    }
  }
  return;
}
#endif
//=========================================================================================
void runLUMatrixTests(Epetra_CrsMatrix * L,  Epetra_MultiVector * bL, Epetra_MultiVector * btL, Epetra_MultiVector * xexactL, 
		      Epetra_CrsMatrix * U,  Epetra_MultiVector * bU, Epetra_MultiVector * btU, Epetra_MultiVector * xexactU, 
		      bool StaticProfile, bool verbose, bool summary) {

  if (L->NoDiagonal()) {
    bL->Update(1.0, *xexactL, 1.0); // Add contribution of a unit diagonal to bL
    btL->Update(1.0, *xexactL, 1.0); // Add contribution of a unit diagonal to btL
  }
  if (U->NoDiagonal()) {
    bU->Update(1.0, *xexactU, 1.0); // Add contribution of a unit diagonal to bU
    btU->Update(1.0, *xexactU, 1.0); // Add contribution of a unit diagonal to btU
  }

  Epetra_MultiVector z(*bL);
  Epetra_MultiVector r(*bL);
  Epetra_SerialDenseVector resvec(bL->NumVectors());

  //Timings
  Epetra_Flops flopcounter;
  L->SetFlopCounter(flopcounter);
  U->SetFlopCounter(flopcounter);
  Epetra_Time timer(L->Comm());
  std::string statdyn =        "dynamic";
  if (StaticProfile) statdyn = "static ";

  for (int j=0; j<4; j++) { // j = 0/2 is notrans, j = 1/3 is trans
    
    bool TransA = (j==1 || j==3);
    std::string contig = "without";
    if (j>1) contig =    "with   ";
    
    if (j==2) {
      L->OptimizeStorage();
      U->OptimizeStorage();
    }

    flopcounter.ResetFlops();
    timer.ResetStartTime();
    
    //10 lower solves
    bool Upper = false;
    bool UnitDiagonal = L->NoDiagonal();  // If no diagonal, then unit must be used
    Epetra_MultiVector * b = TransA ? btL : bL;  // solve with the appropriate b vector
    for( int i = 0; i < 10; ++i )
      L->Solve(Upper, TransA, UnitDiagonal, *b, z); // Solve Lz = bL or L'z = bLt
      
    double elapsed_time = timer.ElapsedTime();
    double total_flops = L->Flops();

    // Compute residual
    r.Update(-1.0, z, 1.0, *xexactL, 0.0); // r = bt - z
    r.Norm2(resvec.Values());

    if (resvec.NormInf()>0.000001) {
      cout << "resvec = " << resvec << endl;
      cout << "z = " << z << endl;
      cout << "xexactL = " << *xexactL << endl;
      cout << "r = " << r << endl;
    }
      
    if (verbose) cout << "ResNorm = " << resvec.NormInf() << ": ";
    double MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Total MFLOPs for 10 " << " Lower solves " << statdyn << " Profile (Trans = " << TransA
		      << ")  and " << contig << " opt storage = " << MFLOPs << " (" << elapsed_time << " s)" <<endl;
    if (summary) {
      if (L->Comm().NumProc()==1) {
	if (TransA) cout << "TransLSv" << statdyn<< "Prof" << contig << "OptStor" << '\t';
	else cout << "NoTransLSv" << statdyn << "Prof" << contig << "OptStor" << '\t';
      }
      cout << MFLOPs << endl;
    }
    flopcounter.ResetFlops();
    timer.ResetStartTime();
    
    //10 upper solves
    Upper = true;
    UnitDiagonal = U->NoDiagonal();  // If no diagonal, then unit must be used
    b = TransA ? btU : bU;  // solve with the appropriate b vector
    for( int i = 0; i < 10; ++i )
      U->Solve(Upper, TransA, UnitDiagonal, *b, z); // Solve Lz = bL or L'z = bLt
      
    elapsed_time = timer.ElapsedTime();
    total_flops = U->Flops();

    // Compute residual
    r.Update(-1.0, z, 1.0, *xexactU, 0.0); // r = bt - z
    r.Norm2(resvec.Values());

    if (resvec.NormInf()>0.001) {
      cout << "U = " << *U << endl;
      //cout << "resvec = " << resvec << endl;
      cout << "z = " << z << endl;
      cout << "xexactU = " << *xexactU << endl;
      //cout << "r = " << r << endl;
      cout << "b = " << *b << endl;
    }

      
    if (verbose) cout << "ResNorm = " << resvec.NormInf() << ": ";
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Total MFLOPs for 10 " << " Upper solves " << statdyn << " Profile (Trans = " << TransA
		      << ")  and " << contig << " opt storage = " << MFLOPs <<endl;
    if (summary) {
      if (L->Comm().NumProc()==1) {
	if (TransA) cout << "TransUSv" << statdyn<< "Prof" << contig << "OptStor" << '\t';
	else cout << "NoTransUSv" << statdyn << "Prof" << contig << "OptStor" << '\t';
      }
      cout << MFLOPs << endl;
    }
  }
  return;
}
