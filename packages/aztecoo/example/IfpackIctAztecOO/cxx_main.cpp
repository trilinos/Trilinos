//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "AztecOO.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"

#include "Epetra_CrsMatrix.h"
#include "Ifpack_CrsIct.h"

  string toString(const int& x) {
     char s[10000];
     sprintf(s, "%d", x);
     return string(s);
}

  string toString(const double& x) {
     char s[10000];
     sprintf(s, "%g", x);
     return string(s);
}


int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  cout << Comm << endl;

  int MyPID = Comm.MyPID();

  bool verbose = false; 
  if (MyPID==0) verbose = true;

  if(argc < 2 && verbose) {
    cerr << "Usage: " << argv[0] 
	 << " HB_filename [level_fill [level_overlap [absolute_threshold [ relative_threshold]]]]" << endl
	 << "where:" << endl
	 << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
	 << "level_fill         - The amount of fill to use for ILU(k) preconditioner (default 0)" << endl
	 << "level_overlap      - The amount of overlap used for overlapping Schwarz subdomains (default 0)" << endl
	 << "absolute_threshold - The minimum value to place on the diagonal prior to factorization (default 0.0)" << endl
	 << "relative_threshold - The relative amount to perturb the diagonal prior to factorization (default 1.0)" << endl << endl
	 << "To specify a non-default value for one of these parameters, you must specify all" << endl
	 << " preceding values but not any subsequent parameters. Example:" << endl
	 << "ifpackHbSerialMsr.exe mymatrix.hb 1  - loads mymatrix.hb, uses level fill of one, all other values are defaults" << endl
	 << endl;
    return(1);

  }

  // Uncomment the next line to debug in mpi mode
  // int tmp; if (MyPID==0) {cout << "Type a character and Enter to continue" <<endl; cin >> tmp; } Comm.Barrier();

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra(argv[1], Comm, readMap, readA, readx, readb, readxexact);

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);

  // Create Exporter to distribute read-in matrix and vectors

  Epetra_Export exporter(*readMap, map);
  Epetra_CrsMatrix A(Copy, map, 0);
  Epetra_Vector x(map);
  Epetra_Vector b(map);
  Epetra_Vector xexact(map);

  Epetra_Time FillTimer(Comm);
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);
  Comm.Barrier();
  double vectorRedistributeTime = FillTimer.ElapsedTime();
  A.Export(*readA, exporter, Add);
  Comm.Barrier();
  double matrixRedistributeTime = FillTimer.ElapsedTime() - vectorRedistributeTime;
  assert(A.FillComplete()==0);    
  Comm.Barrier();
  double fillCompleteTime = FillTimer.ElapsedTime() - matrixRedistributeTime;
  if (Comm.MyPID()==0)	{
    cout << "\n\n****************************************************" << endl;
    cout << "\n Vector redistribute  time (sec) = " << vectorRedistributeTime<< endl;
    cout << "    Matrix redistribute time (sec) = " << matrixRedistributeTime << endl;
    cout << "    Transform to Local  time (sec) = " << fillCompleteTime << endl<< endl;
  }
  Epetra_Vector tmp1(*readMap);
  Epetra_Vector tmp2(map);
  readA->Multiply(false, *readxexact, tmp1);

  A.Multiply(false, xexact, tmp2);
  double residual;
  tmp1.Norm2(&residual);
  if (verbose) cout << "Norm of Ax from file            = " << residual << endl;
  tmp2.Norm2(&residual);
  if (verbose) cout << "Norm of Ax after redistribution = " << residual << endl << endl << endl;

  //cout << "A from file = " << *readA << endl << endl << endl;

  //cout << "A after dist = " << A << endl << endl << endl;

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;

  Comm.Barrier();

  // Construct ICT preconditioner

  double elapsed_time, total_flops, MFLOPs;
  Epetra_Time timer(Comm);

  int Lfil = 0;
  if (argc > 2)  Lfil = atoi(argv[2]);
  if (verbose) cout << "Using Lfil = " << Lfil << endl;
  int Overlap = 0;
  if (argc > 3) Overlap = atoi(argv[3]);
  if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  if (argc > 4) Athresh = atof(argv[4]);
  if (verbose) cout << "Using Absolute Threshold Value of = " << Athresh << endl;

  double Rthresh = 1.0;
  if (argc > 5) Rthresh = atof(argv[5]);
  if (verbose) cout << "Using Relative Threshold Value of = " << Rthresh << endl;

  Ifpack_CrsIct * ICT = 0;

  if (Lfil>-1) {


    Epetra_Flops fact_counter;
  
    elapsed_time = timer.ElapsedTime();
    ICT = new Ifpack_CrsIct(A, 1.0e-6, 100);
    ICT->SetFlopCounter(fact_counter);
    ICT->SetAbsoluteThreshold(Athresh);
    ICT->SetRelativeThreshold(Rthresh);
    //assert(ICT->InitValues()==0);
    int initerr = ICT->InitValues(A);
    if (initerr!=0) cout << Comm << "InitValues error = " << initerr;
    assert(ICT->Factor()==0);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    total_flops = ICT->Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Time to compute preconditioner values = " 
		    << elapsed_time << endl
		    << "MFLOPS for Factorization = " << MFLOPs << endl;
    //cout << *ICT << endl;
  }
  bool transA = false;
  double Condest;
  ICT->Condest(transA, Condest);

  if (verbose) cout << "Condition number estimate for this preconditioner = " << Condest << endl;

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsIct Preconditioner: Lfil = " + toString(Lfil) + 
                                                 " Overlap = " + toString(Overlap) + 
                                                 " Athresh = " + toString(Athresh) + 
                                                 " Rthresh = " + toString(Rthresh); 
  ICT->SetLabel(label.c_str());

  //int Maxiter = 500;
  //double Tolerance = 1.0E-14;

  Epetra_Vector xcomp(map);
  Epetra_Vector resid(map);

  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  xcomp.SetFlopCounter(A);
  b.SetFlopCounter(A);
  resid.SetFlopCounter(A);
  ICT->SetFlopCounter(A);

  if (transA) A.Multiply(transA,xexact,b); // Compute b to be A'*xexact

  elapsed_time = timer.ElapsedTime();
  cout << "Building AztecOO solver" << endl;

  AztecOO solver;
  A.SetUseTranspose(transA); // solve for A transpose if transA is true
  solver.SetUserMatrix(&A);
  solver.SetLHS(&xcomp);
  solver.SetRHS(&b);
  ICT->SetUseTranspose(transA); // solve for A transpose if transA is true
  solver.SetPrecOperator(ICT);


  // Assert symmetric
  // problem->AssertSymmetric();

  // Set Problem Difficulty Level
  //problem->SetPDL(easy);

  //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  //solver.SetAztecOption(AZ_precond, AZ_ls);
  //solver.SetAztecOption(AZ_scaling, 8);
  //solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
  //solver.SetAztecOption(AZ_output, 0);
  //solver.SetAztecOption(AZ_graph_fill, 2);
  //solver.SetAztecOption(AZ_overlap, 1);
  //solver.SetAztecOption(AZ_poly_ord, 9);
  //solver.SetAztecParam(AZ_ilut_fill, 5.0);
  //solver.SetAztecParam(AZ_drop, 0.0);
  //double rthresh = 1.0;
  //cout << "Rel threshold = " << rthresh << endl;
  //solver.SetAztecParam(AZ_rthresh, rthresh);
  //double athresh = 1.0e-2;
  //cout << "Abs threshold = " << athresh << endl;
  //solver.SetAztecParam(AZ_athresh, athresh);

  

  //solver.SetAztecOption(AZ_reorder, 2);

  int Niters = 1200;
  solver.SetAztecOption(AZ_solver, AZ_cg);

  solver.SetAztecOption(AZ_kspace, Niters);
   
  solver.Iterate(Niters, 5.0e-10);

  elapsed_time = timer.ElapsedTime() - elapsed_time;
  total_flops = counter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "Time to compute solution = " 
		    << elapsed_time << endl
		    << "Number of operations in solve = " << total_flops << endl
		    << "MFLOPS for Solve = " << MFLOPs<< endl << endl;

  Epetra_Vector bcomp(map);
  assert(A.Multiply(transA, xcomp, bcomp)==0);
 
 
  assert(resid.Update(1.0, b, -1.0, bcomp, 0.0)==0);

  assert(resid.Norm2(&residual)==0);
  if (Comm.MyPID()==0) 
    cout << "Residual                                                    = " << residual << endl;


  resid.Update(1.0, xcomp, -1.0, xexact, 0.0); // resid = xcomp - xexact

  resid.Norm2(&residual);

  if (verbose) 
    cout << "Norm of the difference between exact and computed solutions = " << residual << endl;

  if (residual>1.0e-5) {
    cout << "Difference between computed and exact solution is large..." << endl      << "Computing norm of A times this difference.  If this norm is small, then matrix is singular"
      << endl;
    assert(A.Multiply(false, resid, bcomp)==0);
    assert(bcomp.Norm2(&residual)==0);
  if (Comm.MyPID()==0)
    cout << "2-norm of A times difference between computed and exact solution  = " << residual << endl;
  
  }


  


  if (ICT!=0) delete ICT;
				       
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
