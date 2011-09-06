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

#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

void BiCGSTAB(Epetra_CrsMatrix &A, bool transA, Epetra_Vector &x, Epetra_Vector &b, 
	      Ifpack_CrsRiluk *M, 
	      int Maxiter, double Tolerance, 
	      double *residual, bool verbose);

  string toString(const int& x) {
     char s[100];
     sprintf(s, "%d", x);
     return string(s);
}

  string toString(const double& x) {
     char s[100];
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

  // Construct ILU preconditioner

  double elapsed_time, total_flops, MFLOPs;
  Epetra_Time timer(Comm);

  int LevelFill = 0;
  if (argc > 2)  LevelFill = atoi(argv[2]);
  if (verbose) cout << "Using Level Fill = " << LevelFill << endl;
  int Overlap = 0;
  if (argc > 3) Overlap = atoi(argv[3]);
  if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  if (argc > 4) Athresh = atof(argv[4]);
  if (verbose) cout << "Using Absolute Threshold Value of = " << Athresh << endl;

  double Rthresh = 1.0;
  if (argc > 5) Rthresh = atof(argv[5]);
  if (verbose) cout << "Using Relative Threshold Value of = " << Rthresh << endl;

  Ifpack_IlukGraph * IlukGraph = 0;
  Ifpack_CrsRiluk * ILUK = 0;

  if (LevelFill>-1) {
    elapsed_time = timer.ElapsedTime();
    IlukGraph = new Ifpack_IlukGraph(A.Graph(), LevelFill, Overlap);
    assert(IlukGraph->ConstructFilledGraph()==0);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    if (verbose) cout << "Time to construct ILUK graph = " << elapsed_time << endl;


    Epetra_Flops fact_counter;
  
    elapsed_time = timer.ElapsedTime();
    ILUK = new Ifpack_CrsRiluk(*IlukGraph);
    ILUK->SetFlopCounter(fact_counter);
    ILUK->SetAbsoluteThreshold(Athresh);
    ILUK->SetRelativeThreshold(Rthresh);
    //assert(ILUK->InitValues()==0);
    int initerr = ILUK->InitValues(A);
    if (initerr!=0) cout << Comm << "InitValues error = " << initerr;
    assert(ILUK->Factor()==0);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    total_flops = ILUK->Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Time to compute preconditioner values = " 
		    << elapsed_time << endl
		    << "MFLOPS for Factorization = " << MFLOPs << endl;
    //cout << *ILUK << endl;
  }
  bool transA = true;
  double Condest;
  ILUK->Condest(transA, Condest);

  if (verbose) cout << "Condition number estimate for this preconditioner = " << Condest << endl;

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsRiluk Preconditioner: LevelFill = " + toString(LevelFill) + 
                                                 " Overlap = " + toString(Overlap) + 
                                                 " Athresh = " + toString(Athresh) + 
                                                 " Rthresh = " + toString(Rthresh); 
  ILUK->SetLabel(label.c_str());

  //int Maxiter = 500;
  //double Tolerance = 1.0E-14;

  Epetra_Vector xcomp(map);
  Epetra_Vector resid(map);

  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  xcomp.SetFlopCounter(A);
  b.SetFlopCounter(A);
  resid.SetFlopCounter(A);
  ILUK->SetFlopCounter(A);

  if (transA) A.Multiply(transA,xexact,b); // Compute b to be A'*xexact

  elapsed_time = timer.ElapsedTime();
  //BiCGSTAB(A, transA, xcomp, b, ILUK, Maxiter, Tolerance, &residual, verbose);
  cout << "Building AztecOO solver" << endl;

  AztecOO solver;
  A.SetUseTranspose(transA); // solve for A transpose if transA is true
  solver.SetUserMatrix(&A);
  solver.SetLHS(&xcomp);
  solver.SetRHS(&b);
  ILUK->SetUseTranspose(transA); // solve for A transpose if transA is true
  solver.SetPrecOperator(ILUK);


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


  


  if (ILUK!=0) delete ILUK;
  if (IlukGraph!=0) delete IlukGraph;
				       
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
void BiCGSTAB(Epetra_CrsMatrix &A, bool transA,
	      Epetra_Vector &x, 
	      Epetra_Vector &b, 
	      Ifpack_CrsRiluk *M, 
	      int Maxiter, 
	      double Tolerance, 
	      double *residual, bool verbose) {

  // Allocate vectors needed for iterations
  Epetra_Vector phat(x.Map()); phat.SetFlopCounter(x);
  Epetra_Vector p(x.Map()); p.SetFlopCounter(x);
  Epetra_Vector shat(x.Map()); shat.SetFlopCounter(x);
  Epetra_Vector s(x.Map()); s.SetFlopCounter(x);
  Epetra_Vector r(x.Map()); r.SetFlopCounter(x);
  Epetra_Vector rtilde(x.Map()); rtilde.Random(); rtilde.SetFlopCounter(x);
  Epetra_Vector v(x.Map()); v.SetFlopCounter(x);
  

  A.Multiply(transA, x, r); // r = A*x

  r.Update(1.0, b, -1.0); // r = b - A*x

  double r_norm, b_norm, scaled_r_norm, rhon, rhonm1 = 1.0;
  double alpha = 1.0, omega = 1.0, sigma;
  double omega_num, omega_den;
  r.Norm2(&r_norm);
  b.Norm2(&b_norm);
  scaled_r_norm = r_norm/b_norm;
  r.Dot(rtilde,&rhon);

  if (verbose) cout << "Initial residual = " << r_norm
		    << " Scaled residual = " << scaled_r_norm << endl;


  for (int i=0; i<Maxiter; i++) { // Main iteration loop   

    double beta = (rhon/rhonm1) * (alpha/omega);
    rhonm1 = rhon;

    /* p    = r + beta*(p - omega*v)       */
    /* phat = M^-1 p                       */
    /* v    = A phat                       */

    double dtemp = - beta*omega;

    p.Update(1.0, r, dtemp, v, beta);
    if (M==0) 
      phat.Scale(1.0, p);
    else
      M->Solve(transA, p, phat);
    A.Multiply(transA, phat, v);

    
    rtilde.Dot(v,&sigma);
    alpha = rhon/sigma;    

    /* s = r - alpha*v                     */
    /* shat = M^-1 s                       */
    /* r = A shat (r is a tmp here for t ) */

    s.Update(-alpha, v, 1.0, r, 0.0);
    if (M==0) 
      shat.Scale(1.0, s);
    else
      M->Solve(transA, s, shat);
    A.Multiply(transA, shat, r);

    r.Dot(s, &omega_num);
    r.Dot(r, &omega_den);
    omega = omega_num / omega_den;

    /* x = x + alpha*phat + omega*shat */
    /* r = s - omega*r */

    x.Update(alpha, phat, omega, shat, 1.0);
    r.Update(1.0, s, -omega); 
    
    r.Norm2(&r_norm);
    scaled_r_norm = r_norm/b_norm;
    r.Dot(rtilde,&rhon);

    if (verbose) cout << "Iter "<< i << " residual = " << r_norm
		      << " Scaled residual = " << scaled_r_norm << endl;

    if (scaled_r_norm < Tolerance) break;
  }
  return;
}
