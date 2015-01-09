/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#include "prototypes.h"
#include "paz_aztec.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_Time.h"
#include "Petra_BlockMap.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_DVBR_Matrix.h"
#include "Petra_RDP_CRS_Matrix.h"
#include "Ifpack_ILUK_Graph.h"
#include "Ifpack_RDP_CRS_RILUK.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }
#define perror1(str,ierr) { fprintf(stderr,"%s %d\n",str,ierr);   exit(-1); }
#define double_quote '"'
void BiCGSTAB(Petra_RDP_CRS_Matrix &A, 
	      Petra_RDP_Vector &x, 
	      Petra_RDP_Vector &b, 
	      Ifpack_RDP_CRS_RILUK *M, 
	      int Maxiter, 
	      double Tolerance, 
	      double *residual, double & FLOPS, bool verbose);

int main(int argc, char *argv[])
{
  int    *update;                  /* vector elements updated on this node. */
  int    *indx;   /* MSR format of real and imag parts */
  int    *bindx;
  int    *bpntr;
  int    *rpntr;
  int    *cpntr;
  int    indexBase = 0;
  double *val;
  double *xguess, *b, *bt, *xexact, *xsolve;
  int    n_nonzeros, n_blk_nonzeros, ierr;
  int    N_update;           /* # of block unknowns updated on this node    */
  int    numLocalEquations;
                                 /* Number scalar equations on this node */
  int    numGlobalEquations, numGlobalBlocks; /* Total number of equations */
  int    numLocalBlocks;
  int    *blockSizes, *numNzBlks, *blkColInds;
  int    *numNz, *ColInds;
  int    N_external, N_blk_eqns;
  int    blk_row, *blk_col_inds;
  int    row,     *col_inds, numEntries;
  double *row_vals;

  double *val_msr;
  int *bindx_msr;
  
  double norm, d ;

  int matrix_type;

  int has_global_indices, option;
  int i, j, m, mp;
  int ione = 1;
  int *proc_config = new int[PAZ_PROC_SIZE];

  double time ;
#ifdef PETRA_MPI
  MPI_Init(&argc,&argv);
#endif

  /* get number of processors and the name of this processor */
 
#ifdef PETRA_MPI
  Petra_Comm& Comm = *new Petra_Comm(MPI_COMM_WORLD);
#else
  Petra_Comm& Comm = *new Petra_Comm();
#endif

  printf("proc %d of %d is alive\n",
      Comm.MyPID(),Comm.NumProc());

  int MyPID = Comm.MyPID();

  bool verbose = false;
  if (MyPID==0) verbose = true;

  
/* Still need Aztec proc info for the HB routines, can switch to Petra later */

#ifdef PETRA_MPI
  PAZ_set_proc_config(proc_config,MPI_COMM_WORLD);
#else
  PAZ_set_proc_config(proc_config,0);
#endif

  Comm.Barrier();

#ifdef VBRMATRIX
  if(argc != 6) 
    perror("error: enter name of data and partition file on command line, followed by levelfill and shift value") ; 
#else
  if(argc != 5) perror("error: enter name of data file on command line, followed by levelfill and shift value") ; 
#endif
  /* Set exact solution to NULL */
  xexact = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

#ifdef VBRMATRIX
  if (Comm.MyPID()==0) 
    read_hb(argv[1], &numGlobalEquations, &n_nonzeros, 
	    &val_msr,  &bindx_msr, &xguess, &b, &bt, &xexact);
  
  create_vbr(argv[2], proc_config, &numGlobalEquations, &numGlobalBlocks,
	     &n_nonzeros, &n_blk_nonzeros, &N_update, &update,
	     bindx_msr, val_msr, &val, &indx, 
	     &rpntr, &cpntr, &bpntr, &bindx);

  if(proc_config[PAZ_node] == 0) 
    {
      free ((void *) val_msr);
      free ((void *) bindx_msr);
      free ((void *) cpntr);
    }
    matrix_type = PAZ_VBR_MATRIX;

  Comm.Barrier();

  distrib_vbr_matrix( proc_config, numGlobalEquations, numGlobalBlocks, 
		      &n_nonzeros, &n_blk_nonzeros,
		      &N_update, &update, 
		      &val, &indx, &rpntr, &cpntr, &bpntr, &bindx, 
		      &xguess, &b, &bt, &xexact);

#else
  if (Comm.MyPID()==0) 
    read_hb(argv[1], &numGlobalEquations, &n_nonzeros,
             &val,  &bindx, &xguess, &b, &bt, &xexact);

  Comm.Barrier();

  distrib_msr_matrix(proc_config, &numGlobalEquations, &n_nonzeros, &N_update,
		  &update, &val, &bindx, &xguess, &b, &bt, &xexact);

#ifdef DEBUG
  for (i = 0; i<N_update; i++)
    if (val[i] == 0.0 ) printf("Zero diagonal at row %d\n",i);
#endif
    matrix_type = PAZ_MSR_MATRIX;
#endif

#ifdef VBRMATRIX
  numLocalEquations = rpntr[N_update];
#else
  numLocalEquations = N_update;
#endif

#ifdef VBRMATRIX

  /********  Make integer data structures needed for this interface *********/

  /* Make blockSizes - number of equations in each block row on this proc */

  numLocalBlocks = N_update;

  blockSizes = new int[numLocalBlocks];
  for (i=0; i<numLocalBlocks; i++)
    blockSizes[i] = rpntr[i+1]-rpntr[i];

  /* Make numNzBlks - number of block entries in each block row */

  numNzBlks = new int[numLocalBlocks];
  for (i=0; i<numLocalBlocks; i++)
    numNzBlks[i] = bpntr[i+1] - bpntr[i];

  /* Make blkColInds - Exactly bindx (just copy pointer) */
  blkColInds = bindx;
  

  Petra_BlockMap& map = *new Petra_BlockMap(numGlobalEquations, numLocalEquations, 
			update, indexBase, Comm, numGlobalBlocks, numLocalBlocks,
			blockSizes);

  Petra_RDP_DVBR_Matrix& A = *new Petra_RDP_DVBR_Matrix(map);
  
  if ((ierr=A.allocate(numNzBlks, blkColInds)))
     perror1("Error in DVBR_Matrix_allocate:",ierr);
 
  /* Add block rows one-at-a-time */

  for (blk_row=0; blk_row<numLocalBlocks; blk_row++)
    {
      row_vals = val + indx[bpntr[blk_row]];
      blk_col_inds = bindx + bpntr[blk_row];
      if ((ierr=A.putBlockRow(update[blk_row], numNzBlks[blk_row], row_vals, 
				      blk_col_inds)))
       {
         printf("Row %d:",update[row]);
         perror1("Error putting block row:",ierr);
       }
    }
  
  if ((ierr=A.FillComplete()))
      perror1("Error in DVBR_Matrix_FillComplete:",ierr);

#else
  /* Make numNzBlks - number of block entries in each block row */

  numNz = new int[numLocalEquations];
  for (i=0; i<numLocalEquations; i++) numNz[i] = bindx[i+1] - bindx[i];

  /* Make ColInds - Exactly bindx, offset by diag (just copy pointer) */
  ColInds = bindx+numLocalEquations+1;

  Petra_Map& map = *new Petra_Map(numGlobalEquations, numLocalEquations, 
			update, 0, Comm);
 
  Comm.Barrier();
  Petra_Time & FillTimer = *new Petra_Time(Comm);
  Petra_RDP_CRS_Matrix& A = *new Petra_RDP_CRS_Matrix(Copy, map, numNz);
  
  /* Add  rows one-at-a-time */

  for (row=0; row<numLocalEquations; row++)
    {
      row_vals = val + bindx[row];
      col_inds = bindx + bindx[row];
      numEntries = bindx[row+1] - bindx[row];
     if ((ierr = A.InsertGlobalValues(update[row], numEntries, row_vals, col_inds)))
       {
         printf("Row %d:",update[row]);
         perror1("Error putting row:",ierr);
       }
     if ((ierr=(A.InsertGlobalValues(update[row], 1, val+row, update+row)<0)))
       perror1("Error putting  diagonal",ierr);
    }
   Comm.Barrier();
   double FillTime = FillTimer.ElapsedTime();
    if ((ierr=A.FillComplete()))    
    perror1("Error in Petra_RDP_Vector_fillComplete",ierr);
   double FillCompleteTime = FillTimer.ElapsedTime() - FillTime;
   if (Comm.MyPID()==0)	{
     cout << "\n\n****************************************************" << endl;
     cout << "\n\nMatrix construction time (sec) = " << FillTime<< endl;
     cout << "    Matrix FillComplete time (sec) = " << FillCompleteTime << endl;
     cout << "    Total construction time (sec) = " << FillTime+FillCompleteTime << endl<< endl;
   }

  
   
#endif
  
#ifdef MULTI_VECTOR

  // Make a second x and b vector that are 2x original x and b; make into a 2-vector.

  int nrhs = 2;
  double **allx = new double*[nrhs];
  double *xexact2 = new double[numLocalEquations];
  for (i = 0;i<numLocalEquations; i++) xexact2[i] = 2.0 * xexact[i];
  allx[0] = xexact; allx[1] = xexact2;

  double **allb = new double*[nrhs];
  double *b2 = new double[numLocalEquations];
  for (i = 0;i<numLocalEquations; i++) b2[i] = 2.0 * b[i];
  allb[0] = b; allb[1] = b2;

  double **allbt = new double*[nrhs];
  double *bt2 = new double[numLocalEquations];
  for (i = 0;i<numLocalEquations; i++) bt2[i] = 2.0 * bt[i];
  allbt[0] = bt; allbt[1] = bt2;

  Petra_RDP_MultiVector& xtrue = *new Petra_RDP_MultiVector(Copy, map, allx, nrhs);
  Petra_RDP_MultiVector& bexact = *new Petra_RDP_MultiVector(Copy, map, allb, nrhs);
  Petra_RDP_MultiVector& btexact = *new Petra_RDP_MultiVector(Copy, map, allbt, nrhs);
  Petra_RDP_MultiVector& bcomp = *new Petra_RDP_MultiVector(map, nrhs);
  Petra_RDP_MultiVector& xcomp = *new Petra_RDP_MultiVector(map, nrhs);
  Petra_RDP_MultiVector& resid = *new Petra_RDP_MultiVector(map, nrhs);

#else
  int nrhs = 1;
  Petra_RDP_Vector& xtrue = *new Petra_RDP_Vector(Copy, map, xexact);
  Petra_RDP_Vector& bexact = *new Petra_RDP_Vector(Copy, map, b);
  Petra_RDP_Vector& btexact = *new Petra_RDP_Vector(Copy, map, bt);
  Petra_RDP_Vector& bcomp = *new Petra_RDP_Vector(map);
  Petra_RDP_Vector& xcomp = *new Petra_RDP_Vector(map);
  Petra_RDP_Vector& resid = *new Petra_RDP_Vector(map);
    
#endif /* MULTI_VECTOR */

  Comm.Barrier();

  // Construct ILU preconditioner

  double elapsed_time, total_flops, MFLOPs;
  Petra_Time & timer = *new Petra_Time(Comm);

  int LevelFill = atoi(argv[argc-3]);
  if (verbose) cout << "Using Level Fill = " << LevelFill << endl;
  double ShiftValue = atof(argv[argc-2]);
  if (verbose) cout << "Using Diagonal Shift Value of = " << ShiftValue << endl;
  int NumThreads = atoi(argv[argc-1]);
  if (verbose) cout << "Using " << NumThreads << " Threads "  << endl;

  Ifpack_RDP_CRS_RILUK * ILUK = 0;
  Ifpack_ILUK_Graph * ILUK_Graph = 0;
  if (LevelFill>-1) {
    elapsed_time = timer.ElapsedTime();
    ILUK_Graph = new Ifpack_ILUK_Graph(A.Graph(), LevelFill);
    assert(ILUK_Graph->ConstructFilledGraph()==0);
    assert(ILUK_Graph->ComputeLevels(NumThreads)==0);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    if (verbose) cout << "Time to construct ILUK graph = " << elapsed_time << endl;
    
    return 0;
    
    elapsed_time = timer.ElapsedTime();
    ILUK = new Ifpack_RDP_CRS_RILUK(A, *ILUK_Graph);
    ILUK->SetShiftValue(ShiftValue);
    assert(ILUK->InitValues()==0);
    assert(ILUK->Factor()==0);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    total_flops = ILUK->Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Time to compute preconditioner values = " 
		      << elapsed_time << endl
		      << "MFLOPS for Factorization = " << MFLOPs << endl;
  }

  int Maxiter = 500;
  double Tolerance = 1.0E-14;
  double * residual = new double[nrhs];


  elapsed_time = timer.ElapsedTime();

  BiCGSTAB(A, xcomp, bexact, ILUK, Maxiter, Tolerance, residual, total_flops, verbose);

  elapsed_time = timer.ElapsedTime() - elapsed_time;
  MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "Time to compute solution = " 
		    << elapsed_time << endl 
		    << "Number of operations in solve = " << total_flops << endl
		    << "MFLOPS for Solve = " << MFLOPs<< endl << endl;
  

  int NumMVs = 100;
  int ii, iii;

  for (ii=0; ii<2; ii++) {
    bool TransA = (ii==1);
    A.ResetFlops();
    elapsed_time = timer.ElapsedTime();
    for (iii=0; iii<NumMVs; iii++)
      if ((ierr=A.Multiply(TransA, xcomp, bcomp))) perror1("Error in matvec",ierr);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    total_flops = A.Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (Comm.MyPID()==0) {
      if (TransA) {
	cout << "\n\n****************************************************" << endl;
	cout << "\n\nResults for transpose multiply with standard storage" << endl;
      }
      else {
	cout << "\n\nMatrix Fill cost = " << (FillTime+FillCompleteTime)/elapsed_time*NumMVs 
	     << " Matrix-Multiplies " << endl;
	cout << "\n\n*******************************************************" << endl;
	cout << "\n\nResults for no transpose multiply with standard storage" << endl;
      }
      
      cout << "\n\nMatrix Fill cost = " << (FillTime+FillCompleteTime)/elapsed_time*NumMVs 
	   << " Matrix-Multiplies " << endl;
      cout << "\n\nTotal FLOPS for standard Storage (" <<NumMVs<< " Multiplies) = " 
	   << total_flops<< endl;
      cout << "    Total time for standard Storage = " << elapsed_time<< endl;
      cout << "    Total MFLOPs for standard matrix multiply = " << MFLOPs << endl<< endl;
    }
    
    // cout << "Vector bcomp" << bcomp << endl;
    
    if (TransA) {
      if ((ierr=resid.Update(-1.0, btexact, 1.0, bcomp, 0.0))) perror1("Error in Update",ierr);}
    else {
      if ((ierr=resid.Update(-1.0, bexact, 1.0, bcomp, 0.0))) perror1("Error in Update",ierr);}
    
    if ((ierr = resid.Norm2(residual))) perror1("Error in Norm2",ierr);
    
    if (Comm.MyPID()==0)
      for (i=0; i< nrhs; i++) printf("Residual[%d]    = %22.16g\n",i,residual[i]);

  }

   // Optimize data layout for memory access

  if ((ierr=A.OptimizeStorage())) perror1("Error in Petra_RDP_CRS_Matrix: OptimizeStorage",ierr);

  for (ii=0; ii<2; ii++) {
    bool TransA = (ii==1);
    A.ResetFlops();
    elapsed_time = timer.ElapsedTime();
    for (iii=0; iii<NumMVs; iii++)
      if ((ierr=A.Multiply(TransA, xcomp, bcomp)))
	perror1("Error in Multiply",ierr);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    total_flops = A.Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (Comm.MyPID()==0) {
      cout << "\n\n****************************************************" << endl;
      if (TransA) cout << "\n\nResults for transpose multiply with optimized storage" << endl;
      else cout << "\n\nResults for no transpose multiply with optimized storage"<< endl;
      cout << "\n\nTotal FLOPS for optimized storage (" <<NumMVs<< " Multiplies) = " 
	   << total_flops<< endl;
      cout << "    Total time for optimized Storage = " << elapsed_time<< endl;
      cout << "    Total MFLOPs for optimized matrix multiply = " << MFLOPs << endl<< endl;
    }
    
    if (TransA) {
      if ((ierr=resid.Update(-1.0, btexact, 1.0, bcomp, 0.0))) perror1("Error in Update",ierr);}
    else {
      if ((ierr=resid.Update(-1.0, bexact, 1.0, bcomp, 0.0))) perror1("Error in Update",ierr); }

    if ((ierr = resid.Norm2(residual))) perror1("Error in Norm2",ierr);

    if (Comm.MyPID()==0)
      for (i=0; i< nrhs; i++) printf("Residual[%d]    = %22.16g\n",i,residual[i]);

  }

  free ((void *) xguess);
  free ((void *) b);
  free ((void *) bt);
  free ((void *) xexact);
  free ((void *) val);
  free ((void *) bindx);
  free ((void *) update);

#ifdef VBRMATRIX
  free ((void *) indx);
  free ((void *) rpntr);
  free ((void *) bpntr);
  delete [] blockSizes;
  delete [] numNzBlks;
#else
  delete [] numNz;

#endif

#ifdef MULTI_VECTOR
  delete [] xexact2;
  delete [] b2;
  delete [] allx;
  delete [] allb;
#endif
  delete &bexact;
  delete &bcomp;
  delete &resid;
  delete &xcomp;
  delete ILUK;
  delete &ILUK_Graph;
  delete &A;
  delete &map;
  delete &timer;
  delete &FillTimer;
  delete &Comm;

  delete proc_config;
				       
#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
void BiCGSTAB(Petra_RDP_CRS_Matrix &A, 
	      Petra_RDP_Vector &x, 
	      Petra_RDP_Vector &b, 
	      Ifpack_RDP_CRS_RILUK *M, 
	      int Maxiter, 
	      double Tolerance, 
	      double *residual, double & FLOPS, bool verbose) {

  // Allocate vectors needed for iterations
  Petra_RDP_Vector& phat = *new Petra_RDP_Vector(x); phat.PutScalar(0.0);
  Petra_RDP_Vector& p = *new Petra_RDP_Vector(x); p.PutScalar(0.0);
  Petra_RDP_Vector& shat = *new Petra_RDP_Vector(x); shat.PutScalar(0.0);
  Petra_RDP_Vector& s = *new Petra_RDP_Vector(x); s.PutScalar(0.0);
  Petra_RDP_Vector& r = *new Petra_RDP_Vector(x); r.PutScalar(0.0);
  Petra_RDP_Vector& rtilde = *new Petra_RDP_Vector(x); rtilde.Random();
  Petra_RDP_Vector& v = *new Petra_RDP_Vector(x); v.PutScalar(0.0);
  

  A.Multiply(false, x, r); // r = A*x

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
      M->LevelSolve(false, p, phat);
    A.Multiply(false, phat, v);

    
    rtilde.Dot(v,&sigma);
    alpha = rhon/sigma;    

    /* s = r - alpha*v                     */
    /* shat = M^-1 s                       */
    /* r = A shat (r is a tmp here for t ) */

    s.Update(-alpha, v, 1.0, r, 0.0);
    if (M==0) 
      shat.Scale(1.0, s);
    else
      M->LevelSolve(false, s, shat);
    A.Multiply(false, shat, r);

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

  FLOPS = phat.Flops()+p.Flops()+shat.Flops()+s.Flops()+r.Flops()+rtilde.Flops()+
          v.Flops()+A.Flops()+x.Flops()+b.Flops();
  if (M!=0) FLOPS += M->Flops();

  delete &phat;
  delete &p;
  delete &shat;
  delete &s;
  delete &r;
  delete &rtilde;
  delete &v;


  return;
}
