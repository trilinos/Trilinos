#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsSingletonFilter.h"

#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

void BiCGSTAB(Epetra_CrsMatrix &A, Epetra_Vector &x, Epetra_Vector &b, 
	      Ifpack_CrsRiluk *M, 
	      int Maxiter, double Tolerance, 
	      double *residual, bool verbose);

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
  bool verbose1 = true;
  if (MyPID==0) verbose = true;

  if(argc < 2 && verbose) {
    cerr << "Usage: " << argv[0] 
	 << " HPC_filename [level_fill [level_overlap [absolute_threshold [ relative_threshold]]]]" << endl
	 << "where:" << endl
	 << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
	 << "level_fill         - The amount of fill to use for ILU(k) preconditioner (default 0)" << endl
	 << "level_overlap      - The amount of overlap used for overlapping Schwarz subdomains (default 0)" << endl
	 << "absolute_threshold - The minimum value to place on the diagonal prior to factorization (default 0.0)" << endl
	 << "relative_threshold - The relative amount to perturb the diagonal prior to factorization (default 1.0)" << endl << endl
	 << "To specify a non-default value for one of these parameters, you must specify all" << endl
	 << " preceding values but not any subsequent parameters. Example:" << endl
	 << "ifpackHpcSerialMsr.exe mymatrix.hpc 1  - loads mymatrix.hpc, uses level fill of one, all other values are defaults" << endl
	 << endl;
    return(1);

  }

  // Uncomment the next three lines to debug in mpi mode
  int tmp;
  if (MyPID==0) cin >> tmp;
  Comm.Barrier();

  Epetra_Map * map;
  Epetra_CrsMatrix * A; 
  Epetra_Vector * x; 
  Epetra_Vector * b;
  Epetra_Vector * xexact;

  Trilinos_Util_ReadHpc2Epetra(argv[1], Comm, map, A, x, b, xexact);

  bool smallProblem = false;
  if (A->RowMap().NumGlobalElements()<100) smallProblem = true;

  if (smallProblem)
    cout << "Original Matrix = " << endl << *A   << endl;

  x->PutScalar(0.0);

  Epetra_LinearProblem FullProblem(A, x, b); 
  
  Epetra_Time ReductionTimer(Comm);
  Epetra_CrsSingletonFilter SingletonFilter(&FullProblem);
  Comm.Barrier();
  double reduceInitTime = ReductionTimer.ElapsedTime();
  SingletonFilter.Analyze();
  Comm.Barrier();
  double reduceAnalyzeTime = ReductionTimer.ElapsedTime() - reduceInitTime;

  if (SingletonFilter.SingletonsDetected())
    cout << "Singletons found" << endl;
  else {
    cout << "Singletons not found" << endl;
    exit(1);
  }
  SingletonFilter.ConstructReducedProblem();
  Comm.Barrier();
  double reduceConstructTime = ReductionTimer.ElapsedTime() - reduceInitTime;

  double totalReduceTime = ReductionTimer.ElapsedTime();

  if (verbose)
    cout << "\n\n****************************************************" << endl
	 << "    Reduction init  time (sec)           = " << reduceInitTime<< endl
	 << "    Reduction Analyze time (sec)         = " << reduceAnalyzeTime << endl
	 << "    Construct Reduced Problem time (sec) = " << reduceConstructTime << endl
	 << "    Reduction Total time (sec)           = " << totalReduceTime << endl<< endl;

  SingletonFilter.Statistics();

  Epetra_LinearProblem * ReducedProblem = SingletonFilter.ReducedProblem();

  Epetra_CrsMatrix * Ap = dynamic_cast<Epetra_CrsMatrix *>(ReducedProblem->GetMatrix());
  Epetra_Vector * bp = (*ReducedProblem->GetRHS())(0);
  Epetra_Vector * xp = (*ReducedProblem->GetLHS())(0);
  

  if (smallProblem)
    cout << " Reduced Matrix = " << endl << *Ap << endl
	 << " LHS before sol = " << endl << *xp << endl
	 << " RHS            = " << endl << *bp << endl;

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
    IlukGraph = new Ifpack_IlukGraph(Ap->Graph(), LevelFill, Overlap);
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
    int initerr = ILUK->InitValues(*Ap);
    if (initerr!=0) {
      cout << endl << Comm << endl << "  InitValues error = " << initerr;
      if (initerr==1) cout << "  Zero diagonal found, warning error only";
      cout << endl << endl;
    }
    assert(ILUK->Factor()==0);
    elapsed_time = timer.ElapsedTime() - elapsed_time;
    total_flops = ILUK->Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Time to compute preconditioner values = " 
		    << elapsed_time << endl
		    << "MFLOPS for Factorization = " << MFLOPs << endl;
    //cout << *ILUK << endl;
  double Condest;
  ILUK->Condest(false, Condest);

  if (verbose) cout << "Condition number estimate for this preconditioner = " << Condest << endl;
  }
  int Maxiter = 50;
  double Tolerance = 1.0E-11;

  Epetra_Flops counter;
  Ap->SetFlopCounter(counter);
  xp->SetFlopCounter(*Ap);
  bp->SetFlopCounter(*Ap);
  if (ILUK!=0) ILUK->SetFlopCounter(*Ap);

  elapsed_time = timer.ElapsedTime();

  double residual;
  BiCGSTAB(*Ap, *xp, *bp, ILUK, Maxiter, Tolerance, &residual, verbose);

  elapsed_time = timer.ElapsedTime() - elapsed_time;
  total_flops = counter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "Time to compute solution = " 
		    << elapsed_time << endl
		    << "Number of operations in solve = " << total_flops << endl
		    << "MFLOPS for Solve = " << MFLOPs<< endl << endl;

  SingletonFilter.ComputeFullSolution();

  if (smallProblem)
  cout << " Reduced LHS after sol = " << endl << *xp << endl
       << " Full    LHS after sol = " << endl << *x << endl
       << " Full  Exact LHS         = " << endl << *xexact << endl;

  Epetra_Vector resid(*x);

  resid.Update(1.0, *x, -1.0, *xexact, 0.0); // resid = xcomp - xexact

  resid.Norm2(&residual);

  if (verbose) 
    cout << "2-norm of difference between computed and exact solution  = " << residual << endl;
    
  if (verbose1 && residual>1.0e-5) {
    if (verbose)
      cout << "Difference between computed and exact solution appears large..." << endl
	   << "Computing norm of A times this difference.  If this norm is small, then matrix is singular"
	   << endl;
    assert(A->Multiply(false, resid, *b)==0);
    assert(b->Norm2(&residual)==0);
    if (verbose) 
      cout << "2-norm of A times difference between computed and exact solution  = " << residual << endl;
    
  }
  


  if (ILUK!=0) delete ILUK;
  if (IlukGraph!=0) delete IlukGraph;
				       
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
void BiCGSTAB(Epetra_CrsMatrix &A, 
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
  Epetra_Vector r(b.Map()); r.SetFlopCounter(x);
  Epetra_Vector rtilde(x.Map()); rtilde.Random(); rtilde.SetFlopCounter(x);
  Epetra_Vector v(x.Map()); v.SetFlopCounter(x);
  

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
      M->Solve(false, p, phat);
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
      M->Solve(false, s, shat);
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

    if (r_norm < Tolerance) break;
  }
  return;
}
