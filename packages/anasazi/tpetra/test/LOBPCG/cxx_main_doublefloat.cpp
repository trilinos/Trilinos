// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This test is for the LOBPCG solver manager on a standard (Ax=xl) symmetric
// eigenvalue problem, with a diagonal preconditioner. The test can be run in either
// float or double.
//
// This driver reads the problem from a Harwell-Boeing (HB) file.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

using namespace Teuchos;
using namespace Anasazi;
using Tpetra::DefaultPlatform;
using Tpetra::Platform;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;
using std::endl;
using std::cout;
using std::string;
using std::setw;
using std::vector;

// TODO: adjust this test to use CrsMatrix<double> regardless of the scalar type for MultiVector
//       this requires the ability to apply those operators to those MultiVectors, which Tpetra
//       currently cannot do.

bool proc_verbose = false, reduce_tol, precond = true, dumpdata = false;
RCP<Map<int> > vmap;
ParameterList pl; 
int dim, nev, blockSize;
double *dvals; 
int *colptr, *rowind;
int MyPID = 0;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar> 
RCP<EigenProblem<Scalar,MultiVector<int,Scalar>,Operator<int,Scalar> > > buildProblem()
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<int,Scalar>         OP;
  typedef MultiVector<int,Scalar>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;
  RCP<CrsMatrix<int,Scalar> > A = rcp(new CrsMatrix<int,Scalar>(*vmap));
  if (MyPID == 0) {
    // HB format is compressed column. CrsMatrix is compressed row.
    const float *dptr = dvals;
    const int *rptr = rowind;
    for (int c=0; c<dim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        A->submitEntry(*rptr-1,c,*dptr);
        A->submitEntry(c,*rptr-1,*dptr);
        ++rptr;
        ++dptr;
      }
    }
  }
  // distribute matrix data to other nodes
  A->fillComplete();
  // Create initial MV
  RCP<MV> X0;
  X0 = rcp( new MV(*vmap,blockSize) );
  MVT::MvRandom( *X0 );
  // Construct a linear problem instance with zero initial MV
  RCP<EigenProblem<Scalar,MV,OP> > problem = rcp( new BasicEigenProblem<Scalar,MV,OP>(A,X0) );
  problem->setNEV(nev);
  // diagonal preconditioner
  if (precond) {
    Array<Scalar> diags(A->getNumMyDiagonals());
    A->getMyDiagCopy(diags());
    typename Array<Scalar>::iterator i;
    for (i=diags.begin(); i != diags.end(); ++i) {
      TEST_FOR_EXCEPTION(*i <= SCT::zero(), std::runtime_error,"Matrix has non-positive diagonal: " << *i);
      *i = SCT::one() / *i;
    }
    RCP<CrsMatrix<int,Scalar> > P = rcp(new CrsMatrix<int,Scalar>(*vmap));
    int gid=vmap->getMinGlobalIndex();
    for (i=diags.begin(); i!=diags.end(); ++i) {
      P->submitEntry(gid,gid,*i);
      ++gid;
    }
    P->fillComplete();
    problem->setPrec(P);
  }
  TEST_FOR_EXCEPT(problem->setProblem() == false);
  return problem;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar>
bool runTest(double ltol, double times[], int &numIters) 
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<int,Scalar>         OP;
  typedef MultiVector<int,Scalar>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;

  const Scalar ONE  = SCT::one();
  pl.set<MT>( "Convergence Tolerance", ltol );         // Relative convergence tolerance requested

  if (MyPID==0) cout << "Testing Scalar == " << typeName(ONE) << endl;

  RCP<EigenProblem<Scalar,MV,OP> > problem;
  Time btimer("Build Timer"), ctimer("Construct Timer"), stimer("Solve Timer");
  ReturnType ret;
  if (MyPID==0) cout << "Building problem..." << endl;
  { 
    TimeMonitor localtimer(btimer);
    problem = buildProblem<Scalar>();
  }
  RCP<SolverManager<Scalar,MV,OP> > solver;
  if (MyPID==0) cout << "Constructing solver..." << endl; 
  {
    TimeMonitor localtimer(ctimer);
    solver = rcp(new LOBPCGSolMgr<Scalar,MV,OP>( problem, rcp(&pl,false) ));
  }
  if (MyPID==0) cout << "Solving problem..." << endl;
  {
    TimeMonitor localtimer(stimer);
    try {
      ret = solver->solve();
    }
    catch (std::exception &e) {
      cout << "Caught exception: " << endl << e.what() << endl;
      ret = Unconverged;
    }
  }
  numIters = solver->getNumIters();
  times[0] = btimer.totalElapsedTime();
  times[1] = ctimer.totalElapsedTime();
  times[2] = stimer.totalElapsedTime();

  bool badRes = false;
  if (ret == Converged) {
    if (proc_verbose) cout << endl;
    if (MyPID==0) cout << "Computing residuals..." << endl;
    //
    // Compute actual residuals.
    //
    Eigensolution<Scalar,MV> solution = problem->getSolution();
    int numsol = solution.numVecs;
    RCP<const OP> A = problem->getA();
    RCP<MV> X = solution.Evecs;
    RCP<MV> R = MVT::MvClone(*X,numsol);
    OPT::Apply(*A,*X,*R); // R = A*X
    vector<Value<Scalar> > lambdas = solution.Evals;
    SerialDenseMatrix<int,Scalar> L(numsol,numsol);
    for (int i=0; i<numsol; ++i) {
      L(i,i) = lambdas[i].realpart;
    }
    MVT::MvTimesMatAddMv(-1.0,*X,L,1.0,*R); // R = A*X - X*L
    vector<MT> resnorms( numsol );
    MVT::MvNorm(*R,resnorms);
    for (int i=0; i<numsol; ++i) {
      if (lambdas[i] != 0.0) {
        resnorms[i] /= SCT::magnitude(lambdas[i]);
      }
    }
    if (proc_verbose) cout << setw(16) << "Lambda" << setw(16) << "Residuals" << endl;
    for ( int i=0; i<numsol; ++i) {
      if (proc_verbose) cout << setw(16) << lambdas[i] << setw(16) << resnorms[i] << endl;
      if (resnorms[i] > ltol) badRes = true;
    }
    if (proc_verbose) cout << endl;
  }
  if (ret != Converged || badRes) return false;
  return true;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) 
{
  GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Platform<int> > platform = DefaultPlatform<int>::getPlatform();
  RCP<Comm<int> > comm = platform->createComm();

  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, debug = false;
  nev = 1;             // number of eigenpairs to find
  int blocksize = 1;   // blocksize used by solver
  int maxiters = -1;   // maximum number of iterations for solver to use
  string filename("bcsstk14.hb");
  double tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("nev",&nev,"Number of eigenpairs to find.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
  cmdp.setOption("reduce-tol","fixed-tol",&reduce_tol,"Require increased accuracy from higher precision scalar types.");
  cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
  cmdp.setOption("dump-data","no-dump-data",&dumpdata,"Dump raw data to data.dat.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }

  MyPID = rank(*comm);
  proc_verbose = ( verbose && (MyPID==0) );
  if (proc_verbose) cout << Anasazi_Version() << endl << endl;

  //
  // Get the data (double) from the HB file and build the Map,Matrix
  //
  int nnz, info;
  nnz = -1;
  if (MyPID == 0) {
    int dim2;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
  }
  else {
    // address uninitialized data warnings
    dvals = NULL;
    colptr = NULL;
    rowind = NULL;
  }
  broadcast(*comm,0,&info);
  broadcast(*comm,0,&nnz);
  broadcast(*comm,0,&dim);
  if (info == 0 || nnz < 0) {
    if (MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // create map
  vmap = rcp(new Map<int>(dim,0,*platform));
  // create the parameter list
  if (maxiters == -1) {
    maxiters = dim/blocksize - 1; // maximum number of iterations to run
  }
  //
  pl.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  pl.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  int verbLevel = Errors + Warnings;
  if (debug) {
    verbLevel += Debug;
  }
  if (verbose) {
    verbLevel += TimingDetails + FinalSummary;
  }
  pl.set( "Verbosity", verbLevel );

  //
  // **********Print out information about problem*******************
  //
  if (MyPID==0) {
    cout << "Filename: " << filename << endl;
    cout << "Dimension of matrix: " << dim << endl;
    cout << "Number of nonzeros: " << nnz << endl;
    cout << "NEV: " << nev << endl;
    cout << "Block size used by solver: " << blocksize << endl;
    cout << "Max number of iterations: " << maxiters << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }

  // run tests for different scalar types
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter);
  /*
#ifdef HAVE_TEUCHOS_QD
  unsigned int old_cw; fpu_fix_start(&old_cw); // fix floating point rounding for intel processors; required for proper implementation
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double ddtime[3]; int dditer; bool ddpass = runTest<dd_real>(ltol,ddtime,dditer);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double qdtime[3]; int qditer; bool qdpass = runTest<qd_real>(ltol,qdtime,qditer);
  fpu_fix_end(&old_cw);                        // restore previous floating point rounding
#endif
  */

  // done with the matrix data now; delete it
  if (MyPID == 0) {
    free( dvals );
    free( colptr );
    free( rowind );
  }

  if (MyPID==0) {
    cout << "Scalar field     Build time     Init time     Solve time     Num Iters     Test Passsed" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << setw(12) << "float"   << "     " << setw(10) <<  ftime[0] << "     " << setw(9) <<  ftime[1] << "     " << setw(10) <<  ftime[2] << "     " << setw(9) <<  fiter << "     " << ( fpass ? "pass" : "fail") << endl;
    cout << setw(12) << "double"  << "     " << setw(10) <<  dtime[0] << "     " << setw(9) <<  dtime[1] << "     " << setw(10) <<  dtime[2] << "     " << setw(9) <<  diter << "     " << ( dpass ? "pass" : "fail") << endl;
    /*
#ifdef HAVE_TEUCHOS_QD
    cout << setw(12) << "dd_real" << "     " << setw(10) << ddtime[0] << "     " << setw(9) << ddtime[1] << "     " << setw(10) << ddtime[2] << "     " << setw(9) << dditer << "     " << (ddpass ? "pass" : "fail") << endl;
    cout << setw(12) << "qd_real" << "     " << setw(10) << qdtime[0] << "     " << setw(9) << qdtime[1] << "     " << setw(10) << qdtime[2] << "     " << setw(9) << qditer << "     " << (qdpass ? "pass" : "fail") << endl;
#endif
    */
  }

  if (dumpdata) {
    std::ofstream fout("data.dat",std::ios_base::app | std::ios_base::out);
    fout << setw(14) << nnz << setw(8) << "float"   << setw(12) <<  ftime[0] << setw(12) <<  ftime[1] << setw(12) <<  ftime[2] << setw(12) <<  fiter << endl;
    fout << setw(14) << nnz << setw(8) << "double"  << setw(12) <<  dtime[0] << setw(12) <<  dtime[1] << setw(12) <<  dtime[2] << setw(12) <<  diter << endl;
    /*
#ifdef HAVE_TEUCHOS_QD
    fout << setw(14) << nnz << setw(8) << "dd_real" << setw(12) << ddtime[0] << setw(12) << ddtime[1] << setw(12) << ddtime[2] << setw(12) << dditer << endl;
    fout << setw(14) << nnz << setw(8) << "qd_real" << setw(12) << qdtime[0] << setw(12) << qdtime[1] << setw(12) << qdtime[2] << setw(12) << qditer << endl;
#endif
    */
    fout.close();
  }

  bool allpass = fpass && dpass;
  /*
#ifdef HAVE_TEUCHOS_QD
  allpass = allpass && ddpass && qdpass;
#endif
  */

  if (!allpass) {
    if (MyPID==0) cout << "\nEnd Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (MyPID==0) cout << "\nEnd Result: TEST PASSED" << endl;
  return 0;
}


