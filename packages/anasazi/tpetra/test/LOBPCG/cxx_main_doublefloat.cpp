// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test is for the LOBPCG solver manager on a standard (Ax=xl) symmetric
// eigenvalue problem, with a diagonal preconditioner. The test can be run in either
// float or double.
//
// This driver reads the problem from a Harwell-Boeing (HB) file.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DiagPrecond.hpp>


using namespace Teuchos;
using namespace Anasazi;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Vector;
using Tpetra::DiagPrecond;
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
RCP<Map<> > vmap;
ParameterList mptestpl;
int rnnzmax;
int dim, mptestnev, blockSize;
double *dvals;
int *colptr, *rowind;
int mptestmypid = 0;
string mptestwhich("LR");

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar>
RCP<Eigenproblem<Scalar, MultiVector<Scalar>, Operator<Scalar> > > buildProblem()
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<Scalar>             OP;
  typedef MultiVector<Scalar>          MV;
  typedef MV::global_ordinal_type      GO;

  RCP<CrsMatrix<Scalar> > A = rcp (new CrsMatrix<Scalar> (*vmap, rnnzmax));
  if (mptestmypid == 0) {
    // HB format is compressed column. CrsMatrix is compressed row.
    const double *dptr = dvals;
    const int *rptr = rowind;
    for (int c=0; c<dim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        // FIXME (mfh 22 Oct 2015) This assumes that Scalar is real.
        A->insertGlobalValues (static_cast<GO> (*rptr-1), tuple<GO> (c), tuple<Scalar> (*dptr));
        if (c != *rptr -1) {
          A->insertGlobalValues (static_cast<GO> (c), tuple<GO> (*rptr-1), tuple<Scalar> (*dptr));
        }
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
  X0->putScalar( SCT::one() );
  // Construct a linear problem instance with zero initial MV
  RCP<Eigenproblem<Scalar,MV,OP> > problem = rcp( new BasicEigenproblem<Scalar,MV,OP>(A,X0) );
  problem->setHermitian(true);
  problem->setNEV(mptestnev);
  // diagonal preconditioner
  if (precond) {
    Vector<Scalar> diags(A->getRowMap());
    A->getLocalDiagCopy(diags);
    for (Teuchos_Ordinal i=0; i<vmap->getNumMyEntries(); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(diags[i] <= SCT::zero(), std::runtime_error,"Matrix is not positive-definite: " << diags[i]);
      diags[i] = SCT::one() / diags[i];
    }
    RCP<Operator<Scalar> > P = rcp (new DiagPrecond<Scalar> (diags));
    problem->setPrec(P);
  }
  TEUCHOS_TEST_FOR_EXCEPT(problem->setProblem() == false);
  return problem;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar>
bool runTest(double ltol, double times[], int &numIters)
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<Scalar>             OP;
  typedef MultiVector<Scalar>          MV;

  const Scalar ONE  = SCT::one();
  mptestpl.set<MT>( "Convergence Tolerance", ltol );         // Relative convergence tolerance requested
  mptestpl.set<MT>( "Locking Tolerance", 0.1*ltol );         // Relative convergence tolerance requested

  if (mptestmypid==0) cout << "\nTesting Scalar == " << typeName(ONE) << endl;

  RCP<Eigenproblem<Scalar,MV,OP> > problem;
  Time btimer("Build Timer"), ctimer("Construct Timer"), stimer("Solve Timer");
  ReturnType ret;
  if (mptestmypid==0) cout << "Building problem..." << endl;
  {
    TimeMonitor localtimer(btimer);
    problem = buildProblem<Scalar>();
  }
  RCP<SolverManager<Scalar,MV,OP> > solver;
  if (mptestmypid==0) cout << "Constructing solver..." << endl;
  {
    TimeMonitor localtimer(ctimer);
    solver = rcp(new LOBPCGSolMgr<Scalar,MV,OP>( problem, mptestpl));
  }
  if (mptestmypid==0) cout << "Solving problem..." << endl;
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
    if (mptestmypid==0) cout << "Computing residuals..." << endl;
    //
    // Compute actual residuals.
    //
    Eigensolution<Scalar,MV> solution = problem->getSolution();
    int numsol = solution.numVecs;
    RCP<const Operator<Scalar> > A = problem->getOperator();
    RCP<MultiVector<Scalar> > X = solution.Evecs;
    RCP<MultiVector<Scalar> > R = rcp(new MultiVector<Scalar>(*X));
    RCP<MultiVector<Scalar> > XL= rcp(new MultiVector<Scalar>(*X));
    A->apply(*X,*R);  // R = A*X
    Array<MT> L( numsol );
    for (int i=0; i<numsol; ++i) {
      L[i] = solution.Evals[i].realpart;
    }
    XL->scale(L);
    R->update(-1.0,*XL,1.0);
    Array<MT> resnorms( numsol );
    R->norm2(resnorms);
    for (int i=0; i<numsol; ++i) {
      if (L[i] != 0.0) {
        resnorms[i] /= SCT::magnitude(L[i]);
      }
    }
    if (proc_verbose) cout << setw(16) << "Lambda" << setw(16) << "Residuals" << endl;
    for ( int i=0; i<numsol; ++i) {
      if (proc_verbose) cout << setw(16) << L[i] << setw(16) << resnorms[i] << endl;
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
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  //
  // Get test parameters from command-line processor
  //
  bool verbose = false, debug = false;
  mptestnev = 1;             // number of eigenpairs to find
  blockSize = 1;   // blocksize used by solver
  int maxiters = -1;   // maximum number of iterations for solver to use
  string filename("bcsstk14.hb");
  double tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("which",&mptestwhich,"Which eigenvalues to compute.");
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("nev",&mptestnev,"Number of eigenpairs to find.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("block-size",&blockSize,"Block size to be used by the CG solver.");
  cmdp.setOption("reduce-tol","fixed-tol",&reduce_tol,"Require increased accuracy from higher precision scalar types.");
  cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
  cmdp.setOption("dump-data","no-dump-data",&dumpdata,"Dump raw data to data.dat.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }

  mptestmypid = rank(*comm);
  proc_verbose = ( verbose && (mptestmypid==0) );
  if (proc_verbose) cout << Anasazi_Version() << endl << endl;

  //
  // Get the data (double) from the HB file and build the Map,Matrix
  //
  int nnz, info;
  nnz = -1;
  if (mptestmypid == 0) {
    int dim2;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
    // find maximum NNZ over all rows
    vector<int> rnnz(dim,0);
    for (int *ri=rowind; ri<rowind+nnz; ++ri) {
      ++rnnz[*ri-1];
    }
    for (int c=0; c<dim; ++c) {
      rnnz[c] += colptr[c+1]-colptr[c];
    }
    rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
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
  broadcast(*comm,0,&rnnzmax);
  if (info == 0 || nnz < 0) {
    if (mptestmypid == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // create map
  vmap = rcp (new Map<> (dim,0,comm));
  // create the parameter list
  if (maxiters == -1) {
    maxiters = dim/blockSize - 1; // maximum number of iterations to run
  }
  //
  mptestpl.set( "Which", mptestwhich );
  mptestpl.set( "Block Size", blockSize );              // Blocksize to be used by iterative solver
  mptestpl.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  int verbLevel = Errors + Warnings;
  if (debug) {
    verbLevel += Debug;
  }
  if (verbose) {
    verbLevel += TimingDetails + FinalSummary;
  }
  mptestpl.set( "Verbosity", verbLevel );

  //
  // **********Print out information about problem*******************
  //
  if (mptestmypid==0) {
    cout << "Filename: " << filename << endl;
    cout << "Dimension of matrix: " << dim << endl;
    cout << "Number of nonzeros: " << nnz << endl;
    cout << "NEV: " << mptestnev << endl;
    cout << "Block size used by solver: " << blockSize << endl;
    cout << "Max number of iterations: " << maxiters << endl;
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }

  // run tests for different scalar types
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter);

  // done with the matrix data now; delete it
  if (mptestmypid == 0) {
    free( dvals );
    free( colptr );
    free( rowind );
  }

  if (mptestmypid==0) {
    cout << endl << endl;
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
    if (mptestmypid==0) cout << "\nEnd Result: TEST FAILED" << endl;
    return -1;
  }
  //
  // Default return value
  //
  if (mptestmypid==0) cout << "\nEnd Result: TEST PASSED" << endl;
  return 0;
}


