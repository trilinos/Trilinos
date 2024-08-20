// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#include <fstream>

#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "MyMemory.h"

#include "ModalAnalysisSolver.h"

#include "BRQMIN.h"
#include "BlockDACG.h"
#include "LOBPCG.h"
#include "LOBPCG_light.h"
#include "KnyazevLOBPCG.h"
#include "ARPACKm3.h"
#include "ModifiedARPACKm3.h"
#include "Davidson.h"
#include "JDPCG.h"

#include "ModalProblem.h"
#include "ModeLaplace1DQ1.h"
#include "ModeLaplace1DQ2.h"
#include "ModeLaplace2DQ1.h"
#include "ModeLaplace2DQ2.h"
#include "ModeLaplace3DQ1.h"
#include "ModeLaplace3DQ2.h"

#include "BlockPCGSolver.h"
#include "AMGOperator.h"
#include "MyIncompleteChol.h"

#include <stdexcept>
#include <Teuchos_Assert.hpp>

const int LOBPCG_CHOL = 1;
const int LOBPCG_LIGHT = 2;
const int LOBPCG_AK_CHOL = 3;
const int ARPACK_ORIG = 5;
const int ARPACK_RESI = 6;
const int GALDAVIDSON = 7;
const int BRQMIN_CHOL = 9;
const int DACG_CHOL = 11;
const int JD_PCG = 13;

const int NO_PREC = 0;
const int AMG_POLYNOMIAL = 1;
const int AMG_GAUSSSEIDEL = 2;
const int INC_CHOL = 3;


int main(int argc, char *argv[]) {

  int i;
  int myPid = 0;
  int numProc = 1;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPid);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Initialize the random number generator
  srand((unsigned) time(0));

  // Initialize memory counters
  initMemCounters();

  int dimension;
  double Lx = 1.0, Ly = 1.0, Lz = 1.0;
  int nX=0, nY=0, nZ=0;
  int casePb = 0;
  int algo = 0;
  int precond = 0, param = 0;
  int numEigen = 0;
  int dimSearch = 0;
  int numBlock = 1;
  int maxIter = 0;
  double tol = 0.0;
  int maxIterCG = 0;
  double tolCG = 0.0;
  int verbose = 0;

  bool paramStop = false;
  for (i=0; i<numProc; ++i) {
    Comm.Barrier();
    if (myPid == i) {
      ifstream fin("control.driver");
      char buff[101];
      TEUCHOS_TEST_FOR_EXCEPTION( !fin, std::runtime_error, "The input file 'control."
			  "driver' could not be opened." );
      fin >> dimension; fin.getline(buff, 100);
      switch (dimension) {
        case 1:
          fin >> Lx; fin.getline(buff, 100);
          fin >> nX; fin.getline(buff, 100);
          fin >> casePb; fin.getline(buff, 100);
          break;
        case 2:
          fin >> Lx >> Ly; fin.getline(buff, 100);
          fin >> nX >> nY; fin.getline(buff, 100);
          fin >> casePb; fin.getline(buff, 100);
          break;
        case 3:
          fin >> Lx >> Ly >> Lz; fin.getline(buff, 100);
          fin >> nX >> nY >> nZ; fin.getline(buff, 100);
          fin >> casePb; fin.getline(buff, 100);
          break;
        default:
          paramStop = true;
      }
      fin >> algo; fin.getline(buff, 100);
      fin >> precond; fin.getline(buff, 100);
      fin >> param; fin.getline(buff, 100);
      fin >> numEigen; fin.getline(buff, 100);
      fin >> dimSearch; fin.getline(buff, 100);
      fin >> numBlock; fin.getline(buff, 100);
      fin >> maxIter; fin.getline(buff, 100);
      fin >> tol; fin.getline(buff, 100);
      fin >> maxIterCG; fin.getline(buff, 100);
      fin >> tolCG; fin.getline(buff, 100);
      fin >> verbose; fin.getline(buff, 100);
    }
    Comm.Barrier();
  }

  // Check the input parameters
  switch (algo) {
    case LOBPCG_CHOL:
    case LOBPCG_LIGHT:
    case LOBPCG_AK_CHOL:
    case 4:
    case ARPACK_ORIG:
    case ARPACK_RESI:
    case GALDAVIDSON:
    case BRQMIN_CHOL:
    case DACG_CHOL:
    case JD_PCG:
      break;
    default:
      paramStop = true;
      break;
  }

  switch (precond) {
    case NO_PREC:
    case AMG_POLYNOMIAL:
    case AMG_GAUSSSEIDEL:
    case INC_CHOL:
      break;
    default:
      paramStop = true;
      break;
  }

  if ((dimension < 1) && (dimension > 3))
    paramStop = true;
  else {
    if ((casePb != 1) && (casePb != 2))
      paramStop = true;
  }

  if (paramStop == true) {
    if (verbose*(myPid==0) > 0) { 
#ifdef EPETRA_MPI
      cerr << "Usage: prun -n NP ./driver" << endl;
#else
      cerr << "Usage: ./driver.serial" << endl;
#endif
      cerr << endl;
      cerr << "Input file: 'control.driver'" << endl;
      cerr << endl;
      cerr << "Example of input file" << endl;
      cerr << endl;
      cerr << "3           Space dimension (1, 2, 3)" << endl;
      cerr << "1.0 2.0 3.0 Lx Ly Lz (Dimension of brick in each direction)" << endl;
      cerr << "10 11 12    nX nY nZ (Number of elements in each direction)" << endl;
      cerr << "1           Case Problem (1=Q1, 2=Q2)" << endl;
      cerr << "1           Algorithm (1 .. 15)" << endl;
      cerr << "1           Preconditioner (0 = None, 1 = AMG Poly, 2 = AMG GS)"
              << endl;
      cerr << "1           Parameter for preconditioner: AMG Poly. Deg."
              << endl;
      cerr << "10          Number of eigenpairs requested" << endl;
      cerr << "5           Dimension of block size" << endl;
      cerr << "4           Number of blocks (ONLY referenced by Davidson (7) and JD-PCG (13))"
              << endl;
      cerr << "1000        Maximum number of iterations for eigensolver" << endl;
      cerr << "1e-06       Tolerance" << endl;
      cerr << "1000        Maximum number of inner iterations (PCG, QMR)" << endl;
      cerr << "1e-07       Tolerance for PCG solve (ONLY referenced in ARPACK)" << endl;
      cerr << "1           Verbose level" << endl;
      cerr << "----------------------------------------------------------------------" << endl;
      cerr << endl;
      cerr << "Eigensolver flags:" << endl;
      cerr << " 1. LOBPCG (RBL & UH version) with Cholesky-based local eigensolver" << endl;
      cerr << " 2. LOBPCG (light version) with Cholesky-based local eigensolver" << endl;
      cerr << " 3. LOBPCG (AK version) with Cholesky-based local eigensolver" << endl;
      cerr << " 5. ARPACK (original version)" << endl;
      cerr << " 6. ARPACK (version with user-provided residual)" << endl;
      cerr << " 7. Generalized Davidson (RBL & UH version)" << endl;
      cerr << " 9. BRQMIN with Cholesky-based local eigensolver" << endl;
      cerr << "11. Block DACG with Cholesky-based local eigensolver" << endl;
      cerr << "13. JDCG (Notay's version of Jacobi-Davidson)" << endl;
      cerr << endl;
      cerr << "The dimension of search space corresponds to" << endl;
      cerr << " * the block size for LOBPCG (1)" << endl;
      cerr << " * NCV for ARPACK (5, 6)" << endl;
      cerr << " * one block size for Generalized Davidson (7)" << endl;
      cerr << " * the block size for BRQMIN (9)" << endl;
      cerr << " * the block size for Block DACG (11)" << endl;
      cerr << " * one block size for Jacobi-Davidson JDCG (13)" << endl;
      cerr << endl;
    }
    // mfh 14 Jan 2011: Replaced "exit(1)" with "return 1" since this
    // is the "main" routine.  If you move this part of the code out
    // of main(), change the line below accordingly.
    return 1;
  }  

  double highMem = currentSize();

  ModalProblem *testCase;
  if (dimension == 1) {
    if (casePb == 1)
      testCase = new ModeLaplace1DQ1(Comm, Lx, nX);
    if (casePb == 2)
      testCase = new ModeLaplace1DQ2(Comm, Lx, nX);
  }
  if (dimension == 2) {
    if (casePb == 1)
      testCase = new ModeLaplace2DQ1(Comm, Lx, nX, Ly, nY);
    if (casePb == 2)
      testCase = new ModeLaplace2DQ2(Comm, Lx, nX, Ly, nY);
  }
  if (dimension == 3) {
    if (casePb == 1)
      testCase = new ModeLaplace3DQ1(Comm, Lx, nX, Ly, nY, Lz, nZ);
    if (casePb == 2)
      testCase = new ModeLaplace3DQ2(Comm, Lx, nX, Ly, nY, Lz, nZ);
  }

  highMem = (highMem < currentSize()) ? currentSize() : highMem;

  // Print some information on the problem to solve
  if (verbose*(myPid==0) > 0) {
    cout << endl;
    testCase->problemInfo();
  }

  // Get the stiffness and mass matrices
  const Epetra_Operator *K = testCase->getStiffness();
  const Epetra_Operator *M = testCase->getMass();

  // Get the map of the equations across processors
  const Epetra_Map Map = K->OperatorDomainMap();
  int localSize = Map.NumMyElements();
  int globalSize = Map.NumGlobalElements();

  // Get the first eigenvalue for the mass matrix
  numEigen = (numEigen > globalSize) ? globalSize : numEigen;

  double *globalWeight = 0;
  double globalMassMin = 0.0;

  if (dynamic_cast<ModeLaplace*>(testCase))
    globalMassMin = dynamic_cast<ModeLaplace*>(testCase)->getFirstMassEigenValue();
  else {
    // Input( Comm, Operator, verbose level, # levels, smoother, parameter, coarse solver)
    AMGOperator precML(Comm, M, 0, 10, 1, 3, 0);
    BlockPCGSolver linMassSolver(Comm, M, 1e-06, 20);
    linMassSolver.setPreconditioner(&precML);
    ARPACKm3 eigMassSolver(Comm, &linMassSolver, 1e-03, globalSize);
    double *lTmp = new double[6];
    Epetra_MultiVector QQ(Map, 6);
    QQ.Random();
    int massNEV = eigMassSolver.solve(1, QQ, lTmp);

    // FIXME (mfh 14 Jan 2011) I'm not sure if std::runtime_error is
    // the right exception to throw.  I'm just replacing exit(1) with
    // an exception, as per Trilinos coding standards.
    TEUCHOS_TEST_FOR_EXCEPTION( massNEV < 1, std::runtime_error, 
			"Error in the computation of smallest eigenvalue for "
			"the mass matrix. Output information from eigensolver"
			" = " << massNEV );
    globalMassMin = lTmp[0];
    delete[] lTmp;
  }

  // Note: The following definition of weight results from the definition
  //       of "Epetra_MultiVector::NormWeighted"
  globalWeight = new double[localSize];
  for (i = 0; i < localSize; ++i) {
    globalWeight[i] = sqrt(globalMassMin/globalSize);
  }

  if (verbose*(myPid==0) > 0) {
    cout << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " Smallest eigenvalue for the mass matrix: " << globalMassMin << endl;
    cout << endl;
  }

  // Define arrays for the eigensolver algorithm
  int qSize;
  switch (algo) {
    case LOBPCG_CHOL:
    case LOBPCG_LIGHT:
      qSize = dimSearch + numEigen;
      break;
    case LOBPCG_AK_CHOL:
      qSize = numEigen;
      dimSearch = numEigen;
      break;
    case ARPACK_ORIG:
    case ARPACK_RESI:
      if (dimSearch <= numEigen)
        dimSearch = 2*numEigen;
      qSize = dimSearch;
      break;
    case GALDAVIDSON:
      numBlock = (numBlock <= 0) ? 1 : numBlock;
      qSize = dimSearch*(numBlock+1);
      break;
    case BRQMIN_CHOL:
      qSize = dimSearch + numEigen;
      break;
    case DACG_CHOL:
      qSize = dimSearch + numEigen;
      break;
    case JD_PCG:
      numBlock = (numBlock <= 0) ? 1 : numBlock;
      qSize = dimSearch*(numBlock+1);
      break;
  }

  double *vQ = new (nothrow) double[qSize*localSize];
  assert(vQ != 0);

  double *lambda = new (nothrow) double[qSize];
  assert(lambda != 0);

  Epetra_MultiVector Q(View, Map, vQ, localSize, qSize);
  Q.Random();

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  char precLabel[100];
  MyIncompleteChol *ICT = 0;
  BlockPCGSolver *opStiffness = new BlockPCGSolver(Comm, K, tolCG, maxIterCG, 
                                                   (verbose) ? verbose-1 : 0);
  AMGOperator *precML = 0;
  switch (precond) {
    case NO_PREC:
      opStiffness->setPreconditioner(0);
      strcpy(precLabel, " No preconditioner");
      break;
    case AMG_POLYNOMIAL:
      // Use AMG preconditioner with polynomial smoother
      if (param < 1)
        param = 2;
      precML = new AMGOperator(Comm, K, verbose, 10, 1, param, 0);
      opStiffness->setPreconditioner(precML);
      strcpy(precLabel, " AMG with polynomial smoother");
      break;
    case AMG_GAUSSSEIDEL:
      // Use AMG preconditioner with Gauss-Seidel smoother
      precML = new AMGOperator(Comm, K, verbose, 10, 2, param, 0);
      opStiffness->setPreconditioner(precML);
      strcpy(precLabel, " AMG with Gauss-Seidel smoother");
      break;
    case INC_CHOL:
      // Use incomplete Cholesky with no-fill in
      Epetra_Operator *KOp = const_cast<Epetra_Operator*>(K);
      if (dynamic_cast<Epetra_CrsMatrix*>(KOp)) {
        ICT = new MyIncompleteChol(Comm, KOp, Epetra_MinDouble, param);
      }
      else {
        cerr << endl;
        cerr << " !!! The incomplete Cholesky factorization can not be done !!!\n";
        cerr << " !!! No preconditioner is used !!!\n";
        cerr << endl;
      }
      opStiffness->setPreconditioner(ICT);
      if (ICT)
        strcpy(precLabel, " Incomplete Cholesky factorization");
      break;
  }

  // Define the GeneralEigenSolver object
  ModalAnalysisSolver *mySolver;

  switch (algo) {
    case LOBPCG_CHOL:
      mySolver = new LOBPCG(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     dimSearch, tol, maxIter, verbose, globalWeight);
      break;
    case LOBPCG_LIGHT:
      mySolver = new LOBPCG_light(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     dimSearch, tol, maxIter, verbose, globalWeight);
      break;
    case LOBPCG_AK_CHOL:
      mySolver = new KnyazevLOBPCG(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     tol, maxIter, verbose, globalWeight);
      break;
    case ARPACK_ORIG:
      mySolver = new ARPACKm3(Comm, opStiffness, M, tol, maxIter, verbose);
      break;
    case ARPACK_RESI:
      mySolver = new ModifiedARPACKm3(Comm, opStiffness, M, tol, maxIter, verbose, 
                                      globalWeight);
      break;
    case GALDAVIDSON:
      mySolver = new Davidson(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     dimSearch, numBlock, tol, maxIter, verbose, globalWeight);
      break;
    case BRQMIN_CHOL:
      mySolver = new BRQMIN(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     dimSearch, tol, maxIter, verbose, globalWeight);
      break;
    case DACG_CHOL:
      mySolver = new BlockDACG(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     dimSearch, tol, maxIter, verbose, globalWeight);
      break;
    case JD_PCG:
      mySolver = new JDPCG(Comm, opStiffness, M, opStiffness->getPreconditioner(),
                     dimSearch, numBlock, tol, maxIter, maxIterCG, verbose, globalWeight);
      break;
  }

  // Solve the eigenproblem
  int knownEV = mySolver->solve(numEigen, Q, lambda);

  // Output information on simulation

  if (knownEV < 0) {
    if (myPid == 0) {
      cerr << endl;
      cerr << " !!! The eigensolver exited with error " << knownEV << " !!! " << endl;
      cerr << endl;
    }
  } // if (knownEV < 0)
  else {

    if (verbose*(myPid == 0) > 1) {
      cout << endl;
      cout << " ---------------------------------------------------- \n";
      cout << "   History of Computation                             \n";
      cout << " ---------------------------------------------------- \n";
      cout << endl; 
      mySolver->algorithmInfo();
      cout << " Preconditioner: " << precLabel << endl;
      cout << endl;
      testCase->problemInfo();
      cout << " Number of computed eigenmodes: " << knownEV << endl;
      mySolver->historyInfo();
    } // if (verbose*(myPid == 0) > 1)

    // Eigenpairs accuracy
    if (knownEV > 0) {
      int nb = (knownEV > numEigen) ? numEigen : knownEV; 
      Epetra_MultiVector copyQ(View, Q, 0, nb);
      if (myPid == 0) {
        cout << endl;
        cout << " ---------------------------------------------------- \n";
        cout << "   Eigenpairs accuracy                                \n";
        cout << " ---------------------------------------------------- \n";
        cout << endl;
        mySolver->algorithmInfo();
        cout << " Preconditioner: " << precLabel << endl;
        cout << endl;
        testCase->problemInfo();
        cout << " Number of computed eigenmodes: " << copyQ.NumVectors() << endl;
        cout.precision(2);
        cout.setf(ios::scientific, ios::floatfield);
        cout << " Tolerance on residuals: " << tol << endl;
        cout << endl;
        cout << " Smallest eigenvalue for the mass matrix: " << globalMassMin << endl;
        cout << endl;
      } // if (myPid == 0)
      testCase->eigenCheck(copyQ, lambda, globalWeight);
    }

    // Print out statistics: memory, time
    if (myPid == 0) {
      cout << endl;
      cout << " ---------------------------------------------------- \n";
      cout << "  Summary of statistics                               \n";
      cout << " ---------------------------------------------------- \n";
      cout << endl;
      mySolver->algorithmInfo();
      cout << " Preconditioner: " << precLabel << endl;
      cout << endl;
      testCase->problemInfo();
      cout << " Number of computed eigenmodes: " << knownEV << endl;
      cout.precision(2);
      cout.setf(ios::scientific, ios::floatfield);
      cout << " Tolerance on residuals: " << tol << endl;
      if ((algo == ARPACK_ORIG) || (algo == ARPACK_RESI)) {
        cout << " Size of Search Space: " << dimSearch << endl;
        cout << endl;
        cout.precision(2);
        cout.setf(ios::scientific, ios::floatfield);
        cout << " Tolerance on PCG solves: " << tolCG << endl;
        cout << endl;
        cout << " Minimum number of PCG iterations per solve: ";
        cout << opStiffness->getMinIter() << endl;
        cout.setf(ios::fixed, ios::floatfield);
        cout << " Average number of PCG iterations per solve: ";
        cout << opStiffness->getAvgIter() << endl;
        cout << " Maximum number of PCG iterations per solve: ";
        cout << opStiffness->getMaxIter() << endl;
      }
      cout << endl;
      if (precond == AMG_POLYNOMIAL) {
        cout << " Number of levels for AMG preconditioner: ";
        cout << precML->getAMG_NLevels() << endl;
        cout << " Polynomial degree of AMG preconditioner: " << param << endl;
        cout << endl;
      }
      if (precond == AMG_GAUSSSEIDEL) {
        cout << " Number of levels for AMG preconditioner: ";
        cout << precML->getAMG_NLevels() << endl;
        cout << endl;
      }
      if ((precond == INC_CHOL) && (ICT)) {
        cout << " Level of fill for incomplete Cholesky factorisation: ";
        cout << param << endl;
        cout << endl;
      }
      cout << " Number of processors: " << numProc << endl;
      cout << endl;
    } // if (myPid == 0)

    // Memory
    double maxHighMem = 0.0;
    Comm.MaxAll(&highMem, &maxHighMem, 1);
    if (myPid == 0) {
      cout << " --- Memory ---\n";
      cout << endl;
      cout.precision(2);
      cout.setf(ios::fixed, ios::floatfield);
      cout << " High water mark in set-up                           = (EST) ";
      cout << maxHighMem << " MB\n";
      cout << endl;
    } // if (myPid == 0)

    testCase->memoryInfo();

    if (myPid == 0) {
      cout.precision(2);
      cout.setf(ios::fixed, ios::floatfield);
      cout << " Memory requested per processor for eigenvectors     = (EST) ";
      cout << Q.GlobalLength()*numEigen*sizeof(double)/1024.0/1024.0/numProc;
      cout << " MB\n";
      cout << endl;
      cout << " Memory requested per processor for working space    = (EST) ";
      cout << Q.GlobalLength()*(Q.NumVectors()-numEigen)*sizeof(double)/1024.0/1024.0/numProc;
      cout << " MB\n";
      cout << endl;
    } // if (myPid == 0)

    mySolver->memoryInfo();

    mySolver->operationInfo();

    mySolver->timeInfo();

  } // if (knownEV < 0)

  // Release all objects
  delete opStiffness;
  delete ICT;
  delete precML;
  delete mySolver;
  
  delete[] lambda;
  delete[] vQ;
  delete[] globalWeight;

  delete testCase;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return 0;

}
