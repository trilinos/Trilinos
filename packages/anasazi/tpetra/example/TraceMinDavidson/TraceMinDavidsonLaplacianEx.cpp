// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example tests TraceMin-Davidson on the problem of finding the Fiedler
//  vector of graph Laplacian (input from a matrix market file)

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for TraceMin-Davidson solver
#include "AnasaziTraceMinDavidsonSolMgr.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide Anasazi with Tpetra adapters
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziOperator.hpp"

// Include header for Tpetra compressed-row storage matrix
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"

// Include headers for reading and writing matrix-market files
#include <MatrixMarket_Tpetra.hpp>

// Include header for sparse matrix operations
#include <TpetraExt_MatrixMatrix_def.hpp>

// Include header for Teuchos serial dense matrix
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Teuchos_ArrayViewDecl.hpp"

#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_ANASAZI_IFPACK2
  #include "Ifpack2_Factory.hpp"
  #include "Ifpack2_Preconditioner.hpp"
#endif

  using Teuchos::RCP;
  using std::cout;
  using std::cin;

  //
  // Specify types used in this example
  //
  using Scalar = double;
  using CrsMatrix = Tpetra::CrsMatrix<Scalar>;
  using Vector = Tpetra::Vector<Scalar>;
  using TMV = Tpetra::MultiVector<Scalar>;
  using TOP = Tpetra::Operator<Scalar>;
  using TMVT = Anasazi::MultiVecTraits<Scalar, TMV>;
  using TOPT = Anasazi::OperatorTraits<Scalar, TMV, TOP>;

void formLaplacian(const RCP<const CrsMatrix>& A, const bool weighted, const bool normalized, RCP<CrsMatrix>& L, RCP<Vector>& auxVec);

int main(int argc, char *argv[]) {
  //
  // Initialize the MPI session
  //
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();

  //
  // Get parameters from command-line processor
  //
  std::string inputFilename("/home/amklinv/matrices/mesh1em6_Laplacian.mtx");
  std::string outputFilename("/home/amklinv/matrices/mesh1em6_Fiedler.mtx");
  Scalar tol = 1e-6;
  int nev = 1;
  int blockSize = 1;
  bool usePrec = false;
  bool useNormalizedLaplacian = false;
  bool useWeightedLaplacian = false;
  bool verbose = true;
  std::string whenToShift = "Always";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("fin",&inputFilename, "Filename for Matrix-Market test matrix.");
  cmdp.setOption("fout",&outputFilename, "Filename for Fiedler vector.");
  cmdp.setOption("tolerance",&tol, "Relative residual used for solver.");
  cmdp.setOption("nev",&nev, "Number of desired eigenpairs.");
  cmdp.setOption("blocksize",&blockSize, "Number of vectors to add to the subspace at each iteration.");
  cmdp.setOption("precondition","no-precondition",&usePrec, "Whether to use a diagonal preconditioner.");
  cmdp.setOption("normalization","no-normalization",&useNormalizedLaplacian, "Whether to normalize the laplacian.");
  cmdp.setOption("weighted","unweighted",&useWeightedLaplacian, "Whether to normalize the laplacian.");
  cmdp.setOption("verbose","quiet",&verbose, "Whether to print a lot of info or a little bit.");
  cmdp.setOption("whenToShift",&whenToShift, "When to perform Ritz shifts. Options: Never, After Trace Levels, Always.");
  if(cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  //
  // Read the matrix from a file
  //
  RCP<const CrsMatrix> fileMat =
    Tpetra::MatrixMarket::Reader<CrsMatrix>::readSparseFile(inputFilename, comm);

  //
  // Form the graph Laplacian and get a const pointer to the data
  //
  RCP<CrsMatrix> L;
  RCP<Vector> auxVec;
  formLaplacian(fileMat, useWeightedLaplacian, useNormalizedLaplacian, L, auxVec);
  RCP<const CrsMatrix> K = L;

  //
  // Compute the norm of the matrix
  //
  Scalar mat_norm = K->getFrobeniusNorm();

  //
  // ************************************
  // Start the block Arnoldi iteration
  // ************************************
  //
  //  Variables used for the Block Arnoldi Method
  //
  int verbosity;
  int numRestartBlocks = 2*nev/blockSize;
  int numBlocks = 10*nev/blockSize;
  if(verbose)
    verbosity = Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
  else
    verbosity = Anasazi::TimingDetails;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );                  // How much information should the solver print?
  MyPL.set( "Saddle Solver Type", "Projected Krylov"); // Use projected minres/gmres to solve the saddle point problem
  MyPL.set( "Block Size", blockSize );                 // Add blockSize vectors to the basis per iteration
  MyPL.set( "Convergence Tolerance", tol*mat_norm );   // How small do the residuals have to be
  MyPL.set( "Relative Convergence Tolerance", false);  // Don't scale residuals by eigenvalues (when checking for convergence)
  MyPL.set( "Use Locking", true);                      // Use deflation
  MyPL.set( "Relative Locking Tolerance", false);      // Don't scale residuals by eigenvalues (when checking whether to lock a vector)
  MyPL.set("Num Restart Blocks", numRestartBlocks);    // When we restart, we start back up with 2*nev blocks
  MyPL.set("Num Blocks", numBlocks);                   // Maximum number of blocks in the subspace
  MyPL.set("When To Shift", whenToShift);
  MyPL.set("Saddle Solver Type", "Block Diagonal Preconditioned Minres");

  //
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  //
  RCP<TMV> ivec = Teuchos::rcp( new TMV(K->getRowMap(), blockSize) );
  TMVT::MvRandom( *ivec );

  //
  // Create the eigenproblem
  //
  RCP<Anasazi::BasicEigenproblem<Scalar,TMV,TOP> > MyProblem =
      Teuchos::rcp( new Anasazi::BasicEigenproblem<Scalar,TMV,TOP>(K, ivec) );

  //
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  //
  MyProblem->setHermitian(true);

  //
  // Set the number of eigenvalues requested
  //
  MyProblem->setNEV( nev );

  if(usePrec)
  {
    #ifdef HAVE_ANASAZI_IFPACK2
    //
    // Construct a diagonal preconditioner
    //
    Ifpack2::Factory factory;
    const std::string precType = "RELAXATION";
    Teuchos::ParameterList PrecPL;
    PrecPL.set( "relaxation: type", "Jacobi");
    RCP<Ifpack2::Preconditioner<Scalar> > Prec = factory.create(precType, K);
    assert(Prec != Teuchos::null);
    Prec->setParameters(PrecPL);
    Prec->initialize();
    Prec->compute();

    //
    // Tell the problem about the preconditioner
    //
    MyProblem->setPrec(Prec);
    #else
    if(myRank == 0)
      cout << "You did not build Trilinos with Ifpack2 preconditioning enabled.  Please either\n1. Reinstall Trilinos with Ifpack2 enabled\n2. Try running this driver again without preconditioning enabled\n";
    return -1;
    #endif
  }

  //
  // Since we are computing the Fiedler vector, we want it to be orthogonal to the null space of the laplacian
  //
  MyProblem->setAuxVecs(auxVec);

  //
  // Inform the eigenproblem that you are finished passing it information
  //
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (myRank == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << std::endl;
    }
    return -1;
  }

  //
  // Initialize the TraceMin-Davidson solver
  //
  Anasazi::Experimental::TraceMinDavidsonSolMgr<Scalar, TMV, TOP> MySolverMgr(MyProblem, MyPL);

  //
  // Solve the problem to the specified tolerances
  //
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && myRank == 0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
  }
  else if (myRank == 0)
    cout << "Anasazi::EigensolverMgr::solve() returned converged." << std::endl;

  //
  // Get the eigenvalues and eigenvectors from the eigenproblem
  //
  Anasazi::Eigensolution<Scalar,TMV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
  RCP<TMV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  //
  // Compute the residual, just as a precaution
  //
  if (numev > 0) {

    Teuchos::SerialDenseMatrix<int,Scalar> T(numev,numev);
    TMV tempvec(K->getRowMap(), TMVT::GetNumberVecs( *evecs ));
    std::vector<Scalar> normR(sol.numVecs);
    TMV Kvec( K->getRowMap(), TMVT::GetNumberVecs( *evecs ) );

    TOPT::Apply( *K, *evecs, Kvec );
    TMVT::MvTransMv( 1.0, Kvec, *evecs, T );
    TMVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, Kvec );
    TMVT::MvNorm( Kvec, normR );

    if (myRank == 0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Error"<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      for (int i=0; i<numev; i++) {
        cout<<std::setw(16)<<T(i,i)
          <<std::setw(16)<<normR[i]/mat_norm
          <<std::endl;
      }
      cout<<"------------------------------------------------------"<<std::endl;
    }
  }

  //
  // Write the Fiedler vector to a file
  //
  if (numev > 0) {
    Tpetra::MatrixMarket::Writer<CrsMatrix>::writeDenseFile(outputFilename,evecs,"","Fiedler vector of "+inputFilename);
  }

  return 0;
}



/* Form the graph Laplacian and also return the null space
 * This function assumes the graph consists of only one strongly connected component
 *
 * If we are forming the weighted-Laplacian,
 * L(i,j) = \sum abs(A(i,j)), if i == j
 * L(i,j) = -abs(A(i,j)), if i ~= j
 *
 * If we are forming the unweighted-Laplacian,
 * L(i,i) = \sum abs(L(i,j)), if i == j
 * L(i,j) = -1, if i ~= j and A(i,j) ~= 0
 *
 * Then we normalize (if desired)
 *
 * This matrix will always be symmetric positive semidefinite with one zero eigenvalue
 *
 * We also want to compute the vector corresponding to the zero eigenvalue.
 * We know the Fiedler vector is orthogonal to it, so we should pass that
 * information to TraceMin via setAuxVecs.
 *
 * If we are not using a normalized Laplacian,
 * auxVec = ones(n,1)
 *
 * Otherwise,
 * auxVec[i] = sqrt(L(i,i))
 */
void formLaplacian(const RCP<const CrsMatrix>& A, const bool weighted, const bool normalized, RCP<CrsMatrix>& L, RCP<Vector>& auxVec)
{
  //
  // Set a few convenient constants
  //
  using SCT = Teuchos::ScalarTraits<Scalar>;
  using LO = Tpetra::Map<>::local_ordinal_type;
  using GO = Tpetra::Map<>::global_ordinal_type;
  Scalar ONE = SCT::one();
  Scalar ZERO = SCT::zero();

  std::vector<GO> diagIndex(static_cast<GO>(1));
  std::vector<LO> diagIndexLcl(static_cast<LO>(1));
  std::vector<Scalar> value(1,ONE);
  Teuchos::ArrayView<const GO> cols(diagIndex);
  Teuchos::ArrayView<const LO> colsLcl(diagIndexLcl);
  Teuchos::ArrayView<const Scalar> vals(value);

  //
  // Get the number of rows of A
  //
  const GO n = A->getGlobalNumRows();

  //
  // Construct the unweighted Laplacian
  // Note: If your graph is not connected, TraceMin-Davidson will not work
  //
  // We assume the matrix is nonsymmetric for generality.
  // If it is in fact symmetric, this step is not necessary
  // L = A + A'
  //
  L = Tpetra::MatrixMatrix::add(ONE,false,*A,ONE,true,*A);
  RCP<const Tpetra::Map<> > rowMap = L->getRowMap();

  // This line tells L that we are going to modify its entries
  L->resumeFill();

  if(weighted)
  {
    using indices_view = typename CrsMatrix::nonconst_global_inds_host_view_type;
    using values_view  = typename CrsMatrix::nonconst_values_host_view_type;
    indices_view colIndices("colIndices");
    values_view values("values");

    // This vector holds the diagonal
    RCP<Vector> diagonal = Teuchos::rcp(new Vector(rowMap));

    //
    // Set the sign of each off-diagonal entry to -
    // Also delete diagonal entries
    // Also compute the sum of off-diagonal entries
    //
    for(GO i=0; i<n; i++)
    {
      // If this process does not own that row of the matrix, do nothing
      // Each process handles its own rows
      if(rowMap->isNodeGlobalElement(i))
      {
        // Figure out how many entries are in the row
        size_t numentries = L->getNumEntriesInGlobalRow(i);
        Kokkos::resize(colIndices,numentries);
        Kokkos::resize(values,numentries);

        // Get a copy of row i
        L->getGlobalRowCopy(i,colIndices,values,numentries);

        for(size_t j=0; j<colIndices.size(); j++)
        {
          // Delete diagonal entries
          if(i == rowMap->getGlobalElement(colIndices[j]))
            values[j] = ZERO;
          // Set the sign of off-diagonal elements
          else
            values[j] = -abs(values[j]);

          // Update the diagonal
          diagonal->sumIntoGlobalValue(i,-values[j]);
        }

        // Reinsert the updated row
        L->replaceGlobalValues(i, colIndices, values);
      }
    }

    // Get a view of the LOCAL diagonal
    Teuchos::ArrayRCP<const Scalar> diagView = diagonal->getData();
    using size_type = Teuchos::ArrayRCP<const Scalar>::size_type;

    // Create the auxiliary vector and set its values
    auxVec = Teuchos::rcp(new Vector(rowMap,false));
    if (normalized) {
      // auxVec[i] = sqrt(L(i,i))
      for (size_type i = 0; i < diagView.size (); ++i) {
        auxVec->replaceLocalValue (static_cast<LO> (i), sqrt (diagView[i]));
      }
    }
    else {
      auxVec->putScalar(ONE); // auxVec[i] = ONE for all i
    }

    // Compute the norm of the vector
    Scalar vecNorm = auxVec->norm2();

    // Normalize the vector
    auxVec->scale(ONE/vecNorm);

    // Finish computing the Laplacian
    if (normalized) {
      // Compute D^{-1/2}
      Vector scaleVec (rowMap, false);
      for(size_type i=0; i<diagView.size(); i++)
      {
        scaleVec.replaceLocalValue(static_cast<LO>(i),ONE/sqrt(diagView[i]));
      }

      // Must be called before scaling
      // Tells the matrix we're (temporarily) done adding things to it
      L->fillComplete();

      // Premultiply L by D^{-1/2}
      L->leftScale(scaleVec);

      // Postmultiply L by D^{-1/2}
      L->rightScale(scaleVec);
      L->resumeFill();

      // Set diagonal entries to 1
      for(GO i=0; i<n; i++)
      {
        diagIndex[0] = i;
        if(rowMap->isNodeGlobalElement(i)) L->replaceGlobalValues(i,cols,vals);
      }
    }
    else
    {
      // Set diagonal entries
      for(size_type i=0; i<diagView.size(); i++)
      {
        diagIndex[0] = rowMap->getGlobalElement(i);
        value[0] = diagView[i];
        L->replaceLocalValues(i,colsLcl,vals);
      }
    }

    // Tell the matrix we're done modifying its entries
    L->fillComplete();
  }
  else {
    //
    // Construct the adjacency matrix by removing the diagonal and setting all off-diagonal entries to 1
    //

    // Here, we ensure that space is reserved for the diagonal entries
    for (GO i = 0; i < n; ++i) {
      diagIndex[0] = i;
      if(rowMap->isNodeGlobalElement(i)) {
        L->replaceGlobalValues(i,cols,vals);
      }
    }

    // Set all entries of the matrix to -1
    L->setAllToScalar(-ONE);

    // Tell the matrix we're done inserting things so that we can scale the entries
    L->fillComplete();

    // Create the auxiliary vector and set its values
    auxVec = Teuchos::rcp( new Vector (rowMap,false) );
    if (normalized) {
      // Compute the degree of each vertex of the graph
      for (GO i = 0; i < n; ++i) {
        // If this process does not own that row of the matrix, do nothing
        // Each process handles its own rows
        if (rowMap->isNodeGlobalElement(i)) {
          Scalar temp;

          // Because we insisted on having a diagonal, the degree is really 1 less than
          // the number of entries in the row
          size_t nnzInRow = L->getNumEntriesInGlobalRow(i) - static_cast<size_t> (1);
          temp = sqrt(nnzInRow);

          // auxVec[i] = sqrt(degree of node i)
          auxVec->replaceGlobalValue(i,temp);
        }
      }
    }
    else {
      auxVec->putScalar(ONE); // auxVec[i] = 1 for all i
    }

    // Compute the norm of the vector
    Scalar vecNorm = auxVec->norm2();

    // Normalize the vector
    auxVec->scale(ONE/vecNorm);

    // Finish computing the Laplacian
    if (normalized) {
      // Compute the degree of each node so we can normalize the Laplacian
      // normalizedL = D^{-1/2} L D^{-1/2}
      Vector scalars(rowMap,false);
      for (GO i = 0; i < n; ++i) {
        // If this process does not own that row of the matrix, do nothing
        // Each process handles its own rows
        if(rowMap->isNodeGlobalElement(i)) {
          Scalar temp;

          // Because we insisted on having a diagonal, the degree is really 1 less than
          // the number of entries in the row
          size_t nnzInRow = L->getNumEntriesInGlobalRow(i) - static_cast<size_t> (1);
          temp = ONE/sqrt(nnzInRow);

          // scalars[i] = 1/sqrt(degree of node i)
          scalars.replaceGlobalValue(i,temp);
        }
      }

      // Premultiply L by D^{-1/2}
      L->leftScale(scalars);

      // Postmultiply L by D^{-1/2}
      L->rightScale(scalars);
      L->resumeFill();

      // Set diagonal entries to 1
      for (GO i = 0; i < n; ++i) {
        diagIndex[0] = i;
        if(rowMap->isNodeGlobalElement(i)) L->replaceGlobalValues(i,cols,vals);
      }
    }
    else {
      L->resumeFill();
      for (GO i = 0; i < n; ++i) {
        // If this process does not own that row of the matrix, do nothing
        // Each process handles its own rows
        if(rowMap->isNodeGlobalElement(i)) {
          // Because we insisted on having a diagonal, the degree is really 1 less than
          // the number of entries in the row
          size_t nnzInRow = L->getNumEntriesInGlobalRow(i) - static_cast<size_t> (1);

          // L[i,i] = degree of node i
          diagIndex[0] = i;
          value[0] = nnzInRow;
          L->replaceGlobalValues(i,cols,vals);
        }
      }
    }

    // Tell the matrix we're done modifying its entries
    L->fillComplete();
  }
}
