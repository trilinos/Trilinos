// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example demonstrates TraceMin-Davidson's ability to converge
//  without the matrix being explicitly stored

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for TraceMin-Davidson solver
#include "AnasaziTraceMinDavidsonSolMgr.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide Anasazi with Tpetra adapters
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziOperator.hpp"

// Include headers for Tpetra
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

// Include header for Teuchos serial dense matrix
#include "Teuchos_SerialDenseMatrix.hpp"

// Include header for Teuchos array views
#include "Teuchos_ArrayViewDecl.hpp"

// Include header for Teuchos parameter lists
#include "Teuchos_ParameterList.hpp"

  //
  // These statements tell the compiler where to look for cout, etc.
  //
  using Teuchos::rcp;
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  //
  // Specify types used in this example
  //
  typedef double                                    Scalar;
  typedef Teuchos::ScalarTraits<Scalar>             TSCT;
  typedef Tpetra::Map<>::local_ordinal_type         LocalOrdinal;
  typedef Tpetra::Map<>::global_ordinal_type        GlobalOrdinal;
  typedef Tpetra::Map<>                             Map;
  typedef Tpetra::MultiVector<Scalar>               TMV;
  typedef Tpetra::Operator<Scalar>                  TOP;
  typedef Anasazi::MultiVecTraits<Scalar, TMV>      TMVT;
  typedef Anasazi::OperatorTraits<Scalar, TMV, TOP> TOPT;
  typedef Tpetra::Import<>                          Import;



//
// Define a class for our user-defined operator.
// In this case, it is the tridiagonal matrix [-1,2,-1].
// You can define it to be whatever you like.
//
// In general, Trilinos does NOT require the user to explicitly deal with MPI.
// If you want to define your own operator though, there's no getting around it.
// Fortunately, Trilinos makes this relatively straightforward with the use of
// Importer objects.  All you have to do is define your initial data distribution
// (which is a block row distribution here), and the data distribution you need
// to perform the operations of your matvec.
//
// For instance, when performing a matvec with a tridiagonal matrix (with a block
// row distribution), each process needs to know the last element owned by the
// previous process and the first element owned by the next process.
//
// If you are only interested in running the code sequentially, you can safely
// ignore everything here regarding maps and importers
//
class MyOp : public virtual TOP {
public:
  //
  // Constructor
  //
  MyOp(const GlobalOrdinal n, const RCP< const Teuchos::Comm<int> > comm)
  {
    //
    // Construct a map for our block row distribution
    //
    opMap_ = rcp( new Map(n, 0, comm) );

    //
    // Get the rank of this process and the number of processes
    // We're going to have to do something special with the first and last processes
    //
    myRank_ = comm->getRank();
    numProcs_ = comm->getSize();

    //
    // Get the local number of rows
    //
    LocalOrdinal nlocal = opMap_->getLocalNumElements();

    //
    // Define the distribution that you need for the matvec.
    // When you define this for your own operator, it is helpful to draw pictures
    // on a sheet of paper to keep track of who needs to receive which elements.
    // Here, each process needs to receive one element from each of its neighbors.
    //

    // All processes but the first will receive one element from the previous process
    if(myRank_ > 0) nlocal++;
    // All processes but the last will receive one element from the next process
    if(myRank_ < numProcs_-1) nlocal++;
    // Construct a list of columns where this process has nonzero elements
    // For our tridiagonal matrix, this is firstRowItOwns-1:lastRowItOwns+1
    std::vector<GlobalOrdinal> indices;
    indices.reserve(nlocal);
    // The first process is a special case...
    if(myRank_ > 0) indices.push_back(opMap_->getMinGlobalIndex()-1);
    for(GlobalOrdinal i=opMap_->getMinGlobalIndex(); i<=opMap_->getMaxGlobalIndex(); i++)
      indices.push_back(i);
    // So is the last process...
    if(myRank_ < numProcs_-1) indices.push_back(opMap_->getMaxGlobalIndex()+1);

    //
    // Wrap our vector in an array view, which is like a pointer
    //
    Teuchos::ArrayView<const GlobalOrdinal> elementList(indices);

    //
    // Make a map for handling the redistribution
    // There will be some redundancies (i.e. some of the elements will be owned by multiple processes)
    //
    GlobalOrdinal numGlobalElements = n + 2*(numProcs_-1);
    redistMap_ = rcp(new Map(numGlobalElements, elementList, 0, comm));

    //
    // Make an Import object that describes how data will be redistributed.
    // It takes a map describing who owns what originally, and a map that describes who you WANT
    // to own what.
    //
    importer_= rcp(new Import(opMap_, redistMap_));
  };



  //
  // These functions are required since we inherit from Tpetra::Operator
  //

  // Destructor
  virtual ~MyOp() {};

  // Returns the maps
  RCP<const Map> getDomainMap() const { return opMap_; };
  RCP<const Map> getRangeMap() const { return opMap_; };

  //
  // Computes Y = alpha Op X + beta Y
  // TraceMin will never use alpha ~= 1 or beta ~= 0,
  // so we have ignored those options for simplicity.
  //
  void apply(const TMV& X, TMV& Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=TSCT::one(), Scalar beta=TSCT::zero()) const
  {
    //
    // Let's make sure alpha is 1 and beta is 0...
    // This will throw an exception if that is not the case.
    //
    TEUCHOS_TEST_FOR_EXCEPTION(alpha != TSCT::one() || beta != TSCT::zero(),std::invalid_argument,
           "MyOp::apply was given alpha != 1 or beta != 0. That's not supposed to happen.");

    //
    // Get the number of local rows
    //
    LocalOrdinal nlocRows = X.getLocalLength();

    //
    // Get the number of vectors
    //
    int numVecs = X.getNumVectors();

    //
    // Make a multivector for holding the redistributed data
    //
    RCP<TMV> redistData = rcp(new TMV(redistMap_, numVecs));

    //
    // Redistribute the data.
    // This will do all the necessary communication for you.
    // All processes now own enough data to do the matvec.
    //
    redistData->doImport(X, *importer_, Tpetra::INSERT);

    //
    // Perform the matvec with the data we now locally own
    //
    // For each column...
    for(int c=0; c<numVecs; c++)
    {
      // Get a view of the desired column
      Teuchos::ArrayRCP<Scalar> colView = redistData->getDataNonConst(c);

      int offset;
      // Y[0,c] = -colView[0] + 2*colView[1] - colView[2] (using local indices)
      if(myRank_ > 0)
      {
        Y.replaceLocalValue(0, c, -colView[0] + 2*colView[1] - colView[2]);
        offset = 0;
      }
      // Y[0,c] = 2*colView[1] - colView[2] (using local indices)
      else
      {
        Y.replaceLocalValue(0, c, 2*colView[0] - colView[1]);
        offset = 1;
      }

      // Y[r,c] = -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset]
      for(LocalOrdinal r=1; r<nlocRows-1; r++)
      {
        Y.replaceLocalValue(r, c, -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset]);
      }

      // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset] - colView[nlocRows+1-offset]
      if(myRank_ < numProcs_-1)
      {
        Y.replaceLocalValue(nlocRows-1, c, -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset] - colView[nlocRows+1-offset]);
      }
      // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
      else
      {
        Y.replaceLocalValue(nlocRows-1, c, -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]);
      }
    }
  }



private:
  RCP<const Map> opMap_, redistMap_;
  RCP<const Import> importer_;
  int myRank_, numProcs_;
};



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
  Scalar tol = 1e-5;
  GlobalOrdinal n = 100;
  int nev = 4;
  int blockSize = 2;
  bool verbose = true;
  std::string saddleSolType = "Projected Krylov";
  std::string which = "SM";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("which",&which, "Which eigenpairs we seek. Options: SM, LM.");
  cmdp.setOption("saddleSolType", &saddleSolType, "Saddle Solver Type. Options: Projected Krylov, Schur Complement, Block Diagonal Preconditioned Minres.");
  if(cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  //
  // Construct the operator.
  // Note that the operator does not have to be an explicitly stored matrix.
  // Here, we are using our user-defined operator.
  //
  RCP<MyOp> K = rcp(new MyOp(n, comm));

  //
  // ************************************
  // Start the TraceMin-Davidson iteration
  // ************************************
  //
  //  Variables used for the TraceMin-Davidson Method
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
  MyPL.set( "Saddle Solver Type", saddleSolType);      // Use projected minres to solve the saddle point problem
  MyPL.set( "Block Size", blockSize );                 // Add blockSize vectors to the basis per iteration
  MyPL.set( "Convergence Tolerance", tol);             // How small do the residuals have to be?
  MyPL.set("Num Restart Blocks", numRestartBlocks);    // When we restart, this is how many blocks we keep
  MyPL.set("Num Blocks", numBlocks);                   // Maximum number of blocks in the subspace
  MyPL.set("Which", which);

  //
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  // We are giving it random entries.
  //
  RCP<TMV> ivec = rcp( new TMV(K->getDomainMap(), numRestartBlocks*blockSize) );
  TMVT::MvRandom( *ivec );

  //
  // Create the eigenproblem
  //
  RCP<Anasazi::BasicEigenproblem<Scalar,TMV,TOP> > MyProblem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,TMV,TOP>(K, ivec) );

  //
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  //
  MyProblem->setHermitian(true);

  //
  // Set the number of eigenvalues requested
  //
  MyProblem->setNEV( nev );

  //
  // Inform the eigenproblem that you are finished passing it information
  //
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (myRank == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
    return -1;
  }

  bool testFailed = false;

  //
  // Initialize the TraceMin-Davidson solver
  //
  Anasazi::Experimental::TraceMinDavidsonSolMgr<Scalar, TMV, TOP> MySolverMgr(MyProblem, MyPL);

  //
  // Solve the problem to the specified tolerances
  //
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if(returnCode != Anasazi::Converged) testFailed = true;

  if (returnCode != Anasazi::Converged && myRank == 0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }
  else if (myRank == 0)
    cout << "Anasazi::EigensolverMgr::solve() returned converged." << endl;

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
    for (int i=0; i < numev; ++i) {
      T(i,i) = evals[i].realpart;
    }
    TMV tempvec (K->getDomainMap (), TMVT::GetNumberVecs (*evecs));
    std::vector<Scalar> normR (numev);
    TMV Kvec (K->getRangeMap (), TMVT::GetNumberVecs (*evecs));

    TOPT::Apply( *K, *evecs, Kvec );
    TMVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, Kvec );
    TMVT::MvNorm( Kvec, normR );

    std::vector<Scalar> true_eigs(numev);
    const double PI  =3.141592653589793238463;
    if (which == "LM") {
      for (int i = n + 1 - numev; i <= n; ++i) {
        Scalar omega = i*PI/(2*n+2);
        true_eigs[i-(n+1-numev)] = 4*sin(omega)*sin(omega);
      }
    }
    else {
      for (int i = 1; i <= numev; ++i) {
        Scalar omega = i*PI/(2*n+2);
        true_eigs[i-1] = 4*sin(omega)*sin(omega);
      }
    }
    for (int i=0; i<numev; ++i) {
      if (normR[i]/T(i,i) > tol) {
        testFailed = true;
      }
      if (std::abs(T(i,i)-true_eigs[i]) > tol) {
        testFailed = true;
      }
    }

    if (myRank == 0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
      cout<<"------------------------------------------------------"<<endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Error"<<endl;
      cout<<"------------------------------------------------------"<<endl;
      for (int i=0; i<numev; i++) {
        cout<<std::setw(16)<<T(i,i)
          <<std::setw(16)<<normR[i]/T(i,i)
          <<endl;
      }
      cout<<"------------------------------------------------------"<<endl;
    }
  }

  if(testFailed) {
    if(myRank == 0)
      cout << "End Result: TEST FAILED" << endl;
    return -1;
  }

  if(myRank == 0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
}
