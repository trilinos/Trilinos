// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
// @HEADER
//  This example demonstrates TraceMin-Davidson's ability to converge
//  without the matrix being explicitly stored

// TODO: Fix the maximum number of iterations

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for TraceMin-Davidson solver
#include "AnasaziTraceMinSolMgr.hpp"

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
  // Instead of constantly typing Tpetra::MultiVector<Scalar, ...>,
  // we can type TMV.
  //
  typedef double                                        Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType  Magnitude;
  typedef Tpetra::Map<>::local_ordinal_type             LocalOrdinal;
  typedef Tpetra::Map<>::global_ordinal_type            GlobalOrdinal;
  typedef Tpetra::Map<>               Map;
  typedef Tpetra::MultiVector<Scalar> TMV;
  typedef Tpetra::Operator<Scalar>    TOP;

//
// Define a class for our user-defined operator.
// In this case, it is the tridiagonal matrix [-1,2,-1].
// You can define it to be whatever you like.
//
// In general, Trilinos does NOT require the user to explicitly deal
// with MPI.  If you want to define your own operator though, there's
// no getting around it.  Fortunately, Trilinos makes this relatively
// straightforward with the use of Tpetra::Import objects.  All you
// have to do is define your initial data distribution (which is a
// block row distribution here), and the data distribution you need to
// perform the operations of your matvec.
//
// For instance, when performing a matvec with a tridiagonal matrix (with a block
// row distribution), each process needs to know the last element owned by the
// previous process and the first element owned by the next process.
//
// If you are only interested in running the code sequentially, you can safely
// ignore everything here regarding maps and importers
//
class MyOp : public virtual TOP {
private:
  typedef Tpetra::Import<> Import;

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
    LocalOrdinal nlocal = opMap_->getNodeNumElements();

    //
    // Define the distribution that you need for the matvec.
    // When you define this for your own operator, it is helpful to draw pictures
    // on a sheet of paper to keep track of who needs to receive which elements.
    // Here, each process needs to receive one element from each of its neighbors.
    //

    // All processes but the first will receive one entry from the
    // previous process.
    if (myRank_ > 0) {
      nlocal++;
    }
    // All processes but the last will receive one entry from the next
    // process.
    if (myRank_ < numProcs_-1) {
      nlocal++;
    }
    // Construct a list of columns where this process has nonzero elements.
    // For our tridiagonal matrix, this is firstRowItOwns-1:lastRowItOwns+1.
    std::vector<GlobalOrdinal> indices;
    indices.reserve (nlocal);
    // The first process is a special case...
    if (myRank_ > 0) {
      indices.push_back (opMap_->getMinGlobalIndex () - 1);
    }
    for (GlobalOrdinal i = opMap_->getMinGlobalIndex ();
         i <= opMap_->getMaxGlobalIndex (); ++i) {
      indices.push_back (i);
    }
    // So is the last process...
    if (myRank_ < numProcs_-1) {
      indices.push_back (opMap_->getMaxGlobalIndex () + 1);
    }

    // Wrap our vector in an array view, which is like a raw pointer.
    Teuchos::ArrayView<const GlobalOrdinal> elementList (indices);

    // Make a Map for handling the redistribution.  There will be some
    // redundancies (i.e., some of the entries will be owned by
    // multiple MPI processes).
    GlobalOrdinal numGlobalElements = n + 2*(numProcs_-1);
    redistMap_ = rcp (new Map (numGlobalElements, elementList, 0, comm));

    // Make an Import object that describes how data will be
    // redistributed.  It takes a Map describing who owns what
    // originally, and a Map that describes who you WANT to own what.
    importer_= rcp (new Import (opMap_, redistMap_));
  }

  //
  // These functions are required since we inherit from Tpetra::Operator
  //

  // Destructor
  virtual ~MyOp() {};

  // Returns the maps
  RCP<const Map> getDomainMap() const { return opMap_; };
  RCP<const Map> getRangeMap() const { return opMap_; };

  // Compute Y = alpha Op X + beta Y.
  //
  // TraceMin will never use alpha ~= 1 or beta ~= 0,
  // so we have ignored those options for simplicity.
  void
  apply (const TMV& X,
         TMV& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const
  {
    typedef Teuchos::ScalarTraits<Scalar> SCT;

    //
    // Let's make sure alpha is 1 and beta is 0...
    // This will throw an exception if that is not the case.
    //
    TEUCHOS_TEST_FOR_EXCEPTION(
      alpha != SCT::one() || beta != SCT::zero(), std::logic_error,
      "MyOp::apply was given alpha != 1 or beta != 0.  It does not currently "
      "implement either of those two cases.");

    // Get the number of local rows
    const LocalOrdinal nlocRows = X.getLocalLength();

    // Get the number of vectors
    const int numVecs = X.getNumVectors();

    // Make a multivector for holding the redistributed data
    RCP<TMV> redistData = rcp(new TMV(redistMap_, numVecs));

    // Redistribute the data.
    // This will do all the necessary communication for you.
    // All processes now own enough data to do the matvec.
    redistData->doImport(X, *importer_, Tpetra::INSERT);

    // Perform the matvec with the data we now locally own
    //
    // For each column...
    for (int c = 0; c < numVecs; ++c) {
      // Get a view of the desired column
      Teuchos::ArrayRCP<Scalar> colView = redistData->getDataNonConst (c);

      int offset;
      // Y[0,c] = -colView[0] + 2*colView[1] - colView[2] (using local indices)
      if (myRank_ > 0) {
        Y.replaceLocalValue (0, c, -colView[0] + 2*colView[1] - colView[2]);
        offset = 0;
      }
      // Y[0,c] = 2*colView[1] - colView[2] (using local indices)
      else {
        Y.replaceLocalValue(0, c, 2*colView[0] - colView[1]);
        offset = 1;
      }

      // Y[r,c] = -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset]
      for(LocalOrdinal r=1; r<nlocRows-1; r++) {
        Y.replaceLocalValue(r, c, -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset]);
      }

      // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset] - colView[nlocRows+1-offset]
      if (myRank_ < numProcs_-1) {
        Y.replaceLocalValue(nlocRows-1, c, -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset] - colView[nlocRows+1-offset]);
      }
      // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
      else {
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
  typedef Anasazi::MultiVecTraits<Scalar, TMV> MVT;
  typedef Anasazi::OperatorTraits<Scalar, TMV, TOP> OPT;

  //
  // Initialize the MPI session
  //
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();

  //
  // Get parameters from command-line processor
  //
  Scalar tol = 1e-5;
  GlobalOrdinal n = 50;
  int nev = 4;
  bool verbose = true;
  std::string saddleSolType = "Projected Krylov";
  std::string which = "SM";
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption("which",&which, "Which eigenpairs we seek. Options: SM, LM.");
  cmdp.setOption ("saddleSolType", &saddleSolType, "Saddle Solver Type. "
                  "Options: Projected Krylov, Schur Complement, Block Diagonal "
                  "Preconditioned Minres.");
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
  MyPL.set( "Block Size", 2*nev );                 // Add blockSize vectors to the basis per iteration
  MyPL.set( "Convergence Tolerance", tol);             // How small do the residuals have to be?
  MyPL.set("Which", which);
  MyPL.set("When To Shift", "Never");

  //
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  // We are giving it random entries.
  //
  RCP<TMV> ivec = rcp( new TMV(K->getDomainMap(), nev) );
  MVT::MvRandom( *ivec );

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
  Anasazi::Experimental::TraceMinSolMgr<Scalar, TMV, TOP> MySolverMgr(MyProblem, MyPL);

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
    for (int i = 0; i < numev; ++i) {
      T(i,i) = evals[i].realpart;
    }
    TMV tempvec(K->getDomainMap(), MVT::GetNumberVecs( *evecs ));
    std::vector<Scalar> normR(numev);
    TMV Kvec( K->getRangeMap(), MVT::GetNumberVecs( *evecs ) );

    OPT::Apply( *K, *evecs, Kvec );
    MVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );

    std::vector<Scalar> true_eigs (numev);
    const double PI = 3.141592653589793238463;
    if (which == "LM") {
      for (size_t i = static_cast<size_t> (n + 1 - numev);
           i <= static_cast<size_t> (n); ++i) {
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
    for (int i = 0; i<numev; ++i) {
      if (normR[i]/T(i,i) > tol) {
        testFailed = true;
        if(myRank == 0)
          cout << "Test is about to fail because "
               << normR[i]/T(i,i) << " > " << tol << endl;
      }

      if (std::abs(T(i,i)-true_eigs[i]) > tol) {
        testFailed = true;
        if (myRank == 0) {
          cout << "Test is about to fail because "
               << std::abs (T(i,i) - true_eigs[i]) << " > " << tol << endl;
        }
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
    cout << myRank << ": TEST FAILED\n";
    if(myRank == 0)
      cout << "End Result: TEST FAILED" << endl;
    return -1;
  }
    cout << myRank << ": TEST PASSED\n";

  if(myRank == 0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
}
