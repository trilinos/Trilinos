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
#include "Tpetra_DefaultPlatform.hpp"
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
  // These statements tell the compiler where to look for cin, cout, etc.
  //
  using Teuchos::rcp;
  using Teuchos::RCP;
  using std::cout;
  using std::cin;
  using std::endl;



  //
  // Specify types used in this example
  // Instead of constantly typing Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>, 
  // we can type MV.
  //                                   
  typedef double                                                       Scalar;
  typedef Teuchos::ScalarTraits<Scalar>                                SCT;
  typedef SCT::magnitudeType                                           Magnitude;
  typedef int                                                          LocalOrdinal;  
  typedef long                                                         GlobalOrdinal;  
  typedef Tpetra::DefaultPlatform::DefaultPlatformType                 Platform; 
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType       Node;     
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>               Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>  MV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>     OP;
  typedef Anasazi::MultiVecTraits<Scalar, MV>                          MVT;
  typedef Anasazi::OperatorTraits<Scalar, MV, OP>                      OPT;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>              Import;



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
class MyOp : public virtual OP {
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
  void apply(const MV& X, MV& Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=SCT::one(), Scalar beta=SCT::zero()) const
  {
    //
    // Let's make sure alpha is 1 and beta is 0...
    // This will throw an exception if that is not the case.
    //
    TEUCHOS_TEST_FOR_EXCEPTION(alpha != SCT::one() || beta != SCT::zero(),std::invalid_argument,
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
    RCP<MV> redistData = rcp(new MV(redistMap_, numVecs));

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
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  // 
  // Get the default communicator and node
  //                                      
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();          
  RCP<Node>                      node = platform.getNode();          
  const int myRank = comm->getRank();  

  //
  // Get parameters from command-line processor
  // 
  Scalar tol = 1e-5;
  GlobalOrdinal n = 100;
  int nev = 1;
  int blockSize = 1;
  bool verbose = true;
  std::string whenToShift = "Always";
  Teuchos::CommandLineProcessor cmdp(false,true); 
  cmdp.setOption("n",&n, "Number of rows of our operator.");
  cmdp.setOption("tolerance",&tol, "Relative residual used for solver.");
  cmdp.setOption("nev",&nev, "Number of desired eigenpairs.");
  cmdp.setOption("blocksize",&blockSize, "Number of vectors to add to the subspace at each iteration.");
  cmdp.setOption("verbose","quiet",&verbose, "Whether to print a lot of info or a little bit.");
  cmdp.setOption("whenToShift",&whenToShift, "When to perform Ritz shifts. Options: Never, After Trace Levels, Always.");
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
  MyPL.set( "Saddle Solver Type", "Projected Krylov"); // Use projected minres to solve the saddle point problem
  MyPL.set( "Block Size", blockSize );                 // Add blockSize vectors to the basis per iteration
  MyPL.set( "Convergence Tolerance", tol);             // How small do the residuals have to be?
  MyPL.set("Num Restart Blocks", numRestartBlocks);    // When we restart, this is how many blocks we keep
  MyPL.set("Num Blocks", numBlocks);                   // Maximum number of blocks in the subspace
  MyPL.set("When To Shift", whenToShift);              // What triggers a Ritz shift?
  
  //
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  // We are giving it random entries.
  //
  RCP<MV> ivec = rcp( new MV(K->getDomainMap(), numRestartBlocks*blockSize) );
  MVT::MvRandom( *ivec );

  //
  // Create the eigenproblem
  //
  RCP<Anasazi::BasicEigenproblem<Scalar,MV,OP> > MyProblem = 
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(K, ivec) );
  
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

  //
  // Initialize the TraceMin-Davidson solver
  //
  Anasazi::Experimental::TraceMinDavidsonSolMgr<Scalar, MV, OP> MySolverMgr(MyProblem, MyPL);
 
  //
  // Solve the problem to the specified tolerances
  //
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && myRank == 0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }
  else if (myRank == 0)
    cout << "Anasazi::EigensolverMgr::solve() returned converged." << endl;
  
  //
  // Get the eigenvalues and eigenvectors from the eigenproblem
  //
  Anasazi::Eigensolution<Scalar,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;
  
  //
  // Compute the residual, just as a precaution
  //
  if (numev > 0) {
    
    Teuchos::SerialDenseMatrix<int,Scalar> T(numev,numev);
    MV tempvec(K->getDomainMap(), MVT::GetNumberVecs( *evecs ));
    std::vector<Scalar> normR(sol.numVecs);
    MV Kvec( K->getRangeMap(), MVT::GetNumberVecs( *evecs ) );

    OPT::Apply( *K, *evecs, Kvec ); 
    MVT::MvTransMv( 1.0, Kvec, *evecs, T );
    MVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );
  
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

  return 0;
}
