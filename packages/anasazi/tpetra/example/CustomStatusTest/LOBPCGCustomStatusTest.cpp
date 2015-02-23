// This example shows how to define a custom operator and status test for Anasazi

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziStatusTestDecl.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "MatrixMarket_Tpetra.hpp"

// These tell the compiler which namespace contains RCP, cout, etc
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayView;
using std::cout;
using std::endl;

// These define convenient aliases
typedef double                                   Scalar;

typedef Tpetra::MultiVector<Scalar>              MV;
typedef Tpetra::Vector<Scalar>                   Vector;
typedef Tpetra::Operator<Scalar>                 OP;

namespace { // (anonymous)

// This class defines an operator A^2 for spectrum folding.  We don't
// explicitly form the matrix A^2; we just implement its action on a
// vector.  We can do this because Anasazi is templated on operators
// rather than matrices.
class FoldOp : public OP {
public:
  typedef Tpetra::Map<> map_type;

  FoldOp (const RCP<const OP> A) { A_ = A; };

  int SetUseTranspose (bool UseTranspose) { return -1; };

  void
  apply (const MV& X, MV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const;

  RCP<const map_type> getDomainMap () const { return A_->getDomainMap (); };
  RCP<const map_type> getRangeMap () const { return A_->getRangeMap (); };

private:
  RCP<const OP> A_;
}; // end of class

void
FoldOp::apply (const MV &X, MV &Y, Teuchos::ETransp mode,
               Scalar alpha, Scalar beta) const
{
  MV Y1 (X.getMap (), X.getNumVectors (), false);
  A_->apply (X, Y1, mode, alpha, beta);
  A_->apply (Y1, Y, mode, alpha, beta);
}

// Define a custom status test for spectrum folding.
//
// Since our operator is defined as A^2, the Anasazi eigensolvers will
// compute the residual as R = A^2 X - X \Sigma^2.  We want to monitor
// the residual R = A X - X \Sigma instead.
class StatusTestFolding : public Anasazi::StatusTest<Scalar,MV,OP> {
public:
  // Constructor
  StatusTestFolding (Scalar tol, int quorum = -1,
                     bool scaled = true,
                     bool throwExceptionOnNan = true,
                     const RCP<const OP>& A = Teuchos::null);

  // Destructor
  virtual ~StatusTestFolding() {};

  // Check whether the test passed or failed
  Anasazi::TestStatus checkStatus (Anasazi::Eigensolver<Scalar,MV,OP>* solver);

  // Return the result of the most recent checkStatus call
  Anasazi::TestStatus getStatus () const { return state_; }

  // Get the indices for the vectors that passed the test
  std::vector<int> whichVecs () const { return ind_; }

  // Get the number of vectors that passed the test
  int howMany () const { return ind_.size (); }

  // Informs the status test that it should reset its internal configuration to the uninitialized state
  void reset () {
    ind_.resize (0);
    state_ = Anasazi::Undefined;
  }

  // Clears the results of the last status test
  void clearStatus () { reset (); };

  // Output formatted description of stopping test to output stream
  std::ostream& print (std::ostream &os, int indent=0) const;

private:
  Anasazi::TestStatus state_;
  Scalar tol_;
  std::vector<int> ind_;
  int quorum_;
  bool scaled_;
  bool throwExceptionOnNaN_;
  RCP<const OP> A_;

  const Scalar ONE;
};

// Constructor for our custom status test
StatusTestFolding::
StatusTestFolding (Scalar tol, int quorum, bool scaled,
                   bool throwExceptionOnNaN, const RCP<const OP>& A)
  : state_ (Anasazi::Undefined),
    tol_ (tol),
    quorum_ (quorum),
    scaled_ (scaled),
    throwExceptionOnNaN_ (throwExceptionOnNaN),
    A_ (A),
    ONE (Teuchos::ScalarTraits<Scalar>::one ())
{}

// Check whether the test passed or failed
Anasazi::TestStatus
StatusTestFolding::checkStatus (Anasazi::Eigensolver<Scalar,MV,OP>* solver)
{
  typedef Anasazi::MultiVecTraits<Scalar, MV> MVT;

  RCP<const MV> X = solver->getRitzVectors ();
  const int numev = X->getNumVectors ();
  std::vector<Scalar> res (numev);
  MV AX (X->getMap (), numev, false);
  Teuchos::SerialDenseMatrix<int, Scalar> T (numev, numev);

  A_->apply (*X, AX);
  MVT::MvTransMv (1.0, AX, *X, T);
  MVT::MvTimesMatAddMv (-1.0, *X, T, 1.0, AX);
  MVT::MvNorm (AX, res);

  // if appropriate, scale the norms by the magnitude of the eigenvalue estimate
  if (scaled_) {
    for (int i = 0; i < numev; ++i) {
      res[i] /= std::abs (T(i,i));
    }
  }

  // test the norms
  ind_.resize (0);
  for (int i = 0; i < numev; ++i) {
    if (res[i] < tol_) {
      ind_.push_back (i);
    }
  }
  const int have = ind_.size ();
  const int need = (quorum_ == -1) ? numev : quorum_;
  state_ = (have >= need) ? Anasazi::Passed : Anasazi::Failed;
  return state_;
}


// Output formatted description of stopping test to output stream
std::ostream&
StatusTestFolding::print (std::ostream& os, int indent) const
{
  std::string ind (indent, ' ');
  os << ind << "- StatusTestFolding: ";
  switch (state_) {
  case Anasazi::Passed:
    os << "Passed\n";
    break;
  case Anasazi::Failed:
    os << "Failed\n";
    break;
  case Anasazi::Undefined:
    os << "Undefined\n";
    break;
  }
  os << ind << "  (Tolerance, WhichNorm,Scaled,Quorum): "
     << "(" << tol_
     << ",RES_2NORM"
     << "," << (scaled_ ? "true" : "false")
     << "," << quorum_
     << ")\n";

  if (state_ != Anasazi::Undefined) {
    os << ind << "  Which vectors: ";
    if (ind_.size () > 0) {
      for (size_t i = 0; i < ind_.size (); ++i) {
        os << ind_[i] << " ";
      }
      os << std::endl;
    }
    else {
      os << "[empty]\n";
    }
  }
  return os;
}
} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  typedef Anasazi::BasicEigenproblem<Scalar,MV,OP> Problem;
  typedef Anasazi::MultiVecTraits<Scalar, MV> MVT;
  typedef Anasazi::OperatorTraits<Scalar, MV, OP> OPT;
  typedef Tpetra::Map<>::node_type Node;
  typedef Tpetra::CrsMatrix<> CrsMatrix;
  typedef Tpetra::MatrixMarket::Reader<CrsMatrix>  Reader;

  //
  // Initialize the MPI session
  //
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);

  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<Node> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank();

  // Read the command line arguments
  std::string fileA ("/u/slotnick_s2/aklinvex/matrices/anderson4.mtx");
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("fileA", &fileA, "Filename for the Matrix-Market stiffness matrix.");
  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Get the matrix
  RCP<const CrsMatrix> A = Reader::readSparseFile (fileA, comm, node);

  // Create the folded operator, K = A^2
  RCP<FoldOp> K = rcp (new FoldOp (A));

  // set parameters for solver
  int blockSize = 4;
  double tol = 1e-5;
  bool scaled = true;
  int nev = 4;

  Teuchos::ParameterList MyPL;
  MyPL.set ("Which", "SR");
  MyPL.set ("Maximum Restarts", 10000);
  MyPL.set ("Maximum Iterations", 1000000);
  MyPL.set ("Block Size", blockSize);
  MyPL.set ("Convergence Tolerance", tol );   // How small do the residuals have to be
  MyPL.set ("Relative Convergence Tolerance", scaled);
  MyPL.set ("Relative Locking Tolerance", scaled);
  MyPL.set ("Verbosity", Anasazi::TimingDetails);

  // Create a MultiVector for an initial subspace to start the solver.
  RCP<MV> ivec = rcp (new MV (A->getRowMap (), blockSize));
  MVT::MvRandom (*ivec);

  // Create the eigenproblem
  RCP<Problem> MyProblem = rcp (new Problem (K, ivec));

  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian (true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV (nev);

  // Tell the problem that you are finished passing it information
  MyProblem->setProblem ();

  // Create the eigensolver and give it your problem and parameters.
  Anasazi::LOBPCGSolMgr<Scalar,MV,OP> solver (MyProblem, MyPL);

  // Give the eigensolver a status test consistent with the spectral transformation.
  RCP<StatusTestFolding> convTest =
    rcp (new StatusTestFolding (tol, nev, scaled, true, A));
  RCP<StatusTestFolding> lockTest =
    rcp (new StatusTestFolding (tol/10., 1, scaled, true, A));
  solver.setGlobalStatusTest (convTest);
  solver.setLockingStatusTest (lockTest);

  // Tell the solver to solve the eigenproblem.
  Anasazi::ReturnType returnCode = solver.solve ();
  if (returnCode != Anasazi::Converged && myRank == 0) {
    cout << "The solve did NOT converge." << endl;
  } else if (myRank == 0) {
    cout << "The solve converged." << endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<Scalar,MV> sol = MyProblem->getSolution ();
  std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  // Compute the residual, just as a precaution
  if (numev > 0) {
    std::vector<Scalar> normR (sol.numVecs);

    MV Avec (A->getRowMap (), MVT::GetNumberVecs (*evecs));
    OPT::Apply (*A, *evecs, Avec);
    Teuchos::SerialDenseMatrix<int,Scalar> T (numev, numev);
    MVT::MvTransMv (1.0, Avec, *evecs, T);

    MVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, Avec);
    MVT::MvNorm (Avec, normR);

    if (myRank == 0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues: "<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Error"<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      for (int i=0; i<numev; i++) {
        cout<<std::setw(16)<<T(i,i)
          <<std::setw(16)<<normR[i]/std::abs(T(i,i))
          <<std::endl;
      }
      cout<<"------------------------------------------------------"<<std::endl;
    }
  }

  return 0;
}
