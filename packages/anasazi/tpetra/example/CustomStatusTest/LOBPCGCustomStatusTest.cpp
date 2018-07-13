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
// This example shows how to define a custom operator and status test for Anasazi

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziStatusTestDecl.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"

// These tell the compiler which namespace contains RCP, cout, etc
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayView;
using std::cout;
using std::endl;

// These define convenient aliases
typedef double                                   Scalar;

typedef Tpetra::MultiVector<Scalar>              TMV;
typedef Tpetra::Vector<Scalar>                   Vector;
typedef Tpetra::Operator<Scalar>                 TOP;

namespace { // (anonymous)

// This class defines an operator A^2 for spectrum folding.  We don't
// explicitly form the matrix A^2; we just implement its action on a
// vector.  We can do this because Anasazi is templated on operators
// rather than matrices.
class FoldOp : public TOP {
public:
  typedef Tpetra::Map<> map_type;

  FoldOp (const RCP<const TOP> A) { A_ = A; };

  int SetUseTranspose (bool UseTranspose) { return -1; };

  void
  apply (const TMV& X, TMV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const;

  RCP<const map_type> getDomainMap () const { return A_->getDomainMap (); };
  RCP<const map_type> getRangeMap () const { return A_->getRangeMap (); };

private:
  RCP<const TOP> A_;
}; // end of class

void
FoldOp::apply (const TMV &X, TMV &Y, Teuchos::ETransp mode,
               Scalar alpha, Scalar beta) const
{
  TMV Y1 (X.getMap (), X.getNumVectors (), false);
  A_->apply (X, Y1, mode, alpha, beta);
  A_->apply (Y1, Y, mode, alpha, beta);
}

// Define a custom status test for spectrum folding.
//
// Since our operator is defined as A^2, the Anasazi eigensolvers will
// compute the residual as R = A^2 X - X \Sigma^2.  We want to monitor
// the residual R = A X - X \Sigma instead.
class StatusTestFolding : public Anasazi::StatusTest<Scalar,TMV,TOP> {
public:
  // Constructor
  StatusTestFolding (Scalar tol, int quorum = -1,
                     bool scaled = true,
                     bool throwExceptionOnNan = true,
                     const RCP<const TOP>& A = Teuchos::null);

  // Destructor
  virtual ~StatusTestFolding() {};

  // Check whether the test passed or failed
  Anasazi::TestStatus checkStatus (Anasazi::Eigensolver<Scalar,TMV,TOP>* solver);

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
  // bool throwExceptionOnNaN_; (unused)
  RCP<const TOP> A_;

  const Scalar ONE;
};

// Constructor for our custom status test
StatusTestFolding::
StatusTestFolding (Scalar tol, int quorum, bool scaled,
                   bool /* throwExceptionOnNaN (unused) */,
                   const RCP<const TOP>& A)
  : state_ (Anasazi::Undefined),
    tol_ (tol),
    quorum_ (quorum),
    scaled_ (scaled),
    /* throwExceptionOnNaN_ (throwExceptionOnNaN), (unused) */
    A_ (A),
    ONE (Teuchos::ScalarTraits<Scalar>::one ())
{}

// Check whether the test passed or failed
Anasazi::TestStatus
StatusTestFolding::checkStatus (Anasazi::Eigensolver<Scalar,TMV,TOP>* solver)
{
  typedef Anasazi::MultiVecTraits<Scalar, TMV> TMVT;

  RCP<const TMV> X = solver->getRitzVectors ();
  const int numev = X->getNumVectors ();
  std::vector<Scalar> res (numev);
  TMV AX (X->getMap (), numev, false);
  Teuchos::SerialDenseMatrix<int, Scalar> T (numev, numev);

  A_->apply (*X, AX);
  TMVT::MvTransMv (1.0, AX, *X, T);
  TMVT::MvTimesMatAddMv (-1.0, *X, T, 1.0, AX);
  TMVT::MvNorm (AX, res);

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
  typedef Anasazi::BasicEigenproblem<Scalar,TMV,TOP> Problem;
  typedef Anasazi::MultiVecTraits<Scalar, TMV> TMVT;
  typedef Anasazi::OperatorTraits<Scalar, TMV, TOP> TOPT;
  typedef Tpetra::CrsMatrix<> CrsMatrix;
  typedef Tpetra::MatrixMarket::Reader<CrsMatrix>  Reader;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank();

  // Read the command line arguments
  std::string fileA ("/u/slotnick_s2/aklinvex/matrices/anderson4.mtx");
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("fileA", &fileA, "Filename for the Matrix-Market stiffness matrix.");
  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Get the matrix
  Teuchos::RCP<Tpetra::Map<>::node_type> node; // for type deduction
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
  RCP<TMV> ivec = rcp (new TMV (A->getRowMap (), blockSize));
  TMVT::MvRandom (*ivec);

  // Create the eigenproblem
  RCP<Problem> MyProblem = rcp (new Problem (K, ivec));

  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian (true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV (nev);

  // Tell the problem that you are finished passing it information
  MyProblem->setProblem ();

  // Create the eigensolver and give it your problem and parameters.
  Anasazi::LOBPCGSolMgr<Scalar,TMV,TOP> solver (MyProblem, MyPL);

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
  Anasazi::Eigensolution<Scalar,TMV> sol = MyProblem->getSolution ();
  std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
  RCP<TMV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  // Compute the residual, just as a precaution
  if (numev > 0) {
    std::vector<Scalar> normR (sol.numVecs);

    TMV Avec (A->getRowMap (), TMVT::GetNumberVecs (*evecs));
    TOPT::Apply (*A, *evecs, Avec);
    Teuchos::SerialDenseMatrix<int,Scalar> T (numev, numev);
    TMVT::MvTransMv (1.0, Avec, *evecs, T);

    TMVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, Avec);
    TMVT::MvNorm (Avec, normR);

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
