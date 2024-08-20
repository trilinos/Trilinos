// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example computes the eigenvalues of smallest magnitude of a quadratic
//  eigenvalue problem
//     M x \lambda^2 + C x \lambda + K x = 0
//  where M and K are symmetric and positive definite and C is skew-symmetric.
//  This is known as a gyroscopic quadratic eigenvalue problem. 
//  This test reads a linearized version of this problem, using the linearization
//    A = [M  0]    B = [0 M]
//        [0 -K]        [M C]
//  and solves it using shift-invert BKS:
//   Op = inv(A) B
//  where inv(A) is provided by an Amesos direct sparse factorization.
//  This test illustrates a few issues:
//     1) A shortcoming in BKSSolMgr which assumed that the dominant eigenpairs
//        were converged, instead of querying the StatusTest to see which pairs
//        are in fact converged. This resulted in unconverged pairs being returned
//        to the caller.
//     2) Some issues sorting the Schur form, presumably due to the Ritz values for
//        the included problem not being well separated.
//  For simplicity (and due to its small size: 536 x 536), this problem is run only in 
//  serial. For a parallel build, only proc 0 does anything.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include <Amesos.h>
#include <EpetraExt_readEpetraLinearSystem.h>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include "AnasaziEpetraAdapter.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using Teuchos::RCP;
using Teuchos::rcp;

// ****************************************************************************
// Class:  AmesosGenOp
// Purpose:  Applies A^{-1}*B*X = Y or A^{-T}*B*X = Y where A is an Amesos solver
//           and B is an Epetra_Operator.  Class supports Epetra_Operator interface.
// Date: Jan 9, 2006
// Author:  Heidi K. Thornquist
// ****************************************************************************

class AmesosGenOp : public virtual Epetra_Operator
{
public:
  // Basic constructor
  AmesosGenOp( Epetra_LinearProblem& problem,
               const RCP<Amesos_BaseSolver>& solver,
               const RCP<Epetra_Operator>& massMtx,
               bool useTranspose = false );
  // Destructor
  ~AmesosGenOp() {};
  
  // Methods for supporting Epetra_Operator interface
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;
  const char* Label() const { return "Amesos direct solver for applying A^{-1}M"; };
  bool UseTranspose() const { return useTranspose_; };
  int SetUseTranspose( bool useTranspose );
  const Epetra_Comm& Comm() const { return solver_->Comm(); };
  const Epetra_Map& OperatorDomainMap() const { return massMtx_->OperatorDomainMap(); };
  const Epetra_Map& OperatorRangeMap() const { return massMtx_->OperatorRangeMap(); };
    
  // Epetra_Operator interface methods that are not supported.
  // Note:  ApplyInverse not defined because M not guaranteed to have an inverse.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const { return -1; };
  bool HasNormInf() const { return false; };
  double NormInf() const { return -1.0; };
   
private:  
  // Default constructor
  AmesosGenOp () {};
 
  // Copy constructor
  AmesosGenOp ( const AmesosGenOp& genOp ) {};

  // Epetra_LinearProblem contained in the Amesos_BaseSolver
  bool useTranspose_;
  RCP<Amesos_BaseSolver> solver_;
  RCP<Epetra_Operator> massMtx_;
  Epetra_LinearProblem* problem_;
  
};

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&cout);

  bool run_me = (mpisess.getRank() == 0);
  bool testFailed = false;

  if (run_me) {
    try {
      // useful typedefs
      typedef double                              ST;
      typedef Teuchos::ScalarTraits<ST>          STT;
      typedef STT::magnitudeType                  MT;
      typedef Teuchos::ScalarTraits<MT>          MTT;
      typedef Epetra_MultiVector                  MV;
      typedef Epetra_Operator                     OP;
      typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
      typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;

      // parse here, so everyone knows about failure
      bool verbose    = false;
      bool debug      = false;
      std::string k_filename = "linearized_qevp_A.hb";
      std::string m_filename = "linearized_qevp_B.hb";
      int blockSize   = 1;
      int numBlocks   = 30;
      int nev         = 20;
      int maxRestarts = 0;
      int extraBlocks = 0;
      int stepSize    = 1;
      int numPrint    = 536;
      MT  tol = 1e-8;
      Teuchos::CommandLineProcessor cmdp(true,true);
      cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
      cmdp.setOption("debug","normal",&debug,"Print debugging information.");
      cmdp.setOption("A-filename",&k_filename,"Filename and path of the stiffness matrix.");
      cmdp.setOption("B-filename",&m_filename,"Filename and path of the mass matrix.");
      cmdp.setOption("extra-blocks",&extraBlocks,"Extra blocks to keep on a restart.");
      cmdp.setOption("block-size",&blockSize,"Block size.");
      cmdp.setOption("num-blocks",&numBlocks,"Number of blocks in Krylov basis.");
      cmdp.setOption("step-size",&stepSize,"Step size.");
      cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
      cmdp.setOption("num-restarts",&maxRestarts,"Maximum number of restarts.");
      cmdp.setOption("tol",&tol,"Convergence tolerance.");
      cmdp.setOption("num-print",&numPrint,"Number of Ritz values to print.");
      // parse() will throw an exception on error
      cmdp.parse(argc,argv);

      // Get the stiffness and mass matrices
      Epetra_SerialComm Comm;
      RCP<Epetra_Map> Map;
      RCP<Epetra_CrsMatrix> A, B;
      EpetraExt::readEpetraLinearSystem( k_filename, Comm, &A, &Map );
      EpetraExt::readEpetraLinearSystem( m_filename, Comm, &B, &Map );

      //
      // *******************************************************
      // Set up Amesos direct solver for inner iteration
      // *******************************************************
      //
      // Create Epetra linear problem class to solve "Kx = b"
      Epetra_LinearProblem AmesosProblem;
      AmesosProblem.SetOperator(A.get());
      // Create Amesos factory and solver for solving "Kx = b" using a direct factorization
      Amesos amesosFactory;
      RCP<Amesos_BaseSolver> AmesosSolver = rcp( amesosFactory.Create( "Klu", AmesosProblem ) );
      // The AmesosGenOp class assumes that the symbolic/numeric factorizations have already
      // been performed on the linear problem.
      AmesosSolver->SymbolicFactorization();
      AmesosSolver->NumericFactorization();

      //
      // ************************************
      // Start the block Arnoldi iteration
      // ************************************
      //
      //  Variables used for the Block Arnoldi Method
      //
      int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
      if (verbose) {
        verbosity += Anasazi::IterationDetails;
      }
      if (debug) {
        verbosity += Anasazi::Debug;
      }
      //
      // Create parameter list to pass into solver
      //
      Teuchos::ParameterList MyPL;
      MyPL.set( "Verbosity", verbosity );
      MyPL.set( "Which", "LM" );
      MyPL.set( "Block Size", blockSize );
      MyPL.set( "Num Blocks", numBlocks );
      MyPL.set( "Maximum Restarts", maxRestarts );
      MyPL.set( "Convergence Tolerance", tol );
      MyPL.set( "Step Size", stepSize );
      MyPL.set( "Extra NEV Blocks", extraBlocks );
      MyPL.set( "Print Number of Ritz Values", numPrint );

      // Create an Epetra_MultiVector for an initial vector to start the solver.
      // Note:  This needs to have the same number of columns as the blocksize.
      RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(A->Map(), blockSize) );
      // MVT::MvRandom( *ivec ); // FINISH: put this back in
      MVT::MvInit(*ivec,1.0);

      // Create the Epetra_Operator for the spectral transformation using the Amesos direct solver.
      RCP<AmesosGenOp> Aop = rcp( new AmesosGenOp(AmesosProblem, AmesosSolver, B) );

      // standard inner product; B is not symmetric positive definite, and Op has no symmetry.
      RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > MyProblem = 
        rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(Aop, ivec) );
      MyProblem->setHermitian(false);
      MyProblem->setNEV( nev );
      // Inform the eigenproblem that you are finished passing it information
      TEUCHOS_TEST_FOR_EXCEPTION( MyProblem->setProblem() != true, 
                          std::runtime_error, "Anasazi::BasicEigenproblem::setProblem() returned with error.");

      // Initialize the Block Arnoldi solver
      Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> MySolverMgr(MyProblem, MyPL);

      // Solve the problem to the specified tolerances or length
      Anasazi::ReturnType returnCode = MySolverMgr.solve();

      // Get the eigenvalues and eigenvectors from the eigenproblem
      Anasazi::Eigensolution<ST,MV> sol = MyProblem->getSolution();
      vector<Anasazi::Value<ST> > evals = sol.Evals;
      RCP<MV> evecs = sol.Evecs;
      vector<int> index = sol.index;
      int numev = sol.numVecs;

      if (returnCode != Anasazi::Converged) {
        cout << "Anasazi::EigensolverMgr::solve() returned unconverged, computing " << numev << " of " << nev << endl;
      }
      else {
        cout << "Anasazi::EigensolverMgr::solve() returned converged, computing " << numev << " of " << nev << endl;
      }

      if (numev > 0) {
        // Compute residuals.
        Teuchos::LAPACK<int,MT> lapack;
        vector<MT> normA(numev);

        // Back-transform the eigenvalues
        for (int i=0; i<numev; ++i) {
          MT mag2 = lapack.LAPY2(evals[i].realpart,evals[i].imagpart);
          mag2 = mag2*mag2;
          evals[i].realpart /=   mag2;
          evals[i].imagpart /= (-mag2);
        }

        // The problem is non-Hermitian.
        vector<int> curind(1);
        vector<MT> resnorm(1), tempnrm(1);
        Teuchos::RCP<const MV> Av_r, Av_i, Bv_r, Bv_i;
        Epetra_MultiVector Aevec(*Map,numev), Bevec(*Map,numev), res(*Map,1);

        // Compute A*evecs, B*evecs
        OPT::Apply( *A, *evecs, Aevec );
        OPT::Apply( *B, *evecs, Bevec );

        for (int i=0; i<numev;) {
          if (index[i]==0) {
            // Get views of the real part of A*evec,B*evec
            curind[0] = i;
            Av_r = MVT::CloneView( Aevec, curind );
            Bv_r = MVT::CloneView( Bevec, curind );

            // Compute set res = lambda*B*evec - A*evec
            // eigenvalue and eigenvector are both real
            MVT::MvAddMv(evals[i].realpart,*Bv_r,-1.0,*Av_r,res);

            // Compute the norm of the residual and increment counter
            MVT::MvNorm( res, resnorm );
            // Scale the residual
            normA[i] = resnorm[0];
            MT mag = MTT::magnitude(evals[i].realpart);
            if ( mag > MTT::one() ) {
              normA[i] /= mag;
            }
            // done with this eigenvector
            i++;
          } else {
            // Get a copy of A*evec, B*evec
            curind[0] = i;
            Av_r = MVT::CloneCopy( Aevec, curind );
            Bv_r = MVT::CloneView( Bevec, curind );
            // Get the imaginary part of A*evec,B*evec
            curind[0] = i+1;
            Av_i = MVT::CloneCopy( Aevec, curind );
            Bv_i = MVT::CloneView( Bevec, curind );
            // grab temp copies of the eigenvalue
            MT l_r = evals[i].realpart,
               l_i = evals[i].imagpart;

            // Compute real part of B*evec*lambda - A*evec
            // B is real
            //   evec =   evec_r +   evec_i*i
            // lambda = lambda_r + lambda_i*i
            // 
            // evec*lambda = evec_r*lambda_r - evec_i*lambda_i + evec_r*lambda_i*i + evec_i*lambda_r*i
            //
            // res = B*evec*lambda - A*evec
            //     = B*(evec_r*lambda_r - evec_i*lambda_i + evec_r*lambda_i*i + evec_i*lambda_r*i) - A*evec_r - A*evec_i*i
            //     = (B*evec_r*lambda_r - B*evec_i*lambda_i - A*evec_r) + (B*evec_r*lambda_i + B*evec_i*lambda_r - A*evec_i)*i
            // real(res) = B*evec_r*lambda_r - B*evec_i*lambda_i - A*evec_r
            // imag(res) = B*evec_r*lambda_i + B*evec_i*lambda_r - A*evec_i

            // compute real part of residual and its norm
            MVT::MvAddMv(l_r,*Bv_r, -l_i,*Bv_i, res);
            MVT::MvAddMv(MTT::one(),res, -MTT::one(),*Av_r, res);
            MVT::MvNorm(res,tempnrm);
            // compute imag part of residual and its norm
            MVT::MvAddMv(l_i,*Bv_r, l_r,*Bv_i, res);
            MVT::MvAddMv(MTT::one(),res, -MTT::one(),*Av_i, res);
            MVT::MvNorm(res,resnorm);

            // Compute the norms and scale by magnitude of eigenvalue
            normA[i] = lapack.LAPY2( tempnrm[0], resnorm[0] );
            MT mag = lapack.LAPY2(evals[i].realpart, evals[i].imagpart);
            if (mag > MTT::one()) {
              normA[i] /= mag;
            }
            normA[i+1] = normA[i];
            // done with this conjugate pair
            i=i+2;
          }
        }

        // Output computed eigenvalues and their direct residuals
        cout.setf(std::ios_base::right, std::ios_base::adjustfield);        
        cout<<endl<< "Actual Residuals"<<endl;
        cout<< std::setw(16) << "Real Part"
          << std::setw(16) << "Imag Part"
          << std::setw(20) << "Direct Residual"<< endl;
        cout<<"-----------------------------------------------------------"<<endl;
        for (int j=0; j<numev; j++) {
          cout<< std::setw(16) << evals[j].realpart 
            << std::setw(16) << evals[j].imagpart 
            << std::setw(20) << normA[j] << endl;
          if (normA[j] > tol) {
            testFailed = true;
          }
        }  
        cout<<"-----------------------------------------------------------"<<endl;
      }
    }
    catch (std::exception &e) {
      cout << "Caught exception: " << endl << e.what() << endl;
      testFailed = true;
    }
  }
  else { // run_me == false 
    cout << "\nNot running on processor " << mpisess.getRank() << ". Serial problem only." << endl;
  }
  if (mpisess.getRank() == 0) {
    if (testFailed) {
      cout << "End Result: TEST FAILED" << endl;
    }
    else {
      cout << "End Result: TEST PASSED" << endl;
    }
  }
  return 0;
}

// ****************************************************************************
// Class:  AmesosGenOp
// Purpose:  Applies A^{-1}*B*X = Y or A^{-T}*B*X = Y where A is an Amesos solver
//           and B is an Epetra_Operator.  Class supports Epetra_Operator interface.
// Date: Jan 9, 2006
// Author:  Heidi K. Thornquist
// ****************************************************************************


AmesosGenOp::AmesosGenOp( Epetra_LinearProblem& problem,
                          const RCP<Amesos_BaseSolver>& solver,
                          const RCP<Epetra_Operator>& massMtx,
                          bool useTranspose )
  : useTranspose_(useTranspose),
    solver_(solver),
    massMtx_(massMtx)
{
  problem_ = const_cast<Epetra_LinearProblem*>( solver->GetProblem() );
  
  if ( solver_->UseTranspose() )
    solver_->SetUseTranspose(!useTranspose); 
  else
    solver_->SetUseTranspose(useTranspose);
  
  if ( massMtx_->UseTranspose() )
    massMtx_->SetUseTranspose(!useTranspose);
  else
    massMtx_->SetUseTranspose(useTranspose);    
}

// Methods for supporting Epetra_Operator interface
int AmesosGenOp::SetUseTranspose(bool useTranspose) 
{ 
  useTranspose_ = useTranspose; 
  if ( solver_->UseTranspose() )
    solver_->SetUseTranspose(!useTranspose); 
  else
    solver_->SetUseTranspose(useTranspose);
  
  if ( massMtx_->UseTranspose() )
    massMtx_->SetUseTranspose(!useTranspose);
  else
    massMtx_->SetUseTranspose(useTranspose);

  return 0;
}

int AmesosGenOp::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const 
{
  if (!useTranspose_) {
    
    // Storage for M*X
    Epetra_MultiVector MX(X.Map(),X.NumVectors());
    
    // Apply M*X
    massMtx_->Apply(X, MX);
    Y.PutScalar(0.0);
    
    // Set the LHS and RHS
    problem_->SetRHS(&MX);
    problem_->SetLHS(&Y);

    // Solve the linear system A*Y = MX
    solver_->Solve();
  }
  else {
    // Storage for A^{-T}*X
    Epetra_MultiVector ATX(X.Map(),X.NumVectors());
    Epetra_MultiVector tmpX = const_cast<Epetra_MultiVector&>(X);
    
    // Set the LHS and RHS
    problem_->SetRHS(&tmpX);
    problem_->SetLHS(&ATX);
    
    // Solve the linear system A^T*Y = X 
    solver_->Solve();
    
    // Apply M*ATX
    massMtx_->Apply(ATX, Y);
  }
  
  return 0;
}
