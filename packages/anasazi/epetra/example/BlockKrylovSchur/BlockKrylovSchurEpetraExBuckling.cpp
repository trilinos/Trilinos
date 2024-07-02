// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example computes the eigenvalues of smallest magnitude of the 
//  discretized 2D Laplacian operator using the block Krylov-Schur method.  
//  This problem shows the construction of an inner-outer iteration using 
//  Amesos as the linear solver within Anasazi.  This operator is 
//  discretized using linear finite elements and constructed as an Epetra 
//  matrix, then passed into the Amesos solver to perform a factorization.
//  The factorization is then used to apply the buckling transform
//  operation to be used within the Krylov decomposition.  The specifics 
//  of the block Krylov-Schur method can be set by the user.

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for block Krylov-Schur solver
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide Anasazi with Epetra adapters
#include "AnasaziEpetraAdapter.hpp"

// Include header for Epetra compressed-row storage matrix and linear problem
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// Include header for Amesos solver and solver interface for Epetra_Operator
#include "Amesos.h"

// Include header for Teuchos serial dense matrix
#include "Teuchos_SerialDenseMatrix.hpp"

// Include header for the problem definition
#include "ModeLaplace2DQ2.h"

// Include EpetraExt MatrixMatrix helpers.
#include "EpetraExt_MatrixMatrix.h"

// Include selected communicator class and map required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

// ****************************************************************************
// Class:  AmesosBucklingOp
// Purpose:  Applies (K - sigma*M)^{-1}*K*X = Y where (K - sigma*M) is an Amesos solver
//           and K is an Epetra_Operator.  Class supports Epetra_Operator interface.
// Date: July 9, 2009
// Author:  Heidi K. Thornquist
// ****************************************************************************

class AmesosBucklingOp : public virtual Epetra_Operator
{
public:
  // Basic constructor
  AmesosBucklingOp( Epetra_LinearProblem& problem,
		    const Teuchos::RCP<Amesos_BaseSolver>& solver,
		    const Teuchos::RCP<Epetra_Operator>& stiffMtx
		    );

  // Destructor
  ~AmesosBucklingOp() {};
  
  // Methods for supporting Epetra_Operator interface
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;
  const char* Label() const { return "Amesos direct solver for applying (K - sigma*M)^{-1}K"; };
  bool UseTranspose() const { return false; };
  int SetUseTranspose( bool useTranspose ) { return -1; };
  const Epetra_Comm& Comm() const { return solver_->Comm(); };
  const Epetra_Map& OperatorDomainMap() const { return stiffMtx_->OperatorDomainMap(); };
  const Epetra_Map& OperatorRangeMap() const { return stiffMtx_->OperatorRangeMap(); };
    
  // Epetra_Operator interface methods that are not supported.
  // Note:  ApplyInverse not defined because K is not guaranteed to have an inverse.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const { return -1; };
  bool HasNormInf() const { return false; };
  double NormInf() const { return -1.0; };
   
private:  
  // Default constructor
  AmesosBucklingOp () {};
 
  // Copy constructor
  AmesosBucklingOp ( const AmesosBucklingOp& genOp ) {};

  // Epetra_LinearProblem contained in the Amesos_BaseSolver
  Teuchos::RCP<Amesos_BaseSolver> solver_;
  Teuchos::RCP<Epetra_Operator> stiffMtx_;
  Epetra_LinearProblem* problem_;
  
};

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int main(int argc, char *argv[]) {
  int i;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  // Number of dimension of the domain
  int space_dim = 2;
  
  // Size of each of the dimensions of the domain
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  brick_dim[1] = 1.0;
  
  // Number of elements in each of the dimensions of the domain
  std::vector<int> elements( space_dim );
  elements[0] = 10;
  elements[1] = 10;
  
  // Create problem
  Teuchos::RCP<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );
  
  // Get the stiffness and mass matrices
  Teuchos::RCP<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RCP<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
	
  //
  // *******************************************************
  // Set up Amesos direct solver for inner iteration
  // *******************************************************
  //

  // Create the shifted system K - sigma * M.
  // For the buckling transformation, this shift must be nonzero.

  double sigma = 1.0;
  Epetra_CrsMatrix Kcopy( *K );

  int addErr = EpetraExt::MatrixMatrix::Add( *M, false, -sigma, Kcopy, 1.0 );
  if (addErr != 0) {
    if (MyPID == 0) {
      std::cout << "EpetraExt::MatrixMatrix::Add returned with error: " << addErr << std::endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Create Epetra linear problem class to solve "x = b"
  Epetra_LinearProblem AmesosProblem;
  AmesosProblem.SetOperator(&Kcopy);
  
  // Create Amesos factory and solver for solving "(K - sigma*M)x = b" using a direct factorization
  Amesos amesosFactory;
  Teuchos::RCP<Amesos_BaseSolver> AmesosSolver = 
    Teuchos::rcp( amesosFactory.Create( "Klu", AmesosProblem ) );

  // The AmesosBucklingOp class assumes that the symbolic/numeric factorizations have already
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
  int nev = 10;
  int blockSize = 3;  
  int numBlocks = 3*nev/blockSize;
  int maxRestarts = 5;
  //int step = 5;
  double tol = 1.0e-8;
  std::string which = "LM";
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  //MyPL.set( "Step Size", step );
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
  
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->Map(), blockSize) );
  MVT::MvRandom( *ivec );
  
  // Create the Epetra_Operator for the buckling transformation using the Amesos direct solver.
  Teuchos::RCP<AmesosBucklingOp> BucklingOp 
    = Teuchos::rcp( new AmesosBucklingOp(AmesosProblem, AmesosSolver, K) );
  
  Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(BucklingOp, K, ivec) );
  
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian(true);
  
  // Set the number of eigenvalues requested 
  MyProblem->setNEV( nev );
  
  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (MyPID == 0) {
      std::cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << std::endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);
  
  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0) {
    std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;
  
  if (numev > 0) {

    // Undo buckling transformation; computed eigenvalues are real
    std::vector<double> compEvals(numev);
    for (i=0; i<numev; ++i) {
      compEvals[i] = sigma*evals[i].realpart/(evals[i].realpart-1.0);
    }
    
    // Remember, eigenvectors are constructed K-orthogonal to preserve symmetry,
    // so numerator of the Rayleigh quotient is 1.0.
    Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
    Epetra_MultiVector tempvec(M->Map(), MVT::GetNumberVecs( *evecs ));
    OPT::Apply( *M, *evecs, tempvec );
    MVT::MvTransMv( 1.0, tempvec, *evecs, dmatr );
    
    if (MyPID==0) {
      double rq_eval = 0.0;
      std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      std::cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<std::endl;
      std::cout<<"------------------------------------------------------"<<std::endl;
      std::cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Rayleigh Error"<<std::endl;
      std::cout<<"------------------------------------------------------"<<std::endl;
      for (i=0; i<numev; i++) {
        rq_eval = 1.0 / dmatr(i,i);
        std::cout<<std::setw(16)<<rq_eval
          <<std::setw(16)<<Teuchos::ScalarTraits<double>::magnitude(rq_eval-compEvals[i])
          <<std::endl;
      }
      std::cout<<"------------------------------------------------------"<<std::endl;
    }
    
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

// ****************************************************************************
// Class:  AmesosBucklingOp
// Purpose:  Applies (A - sigma*M)^{-1}*A*X = Y where A - sigma*M is an Amesos solver
//           and A is an Epetra_Operator.  Class supports Epetra_Operator interface.
// Date: July 9, 2009
// Author:  Heidi K. Thornquist
// ****************************************************************************


AmesosBucklingOp::AmesosBucklingOp( Epetra_LinearProblem& problem,
				    const Teuchos::RCP<Amesos_BaseSolver>& solver,
				    const Teuchos::RCP<Epetra_Operator>& stiffMtx
				    )
  : solver_(solver),
    stiffMtx_(stiffMtx)
{
  problem_ = const_cast<Epetra_LinearProblem*>( solver->GetProblem() );
}

int AmesosBucklingOp::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const 
{
    
  // Storage for A*X
  Epetra_MultiVector AX(X.Map(),X.NumVectors());
    
  // Apply A*X
  stiffMtx_->Apply(X, AX);
  Y.PutScalar(0.0);
    
  // Set the LHS and RHS
  problem_->SetRHS(&AX);
  problem_->SetLHS(&Y);

  // Solve the linear system (A-sigma*M)*Y = AX
  solver_->Solve();
  
  return 0;
}
