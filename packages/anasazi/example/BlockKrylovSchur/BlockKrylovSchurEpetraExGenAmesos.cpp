//  This example computes the eigenvalues of smallest magnitude of the 
//  discretized 2D Laplacian operator using the block Krylov-Schur method.  
//  This problem shows the construction of an inner-outer iteration using 
//  Amesos as the linear solver within Anasazi.  This operator is 
//  discretized using linear finite elements and constructed as an Epetra 
//  matrix, then passed into the Amesos solver to perform a factorization.
//  The factorization is then used to apply the shift-invert
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

// Include selected communicator class and map required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

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
               const Teuchos::RCP<Amesos_BaseSolver>& solver,
               const Teuchos::RCP<Epetra_Operator>& massMtx,
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
  Teuchos::RCP<Amesos_BaseSolver> solver_;
  Teuchos::RCP<Epetra_Operator> massMtx_;
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

  // Create Epetra linear problem class to solve "Kx = b"
  Epetra_LinearProblem AmesosProblem;
  AmesosProblem.SetOperator(K.get());
  
  // Create Amesos factory and solver for solving "Kx = b" using a direct factorization
  Amesos amesosFactory;
  Teuchos::RCP<Amesos_BaseSolver> AmesosSolver = 
    Teuchos::rcp( amesosFactory.Create( "Klu", AmesosProblem ) );

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
  
  // Create the Epetra_Operator for the spectral transformation using the Amesos direct solver.
  Teuchos::RCP<AmesosGenOp> Aop 
    = Teuchos::rcp( new AmesosGenOp(AmesosProblem, AmesosSolver, M) );
  
  Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Aop, M, ivec) );
  
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian(true);
  
  // Set the number of eigenvalues requested 
  MyProblem->setNEV( nev );
  
  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
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
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;
  
  if (numev > 0) {
    
    Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
    Epetra_MultiVector tempvec(K->Map(), MVT::GetNumberVecs( *evecs ));
    OPT::Apply( *K, *evecs, tempvec );
    MVT::MvTransMv( 1.0, tempvec, *evecs, dmatr );
    
    if (MyPID==0) {
      double compeval = 0.0;
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
      cout<<"------------------------------------------------------"<<endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Rayleigh Error"<<endl;
      cout<<"------------------------------------------------------"<<endl;
      for (i=0; i<numev; i++) {
        compeval = dmatr(i,i);
        cout<<std::setw(16)<<compeval
          <<std::setw(16)<<Teuchos::ScalarTraits<double>::magnitude(compeval-1.0/evals[i].realpart)
          <<endl;
      }
      cout<<"------------------------------------------------------"<<endl;
    }
    
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

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
                          const Teuchos::RCP<Amesos_BaseSolver>& solver,
                          const Teuchos::RCP<Epetra_Operator>& massMtx,
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
