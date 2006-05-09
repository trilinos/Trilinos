// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
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
#include "AnasaziBlockKrylovSchur.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide basic sorting utility required by block Krylov-Schur method
#include "AnasaziBasicSort.hpp"

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
	       const Teuchos::RefCountPtr<Amesos_BaseSolver>& solver,
	       const Teuchos::RefCountPtr<Epetra_Operator>& massMtx,
	       bool useTranspose = false );
  // Destructor
  ~AmesosGenOp() {};
  
  // Methods for supporting Epetra_Operator interface
  int SetUseTranpose(bool useTranspose);
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
  Teuchos::RefCountPtr<Amesos_BaseSolver> solver_;
  Teuchos::RefCountPtr<Epetra_Operator> massMtx_;
  Epetra_LinearProblem* problem_;
  
};

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int main(int argc, char *argv[]) {
  int i, info;
  Anasazi::ReturnType returnCode = Anasazi::Ok;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
    Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
  MyOM->SetVerbosity( Anasazi::FinalSummary );  

  std::string which;
  if (argc > 1) {
    which = argv[1];
  }
  else {
    which = "LM";
  }
  if ( which != "SM" && which != "LM" && which != "SR" && which != "LR" ) {
    if (MyOM->doPrint()) {
      std::cout << "Usage: " << argv[0] << " [sort string]" << endl
        << "where:" << endl
        << "sort string       - SM | LM | SR | LR" << endl << endl;
    }
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

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
  Teuchos::RefCountPtr<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );
  
  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  //
  // *******************************************************
  // Set up Belos GMRES operator for inner iteration
  // *******************************************************
  //

  // Create Epetra linear problem class to solve "Kx = b"
  Epetra_LinearProblem AmesosProblem;
  AmesosProblem.SetOperator(K.get());
  
  // Create Amesos factory and solver for solving "Kx = b" using a direct factorization
  Amesos amesosFactory;
  Teuchos::RefCountPtr<Amesos_BaseSolver> AmesosSolver = 
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
  int maxBlocks = 3*nev/blockSize;
  int maxRestarts = 5;
  //int step = 5;
  double tol = 1.0e-8;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Max Restarts", maxRestarts );
  MyPL.set( "Tol", tol );
  //MyPL.set( "Step Size", step );
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;
 
  // Create an Anasazi::EpetraMultiVec for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = 
    Teuchos::rcp( new Epetra_MultiVector(K->Map(), blockSize) );
  MVT::MvRandom(*ivec);
  
  // Create the Epetra_Operator for the spectral transformation using the Amesos direct solver.
  Teuchos::RefCountPtr<AmesosGenOp> Aop = Teuchos::rcp( new AmesosGenOp(AmesosProblem,
									AmesosSolver, M) );	
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Aop, M, ivec) );
  
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->SetSymmetric(true);
  
  // Set the number of eigenvalues requested 
  MyProblem->SetNEV( nev );
  
  // Inform the eigenproblem that you are finishing passing it information
  info = MyProblem->SetProblem();
  if (info) {
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;
  }

  // Create a sorting manager to handle the sorting of eigenvalues in the solver
  Teuchos::RefCountPtr<Anasazi::BasicSort<double,MV,OP> > MySort = 
    Teuchos::rcp( new Anasazi::BasicSort<double,MV,OP>(which) );
  
  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchur<double,MV,OP> MySolver(MyProblem, MySort, MyOM, MyPL);
  
  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();

  // Check that the solver returned OK, if not exit example
  if (returnCode != Anasazi::Ok) {
#ifdef EPETRA_MPI
        MPI_Finalize();
#endif
    return -1;
  }
  
  // Obtain eigenvectors directly
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals(); 

  // Retrieve real and imaginary parts of the eigenvectors
  // The size of the eigenvector storage is nev.
  // The real part of the eigenvectors is stored in the first nev vectors.
  // The imaginary part of the eigenvectors is stored in the second nev vectors.
  Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();

  Teuchos::SerialDenseMatrix<int,double> dmat(nev,nev);
  Epetra_MultiVector tempvec(K->Map(), evecs->NumVectors());	
  K->Apply( *evecs, tempvec );
  MVT::MvTransMv( 1.0, *evecs, tempvec, dmat );

  if (MyOM->doPrint()) {
    double compeval = 0.0;
    cout.setf(ios_base::right, ios_base::adjustfield);
    cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
    cout<<"------------------------------------------------------"<<endl;
    cout<<std::setw(16)<<"Real Part" 
	<<std::setw(16)<<"Rayleigh Error"<<endl;
    cout<<"------------------------------------------------------"<<endl;
    for (i=0; i<nev; i++) {
      compeval = dmat(i,i);
      cout <<std::setw(16)<<compeval
	   <<std::setw(16)<<Teuchos::ScalarTraits<double>::magnitude(compeval-1.0/(*evals)[i])
	   <<endl;
    }
    cout<<"------------------------------------------------------"<<endl;
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
			  const Teuchos::RefCountPtr<Amesos_BaseSolver>& solver,
			  const Teuchos::RefCountPtr<Epetra_Operator>& massMtx,
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
