// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_AnasaziPOD.h"
#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziSpecializedEpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "AnasaziGlobalComm.hpp"
#else
#include "Epetra_SerialComm.h"
#endif

namespace RBGen {

  AnasaziPOD::AnasaziPOD() :
    isInitialized_( false ),
    isInner_( true ),
    basis_size_( 16 ),
    comp_time_( 0.0 )
  {
  }

  void AnasaziPOD::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                               const Teuchos::RCP< const Epetra_MultiVector >& ss,
                               const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio )
  {

   // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList& rbmethod_params = params->sublist( "Reduced Basis Method" );

    if ( rbmethod_params.get("Basis Size", 16) < ss->NumVectors() ) {
      basis_size_ = rbmethod_params.get("Basis Size", 16);
    }
    else {
      basis_size_ = ss->NumVectors();
    }
    // Get the inner / outer product form of the operator
    isInner_ = ( rbmethod_params.get("Anasazi POD Operator Form","Inner")=="Inner"? true : false );

    // See if there is a matrix to be used for an inner-product form
    if (rbmethod_params.isParameter( "POD Operator Weighting Matrix" )) {
      std::string matFile = Teuchos::getParameter<std::string>( rbmethod_params, "POD Operator Weighting Matrix" );
      std::vector<std::string> filename(1,matFile);
      op_ = fileio->Read( filename );
    }

    // See if there is a matrix to be used in the orthogonal basis construction
    if (rbmethod_params.isParameter( "Inner Product Matrix" )) {
      std::string matFile = Teuchos::getParameter<std::string>( rbmethod_params, "Inner Product Matrix" );
      std::vector<std::string> filename(1,matFile);
      inner_prod_op_ = fileio->Read( filename );
    }

    // Resize the singular value vector
    sv_.resize( basis_size_ );

    // Set the snapshot set
    ss_ = ss;
  }

  void AnasaziPOD::computeBasis()
  {
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
    Epetra_SerialComm comm;
#endif

    int MyPID = comm.MyPID();
    //
    //  Variables used for the Block Krylov Method
    //
    int step = 5;
    int num_vecs = ss_->NumVectors();
    int vec_length = ss_->GlobalLength();
    //
    //  If the user is requesting more basis vectors than there are snapshots,
    //  compute the basis vectors using an outer product formulation.
    //
    if (basis_size_ > num_vecs) {
       step = 1;
       basis_size_ = num_vecs;
       isInner_ = false;
    }
    Epetra_Time timer( comm );
    int i, blockSize = 1;
    int nev = basis_size_;
    int maxBlocks = 3*basis_size_;
    int maxRestarts = 300;
    int verbosity = Anasazi::Warnings + Anasazi::Errors;
    double tol = 1e-14;
    std::string which="LM";
    //
    // If the user is requesting a large portion of basis vectors, reduce the
    // maximum number of blocks for BKS.
    //
    if (isInner_) {
      if ( maxBlocks > num_vecs )
        maxBlocks = num_vecs-1;
    }
    else {
      if ( maxBlocks > vec_length )
        maxBlocks = vec_length-1;
    }
    //
    // Create parameter list to pass into solver
    //
    Teuchos::ParameterList MyPL;
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Num Blocks", maxBlocks );
    MyPL.set( "Maximum Restarts", maxRestarts );
    MyPL.set( "Which", which );
    MyPL.set( "Step Size", step );
    MyPL.set( "Convergence Tolerance", tol );
    //
    //  Typedefs for Anasazi solvers
    //
    typedef Anasazi::MultiVec<double> MV;
    typedef Anasazi::Operator<double> OP;
    //
    // Create the operator.
    //
    Teuchos::RCP<OP> Amat;
    if (op_ != Teuchos::null) {
      // Call the constructor for the A^T*W*A operator
      Amat = Teuchos::rcp( new Anasazi::EpetraWSymMVOp( ss_, op_ ) );
      isInner_ = true;  // Only inner products are supported at this time.
    }
    else {
      // Call the constructor for the (A^T*A) operator
      Amat = Teuchos::rcp( new Anasazi::EpetraSymMVOp( ss_, !isInner_ ) );
    }
    //
    // Create a map for the columns of the snapshot set.
    //
    Epetra_LocalMap localMap(num_vecs, 0, comm);
    //
    // Create the initial vector and randomize it.
    //
    Teuchos::RCP<MV> ivec;
    if (isInner_) {
      if (inner_prod_op_ != Teuchos::null)
        ivec=Teuchos::rcp( new Anasazi::EpetraOpMultiVec( inner_prod_op_, localMap, blockSize ) );
      else
        ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec( localMap, blockSize ) );
    }
    else {
      if (inner_prod_op_ != Teuchos::null)
        ivec = Teuchos::rcp( new Anasazi::EpetraOpMultiVec( inner_prod_op_, ss_->Map(), blockSize ) );
      else
        ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec( ss_->Map(), blockSize ) );
    }
    ivec->MvRandom();
    //

    // Create the eigenproblem
    Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem =
      Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Amat, ivec) );

    // Inform the eigenproblem that the operator A is symmetric
    MyProblem->setHermitian( true );

    // Set the number of eigenvalues requested
    MyProblem->setNEV( nev );

    // Inform the eigenproblem that you are finishing passing it information
    bool boolret = MyProblem->setProblem();
    if (boolret != true) {
      if (MyPID == 0) {
        std::cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << std::endl;
      }
    }

    // Initialize the Block Arnoldi solver
    Anasazi::BlockKrylovSchurSolMgr<double,MV,OP> MySolverMgr(MyProblem, MyPL);

    timer.ResetStartTime();

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr.solve();
    if (returnCode != Anasazi::Converged && MyPID==0) {
      std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
    }

    comp_time_ = timer.ElapsedTime();

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    int numev = sol.numVecs;

    if (numev < nev && MyPID==0) {
      std::cout << "RBGen::AnasaziPOD will return " << numev << " singular vectors." << std::endl;
    }

    // Resize the singular value vector to the number of computed singular pairs.
    sv_.resize( numev );

    if (numev > 0) {
      // Retrieve eigenvectors
      std::vector<int> index(numev);
      for (i=0; i<numev; i++) { index[i] = i; }
      Teuchos::RCP<const Anasazi::EpetraMultiVec> evecs = Teuchos::rcp_dynamic_cast<const Anasazi::EpetraMultiVec>(
        Anasazi::MultiVecTraits<double,MV>::CloneView(*sol.Evecs, index) );
      // const Anasazi::EpetraMultiVec* evecs = dynamic_cast<const Anasazi::EpetraMultiVec *>(sol.Evecs->CloneView( index ));
      //
      // Compute singular values which are the square root of the eigenvalues
      //
      for (i=0; i<numev; i++) {
        sv_[i] = Teuchos::ScalarTraits<double>::squareroot( Teuchos::ScalarTraits<double>::magnitude(evals[i].realpart) );
      }
      //
      // If we are using inner product formulation,
      // then we must compute left singular vectors:
      //             u = Av/sigma
      //
      std::vector<double> tempnrm( numev );
      if (isInner_) {
        basis_ = Teuchos::rcp( new Epetra_MultiVector(ss_->Map(), numev) );
        Epetra_MultiVector AV( ss_->Map(), numev );

        /* A*V */
        AV.Multiply( 'N', 'N', 1.0, *ss_, *evecs, 0.0 );
        AV.Norm2( &tempnrm[0] );

        /* U = A*V(i)/S(i) */
        Epetra_LocalMap localMap2( numev, 0, ss_->Map().Comm() );
        Epetra_MultiVector S( localMap2, numev );
        for( i=0; i<numev; i++ ) { S[i][i] = 1.0/tempnrm[i]; }
        basis_->Multiply( 'N', 'N', 1.0, AV, S, 0.0 );

        /* Compute direct residuals : || Av - sigma*u ||
           for (i=0; i<nev; i++) { S[i][i] = sv_[i]; }
           info = AV.Multiply( 'N', 'N', -1.0, *basis_, S, 1.0 );
           AV.Norm2( &tempnrm[0] );
           if (MyPID == 0) {
           std::cout<<"Singular Value"<<"\t\t"<<"Direct Residual"<<std::endl;
           std::cout<<"------------------------------------------------------"<<std::endl;
           for (i=0; i<nev; i++) {
           std::cout<< sv_[i] << "\t\t\t" << tempnrm[i] << std::endl;
           }
           std::cout<<"------------------------------------------------------"<<std::endl;
           }
        */
      } else {
        basis_ = Teuchos::rcp( new Epetra_MultiVector( Copy, *evecs, 0, numev ) );
        Epetra_MultiVector ATU( localMap, numev );
        Epetra_MultiVector V( localMap, numev );

        /* A^T*U */
        ATU.Multiply( 'T', 'N', 1.0, *ss_, *evecs, 0.0 );
        ATU.Norm2( &tempnrm[0] );

        /* V = A^T*U(i)/S(i) */
        Epetra_LocalMap localMap2( numev, 0, ss_->Map().Comm() );
        Epetra_MultiVector S( localMap2, numev );
        for( i=0; i<numev; i++ ) { S[i][i] = 1.0/tempnrm[i]; }
        V.Multiply( 'N', 'N', 1.0, ATU, S, 0.0 );

        /* Compute direct residuals : || (A^T)u - sigma*v ||
           for (i=0; i<nev; i++) { S[i][i] = sv_[i]; }
           info = ATU.Multiply( 'N', 'N', -1.0, V, S, 1.0 );
           ATU.Norm2( tempnrm );
           if (comm.MyPID() == 0) {
           std::cout<<"Singular Value"<<"\t\t"<<"Direct Residual"<<std::endl;
           std::cout<<"------------------------------------------------------"<<std::endl;
           for (i=0; i<nev; i++) {
           std::cout<< sv_[i] << "\t\t\t" << tempnrm[i] << std::endl;
           }
           std::cout<<"------------------------------------------------------"<<std::endl;
           }
        */
      }
    }
  }
} // end of RBGen namespace




