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
//  This example computes the eigenvalues of smallest magnitude of the discretized 1D Laplacian
//  equation using the block Implicitly-Restarted Arnoldi method.  This problem shows the 
//  construction of an inner-outer iteration using Belos as the linear solver within Anasazi.  
//  An Ifpack preconditioner is constructed to precondition the linear solver.  This operator 
//  is discretized using finite elements and constructed as an Epetra matrix, then passed into 
//  the AnasaziPetraMat to be used in the construction of the Krylov decomposition.  The 
//  specifics of the block Arnoldi method can be set by the user.

#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Ifpack_CrsIct.h"
#include "Epetra_CrsMatrix.h"

#include "BelosPetraInterface.hpp"
#include "BelosEpetraOperator.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestCombo.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
	int i, info;

#ifdef EPETRA_MPI

	// Initialize MPI
	MPI_Init(&argc,&argv);

#endif

#ifdef EPETRA_MPI
	Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm Comm;
#endif

	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

	bool verbose = (MyPID==0);

        //  Dimension of the matrix
        int NumGlobalElements = 1000;

        // Construct a Map that puts approximately the same number of
        // equations on each processor.

        Epetra_Map Map(NumGlobalElements, 0, Comm);

        // Get update list and number of local equations from newly created Map.
        
        int NumMyElements = Map.NumMyElements();

	std::vector<int> MyGlobalElements(NumMyElements);
        Map.MyGlobalElements(&MyGlobalElements[0]);

        // Create an integer vector NumNz that is used to build the Petra Matrix.
        // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
        // on this processor
	std::vector<int> NumNz(NumMyElements);

        // We are building two tridiagonal matrices
        // So we need 2 off-diagonal terms (except for the first and last equation)

        for (i=0; i<NumMyElements; i++) {
                if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
                        NumNz[i] = 2;
                else
                        NumNz[i] = 3;  
        }

        // Create both the stiffness and mass Epetra_Matrix        
        Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]));
        Teuchos::RefCountPtr<Epetra_CrsMatrix> B = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]));

        const double one = 1.0;
        std::vector<double> ValuesA(2);
        std::vector<double> ValuesB(2);
	std::vector<int> Indices(2);

	// Set values of stiffness matrix.
        double h = one /(NumGlobalElements + one);
        ValuesA[0] = -one/h; ValuesA[1] = -one/h;
        double diagA = 2.0/h;

	// Set values of mass matrix.
        h = one /(6.0*(NumGlobalElements + one));
        ValuesB[0] = one/h; ValuesB[1] = one/h;
        double diagB = 4.0/h;
	int NumEntries;

        for (i=0; i<NumMyElements; i++)
        {
                if (MyGlobalElements[i]==0)
                {
                        Indices[0] = 1;
                        NumEntries = 1;
                        info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &ValuesA[1], &Indices[0]);
			assert( info==0 );
                        info = B->InsertGlobalValues(MyGlobalElements[i], NumEntries, &ValuesB[1], &Indices[0]);
			assert( info==0 );
                }
                else if (MyGlobalElements[i] == NumGlobalElements-1)
                {
                        Indices[0] = NumGlobalElements-2;
                        NumEntries = 1;
                        info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &ValuesA[0], &Indices[0]);
			assert( info==0 );
                        info = B->InsertGlobalValues(MyGlobalElements[i], NumEntries, &ValuesB[0], &Indices[0]);
			assert( info==0 );
                }
                else
                {
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
                        NumEntries = 2;
                        info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &ValuesA[0], &Indices[0]);
			assert( info==0 );
                        info = B->InsertGlobalValues(MyGlobalElements[i], NumEntries, &ValuesB[0], &Indices[0]);
			assert( info==0 );
                }
                // Put in the diagonal entry
                info = A->InsertGlobalValues(MyGlobalElements[i], 1, &diagA, &MyGlobalElements[i]);
		assert( info==0 );
                info = B->InsertGlobalValues(MyGlobalElements[i], 1, &diagB, &MyGlobalElements[i]);
		assert( info==0 );
        }
         
        // Finish up
        info = A->TransformToLocal();
	assert( info==0 );
	info = A->OptimizeStorage();
	assert( info==0 );
        A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
        info = B->TransformToLocal();
	assert( info==0 );
	info = B->OptimizeStorage();
	assert( info==0 );
        B->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
        //*****Select the Preconditioner*****
        //
	if (verbose) cout << endl << endl;
	if (verbose) cout << "Constructing ICT preconditioner" << endl;
        int Lfill = 0;
        if (argc > 2) Lfill = atoi(argv[2]);
        if (verbose) cout << "Using Lfill = " << Lfill << endl;
        int Overlap = 0;
        if (argc > 3) Overlap = atoi(argv[3]);
        if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
        double Athresh = 0.0;
        if (argc > 4) Athresh = atof(argv[4]);
        if (verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
        double Rthresh = 1.0;
        if (argc >5) Rthresh = atof(argv[5]);
        if (verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
        double dropTol = 1.0e-6;
        //
        Teuchos::RefCountPtr<Ifpack_CrsIct> ICT;
        //
        if (Lfill > -1) {
                ICT = Teuchos::rcp( new Ifpack_CrsIct(*A, dropTol, Lfill) );
                ICT->SetAbsoluteThreshold(Athresh);
                ICT->SetRelativeThreshold(Rthresh);
                int initerr = ICT->InitValues(*A);
                if (initerr != 0) cout << "InitValues error = " << initerr;
                info = ICT->Factor();
		assert( info==0 );
        } 
        //
        bool transA = false;
        double Cond_Est;
        ICT->Condest(transA, Cond_Est);
        if (verbose) {
                cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
                cout << endl;
        } 
	//
	//*******************************************************
	// Set up Belos Block GMRES operator for inner iteration
	//*******************************************************
	//
	int blockSize = 3;  // block size used by linear solver and eigensolver [ not required to be the same ]
        int maxits = NumGlobalElements/blockSize - 1; // maximum number of iterations to run
        double btol = 1.0e-7;  // relative residual tolerance
        //
        // Create the Belos::LinearProblemManager
        //
	Belos::PetraMat<double> BelosMat(A.get());
        Belos::PetraPrec<double> BelosPrec(ICT.get());
	Belos::LinearProblemManager<double> My_LP;
	My_LP.SetOperator( &BelosMat );
	My_LP.SetLeftPrec( &BelosPrec );
	My_LP.SetBlockSize( blockSize );
	//
	// Create the Belos::StatusTest
	//
	Belos::StatusTestMaxIters<double> test1( maxits );
	Belos::StatusTestResNorm<double> test2( btol );
	Belos::StatusTestCombo<double> My_Test( Belos::StatusTestCombo<double>::OR, test1, test2 );
	//
	// Create the Belos::OutputManager
	//
	Belos::OutputManager<double> My_OM( MyPID );
	//My_OM.SetVerbosity( 2 );
	//
	// Create the ParameterList for the Belos Operator
	// 
	Teuchos::ParameterList My_List;
	My_List.set( "Solver", "BlockCG" );
	My_List.set( "MaxIters", maxits );
	//
	// Create the Belos::EpetraOperator
	//
	Teuchos::RefCountPtr<Belos::EpetraOperator<double> > BelosOp = 
	  Teuchos::rcp( new Belos::EpetraOperator<double>(My_LP, My_Test, My_OM, My_List ));
	//
	//************************************
	// Start the block Arnoldi iteration
	//************************************
	//
	//  Variables used for the Block Arnoldi Method
	//
	int nev = 10;
	int maxBlocks = 20;
	int maxRestarts = 5;
	int step = 5;
	double tol = 1.0e-6;
	string which="LM";
	//
	// Create parameter list to pass into solver
	//
	Teuchos::ParameterList MyPL;
	MyPL.set( "Block Size", blockSize );
	MyPL.set( "Max Blocks", maxBlocks );
	MyPL.set( "Max Restarts", maxRestarts );
	MyPL.set( "Tol", tol );
	MyPL.set( "Step Size", step );

	typedef Anasazi::MultiVec<double> MV;
	typedef Anasazi::Operator<double> OP;

	// Create a AnasaziEpetraMultiVec for an initial vector to start the solver.
	// Note:  This needs to have the same number of columns as the blocksize.
	Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> ivec = 
	  Teuchos::rcp( new Anasazi::EpetraMultiVec(Map, blockSize) );
	ivec->MvRandom();
    
	// Call the ctor that calls the petra ctor for a matrix
	Teuchos::RefCountPtr<Anasazi::EpetraOp> Amat = Teuchos::rcp( new Anasazi::EpetraOp(A) );
	Teuchos::RefCountPtr<Anasazi::EpetraOp> Bmat = Teuchos::rcp( new Anasazi::EpetraOp(B) );
	Teuchos::RefCountPtr<Anasazi::EpetraGenOp> Aop = Teuchos::rcp( new Anasazi::EpetraGenOp(BelosOp, B) );	

	Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
	  Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Aop, Bmat, ivec) );

	// Inform the eigenproblem that the matrix pencil (A,B) is symmetric
	MyProblem->SetSymmetric(true);

	// Set the number of eigenvalues requested 
	MyProblem->SetNEV( nev );

        // Inform the eigenproblem that you are finishing passing it information
        info = MyProblem->SetProblem();
	if (info)
	  cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;

        // Create a sorting manager to handle the sorting of eigenvalues in the solver
	Teuchos::RefCountPtr<Anasazi::BasicSort<double,MV,OP> > MySort = 
	  Teuchos::rcp( new Anasazi::BasicSort<double,MV,OP>(which) );
	
        // Create an output manager to handle the I/O from the solver
        Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
	  Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
        MyOM->SetVerbosity( Anasazi::FinalSummary );

	// Initialize the Block Arnoldi solver
	Anasazi::BlockKrylovSchur<double,MV,OP> MySolver(MyProblem, MySort, MyOM, MyPL);
	
	// solve the problem to the specified tolerances or length
	MySolver.solve();
	
	// obtain eigenvectors directly
	Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals(); 

	// retrieve real and imaginary parts of the eigenvectors
	// The size of the eigenvector storage is nev.
	// The real part of the eigenvectors is stored in the first nev vectors.
	// The imaginary part of the eigenvectors is stored in the second nev vectors.
	Anasazi::EpetraMultiVec* evecr = dynamic_cast<Anasazi::EpetraMultiVec*>(MyProblem->GetEvecs()->CloneCopy());

	Teuchos::SerialDenseMatrix<int,double> dmatr(nev,nev);
	Anasazi::EpetraMultiVec tempvec(Map, evecr->GetNumberVecs());	
	A->Apply( *evecr, tempvec );
	tempvec.MvTransMv( 1.0, *evecr, dmatr );

	double compeval = 0.0;
	cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
	cout<<"Real Part \t Rayleigh Error"<<endl;
	for (i=0; i<nev; i++) {
		compeval = dmatr(i,i);
		cout<<compeval<<"\t"<<Teuchos::ScalarTraits<double>::magnitude(compeval-one/(*evals)[i])<<endl;
	}

	return 0;
}
