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

#include "AnasaziPetraInterface.hpp"
#include "BelosPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziConfigDefs.hpp"
#include "Ifpack_CrsIct.h"
#include "Epetra_CrsMatrix.h"

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
	int i;

#ifdef EPETRA_MPI

	// Initialize MPI
	MPI_Init(&argc,&argv);

#endif

#ifdef EPETRA_MPI
	Epetra_MpiComm & Comm = *new Epetra_MpiComm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm & Comm = *new Epetra_SerialComm();
#endif

	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

	bool verbose = (MyPID==0);

        //  Dimension of the matrix
        int NumGlobalElements = 1000;

        // Construct a Map that puts approximately the same number of
        // equations on each processor.

        Epetra_Map& Map = *new Epetra_Map(NumGlobalElements, 0, Comm);

        // Get update list and number of local equations from newly created Map.
        
        int NumMyElements = Map.NumMyElements();

        int * MyGlobalElements = new int[NumMyElements];
        Map.MyGlobalElements(MyGlobalElements);

        // Create an integer vector NumNz that is used to build the Petra Matrix.
        // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
        // on this processor
        int * NumNz = new int[NumMyElements];

        // We are building two tridiagonal matrices
        // So we need 2 off-diagonal terms (except for the first and last equation)

        for (i=0; i<NumMyElements; i++) {
                if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
                        NumNz[i] = 2;
                else
                        NumNz[i] = 3;  
        }

        // Create both the stiffness and mass Epetra_Matrix        
        Epetra_CrsMatrix& A = *new Epetra_CrsMatrix(Copy, Map, NumNz);
        Epetra_CrsMatrix& B = *new Epetra_CrsMatrix(Copy, Map, NumNz);

        const double one = 1.0;
        double *ValuesA = new double[2];
        double *ValuesB = new double[2];

	// Set values of stiffness matrix.
        double h = one /(NumGlobalElements + one);
        ValuesA[0] = -one/h; ValuesA[1] = -one/h;
        double diagA = 2.0/h;

	// Set values of mass matrix.
        h = one /(6.0*(NumGlobalElements + one));
        ValuesB[0] = one/h; ValuesB[1] = one/h;
        double diagB = 4.0/h;
        int *Indices = new int[2];
        int NumEntries;

        for (i=0; i<NumMyElements; i++)
        {
                if (MyGlobalElements[i]==0)
                {
                        Indices[0] = 1;
                        NumEntries = 1;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesA+1, Indices)==0);
                        assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesB+1, Indices)==0);
                }
                else if (MyGlobalElements[i] == NumGlobalElements-1)
                {
                        Indices[0] = NumGlobalElements-2;
                        NumEntries = 1;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesA, Indices)==0);
                        assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesB, Indices)==0);
                }
                else
                {
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
                        NumEntries = 2;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesA, Indices)==0);
                        assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesB, Indices)==0);
                }
                // Put in the diagonal entry
                assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &diagA, MyGlobalElements+i)==0);
                assert(B.InsertGlobalValues(MyGlobalElements[i], 1, &diagB, MyGlobalElements+i)==0);
        }
         
        // Finish up
        assert(A.TransformToLocal()==0);
	assert(A.OptimizeStorage()==0);
        A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
        assert(B.TransformToLocal()==0);
	assert(B.OptimizeStorage()==0);
        B.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks


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
        Ifpack_CrsIct* ICT = 0;
        //
        if (Lfill > -1) {
                ICT = new Ifpack_CrsIct(A, dropTol, Lfill);

                ICT->SetAbsoluteThreshold(Athresh);
                ICT->SetRelativeThreshold(Rthresh);
                int initerr = ICT->InitValues(A);
                if (initerr != 0) cout << "InitValues error = " << initerr;
                assert(ICT->Factor() == 0);
        } 
        //
        bool transA = false;
        double Cond_Est;
        ICT->Condest(transA, Cond_Est);
        if (verbose) {
                cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
                cout << endl;
        } 
        Epetra_Operator& prec = dynamic_cast<Epetra_Operator&>(*ICT);
	//
	//*******************************************************
	// Set up Belos Block GMRES operator for inner iteration
	//*******************************************************
	//
	int block = 3;  // blocksize used by linear solver and eigensolver [ not required to be the same ]
        int maxits = NumGlobalElements/block - 1; // maximum number of iterations to run
        double btol = 1.0e-7;  // relative residual tolerance
        //
        // Create the Belos::LinearProblemManager
        //
	Belos::PetraMat<double> BelosMat(&A);
        Belos::PetraPrec<double> BelosPrec(&prec);
	Belos::LinearProblemManager<double> My_LP;
	My_LP.SetOperator( &BelosMat );
	My_LP.SetLeftPrec( &BelosPrec );
	My_LP.SetBlockSize( block );
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
	Belos::EpetraOperator<double> BelosOp( My_LP, My_Test, My_OM, My_List );
	//
	//************************************
	// Start the block Arnoldi iteration
	//************************************
	//
	//  Variables used for the Block Arnoldi Method
	//
	int length = 20;
	int nev = 10;
	double tol = 1.0e-6;
	string which="LM";
	int step = 5;
	int restarts = 5;

	// create a PetraAnasaziVec. Note that the decision to make a view or
	// or copy is determined by the petra constructor called by AnasaziPetraVec.
	// This is possible because I pass in arguements needed by petra.

	Anasazi::PetraVec<double> ivec(Map, block);
	ivec.MvRandom();
    
	// call the ctor that calls the petra ctor for a matrix

	Anasazi::PetraMat<double> Amat(A);
	Anasazi::PetraMat<double> Bmat(B);
	Anasazi::PetraGenOp<double> Aop(BelosOp, B);	
	Anasazi::Eigenproblem<double> MyProblem(&Amat, &Bmat, &Aop, &ivec);

	// initialize the Block Arnoldi solver
	Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, tol, nev, length, block, 
						which, step, restarts);
	
	// inform the solver that the problem is symmetric
	MyBlockArnoldi.setSymmetric(true);
	//MyBlockArnoldi.setDebugLevel(3);
	
#ifdef UNIX
	Epetra_Time & timer = *new Epetra_Time(Comm);
#endif

	// iterate a few steps (if you wish)
	//MyBlockArnoldi.iterate(5);

	// solve the problem to the specified tolerances or length
	MyBlockArnoldi.solve();

#ifdef UNIX
	double elapsed_time = timer.ElapsedTime();
	double total_flops = A.Flops();
	double MFLOPs = total_flops/elapsed_time/1000000.0;
#endif

	// obtain results directly
	double* resids = MyBlockArnoldi.getResiduals();
	double* evalr = MyBlockArnoldi.getEvals(); 

	// retrieve real and imaginary parts of the eigenvectors
	Anasazi::PetraVec<double> evecr(Map, nev);
	MyBlockArnoldi.getEvecs( evecr );

	Teuchos::SerialDenseMatrix<int,double> dmatr(nev,nev);
	MyProblem.AInProd( one, evecr, evecr, dmatr );
	double compeval;

	if (verbose) {
	  // output results to screen
	  MyBlockArnoldi.currentStatus();
	  
	  cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
	  cout<<"Real Part \t Rayleigh Error"<<endl;
	  for (i=0; i<nev; i++) {
	    compeval = dmatr(i,i);
	    cout<<compeval<<"\t"<<abs(compeval-one/evalr[i])<<endl;
	  }

	}

#ifdef UNIX
	if (verbose)
		cout << "\n\nTotal MFLOPs for Arnoldi = " << MFLOPs << " Elapsed Time = "<<  elapsed_time << endl;
#endif


	// Release all objects
        if (resids) delete [] resids;
	if (evalr) delete [] evalr;
	if (ICT) delete ICT;

	delete [] NumNz;
	delete [] MyGlobalElements;
	delete [] Indices;
	delete [] ValuesA;
	delete [] ValuesB;

	delete &A, &B;
	delete &Map;
	delete &Comm;

	return 0;
}
