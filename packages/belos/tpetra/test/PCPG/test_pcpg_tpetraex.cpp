// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Purpose
// The example tests the successive right-hand sides capabilities of MueLu
// and Belos on a heat flow u_t = u_xx problem.
//
// A sequence of linear systems with the same coefficient matrix and
// different right-hand sides is solved.  A seed space is generated dynamically,
// and a deflated linear system is solved.  After each solves, the first
// few Krylov vectors are saved, and used to reduce the number of iterations
// for later solves.
// The optimal numbers of vectors to deflate and save are not known.
// Presently, the maximum number of vectors to deflate (seed space dimension)
// and to save are user paraemters.
// The seed space dimension is less than or equal to total number of vectors saved.
// The difference between the seed space dimension and the total number of vectors,
// is the number of vectors used to update the seed space after each solve.
// I guess that a seed space whose dimension is a small fraction of the total space
// will be best.
//
// maxSave=1 and maxDeflate=0 uses no recycling (not tested ).
//
// TODO: Instrument with timers, so that we can tell what is going on besides
//       by counting the numbers of iterations.

// Adapted from test_pcpg_epetraex.cpp by David M. Day (with original comments)

// All preconditioning has been commented out

// Belos
#include <BelosPCPGSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_Map_fwd.hpp>
#include <Tpetra_Vector_fwd.hpp>
#include <Tpetra_CrsMatrix_fwd.hpp>
#include <TpetraExt_MatrixMatrix.hpp>

// MueLu
// #include <MueLu_TpetraOperator.hpp>
// #include <MueLu_CreateTpetraPreconditioner.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_CommandLineProcessor.hpp>


template<typename ScalarType>
int run(int argc, char *argv[]) {
    using ST = typename Tpetra::Vector<ScalarType>::scalar_type;
    using LO = typename Tpetra::Vector<>::local_ordinal_type;
    using GO = typename Tpetra::Vector<>::global_ordinal_type;
    using NT = typename Tpetra::Vector<>::node_type;

    using SCT = typename Teuchos::ScalarTraits<ST>;
    using MT  = typename SCT::magnitudeType;
    using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
    using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
    using MVT = typename Belos::MultiVecTraits<ST,MV>;
    using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

    using tcrsmatrix_t   = Tpetra::CrsMatrix<ST,LO,GO,NT>;
    using tmap_t         = Tpetra::Map<LO,GO,NT>;
    using tvector_t      = Tpetra::Vector<ST,LO,GO,NT>;
    using tmultivector_t = Tpetra::MultiVector<ST,LO,GO,NT>;

    using toperator_t  = Tpetra::Operator<ST,LO,GO,NT>;
    // using mtoperator_t = MueLu::TpetraOperator<ST,LO,GO,NT>;

    using starray_t = Teuchos::Array<ST>;
    using goarray_t = Teuchos::Array<GO>;

    using Teuchos::ParameterList; // all of this may look fine but do not be fooled ...
    using Teuchos::RCP;           // it is not so clear what any of this does
    using Teuchos::rcp;
    using Teuchos::Comm;
    using Teuchos::rcp_dynamic_cast;

    Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
    const auto comm = Tpetra::getDefaultComm();
    const int MyPID = comm->getRank();
    const int numProc = comm->getSize();

    // Laplace's equation, homogenous Dirichlet boundary counditions, [0,1]^2
    // regular mesh, Q1 finite elements
    bool success = false;
    bool verbose = false;

    try {
        bool proc_verbose = false;
        int frequency = -1;        // frequency of status test output.
        int blocksize = 1;         // blocksize, PCPGIter
        int numrhs = 1;            // number of right-hand sides to solve for
        int maxiters = 30;         // maximum number of iterations allowed per linear system

        int maxDeflate = 4; // maximum number of vectors deflated from the linear system;
        // There is no overhead cost assoc with changing maxDeflate between solves
        int maxSave = 8;    // maximum number of vectors saved from current and previous .");
        // If maxSave changes between solves, then re-initialize (setSize).

        // Hypothesis: seed vectors are conjugate.
        // Initial versions allowed users to supply a seed space et cetera, but no longer.

        // The documentation is suitable for certain tasks, like defining a modules grammar,
        std::string ortho("ICGS"); // The Belos documentation obscures the fact that
        // IMGS is Iterated Modified Gram Schmidt,
        // ICGS is Iterated Classical Gram Schmidt, and
        // DKGS is another Iterated Classical Gram Schmidt.
        // Mathematical issues, such as the difference between ICGS and DKGS, are not documented at all.
        // UH tells me that Anasazi::SVQBOrthoManager is available;  I need it for Belos
        MT tol = 1.0e-8;           // relative residual tolerance
        // How do command line parsers work?
        Teuchos::CommandLineProcessor cmdp(false,true);

        cmdp.setOption("verbose","quiet",&verbose,"Print messages and results");
        cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters)");
        cmdp.setOption("tol",&tol,"Relative residual tolerance used by PCPG solver");
        cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for");
        cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size)");
        cmdp.setOption("num-deflate",&maxDeflate,"Number of vectors deflated from the linear system");
        cmdp.setOption("num-save",&maxSave,"Number of vectors saved from old Krylov subspaces");
        cmdp.setOption("ortho-type",&ortho,"Orthogonalization type, either DGKS, ICGS or IMGS");

        if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
            return -1;
        }
        if (!verbose)
            frequency = -1;  // reset frequency if test is not verbose

        ////////////////////////////////////////////////////
        //                Form the problem                //
        ////////////////////////////////////////////////////

        int num_time_step = 4;
        GO numElePerDirection = 14 * numProc;
        size_t numNodes = (numElePerDirection - 1)*(numElePerDirection - 1);
        GO base = 0;

        //By the way, either matrix has (3*numElePerDirection - 2)^2 nonzeros.
        RCP<tmap_t> Map         = rcp(new tmap_t(numNodes, base, comm));
        RCP<tcrsmatrix_t> Stiff = rcp(new tcrsmatrix_t(Map, numNodes));
        RCP<tcrsmatrix_t> Mass  = rcp(new tcrsmatrix_t(Map, numNodes));
        RCP<tvector_t> vecLHS   = rcp(new tvector_t(Map));
        RCP<tvector_t> vecRHS   = rcp(new tvector_t(Map));
        RCP<tmultivector_t> LHS, RHS;

        ST ko = 8.0 / 3.0, k1 = -1.0 / 3.0;
        starray_t k_arr(2);
        k_arr[0] = ko;
        k_arr[1] = k1;

        ST h = 1.0 / static_cast<ST>(numElePerDirection);  // x=(iX,iY)h

        ST mo = h*h*4.0/9.0, m1 = h*h/9.0, m2 = h*h/36.0;
        starray_t m_arr(3);
        m_arr[0] = mo;
        m_arr[1] = m1;
        m_arr[2] = m2;

        ST pi = 4.0*atan(1.0), valueLHS;
        GO lid, node, iX, iY;

        goarray_t pos_arr(1);

        for (lid = Map->getMinLocalIndex(); lid <= Map->getMaxLocalIndex(); lid++) {

            node = Map->getGlobalElement(lid);
            iX  = node  % (numElePerDirection-1);
            iY  = ( node - iX )/(numElePerDirection-1);
            pos_arr[0] = node;
            Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(0,1)); // global row ID, global col ID, value
            Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(0,1)); // init guess violates hom Dir bc
            valueLHS = sin( pi*h*((ST) iX+1) )*cos( 2.0 * pi*h*((ST) iY+1) );
            vecLHS->replaceGlobalValue(node, valueLHS);

            if (iY > 0) {
                pos_arr[0] = iX + (iY-1)*(numElePerDirection-1);
                Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); //North
                Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(1,1));
            }

            if (iY < numElePerDirection-2) {
                pos_arr[0] = iX + (iY+1)*(numElePerDirection-1);
                Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); //South
                Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(1,1));
            }

            if (iX > 0) {
                pos_arr[0] = iX-1 + iY*(numElePerDirection-1);
                Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); // West
                Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(1,1));
                if (iY > 0) {
                    pos_arr[0] = iX-1 + (iY-1)*(numElePerDirection-1);
                    Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); // North West
                    Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(2,1));
                }
                if (iY < numElePerDirection-2) {
                    pos_arr[0] = iX-1 + (iY+1)*(numElePerDirection-1);
                    Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); // South West
                    Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(2,1));
                }
            }

            if (iX < numElePerDirection - 2) {
                pos_arr[0] = iX+1 + iY*(numElePerDirection-1);
                Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); // East
                Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(1,1));
                if (iY > 0) {
                    pos_arr[0] = iX+1 + (iY-1)*(numElePerDirection-1);
                    Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); // North East
                    Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(2,1));
                }
                if (iY < numElePerDirection-2) {
                    pos_arr[0] = iX+1 + (iY+1)*(numElePerDirection-1);
                    Stiff->insertGlobalValues(node, pos_arr.view(0,1), k_arr.view(1,1)); // South East
                    Mass->insertGlobalValues(node, pos_arr.view(0,1), m_arr.view(2,1));
                }
            }
        }

        Stiff->fillComplete();
        Mass->fillComplete();

        ST one = 1.0, hdt = .00005; // A = Mass+Stiff*dt/2
        RCP<tcrsmatrix_t> A = Tpetra::MatrixMatrix::add(one, false, *Mass, hdt, false, *Stiff);

        hdt = -hdt; // B = Mass-Stiff*dt/2
        RCP<tcrsmatrix_t> B = Tpetra::MatrixMatrix::add(one, false, *Mass, hdt, false, *Stiff);

        B->apply(*vecLHS, *vecRHS); // rhs_new := B*lhs_old,

        proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

        LHS = Teuchos::rcp_implicit_cast<tmultivector_t>(vecLHS);
        RHS = Teuchos::rcp_implicit_cast<tmultivector_t>(vecRHS);

        ////////////////////////////////////////////////////
        //            Construct Preconditioner            //
        ////////////////////////////////////////////////////

//         ParameterList MueLuList; // Set MueLuList for Smoothed Aggregation

//         MueLuList.set("smoother: type", "CHEBYSHEV");
//         MueLuList.set("smoother: pre or post", "both"); // both pre- and post-smoothing

// #ifdef HAVE_MUELU_AMESOS2
//         MueLuList.set("coarse: type", "KLU2");
// #else
//         MueLuList.set("coarse: type", "none")
// #endif
//         RCP<toperator_t> A_op = A;
//         RCP<mtoperator_t> Prec = MueLu::CreateTpetraPreconditioner(A_op, MueLuList);

        ///////////////////////////////////////////////////
        //             Create Parameter List             //
        ///////////////////////////////////////////////////

        const size_t NumGlobalElements = RHS->getGlobalLength();

        if (maxiters == -1)
            maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run

        ParameterList belosList;
        belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
        belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
        belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
        belosList.set( "Num Deflated Blocks", maxDeflate );    // Number of vectors in seed space
        belosList.set( "Num Saved Blocks", maxSave );          // Number of vectors saved from old spaces
        belosList.set( "Orthogonalization", ortho );           // Orthogonalization type

        if (numrhs > 1) {
        belosList.set( "Show Maximum Residual Norm Only", true );  // although numrhs = 1.
        }
        if (verbose) {
            belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
            Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
            if (frequency > 0)
                belosList.set( "Output Frequency", frequency );
        }
        else
            belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::FinalSummary );

        ///////////////////////////////////////////////////
        //  Construct /*Preconditioned*/ Linear Problem  //
        ///////////////////////////////////////////////////

        RCP<Belos::LinearProblem<ST,MV,OP> > problem
            = rcp( new Belos::LinearProblem<ST,MV,OP>( A, LHS, RHS ) );

        // problem->setLeftPrec( Prec ); // for Preconditioned Problem

        bool set = problem->setProblem();
        if (set == false) {
            if (proc_verbose) {
                std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
            }
            return -1;
        }

        // Create an iterative solver manager.
        RCP< Belos::SolverManager<ST,MV,OP> > solver
        = rcp( new Belos::PCPGSolMgr<ST,MV,OP>(problem, rcp(&belosList,false)) );

        ////////////////////////////////////////////////////
        //                  Iterate PCPG                  //
        ////////////////////////////////////////////////////

        if (proc_verbose) {
            std::cout << std::endl << std::endl;
            std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
            std::cout << "Number of right-hand sides: " << numrhs << std::endl;
            std::cout << "Block size used by solver: " << blocksize << std::endl;
            std::cout << "Maximum number of iterations allowed: " << maxiters << std::endl;
            std::cout << "Relative residual tolerance: " << tol << std::endl;
            std::cout << std::endl;
        }
        bool badRes;
        for( int time_step = 0; time_step < num_time_step; time_step++){
        if (time_step) {
            B->apply(*LHS, *RHS); // rhs_new := B*lhs_old,
            set = problem->setProblem(LHS, RHS);
            if (set == false) {
                if (proc_verbose)
                    std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
                return -1;
            }
        } // if time_step
        std::vector<ST> rhs_norm(numrhs);
        MVT::MvNorm(*RHS, rhs_norm);
        std::cout << "\t\t\t\tRHS norm is ... " << rhs_norm[0] << std::endl;

        // Perform solve

        Belos::ReturnType ret = solver->solve();

        // Compute actual residuals.

        badRes = false;
        std::vector<ST> actual_resids(numrhs);
        tmultivector_t resid(Map, numrhs);
        OPT::Apply( *A, *LHS, resid );
        MVT::MvAddMv( -1.0, resid, 1.0, *RHS, resid );
        MVT::MvNorm( resid, actual_resids );
        MVT::MvNorm( *RHS, rhs_norm );
        std::cout << "\t\t\t\tRHS norm is ... " << rhs_norm[0] << std::endl;

        if (proc_verbose) {
            std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
            for ( int i=0; i<numrhs; i++) {
                ST actRes = actual_resids[i]/rhs_norm[i];
                std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
                if (actRes > tol) badRes = true;
            }
        }

        success = ret==Belos::Converged && !badRes;
        if (!success)
            break;
        } // for time_step

        if (success) {
            if (proc_verbose)
                std::cout << "End Result: TEST PASSED" << std::endl;
        } else {
            if (proc_verbose)
                std::cout << "End Result: TEST FAILED" << std::endl;
        }
    } // try block
    TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

    return (success ? EXIT_SUCCESS : EXIT_FAILURE);
} // run

int main(int argc, char *argv[]) {
    // run with different ST
    run<double>(argc, argv);
    // run<float>(argc, argv); // FAILS -- will need to change tolerance
}
