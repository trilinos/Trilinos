#include <iostream>
// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra includes
#include <Epetra_Vector.h>


// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_FROSchHelpers.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// MueLu includes
#include <Thyra_FROSch_TwoLevelPreconditionerFactory_def.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
// Galeri
#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
#include <Galeri_Utils.h>




#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
typedef EpetraNode NO;


using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace Belos;


int main(int argc, char *argv[])
{
    
    using Teuchos::RCP; // reference count pointers
    using Teuchos::rcp;
    using Teuchos::TimeMonitor;
    
    // TUTORIALSPLIT ===========================================================
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    Teuchos::RCP<Teuchos::FancyOStream> fancy = fancyOStream(Teuchos::rcpFromRef(std::cout));

        bool success = false;
        RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
        int MyPID   = comm->getRank();
        int NumProc = comm->getSize();
        
        const Teuchos::RCP<Epetra_Comm> epComm = Teuchos::rcp_const_cast<Epetra_Comm>(Xpetra::toEpetra(comm));
        
        // TUTORIALSPLIT ===========================================================
        // ================================
        // Convenient definitions
        // ================================
        //SC zero = Teuchos::ScalarTraits<SC>::zero();
        SC one = Teuchos::ScalarTraits<SC>::one();
        
        // Instead of checking each time for rank, create a rank 0 stream
        Teuchos::FancyOStream& fancyout = *fancy;
        fancyout.setOutputToRootOnly(0);
        
        
        
        // ================================
        // Parameters initialization
        // ================================
        Teuchos::CommandLineProcessor clp(false);
        GO nx                    = 100;           clp.setOption("nx",                       &nx, "mesh size in x direction");
        GO ny                    = 100;           clp.setOption("ny",                       &ny, "mesh size in y direction");
        int mgridSweeps          = 1;             clp.setOption("mgridSweeps",     &mgridSweeps, "number of multigrid sweeps within Multigrid solver.");
        std::string printTimings = "no";          clp.setOption("timings",        &printTimings, "print timings to screen [yes/no]");
        double tol               = 1e-12;         clp.setOption("tol",                     &tol, "solver convergence tolerance");
        int importOldData = 0;                    clp.setOption("importOldData", &importOldData, "import map and matrix from previous run (highly experimental).");
        
        std::string xmlFileName = "stratimikos_ParameterList.xml"; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'stratimikos_ParameterList.xml'.");
        
        switch (clp.parse(argc,argv)) {
            case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
            case Teuchos::CommandLineProcessor::PARSE_ERROR:
            case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
            case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
        }
        
        
        
        Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
        // ================================
        // Problem construction
        // ================================
        RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time"))), tm;
        
        comm->barrier();
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));
        
        Teuchos::ParameterList GaleriList;
        GaleriList.set("nx", nx);
        GaleriList.set("ny", ny);
        GaleriList.set("mx", epComm->NumProc());
        GaleriList.set("my", 1);
        GaleriList.set("lx", 1.0); // length of x-axis
        GaleriList.set("ly", 1.0); // length of y-axis
        
        Teuchos::RCP<Epetra_Map> epMap = Teuchos::null;
        Teuchos::RCP<Epetra_MultiVector> epCoord = Teuchos::null;
        Teuchos::RCP<Epetra_CrsMatrix> epA = Teuchos::null;
        
    
            // TUTORIALSPLIT ===========================================================
            // create map
            epMap = Teuchos::rcp(Galeri::CreateMap("Cartesian2D", *epComm, GaleriList));
            
            // create coordinates
            epCoord = Teuchos::rcp(Galeri::CreateCartesianCoordinates("2D", epMap.get(), GaleriList));
            
            // create matrix
            epA = Teuchos::rcp(Galeri::CreateCrsMatrix("Laplace2D", epMap.get(), GaleriList));
            
            double hx = 1./(nx-1); double hy = 1./(ny-1);
            epA->Scale(1./(hx*hy));
    
    
        // TUTORIALSPLIT ===========================================================
        // Epetra -> Xpetra
        Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO > > exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int,NO>(epA));
    Teuchos::RCP<CrsMatrixWrap<SC,LO,GO> > exAWrap = Teuchos::rcp(new CrsMatrixWrap<SC,LO,GO> (exA));
        
        RCP<Matrix<SC,LO,GO,NO> > A = Teuchos::rcp_dynamic_cast<Matrix<SC,LO,GO,NO>>(exAWrap);
        A->SetFixedBlockSize(1);
        
        // TUTORIALSPLIT ===========================================================
        // set rhs and solution vector
        RCP<Epetra_Vector> B = Teuchos::rcp(new Epetra_Vector(*epMap));
        RCP<Epetra_Vector> X = Teuchos::rcp(new Epetra_Vector(*epMap));
        B->PutScalar(1.0);
        X->PutScalar(0.0);
        
        // Epetra -> Xpetra
        RCP<Vector<SC,LO,GO,NO> > xB = Teuchos::rcp(new Xpetra::EpetraVectorT<int,NO>(B));
        RCP<Vector<SC,LO,GO,NO> > xX = Teuchos::rcp(new Xpetra::EpetraVectorT<int,NO>(X));
        RCP<MultiVector<SC,LO,GO,NO> > coords = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int,NO>(epCoord));
        
        xX->setSeed(100);
        xX->randomize();
        
        // TUTORIALSPLIT ===========================================================
        // build null space vector
        RCP<const Map<LO,GO,NO> > map = A->getRowMap();
    
        // TUTORIALSPLIT ===========================================================
        comm->barrier();
        tm = Teuchos::null;
        
        fancyout << "========================================================\nGaleri complete.\n========================================================" << std::endl;
        
        // ================================
        //
        // Build Thyra linear algebra objects
        //
        
        RCP<const Thyra::LinearOpBase<SC> > thyraA = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(exAWrap->getCrsMatrix());
        
        RCP<      Thyra::VectorBase<SC> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<SC> >(Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraVector(xX));
        RCP<const Thyra::VectorBase<SC> >thyraB = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraVector(xB);
    
    //
    //Stratimikos
    //
    
    
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
    Stratimikos::enableFROSchTwoLevel<LO,GO,NO>(linearSolverBuilder);

    comm->barrier();
    comm->barrier();
    comm->barrier();

    fancyout <<"Stratimikos enable FROSCh\n";

    // Register FROSCH TwoLevelPreconditioner as a Stratimikos preconditioner strategy.
    linearSolverBuilder.setParameterList(paramList);              // Setup solver parameters using a Stratimikos parameter list.
    comm->barrier();
    comm->barrier();
    comm->barrier();

    fancyout <<"Stratimikos Op\n";

    // Build a new "solver factory" according to the previously specified parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    comm->barrier();
    comm->barrier();
    comm->barrier();

    fancyout <<"Solver Factory\n";

    //Build Thyra operator corresponding to A^-1 computed using stratimikos
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > thyraInversA = Thyra::linearOpWithSolve(*solverFactory,thyraA);
    comm->barrier();
    comm->barrier();
    comm->barrier();

    fancyout <<"Thyra Op A-1\n";
    //Solve Ax = b
    
    Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*thyraInversA,Thyra::NOTRANS,*thyraB,thyraX.ptr());
    comm->barrier();
    comm->barrier();
    comm->barrier();
    
    fancyout <<"Status\n";
    
    success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);


    MPI_Finalize();

    return(success ? EXIT_SUCCESS:EXIT_FAILURE);
}


