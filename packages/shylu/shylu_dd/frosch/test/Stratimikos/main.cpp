
#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>

#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
//#include <BelosPseudoBlockGmresSolMgr.hpp>

#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_RGDSWPreconditioner_def.hpp>

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


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode; // Hier Default verwenden???
typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    {
        
        Epetra_MpiComm CommWorld(MPI_COMM_WORLD);
        
        CommandLineProcessor My_CLP;
        
        Teuchos::RCP<Teuchos::FancyOStream> fancy = fancyOStream(Teuchos::rcpFromRef(std::cout));
        
        Teuchos::FancyOStream& fancyout = *fancy;
        fancyout.setOutputToRootOnly(0);
        
        int M = 4;
        My_CLP.setOption("M",&M,"H / h.");
        int Dimension = 3;
        My_CLP.setOption("DIM",&Dimension,"Dimension.");
        int Overlap = 1;
        My_CLP.setOption("OL",&Overlap,"Overlap.");
        bool Reduced = false;
        My_CLP.setOption("RGDSW","GDSW",&Reduced,"Using the reduced coarse space.");
        
        My_CLP.recogniseAllOptions(true);
        My_CLP.throwExceptions(false);
        CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
        if(parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
            MPI_Finalize();
            return 0;
        }
        
        int N;
        MPI_Comm COMM;
        int color=1;
        //bool onFirstLevelComm=false;
        if (Dimension == 2) {
            N = (int) (pow(CommWorld.NumProc(),1/2.) + 100*numeric_limits<double>::epsilon()); // 1/H
            if (CommWorld.MyPID()<N*N) {
                color=0;
            }
        } else if (Dimension == 3) {
            N = (int) (pow(CommWorld.NumProc(),1/3.) + 100*numeric_limits<double>::epsilon()); // 1/H
            if (CommWorld.MyPID()<N*N*N) {
                color=0;
            }
        } else {
            assert(0!=0);
        }
        
        MPI_Comm_split(CommWorld.Comm(),color,CommWorld.MyPID(),&COMM);
        RCP<Epetra_MpiComm> Comm(new Epetra_MpiComm(COMM));
        
        if (color==0) {
            
            RCP<ParameterList> parameterList;
            if (!Reduced) {
                parameterList = getParametersFromXmlFile("stratimikos_ParameterList.xml");
            } else {
                parameterList = getParametersFromXmlFile("stratimikos_ParameterList.xml");
            }
            /*if (Comm->MyPID()==0) {
                cout << "--------------------------------------------------------------------------------\nPARAMETERS:" << endl;
                parameterList->print(cout);
                cout << "--------------------------------------------------------------------------------\n\n";
            }*/
            
            if (Comm->MyPID()==0) cout << "ASSEMBLY...";
            
            ParameterList GalerList;
            GalerList.set("nx", N*M);
            GalerList.set("ny", N*M);
            GalerList.set("nz", N*M);
            GalerList.set("mx", N);
            GalerList.set("my", N);
            GalerList.set("mz", N);
            
            RCP<Epetra_Map> UniqueMap;
            RCP<Epetra_Map> RepeatedMap;
            RCP<Epetra_CrsMatrix> K;
            
            if (Dimension==2) {
                UniqueMap.reset(Galeri::CreateMap("Cartesian2D", *Comm, GalerList));
                K.reset(Galeri::CreateCrsMatrix("Laplace2D", UniqueMap.get(), GalerList));
            } else if (Dimension==3) {
                UniqueMap.reset(Galeri::CreateMap("Cartesian3D", *Comm, GalerList));
                K.reset(Galeri::CreateCrsMatrix("Laplace3D", UniqueMap.get(), GalerList));
            }
            
            EpetraCrsMatrixT<GO,NO> xK(K);
            RCP<CrsMatrix<SC,LO,GO,NO> > xCrsMat = rcpFromRef(xK);
            RCP<Matrix<SC,LO,GO,NO> > xMat = rcp(new CrsMatrixWrap<SC,LO,GO,NO>(xCrsMat));
            xMat->describe(*fancy,Teuchos::VERB_EXTREME);
            if (Comm->MyPID()==0) cout << "done" <<std::endl;
            
            RCP<Epetra_MultiVector> solution(new Epetra_MultiVector(*UniqueMap,1));
            EpetraMultiVectorT<GO,NO> eSolution(solution);
            RCP<MultiVector<SC,LO,GO,NO> > xSolution = rcpFromRef(eSolution);

            
            
            RCP<Epetra_MultiVector> rightHandSide(new Epetra_MultiVector(*UniqueMap,1));
            EpetraMultiVectorT<GO,NO> eRightHandSide(rightHandSide);
            RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = rcpFromRef(eRightHandSide);
            
            xSolution->putScalar(0.0);
            xRightHandSide->putScalar(1.0);
            
            RCP<const Map<LO,GO,NO> > map = xMat->getRowMap();
            Comm->Barrier();
            
            fancyout<<"-------Create Null Space----------\n";
            
            RCP<const Thyra::LinearOpBase<SC> > thyraA = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(xCrsMat
                                                                                                  );
            RCP<Thyra::MultiVectorBase<SC> >thyraX = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<SC> >(Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
            RCP<const Thyra::MultiVectorBase<SC> >thyraB = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);
            
            Comm->Barrier();
            
            fancyout<<"-------To Thyra----------\n";
            //Stratimikos
            Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
            Stratimikos::enableFROSchTwoLevel<LO,GO,NO>(linearSolverBuilder);
            
            Comm->Barrier();
            
            fancyout<<"-------Stratimikos enable FROSch----------\n";
            
            linearSolverBuilder.setParameterList(parameterList);
            
            
            RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
            
            Comm->Barrier();
            
            fancyout<<"-------Thyra_LinearOpWithSolveFactory----------\n";

            Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > thyraInversA = Thyra::linearOpWithSolve(*solverFactory,thyraA);
            
            Comm->Barrier();
            
            fancyout<<"-------Thyra Invers----------\n";
            Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*thyraInversA,Thyra::NOTRANS,*thyraB,thyraX.ptr());
            Comm->Barrier();
            
            fancyout <<"Status\n";
            bool success = false;

            success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
            

        }
        MPI_Comm_free(&COMM);
        
    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);
}

