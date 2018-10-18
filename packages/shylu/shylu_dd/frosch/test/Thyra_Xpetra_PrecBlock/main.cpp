#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_HDF5.h"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>

#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>

#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_RGDSWPreconditioner_def.hpp>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_VectorSpaceBase_def.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>
#include <Thyra_BelosLinearOpWithSolve_def.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory_def.hpp>

// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include "Stratimikos_FROSchXpetra.hpp"

// Xpetra include
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"

#include "Thyra_FROSchXpetraTwoLevelBlockPrec_def.hpp"
#include "EpetraExt_HDF5.h"

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
using namespace EpetraExt;

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    {        
        Epetra_MpiComm CommWorld(MPI_COMM_WORLD);
        
        CommandLineProcessor My_CLP;
        
        Teuchos::RCP<Teuchos::FancyOStream>
        
        out = Teuchos::VerboseObjectBase::getDefaultOStream();
        
        string xmlFile = "xpetra_ParameterList.xml";
        My_CLP.setOption("List",&xmlFile,"File name of the parameter list.");

        My_CLP.recogniseAllOptions(true);
        My_CLP.throwExceptions(false);
        CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
        if(parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
            MPI_Finalize();
            return 0;
        }
        
        MPI_Comm COMM;
        int color=0;
        //bool onFirstLevelComm=false;
        
        
        MPI_Comm_split(CommWorld.Comm(),color,CommWorld.MyPID(),&COMM);
        RCP<Epetra_MpiComm> Comm(new Epetra_MpiComm(COMM));
        RCP<const Teuchos::Comm<int> > TeuchosComm = rcp(new MpiComm<int> (COMM));
        if (color==0) {    
        
            ///////////////////
            // ParameterList //
            ///////////////////
            RCP<ParameterList> parameterList = getParametersFromXmlFile(xmlFile);
            if (Comm->MyPID()==0) {
                cout << "--------------------------------------------------------------------------------\nPARAMETERS:" << endl;
                parameterList->print(cout);
                cout << "--------------------------------------------------------------------------------\n\n";
            }
            
            unsigned Dimension = 2;
            RCP<HDF5> hDF5IO(new HDF5(*Comm));
            hDF5IO->Open("stokes.h5");
            string groupNameMatrix = "Matrix";
            string groupNameRepeatedMapVelo =  "RepeatedMapVelocity";
            string groupNameRepeatedMapPress =  "RepeatedMapPressure";
            string groupNameRHS = "RHS";
            
            
            ///////////////////
            // Repeated Maps //
            ///////////////////
            Epetra_Map *repeatedMapEpetraVelo;
            hDF5IO->Read(groupNameRepeatedMapVelo,repeatedMapEpetraVelo);
            RCP<Map<LO,GO,NO> > repeatedMapVelo = ConvertToXpetra<LO,GO,NO>(UseTpetra,*repeatedMapEpetraVelo,TeuchosComm);
            
            Epetra_Map *repeatedMapEpetraPress;
            hDF5IO->Read(groupNameRepeatedMapPress,repeatedMapEpetraPress);
            GO offsetVelocityMap = repeatedMapVelo->getMaxAllGlobalIndex()+1;
            Array<GO> elementList(repeatedMapEpetraPress->NumMyElements());
            for (unsigned i=0; i<elementList.size(); i++) {
                elementList[i] = repeatedMapEpetraPress->GID(i) + offsetVelocityMap;
            }
            RCP<Map<LO,GO,NO> > repeatedMapPress = Xpetra::MapFactory<LO,GO,NO>::Build(UseTpetra,-1,elementList,0,TeuchosComm);
            
            Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > repeatedMapsVector(2);
            
            repeatedMapsVector[0] = repeatedMapVelo;
            repeatedMapsVector[1] = repeatedMapPress;
            
            RCP<Map<LO,GO,NO> > repeatedMap = MergeMaps(repeatedMapsVector);
            
            ////////////////
            // Unique Map //
            ////////////////
            RCP<Map<LO,GO,NO> > uniqueMap = BuildUniqueMap<LO,GO,NO>(repeatedMap);
            
            
            /////////
            // RHS //
            /////////
            Epetra_MultiVector *rhsEpetra;
            hDF5IO->Read(groupNameRHS,rhsEpetra);
            
            RCP<MultiVector<SC,LO,GO,NO> > rhsTmp = ConvertToXpetra<SC,LO,GO,NO>(UseTpetra,*rhsEpetra,TeuchosComm);
            RCP<MultiVector<SC,LO,GO,NO> > rhs = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);
            RCP<Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(rhsTmp->getMap(),uniqueMap);
            rhs->doImport(*rhsTmp,*scatter,Xpetra::ADD);
            
            ////////////
            // Matrix //
            ////////////
            Epetra_CrsMatrix *matrixEpetra;
            hDF5IO->Read(groupNameMatrix,matrixEpetra);
            RCP<Matrix<SC,LO,GO,NO> > matrixTmp = ConvertToXpetra<SC,LO,GO,NO>(UseTpetra,*matrixEpetra,TeuchosComm);
            RCP<Matrix<SC,LO,GO,NO> > matrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,matrixTmp->getGlobalMaxNumRowEntries());
            
            matrix->doImport(*matrixTmp,*scatter,Xpetra::ADD);
            matrix->fillComplete();
           
            RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);
            RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);
            
            xSolution->putScalar(0.0);
            xRightHandSide->putScalar(1.0);
            
            Teuchos::RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > tmpCrsWrap = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >(matrix);
                        
            RCP<const Thyra::LinearOpBase<SC> > K_thyra = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(tmpCrsWrap->getCrsMatrix());
            
            RCP<Thyra::MultiVectorBase<SC> >thyraX =
            Teuchos::rcp_const_cast<Thyra::MultiVectorBase<SC> >(Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
            RCP<const Thyra::MultiVectorBase<SC> >thyraB = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------Stratimikos LinearSolverBuilder-----------\n";
           
            //-----------Set Coordinates and RepMap in ParameterList--------------------------
            RCP<ParameterList> plList =  sublist(parameterList,"Preconditioner Types");
//            sublist(plList,"TwoLevelPreconditioner")->set("Coordinates",Coord);
            
            sublist(plList,"FROSchBlock")->set("RepeatedMap1",repeatedMapsVector[0]);
            sublist(plList,"FROSchBlock")->set("RepeatedMap2",repeatedMapsVector[1]);
            sublist(plList,"FROSchBlock")->set("DofsPerNode1",(int)Dimension);
            sublist(plList,"FROSchBlock")->set("DofsPerNode2",1);
            sublist(plList,"FROSchBlock")->set("Ordering1",FROSch::NodeWise);
            sublist(plList,"FROSchBlock")->set("Ordering2",FROSch::NodeWise);
            
            Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
            Stratimikos::enableFROSch<LO,GO,NO>(linearSolverBuilder,"FROSchBlock");
            
            linearSolverBuilder.setParameterList(parameterList);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------Thyra CreatePrec-----------\n";
            
            RCP<Thyra::PreconditionerFactoryBase<SC> > pfbFactory = linearSolverBuilder.createPreconditioningStrategy("");
            
            RCP<Thyra::PreconditionerBase<SC> > ThyraPrec = Thyra::prec(*pfbFactory,K_thyra);

            RCP<const Thyra::LinearOpBase<SC> > LinearPrecOp = ThyraPrec->getUnspecifiedPrecOp();
            assert(nonnull(LinearPrecOp));
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            if (Comm->MyPID()==0) cout << "----------------Apply Prec-----------\n";
            
            LinearPrecOp->apply(Thyra::NOTRANS,*thyraB,thyraX.ptr(),1.0,0.0);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
        }
        MPI_Comm_free(&COMM);
    
    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);

}



