#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
#include "Galeri_Utils.h"


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
#include <Stratimikos_FROSchXpetra.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_FROSchXpetraFactory_def.hpp"


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
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
        
        Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
        
        int M = 4;
        My_CLP.setOption("M",&M,"H / h.");
        int Dimension = 2;
        My_CLP.setOption("Dim",&Dimension,"Dimension.");
        int Overlap = 1;
        My_CLP.setOption("Overlap",&Overlap,"Overlap.");
        int DofsPerNode = 1;
        My_CLP.setOption("DPN",&DofsPerNode,"Dofs per node.");
        int DOFOrdering = 0;
        My_CLP.setOption("Ordering",&DOFOrdering,"Dofs ordering (NodeWise=0,DimensionWise=1,Custom=2).");
        string xmlFile = "xpetra_ParameterList.xml";
        My_CLP.setOption("List",&xmlFile,"File name of the parameter list.");
        
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
        RCP<const Teuchos::Comm<int> > TeuchosComm = rcp(new MpiComm<int> (COMM));
        if (color==0) {
            
            RCP<ParameterList> parameterList = getParametersFromXmlFile(xmlFile);
            
             /*if(Comm->MyPID()==0) {
                cout << "--------------------------------------------------------------------------------\nPARAMETERS:" << endl;
                parameterList->print(cout);
                cout << "--------------------------------------------------------------------------------\n\n";
            }*/
            
            if (Comm->MyPID()==0) cout << "----------------ASSEMBLY-----------\n";
            
            ParameterList GalerList;
            GalerList.set("nx", N*M);
            GalerList.set("ny", N*M);
            GalerList.set("nz", N*M);
            GalerList.set("mx", N);
            GalerList.set("my", N);
            GalerList.set("mz", N);
            
            RCP<Epetra_Map> UniqueMapEpetra;
            RCP<Epetra_CrsMatrix> KEpetra;
            Teuchos::RCP<Epetra_MultiVector> epCoord = Teuchos::null;
            
            if (Dimension==2) {
                UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian2D", *Comm, GalerList));
                KEpetra.reset(Galeri::CreateCrsMatrix("Laplace2D", UniqueMapEpetra.get(), GalerList));
                epCoord = Teuchos::rcp(Galeri::CreateCartesianCoordinates("2D", UniqueMapEpetra.get(), GalerList));

            } else if (Dimension==3) {
                UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian3D", *Comm, GalerList));
                KEpetra.reset(Galeri::CreateCrsMatrix("Laplace3D", UniqueMapEpetra.get(), GalerList));
                epCoord = Teuchos::rcp(Galeri::CreateCartesianCoordinates("3D", UniqueMapEpetra.get(), GalerList));
            }
            
            RCP<Map<LO,GO,NO> > UniqueMap;
            RCP<Matrix<SC,LO,GO,NO> > K;
            
            if (DOFOrdering == 0) {
                Array<GO> uniqueMapArray(DofsPerNode*UniqueMapEpetra->NumMyElements());
                for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                    for (LO j=0; j<DofsPerNode; j++) {
                        uniqueMapArray[DofsPerNode*i+j] = DofsPerNode*UniqueMapEpetra->GID(i)+j;
                    }
                }
                
                UniqueMap = MapFactory<LO,GO,NO>::Build(UseEpetra,-1,uniqueMapArray(),0,TeuchosComm);
                K = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
                for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                    LO numEntries;
                    GO* indices;
                    SC* values;
                    KEpetra->ExtractMyRowView(i,numEntries,values,indices);
                    
                    for (LO j=0; j<DofsPerNode; j++) {
                        Array<GO> indicesArray(numEntries);
                        ArrayView<SC> valuesArrayView(values,numEntries);
                        for (LO k=0; k<numEntries; k++) {
                            indicesArray[k] = DofsPerNode*KEpetra->ColMap().GID(indices[k])+j;
                        }
                        K->insertGlobalValues(DofsPerNode*KEpetra->RowMap().GID(i)+j,indicesArray(),valuesArrayView);
                    }
                }
                K->fillComplete();
            } else if (DOFOrdering == 1) {
                Array<GO> uniqueMapArray(DofsPerNode*UniqueMapEpetra->NumMyElements());
                for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                    for (LO j=0; j<DofsPerNode; j++) {
                        uniqueMapArray[i+UniqueMapEpetra->NumMyElements()*j] = UniqueMapEpetra->GID(i)+(UniqueMapEpetra->MaxAllGID()+1)*j;
                    }
                }
                
                UniqueMap = MapFactory<LO,GO,NO>::Build(UseEpetra,-1,uniqueMapArray(),0,TeuchosComm);
                K = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
                for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                    LO numEntries;
                    GO* indices;
                    SC* values;
                    KEpetra->ExtractMyRowView(i,numEntries,values,indices);
                    
                    for (LO j=0; j<DofsPerNode; j++) {
                        Array<GO> indicesArray(numEntries);
                        ArrayView<SC> valuesArrayView(values,numEntries);
                        for (LO k=0; k<numEntries; k++) {
                            indicesArray[k] = KEpetra->ColMap().GID(indices[k])+(KEpetra->ColMap().MaxAllGID()+1)*j;
                        }
                        K->insertGlobalValues(UniqueMapEpetra->GID(i)+(UniqueMapEpetra->MaxAllGID()+1)*j,indicesArray(),valuesArrayView);
                    }
                }
                K->fillComplete();
            } else if (DOFOrdering == 2) {
                assert(0!=0); // TODO: Andere Sortierung implementieren
            } else {
                assert(0!=0);
            }
            
           
            RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
            RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
            
            xSolution->putScalar(0.0);
            xRightHandSide->putScalar(1.0);
            
            
            Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO > > exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int,NO>(KEpetra));
            Teuchos::RCP<CrsMatrixWrap<SC,LO,GO> > exAWrap = Teuchos::rcp(new CrsMatrixWrap<SC,LO,GO> (exA));
            
            RCP<Xpetra::MultiVector<SC,LO,GO,NO> >Coord = ConvertToXpetra<SC,LO,GO,NO>(UseEpetra,*epCoord,TeuchosComm);
            
            RCP<Xpetra::Map<LO,GO,NO> > RepMapX = FROSch::BuildRepeatedMap<SC,LO,GO,NO>(K);
      
            
            
            RCP<const Thyra::LinearOpBase<SC> > K_thyra = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(exAWrap->getCrsMatrix());
            RCP<Thyra::MultiVectorBase<SC> >thyraX =
            Teuchos::rcp_const_cast<Thyra::MultiVectorBase<SC> >(Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
            RCP<const Thyra::MultiVectorBase<SC> >thyraB = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------Stratimikos LinearSolverBuilder-----------\n";
           
            //-----------Set Coordinates and RepMap in ParameterList--------------------------
           // RCP<ParameterList> plList =  sublist(parameterList,"Preconditioner Types");
           // sublist(plList,"TwoLevelPreconditioner")->set("Coordinates",Coord);
           // sublist(plList,"TwoLevelPreconditioner")->set("RepeatedMap",RepMapX);

            Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
            Stratimikos::enableFROSch<LO,GO,NO>(linearSolverBuilder);
            
            linearSolverBuilder.setParameterList(parameterList);

            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------Thyra PrepForSolve-----------\n";
        
            
            
            RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory =
            linearSolverBuilder.createLinearSolveStrategy("");
            
            lowsFactory->setOStream(out);
            lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            if (Comm->MyPID()==0) cout << "----------------Thyra LinearOpWithSolve-----------\n";

            RCP<Thyra::LinearOpWithSolveBase<SC> > lows =
                 Thyra::linearOpWithSolve(*lowsFactory, K_thyra);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
  
            if (Comm->MyPID()==0) cout << "----------------Solve-----------\n";
            Thyra::SolveStatus<double> status =
            Thyra::solve<double>(*lows, Thyra::NOTRANS, *thyraB, thyraX.ptr());
            
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            //Check Solution
            const RCP<const Thyra::VectorSpaceBase<double> > XY_domain = thyraX->domain();
            
            Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > DomainM = K->getDomainMap();
            
            Teuchos::RCP<const Xpetra::Map<LO,GO,NO> >RangeM = K->getRangeMap();
            
            RCP<const Xpetra::EpetraMapT<GO,NO> > eDomainM = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GO,NO> >(DomainM);
            
            const Epetra_Map epetraDomain = eDomainM->getEpetra_Map();
            
            RCP<const Xpetra::EpetraMapT<GO,NO> > eRangeM = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GO,NO> >(RangeM);
            
            const Epetra_Map epetraRange = eRangeM->getEpetra_Map();
            
            RCP<const Epetra_MultiVector> X;
            
            X = Thyra::get_Epetra_MultiVector(epetraDomain, thyraX );
            std::cout<<*X;
            
        }
        MPI_Comm_free(&COMM);
    
    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);

}



