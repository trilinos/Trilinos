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

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
//#include <Thyra_BelosLinearOpWithSolveFactory.hpp>
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
//#include <Thyra_BelosLinearOpWithSolve_def.hpp>
//#include <Thyra_BelosLinearOpWithSolveFactory_def.hpp>

// Stratimikos includes
//#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_FROSchXpetra.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_FROSchFactory_def.hpp"
#include <FROSch_Tools_def.hpp>

#include "EpetraExt_RowMatrixOut.h"


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType DefaultNode;
typedef DefaultNode NO;


using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;
using namespace Thyra;

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
        My_CLP.setOption("DIM",&Dimension,"Dimension.");
        int Overlap = 0;
        My_CLP.setOption("O",&Overlap,"Overlap.");
        int NumberOfBlocks = 1;
        My_CLP.setOption("NB",&NumberOfBlocks,"Number of blocks.");
        int DofsPerNode = 1;
        My_CLP.setOption("DPN",&DofsPerNode,"Dofs per node.");
        int DOFOrdering = 0;
        My_CLP.setOption("ORD",&DOFOrdering,"Dofs ordering (NodeWise=0, DimensionWise=1, Custom=2).");
        string xmlFile = "ParameterList.xml";
        My_CLP.setOption("PLIST",&xmlFile,"File name of the parameter list.");
        bool useepetra = false;
        My_CLP.setOption("USEEPETRA","USETPETRA",&useepetra,"Use Epetra infrastructure for the linear algebra.");
        
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
            assert(false);
        }
        
        UnderlyingLib xpetraLib = UseTpetra;
        if (useepetra) {
            xpetraLib = UseEpetra;
        } else {
            xpetraLib = UseTpetra;
        }
        
        MPI_Comm_split(CommWorld.Comm(),color,CommWorld.MyPID(),&COMM);
        RCP<Epetra_MpiComm> Comm(new Epetra_MpiComm(COMM));
        RCP<const Teuchos::Comm<int> > TeuchosComm = rcp(new MpiComm<int> (COMM));
        if (color==0) {
            
            RCP<ParameterList> parameterList = getParametersFromXmlFile(xmlFile);
            
            ArrayRCP<RCP<Matrix<SC,LO,GO,NO> > > K(NumberOfBlocks);
            ArrayRCP<RCP<Map<LO,GO,NO> > > RepeatedMaps(NumberOfBlocks);
            ArrayRCP<RCP<MultiVector<SC,LO,GO,NO> > > Coordinates(NumberOfBlocks);
            ArrayRCP<UN> dofsPerNodeVector(NumberOfBlocks);
            
            for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                Comm->Barrier();
                if (Comm->MyPID()==0) cout << "###################\n# Assembly Block " << block << " #\n###################\n" << endl;
                
                dofsPerNodeVector[block] = max(DofsPerNode-block,(UN) 1);
                
                ParameterList GalerList;
                GalerList.set("nx", int(N*(block+1)*M));
                GalerList.set("ny", int(N*(block+1)*M));
                GalerList.set("nz", int(N*(block+1)*M));
                GalerList.set("mx", int(N));
                GalerList.set("my", int(N));
                GalerList.set("mz", int(N));
                
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
                
                Coordinates[block] = ConvertToXpetra<SC,LO,GO,NO>(xpetraLib,*epCoord,TeuchosComm);
                
                RCP<Map<LO,GO,NO> > UniqueMap;
                
                if (DOFOrdering == 0) {
                    Array<GO> uniqueMapArray(dofsPerNodeVector[block]*UniqueMapEpetra->NumMyElements());
                    for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                        for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                            uniqueMapArray[dofsPerNodeVector[block]*i+j] = dofsPerNodeVector[block]*UniqueMapEpetra->GID(i)+j;
                        }
                    }
                    
                    UniqueMap = MapFactory<LO,GO,NO>::Build(xpetraLib,-1,uniqueMapArray(),0,TeuchosComm);
                    K[block] = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
                    for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                        LO numEntries;
                        GO* indices;
                        SC* values;
                        KEpetra->ExtractMyRowView(i,numEntries,values,indices);
                        
                        for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                            Array<GO> indicesArray(numEntries);
                            ArrayView<SC> valuesArrayView(values,numEntries);
                            for (LO k=0; k<numEntries; k++) {
                                indicesArray[k] = dofsPerNodeVector[block]*KEpetra->ColMap().GID(indices[k])+j;
                            }
                            K[block]->insertGlobalValues(dofsPerNodeVector[block]*KEpetra->RowMap().GID(i)+j,indicesArray(),valuesArrayView);
                        }
                    }
                    K[block]->fillComplete();
                } else if (DOFOrdering == 1) {
                    Array<GO> uniqueMapArray(dofsPerNodeVector[block]*UniqueMapEpetra->NumMyElements());
                    for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                        for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                            uniqueMapArray[i+UniqueMapEpetra->NumMyElements()*j] = UniqueMapEpetra->GID(i)+(UniqueMapEpetra->MaxAllGID()+1)*j;
                        }
                    }
                    
                    UniqueMap = MapFactory<LO,GO,NO>::Build(xpetraLib,-1,uniqueMapArray(),0,TeuchosComm);
                    K[block] = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
                    for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                        LO numEntries;
                        GO* indices;
                        SC* values;
                        KEpetra->ExtractMyRowView(i,numEntries,values,indices);
                        
                        for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                            Array<GO> indicesArray(numEntries);
                            ArrayView<SC> valuesArrayView(values,numEntries);
                            for (LO k=0; k<numEntries; k++) {
                                indicesArray[k] = KEpetra->ColMap().GID(indices[k])+(KEpetra->ColMap().MaxAllGID()+1)*j;
                            }
                            K[block]->insertGlobalValues(UniqueMapEpetra->GID(i)+(UniqueMapEpetra->MaxAllGID()+1)*j,indicesArray(),valuesArrayView);
                        }
                    }
                    K[block]->fillComplete();
                } else if (DOFOrdering == 2) {
                    assert(false); // TODO: Andere Sortierung implementieren
                } else {
                    assert(false);
                }
                
                RepeatedMaps[block] = FROSch::BuildRepeatedMap<SC,LO,GO,NO>(K[block]); //Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); RepeatedMaps[block]->describe(*fancy,Teuchos::VERB_EXTREME);
            }
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "##############################\n# Assembly Monolythic System #\n##############################\n" << endl;
            
            RCP<Matrix<SC,LO,GO,NO> > KMonolithic;
            if (NumberOfBlocks>1) {
                
                Array<GO> uniqueMapArray(0);
                GO tmpOffset = 0;
                for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                    Teuchos::ArrayView<const GO> tmpGIDs = K[block]->getMap()->getNodeElementList();
                    for (LO i=0; i<tmpGIDs.size(); i++) {
                        uniqueMapArray.push_back(tmpGIDs[i]+tmpOffset);
                    }
                    tmpOffset += K[block]->getMap()->getMaxAllGlobalIndex();
                }
                RCP<Map<LO,GO,NO> > UniqueMapMonolithic = MapFactory<LO,GO,NO>::Build(xpetraLib,-1,uniqueMapArray(),0,TeuchosComm);
                
                tmpOffset = 0;
                KMonolithic = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMapMonolithic,K[0]->getGlobalMaxNumRowEntries());
                for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                    for (LO i=0; i<(LO) K[block]->getNodeNumRows(); i++) {
                        ArrayView<const LO> indices;
                        ArrayView<const SC> values;
                        K[block]->getLocalRowView(i,indices,values);
                        Array<GO> indicesGlobal(indices.size());
                        for (UN j=0; j<indices.size(); j++) {
                            indicesGlobal[j] = K[block]->getColMap()->getGlobalElement(indices[j])+tmpOffset;
                        }
                        KMonolithic->insertGlobalValues(K[block]->getMap()->getGlobalElement(i)+tmpOffset,indicesGlobal(),values);
                    }
                    tmpOffset += K[block]->getMap()->getMaxAllGlobalIndex();
                }
                KMonolithic->fillComplete();
            } else if (NumberOfBlocks==1) {
                KMonolithic = K[0];
            } else {
                assert(false);
            }

            RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(KMonolithic->getMap(),1);
            RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(KMonolithic->getMap(),1);
            
            xSolution->putScalar(0.0);
            xRightHandSide->putScalar(1.0);
      
            CrsMatrixWrap<SC,LO,GO,NO>& crsWrapK = dynamic_cast<CrsMatrixWrap<SC,LO,GO,NO>&>(*KMonolithic);
            RCP<const LinearOpBase<SC> > K_thyra = ThyraUtils<SC,LO,GO,NO>::toThyra(crsWrapK.getCrsMatrix());
            RCP<MultiVectorBase<SC> >thyraX = Teuchos::rcp_const_cast<MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
            RCP<const MultiVectorBase<SC> >thyraB = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);

            //-----------Set Coordinates and RepMap in ParameterList--------------------------
            RCP<ParameterList> plList =  sublist(parameterList,"Preconditioner Types");
            if (NumberOfBlocks>1) {
                sublist(plList,"FROSch")->set("Repeated Map Vector",RepeatedMaps);
                
                sublist(plList,"FROSch")->set("Dimension",Dimension);
                sublist(plList,"FROSch")->set("Overlap",Overlap);
                
                ArrayRCP<DofOrdering> dofOrderings(NumberOfBlocks);
                if (DOFOrdering == 0) {
                    for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                        dofOrderings[block] = NodeWise;
                    }
                } else if (DOFOrdering == 1) {
                    for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                        dofOrderings[block] = DimensionWise;
                    }
                } else {
                    assert(false);
                }
                
                sublist(plList,"FROSch")->set("DofOrdering Vector",dofOrderings);
                sublist(plList,"FROSch")->set("DofsPerNode Vector",dofsPerNodeVector);
            } else if (NumberOfBlocks==1) {
                sublist(plList,"FROSch")->set("Repeated Map",RepeatedMaps[0]);
                // sublist(plList,"FROSch")->set("Coordinates List",Coordinates[0]); // Does not work yet...
                
                sublist(plList,"FROSch")->set("Dimension",Dimension);
                sublist(plList,"FROSch")->set("Overlap",Overlap);
                
                string DofOrderingString;
                if (DOFOrdering == 0) {
                    DofOrderingString = "NodeWise";
                } else if (DOFOrdering == 1) {
                    DofOrderingString = "DimensionWise";
                } else {
                    assert(false);
                }
                sublist(plList,"FROSch")->set("DofOrdering",DofOrderingString);
                sublist(plList,"FROSch")->set("DofsPerNode",DofsPerNode);
            } else {
                assert(false);
            }
            
            if(Comm->MyPID()==0) {
                cout << "##################\n# Parameter List #\n##################" << endl;
                parameterList->print(cout);
                cout << endl;
            }
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "###################################\n# Stratimikos LinearSolverBuilder #\n###################################\n" << endl;
            Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
            Stratimikos::enableFROSch<LO,GO,NO>(linearSolverBuilder);
            linearSolverBuilder.setParameterList(parameterList);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "######################\n# Thyra PrepForSolve #\n######################\n" << endl;
            
            RCP<LinearOpWithSolveFactoryBase<SC> > lowsFactory =
            linearSolverBuilder.createLinearSolveStrategy("");
            
            lowsFactory->setOStream(out);
            lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "###########################\n# Thyra LinearOpWithSolve #\n###########################" << endl;

            RCP<LinearOpWithSolveBase<SC> > lows =
                 linearOpWithSolve(*lowsFactory, K_thyra);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "\n#########\n# Solve #\n#########" << endl;
            SolveStatus<double> status =
            solve<double>(*lows, Thyra::NOTRANS, *thyraB, thyraX.ptr());
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "\n#############\n# Finished! #\n#############" << endl;
            
        }
        MPI_Comm_free(&COMM);
    
    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);

}



