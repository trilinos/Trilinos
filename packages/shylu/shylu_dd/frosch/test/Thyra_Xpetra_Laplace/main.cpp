// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <ShyLU_DDFROSch_config.h>

#include <mpi.h>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

// Galeri::Xpetra
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Thyra_EpetraLinearOp.hpp>
#endif
#include <Thyra_VectorSpaceBase_def.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>

// Stratimikos includes
#include <Stratimikos_FROSch_def.hpp>

#include <Tpetra_Core.hpp>

// Xpetra include
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_FROSchFactory_def.hpp"
#include <FROSch_Tools_def.hpp>


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Thyra;

int main(int argc, char *argv[])
{
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();

    CommandLineProcessor My_CLP;

    RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

    int M = 3;
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
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Thyra Laplace Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    int N = 0;
    int color=1;
    if (Dimension == 2) {
        N = (int) (pow(CommWorld->getSize(),1/2.) + 100*numeric_limits<double>::epsilon()); // 1/H
        if (CommWorld->getRank()<N*N) {
            color=0;
        }
    } else if (Dimension == 3) {
        N = (int) (pow(CommWorld->getSize(),1/3.) + 100*numeric_limits<double>::epsilon()); // 1/H
        if (CommWorld->getRank()<N*N*N) {
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

    RCP<const Comm<int> > Comm = CommWorld->split(color,CommWorld->getRank());

    if (color==0) {

        RCP<ParameterList> parameterList = getParametersFromXmlFile(xmlFile);

        ArrayRCP<RCP<Matrix<SC,LO,GO,NO> > > K(NumberOfBlocks);
        ArrayRCP<RCP<Map<LO,GO,NO> > > RepeatedMaps(NumberOfBlocks);
        ArrayRCP<RCP<MultiVector<SC,LO,GO,NO> > > Coordinates(NumberOfBlocks);
        ArrayRCP<UN> dofsPerNodeVector(NumberOfBlocks);

        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        for (UN block=0; block<(UN) NumberOfBlocks; block++) {
            Comm->barrier(); if (Comm->getRank()==0) cout << "###################\n# Assembly Block " << block << " #\n###################\n" << endl;

            dofsPerNodeVector[block] = (UN) max(int(DofsPerNode-block),1);

            ParameterList GaleriList;
            GaleriList.set("nx", GO(N*(M+block)));
            GaleriList.set("ny", GO(N*(M+block)));
            GaleriList.set("nz", GO(N*(M+block)));
            GaleriList.set("mx", GO(N));
            GaleriList.set("my", GO(N));
            GaleriList.set("mz", GO(N));

            RCP<const Map<LO,GO,NO> > UniqueMapTmp;
            RCP<MultiVector<SC,LO,GO,NO> > CoordinatesTmp;
            RCP<Matrix<SC,LO,GO,NO> > KTmp;
            if (Dimension==2) {
                UniqueMapTmp = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian2D",Comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
                CoordinatesTmp = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("2D",UniqueMapTmp,GaleriList);
                RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace2D",UniqueMapTmp,GaleriList);
                KTmp = Problem->BuildMatrix();
            } else if (Dimension==3) {
                UniqueMapTmp = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",Comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
                CoordinatesTmp = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("3D",UniqueMapTmp,GaleriList);
                RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace3D",UniqueMapTmp,GaleriList);
                KTmp = Problem->BuildMatrix();
            }

            RCP<Map<LO,GO,NO> > UniqueMap;

            if (DOFOrdering == 0) {
                Array<GO> uniqueMapArray(dofsPerNodeVector[block]*UniqueMapTmp->getLocalNumElements());
                for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
                    for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                        uniqueMapArray[dofsPerNodeVector[block]*i+j] = dofsPerNodeVector[block]*UniqueMapTmp->getGlobalElement(i)+j;
                    }
                }

                UniqueMap = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,uniqueMapArray(),0,Comm);
                K[block] = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KTmp->getGlobalMaxNumRowEntries());
                for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
                    ArrayView<const LO> indices;
                    ArrayView<const SC> values;
                    KTmp->getLocalRowView(i,indices,values);

                    for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                        Array<GO> indicesArray(indices.size());
                        for (LO k=0; k<indices.size(); k++) {
                            indicesArray[k] = dofsPerNodeVector[block]*KTmp->getColMap()->getGlobalElement(indices[k])+j;
                        }
                        K[block]->insertGlobalValues(dofsPerNodeVector[block]*KTmp->getRowMap()->getGlobalElement(i)+j,indicesArray(),values);
                    }
                }
                K[block]->fillComplete();
            } else if (DOFOrdering == 1) {
                Array<GO> uniqueMapArray(dofsPerNodeVector[block]*UniqueMapTmp->getLocalNumElements());
                for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
                    for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                        uniqueMapArray[i+UniqueMapTmp->getLocalNumElements()*j] = UniqueMapTmp->getGlobalElement(i)+(UniqueMapTmp->getMaxAllGlobalIndex()+1)*j;
                    }
                }

                UniqueMap = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,uniqueMapArray(),0,Comm);
                K[block] = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KTmp->getGlobalMaxNumRowEntries());
                for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
                    ArrayView<const LO> indices;
                    ArrayView<const SC> values;
                    KTmp->getLocalRowView(i,indices,values);

                    for (UN j=0; j<dofsPerNodeVector[block]; j++) {
                        Array<GO> indicesArray(indices.size());
                        for (LO k=0; k<indices.size(); k++) {
                            indicesArray[k] = KTmp->getColMap()->getGlobalElement(indices[k])+(KTmp->getColMap()->getMaxAllGlobalIndex()+1)*j;
                        }
                        K[block]->insertGlobalValues(UniqueMapTmp->getGlobalElement(i)+(UniqueMapTmp->getMaxAllGlobalIndex()+1)*j,indicesArray(),values);
                    }
                }
                K[block]->fillComplete();
            } else if (DOFOrdering == 2) {
                assert(false); // TODO: Andere Sortierung implementieren
            } else {
                assert(false);
            }

            RepeatedMaps[block] = BuildRepeatedMapNonConst<LO,GO,NO>(K[block]->getCrsGraph()); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); RepeatedMaps[block]->describe(*fancy,VERB_EXTREME);
        }

        Comm->barrier(); if (Comm->getRank()==0) cout << "##############################\n# Assembly Monolithic System #\n##############################\n" << endl;

        RCP<Matrix<SC,LO,GO,NO> > KMonolithic;
        if (NumberOfBlocks>1) {

            Array<GO> uniqueMapArray(0);
            GO tmpOffset = 0;
            for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                ArrayView<const GO> tmpgetGlobalElements = K[block]->getMap()->getLocalElementList();
                for (LO i=0; i<tmpgetGlobalElements.size(); i++) {
                    uniqueMapArray.push_back(tmpgetGlobalElements[i]+tmpOffset);
                }
                tmpOffset += K[block]->getMap()->getMaxAllGlobalIndex()+1;
            }
            RCP<Map<LO,GO,NO> > UniqueMapMonolithic = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,uniqueMapArray(),0,Comm);

            tmpOffset = 0;
            KMonolithic = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMapMonolithic,K[0]->getGlobalMaxNumRowEntries());
            for (UN block=0; block<(UN) NumberOfBlocks; block++) {
                for (LO i=0; i<(LO) K[block]->getLocalNumRows(); i++) {
                    ArrayView<const LO> indices;
                    ArrayView<const SC> values;
                    K[block]->getLocalRowView(i,indices,values);
                    Array<GO> indicesGlobal(indices.size());
                    for (UN j=0; j<indices.size(); j++) {
                        indicesGlobal[j] = K[block]->getColMap()->getGlobalElement(indices[j])+tmpOffset;
                    }
                    KMonolithic->insertGlobalValues(K[block]->getMap()->getGlobalElement(i)+tmpOffset,indicesGlobal(),values);
                }
                tmpOffset += K[block]->getMap()->getMaxAllGlobalIndex()+1;
            }
            KMonolithic->fillComplete();
        } else if (NumberOfBlocks==1) {
            KMonolithic = K[0];
        } else {
            assert(false);
        }

        RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(KMonolithic->getMap(),1);
        RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(KMonolithic->getMap(),1);

        xSolution->putScalar(ScalarTraits<SC>::zero());
        xRightHandSide->putScalar(ScalarTraits<SC>::one());

        CrsMatrixWrap<SC,LO,GO,NO>& crsWrapK = dynamic_cast<CrsMatrixWrap<SC,LO,GO,NO>&>(*KMonolithic);
        RCP<const LinearOpBase<SC> > K_thyra = ThyraUtils<SC,LO,GO,NO>::toThyra(crsWrapK.getCrsMatrix());
        RCP<MultiVectorBase<SC> >thyraX = rcp_const_cast<MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
        RCP<const MultiVectorBase<SC> >thyraB = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);

        //-----------Set Coordinates and RepMap in ParameterList--------------------------
        RCP<ParameterList> plList =  sublist(parameterList,"Preconditioner Types");
        sublist(plList,"FROSch")->set("Dimension",Dimension);
        sublist(plList,"FROSch")->set("Overlap",Overlap);
        if (NumberOfBlocks>1) {
            sublist(plList,"FROSch")->set("Repeated Map Vector",RepeatedMaps);

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

        Comm->barrier();
        if (Comm->getRank()==0) {
            cout << "##################\n# Parameter List #\n##################" << endl;
            parameterList->print(cout);
            cout << endl;
        }

        Comm->barrier(); if (Comm->getRank()==0) cout << "###################################\n# Stratimikos LinearSolverBuilder #\n###################################\n" << endl;
        Stratimikos::LinearSolverBuilder<SC> linearSolverBuilder;
        Stratimikos::enableFROSch<LO,GO,NO>(linearSolverBuilder);
        linearSolverBuilder.setParameterList(parameterList);

        Comm->barrier(); if (Comm->getRank()==0) cout << "######################\n# Thyra PrepForSolve #\n######################\n" << endl;

        RCP<LinearOpWithSolveFactoryBase<SC> > lowsFactory =
        linearSolverBuilder.createLinearSolveStrategy("");

        lowsFactory->setOStream(out);
        lowsFactory->setVerbLevel(VERB_HIGH);

        Comm->barrier(); if (Comm->getRank()==0) cout << "###########################\n# Thyra LinearOpWithSolve #\n###########################" << endl;

        RCP<LinearOpWithSolveBase<SC> > lows =
        linearOpWithSolve(*lowsFactory, K_thyra);

        Comm->barrier(); if (Comm->getRank()==0) cout << "\n#########\n# Solve #\n#########" << endl;
        SolveStatus<SC> status =
        solve<SC>(*lows, Thyra::NOTRANS, *thyraB, thyraX.ptr());

        Comm->barrier(); if (Comm->getRank()==0) cout << "\n#############\n# Finished! #\n#############" << endl;
    }

    CommWorld->barrier();
    stackedTimer->stop("Thyra Laplace Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    return(EXIT_SUCCESS);

}
