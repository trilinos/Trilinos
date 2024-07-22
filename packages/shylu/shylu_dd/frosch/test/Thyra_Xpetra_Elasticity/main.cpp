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
    string xmlFile = "ParameterList.xml";
    My_CLP.setOption("PLIST",&xmlFile,"File name of the parameter list.");
    bool useepetra = false;
    My_CLP.setOption("USEEPETRA","USETPETRA",&useepetra,"Use Epetra infrastructure for the linear algebra.");
    bool useGeoMap = false;
    My_CLP.setOption("useGeoMap","useAlgMap",&useGeoMap,"Use Geometric Map");
    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Thyra Elasticity Test"));
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

        Comm->barrier();
        if (Comm->getRank()==0) {
            cout << "##################\n# Parameter List #\n##################" << endl;
            parameterList->print(cout);
            cout << endl;
        }

        Comm->barrier(); if (Comm->getRank()==0) cout << "##############################\n# Assembly Laplacian #\n##############################\n" << endl;

        ParameterList GaleriList;
        GaleriList.set("nx", GO(N*M));
        GaleriList.set("ny", GO(N*M));
        GaleriList.set("nz", GO(N*M));
        GaleriList.set("mx", GO(N));
        GaleriList.set("my", GO(N));
        GaleriList.set("mz", GO(N));

        RCP<const Map<LO,GO,NO> > UniqueNodeMap;
        RCP<const Map<LO,GO,NO> > UniqueMap;
        RCP<MultiVector<SC,LO,GO,NO> > Coordinates;
        RCP<Matrix<SC,LO,GO,NO> > K;
        if (Dimension==2) {
            UniqueNodeMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian2D",Comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
            UniqueMap = Xpetra::MapFactory<LO,GO,NO>::Build(UniqueNodeMap,2);
            Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("2D",UniqueMap,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Elasticity2D",UniqueMap,GaleriList);
            K = Problem->BuildMatrix();
        } else if (Dimension==3) {
            UniqueNodeMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",Comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
            UniqueMap = Xpetra::MapFactory<LO,GO,NO>::Build(UniqueNodeMap,3);
            Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("3D",UniqueMap,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Elasticity3D",UniqueMap,GaleriList);
            K = Problem->BuildMatrix();
        }


        RCP<Map<LO,GO,NO> > FullRepeatedMap;
        RCP<Map<LO,GO,NO> > RepeatedMap;
        RCP<const Map<LO,GO,NO> > FullRepeatedMapNode;
        if (useGeoMap) {
            if (Dimension == 2) {
                FullRepeatedMap = BuildRepeatedMapGaleriStruct2D<SC,LO,GO,NO>(K,M,Dimension);
                RepeatedMap = FullRepeatedMap;
            } else if (Dimension == 3) {
                FullRepeatedMapNode = BuildRepeatedMapGaleriStruct3D<SC,LO,GO,NO>(K->getMap(),M,Dimension);
                FullRepeatedMap = BuildMapFromNodeMap(FullRepeatedMapNode,Dimension,NodeWise);
                //FullRepeatedMapNode->describe(*fancy,Teuchos::VERB_EXTREME);
                RepeatedMap = FullRepeatedMap;
            }
        } else {
            RepeatedMap = BuildRepeatedMapNonConst<LO,GO,NO>(K->getCrsGraph());
        }


        RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
        RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);

        xSolution->putScalar(ScalarTraits<SC>::zero());
        xRightHandSide->putScalar(ScalarTraits<SC>::one());

        CrsMatrixWrap<SC,LO,GO,NO>& crsWrapK = dynamic_cast<CrsMatrixWrap<SC,LO,GO,NO>&>(*K);
        RCP<const LinearOpBase<SC> > K_thyra = ThyraUtils<SC,LO,GO,NO>::toThyra(crsWrapK.getCrsMatrix());
        RCP<MultiVectorBase<SC> >thyraX = rcp_const_cast<MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
        RCP<const MultiVectorBase<SC> >thyraB = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);

        //-----------Set Coordinates and RepMap in ParameterList--------------------------
        RCP<ParameterList> plList = sublist(parameterList,"Preconditioner Types");
        sublist(plList,"FROSch")->set("Dimension",Dimension);
        sublist(plList,"FROSch")->set("Overlap",Overlap);
        sublist(plList,"FROSch")->set("DofOrdering","NodeWise");
        sublist(plList,"FROSch")->set("DofsPerNode",Dimension);

        sublist(plList,"FROSch")->set("Repeated Map",RepeatedMap);
        sublist(plList,"FROSch")->set("Coordinates List",Coordinates);

        Comm->barrier();
        if (Comm->getRank()==0) {
            cout << "##################\n# Parameter List #\n##################" << endl;
            parameterList->print(cout);
            cout << endl;
        }

        Comm->barrier(); if (Comm->getRank()==0) cout << "###################################\n# Stratimikos LinearSolverBuilder #\n###################################\n" << endl;
        Stratimikos::LinearSolverBuilder<SC> linearSolverBuilder;
        Stratimikos::enableFROSch<SC,LO,GO,NO>(linearSolverBuilder);
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
    stackedTimer->stop("Thyra Elasticity Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);
    std::string watchrProblemName = std::string("FROSch Elasticity ") + std::to_string(Comm->getSize()) + " ranks";
    auto xmlOut = stackedTimer->reportWatchrXML(watchrProblemName, Comm);
    if (xmlOut.length())
      std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';

    return(EXIT_SUCCESS);

}
