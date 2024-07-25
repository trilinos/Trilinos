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

    int M = 4;
    My_CLP.setOption("M",&M,"H / h.");
    int Dimension = 2;
    My_CLP.setOption("DIM",&Dimension,"Dimension.");
    int Iterations = 2;
    My_CLP.setOption("ITS",&Iterations,"Iterations.");
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
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Gathering Test"));
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

        RCP<const Map<LO,GO,NO> > UniqueMap;
        RCP<MultiVector<SC,LO,GO,NO> > Coordinates;
        RCP<Matrix<SC,LO,GO,NO> > K;
        if (Dimension==2) {
            UniqueMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian2D",Comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
            Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("2D",UniqueMap,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace2D",UniqueMap,GaleriList);
            K = Problem->BuildMatrix();
        } else if (Dimension==3) {
            UniqueMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",Comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
            Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("3D",UniqueMap,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace3D",UniqueMap,GaleriList);
            K = Problem->BuildMatrix();
        }
        RCP<const Map<LO,GO,NO> > RepeatedMap = BuildRepeatedMapNonConst<LO,GO,NO>(K->getCrsGraph());

        Comm->barrier(); if (Comm->getRank()==0) cout << "##############################\n# Set Up Coarse Operator #\n##############################\n" << endl;

        RCP<GDSWCoarseOperator<SC,LO,GO,NO> > CoarseOperator(new GDSWCoarseOperator<SC,LO,GO,NO>(K,parameterList));
        CoarseOperator->initialize(Dimension,RepeatedMap);
        CoarseOperator->compute();

        RCP<MultiVector<SC,LO,GO,NO> > x = MultiVectorFactory<SC,LO,GO,NO>::Build(K->getMap(),1);
        RCP<MultiVector<SC,LO,GO,NO> > y = MultiVectorFactory<SC,LO,GO,NO>::Build(K->getMap(),1);
        x->putScalar(ScalarTraits<SC>::one());
        y->putScalar(ScalarTraits<SC>::zero());

        Comm->barrier(); if (Comm->getRank()==0) cout << "##############################\n# Apply Coarse Operator #\n##############################\n" << endl;
        for (int i=0; i<Iterations; i++) {
            CoarseOperator->apply(*x,*y,true);
            *x = *y;
        }

        Comm->barrier();
        SC norm = y->getVector(0)->norm2();
        if (Comm->getRank()==0) cout << "Norm2: " << norm << endl;
    }

    CommWorld->barrier();
    stackedTimer->stop("Gathering Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    return(EXIT_SUCCESS);

}
