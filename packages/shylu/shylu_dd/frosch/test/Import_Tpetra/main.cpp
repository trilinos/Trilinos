// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>

#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif

#include "FROSch_Tools_def.hpp"


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;

using namespace std;

int main(int argc, char *argv[])
{
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();

    CommandLineProcessor My_CLP;

    RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

    bool useepetra = true;
    My_CLP.setOption("USEEPETRA","USETPETRA",&useepetra,"Use Epetra infrastructure for the linear algebra.");

    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Import Tpetra Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    UnderlyingLib xpetraLib = UseTpetra;
    if (useepetra) {
        xpetraLib = UseEpetra;
    } else {
        xpetraLib = UseTpetra;
    }

    assert(CommWorld->getSize()==8);

    CommWorld->barrier(); if (CommWorld->getRank()==0) cout << "#############\n# Assembly #\n#############\n" << endl;

    ParameterList GaleriList;
    GaleriList.set("nx",GO(8));
    GaleriList.set("ny",GO(8));
    GaleriList.set("nz",GO(8));
    GaleriList.set("mx",GO(2));
    GaleriList.set("my",GO(2));
    GaleriList.set("mz",GO(2));

    RCP<const Map<LO,GO,NO> > uniqueMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",CommWorld,GaleriList);
    RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace3D",uniqueMap,GaleriList);
    RCP<Matrix<SC,LO,GO,NO> > K = Problem->BuildMatrix();

    RCP<const Matrix<SC,LO,GO,NO> > tmpMatrix;
    RCP<const Map<LO,GO,NO> > overlappingMap = MapFactory<LO,GO,NO>::Build(uniqueMap,1);
    FROSch::ExtendOverlapByOneLayer<SC,LO,GO,NO>(K.getConst(),overlappingMap,tmpMatrix,overlappingMap);

    K = MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,tmpMatrix->getGlobalMaxNumRowEntries());

    CommWorld->barrier(); if (CommWorld->getRank()==0) cout << "#############\n# Performing Import #\n#############\n" << endl;
    RCP<Export<LO,GO,NO> > gather = ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);
    K->doImport(*tmpMatrix,*gather,ADD);

    CommWorld->barrier(); if (CommWorld->getRank()==0) cout << "\n#############\n# Finished! #\n#############" << endl;

    CommWorld->barrier();
    stackedTimer->stop("Import Tpetra Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    return(EXIT_SUCCESS);
}
