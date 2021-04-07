//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#include "mpi.h"

#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

#include <Teuchos_DefaultSerialComm.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include <Xpetra_DefaultPlatform.hpp>

#include <FROSch_GDSWInterfacePartitionOfUnity_def.hpp>

#include <FROSch_Tools_decl.hpp>


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = KokkosClassic::DefaultNode::DefaultNodeType;

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;
    using namespace FROSch;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();

    CommandLineProcessor My_CLP;

    RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

    int M = 4;
    My_CLP.setOption("M",&M,"H / h.");
    int Dimension = 3;
    My_CLP.setOption("DIM",&Dimension,"Dimension.");
    bool useepetra = true;
    My_CLP.setOption("USEEPETRA","USETPETRA",&useepetra,"Use Epetra infrastructure for the linear algebra.");

    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Interface Partition of Unity Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    int N=4;
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

    RCP<const Comm<int> > comm = CommWorld->split(color,CommWorld->getRank());

    if (color==0) {

        comm->barrier(); if (comm->getRank()==0) cout << "#############\n# Assembly #\n#############\n" << endl;

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
            UniqueMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian2D",comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
            Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("2D",UniqueMap,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace2D",UniqueMap,GaleriList);
            K = Problem->BuildMatrix();
        } else if (Dimension==3) {
            UniqueMap = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",comm,GaleriList); // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); nodeMap->describe(*fancy,VERB_EXTREME);
            Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map<LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("3D",UniqueMap,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace3D",UniqueMap,GaleriList);
            K = Problem->BuildMatrix();
        }

        comm->barrier(); if (comm->getRank()==0) cout << "#############\n# Constructing Repeated Map #\n#############\n" << endl;
        RCP<const Map<LO,GO,NO> > RepeatedMap = BuildRepeatedMap<LO,GO,NO>(K->getCrsGraph());

        RCP<const Map<LO,GO,NO> > RepeatedNodesMap;
        ArrayRCP<RCP<const Map<LO,GO,NO> > > RepeatedDofMaps;
        BuildDofMaps(RepeatedMap,1,NodeWise,RepeatedNodesMap,RepeatedDofMaps);

        comm->barrier(); if (comm->getRank()==0) cout << "#############\n# Constructing Interface Partition of Unity #\n#############\n" << endl;
        RCP<const Comm<int> > SerialComm = createSerialComm<int>();

        RCP<ParameterList> parameterList = getParametersFromXmlFile("ParametersIPOU.xml");
        RCP<InterfacePartitionOfUnity<SC,LO,GO,NO> > IPOU(new GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>(RepeatedMap->getComm(),SerialComm,Dimension,1,RepeatedNodesMap,RepeatedDofMaps,sublist(parameterList,"GDSW"),All,UN(1)));
        IPOU->sortInterface(K);
        IPOU->computePartitionOfUnity();

        comm->barrier(); if (comm->getRank()==0) cout << "\n#############\n# Finished! #\n#############" << endl;
    }

    CommWorld->barrier();
    stackedTimer->stop("Interface Partition of Unity Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    return(EXIT_SUCCESS);
}
