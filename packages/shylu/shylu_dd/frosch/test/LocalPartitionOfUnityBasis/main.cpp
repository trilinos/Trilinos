// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//#ifdef HAVE_MPI
#include "mpi.h"

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "FROSch_LocalPartitionOfUnityBasis_def.hpp"


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;

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
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Local Partition of Unity Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    UnderlyingLib xpetraLib = UseTpetra;
    if (useepetra) {
        xpetraLib = UseEpetra;
    } else {
        xpetraLib = UseTpetra;
    }

    RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));

    Array<GO> nodes(9);
    Array<GO> dofs(27);
    Array<GO> dofs1(9);
    Array<GO> dofs2(9);
    Array<GO> dofs3(9);

    for (UN i=0; i<9; i++) {
        nodes[i] = i;
        dofs[3*i] = 3*i;
        dofs[3*i+1] = 3*i+1;
        dofs[3*i+2] = 3*i+2;
        dofs1[i] = 3*i;
        dofs2[i] = 3*i+1;
        dofs3[i] = 3*i+2;
    }

    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    RCP<Map<LO,GO,NO> > nodesMap = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,nodes(),0,SerialComm);
    RCP<Map<LO,GO,NO> > dofsMap = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,dofs(),0,SerialComm);
    ArrayRCP<RCP<Map<LO,GO,NO> > > dofsMaps (3);
    dofsMaps[0] = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,dofs1(),0,SerialComm);
    dofsMaps[1] = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,dofs2(),0,SerialComm);
    dofsMaps[2] = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,dofs3(),0,SerialComm);

    ArrayRCP<RCP<const MultiVector<SC,LO,GO,NO> > > partitionOfUnity(2);

    RCP<MultiVector<SC,LO,GO,NO> > tmpVec = MultiVectorFactory<SC,LO,GO,NO>::Build(dofsMap,2);
    tmpVec->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[0]),0,ScalarTraits<SC>::one());
    tmpVec->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[0]),0,ScalarTraits<SC>::one());
    tmpVec->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[0]),0,ScalarTraits<SC>::one());

    tmpVec->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[8]),1,ScalarTraits<SC>::one());
    tmpVec->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[8]),1,ScalarTraits<SC>::one());
    tmpVec->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[8]),1,ScalarTraits<SC>::one());
    partitionOfUnity[0] = tmpVec;

    tmpVec = MultiVectorFactory<SC,LO,GO,NO>::Build(dofsMap,1);
    for (UN i=1; i<8; i++) {
        tmpVec->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),0,ScalarTraits<SC>::one());
        tmpVec->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),0,ScalarTraits<SC>::one());
        tmpVec->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),0,ScalarTraits<SC>::one());
    }
    partitionOfUnity[1] = tmpVec;

    RCP<MultiVector<SC,LO,GO,NO> > nullspace = MultiVectorFactory<SC,LO,GO,NO>::Build(dofsMap,6);
    for (UN i=0; i<9; i++) {
        nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),0,ScalarTraits<SC>::one());
        nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),1,ScalarTraits<SC>::one());
        nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),2,ScalarTraits<SC>::one());

        nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),3,ScalarTraits<SC>::one()*double(i));
        nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),3,-ScalarTraits<SC>::one()*double(i));
        nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),3,ScalarTraits<SC>::zero());

        nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),4,-ScalarTraits<SC>::one()*double(i));
        nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),4,ScalarTraits<SC>::zero());
        nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),4,ScalarTraits<SC>::one()*double(i));

        nullspace->replaceLocalValue(dofsMaps[0]->getGlobalElement(nodes[i]),5,ScalarTraits<SC>::zero());
        nullspace->replaceLocalValue(dofsMaps[1]->getGlobalElement(nodes[i]),5,ScalarTraits<SC>::one()*double(i));
        nullspace->replaceLocalValue(dofsMaps[2]->getGlobalElement(nodes[i]),5,-ScalarTraits<SC>::one()*double(i));
    }

    Array<GO> partitionOfUnityMapVec1(2);
    partitionOfUnityMapVec1[0] = 0;
    partitionOfUnityMapVec1[1] = 1;
    Array<GO> partitionOfUnityMapVec2(1);
    partitionOfUnityMapVec1[0] = 0;

    ArrayRCP<RCP<const Map<LO,GO,NO> > > partitionOfUnityMaps(2);
    partitionOfUnityMaps[0] = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,partitionOfUnityMapVec1(),0,SerialComm);
    partitionOfUnityMaps[1] = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,partitionOfUnityMapVec2(),0,SerialComm);

    RCP<ParameterList> parameterList = getParametersFromXmlFile("Parameters.xml");

    LocalPartitionOfUnityBasis<SC,LO,GO,NO> TestBasis(SerialComm,SerialComm,3,parameterList,nullspace,partitionOfUnity,partitionOfUnityMaps);
    TestBasis.buildLocalPartitionOfUnityBasis();

    RCP<CoarseSpace<SC,LO,GO,NO> > coarseSpace = TestBasis.getLocalPartitionOfUnitySpace();

    for (UN i=0; i<coarseSpace->getAssembledBasis()->getLocalLength(); i++) {
        for (UN j=0; j<coarseSpace->getAssembledBasis()->getNumVectors(); j++) {
            std::cout << coarseSpace->getAssembledBasis()->getData(j)[i] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    CommWorld->barrier();
    stackedTimer->stop("Local Partition of Unity Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    return(EXIT_SUCCESS);
}
