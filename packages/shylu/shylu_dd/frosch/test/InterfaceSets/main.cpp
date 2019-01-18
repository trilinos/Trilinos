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

#include "Teuchos_CommandLineProcessor.hpp"

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include <Xpetra_DefaultPlatform.hpp>

#include <FROSch_DDInterface_def.hpp>

#include <FROSch_Tools_decl.hpp>

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;
    using namespace FROSch;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

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
    if(parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    int N;
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

        Comm->barrier(); if (Comm->getRank()==0) cout << "#############\n# Assembly #\n#############\n" << endl;

        ParameterList GaleriList;
        GaleriList.set("nx", int(N*M));
        GaleriList.set("ny", int(N*M));
        GaleriList.set("nz", int(N*M));
        GaleriList.set("mx", int(N));
        GaleriList.set("my", int(N));
        GaleriList.set("mz", int(N));

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

        Comm->barrier(); if (Comm->getRank()==0) cout << "#############\n# Constructing Repeated Map #\n#############\n" << endl;
        RCP<Map<LO,GO,NO> > RepeatedMap = BuildRepeatedMap<SC,LO,GO,NO>(K);

        Comm->barrier(); if (Comm->getRank()==0) cout << "#############\n# Identification of Interface Sets #\n#############\n" << endl;
        RCP<EntitySet<SC,LO,GO,NO> > vertices,shortEdges,straightEdges,edges,faces,interface,interior;
        //RCP<Map<LO,GO,NO> > verticesMap,shortEdgesMap,straightEdgesMap,edgesMap,facesMap;

        DDInterface<SC,LO,GO,NO> dDInterface(Dimension,1,RepeatedMap);
        dDInterface.divideUnconnectedEntities(K);
        dDInterface.sortEntities();

        ////////////////////////////////
        // Build Processor Map Coarse //
        ////////////////////////////////
        ArrayRCP<RCP<Map<LO,GO,NO> > > MapVector(5);

        vertices = dDInterface.getVertices();
        vertices->buildEntityMap(RepeatedMap);

        shortEdges = dDInterface.getShortEdges();
        shortEdges->buildEntityMap(RepeatedMap);

        straightEdges = dDInterface.getStraightEdges();
        straightEdges->buildEntityMap(RepeatedMap);

        edges = dDInterface.getEdges();
        edges->buildEntityMap(RepeatedMap);

        faces = dDInterface.getFaces();
        faces->buildEntityMap(RepeatedMap);

        // Vertices
        LO ii=0;
        for (UN i=0; i<1; i++) {
            MapVector[ii] = vertices->getEntityMap();
            ii++;
        }
        // ShortEdges
        for (UN i=0; i<1; i++) {
            MapVector[ii] = shortEdges->getEntityMap();
            ii++;
        }
        // StraightEdges
        for (UN i=0; i<1; i++) {
            MapVector[ii] = straightEdges->getEntityMap();
            ii++;
        }
        // Edges
        for (UN i=0; i<1; i++) {
            MapVector[ii] = edges->getEntityMap();
            ii++;
        }
        // Faces
        for (UN i=0; i<1; i++) {
            MapVector[ii] = faces->getEntityMap();
            ii++;
        }

        vector<LO> NumEntitiesGlobal(5);
        NumEntitiesGlobal[0] = vertices->getEntityMap()->getMaxAllGlobalIndex();
        if (vertices->getEntityMap()->lib()==Xpetra::UseEpetra || vertices->getEntityMap()->getGlobalNumElements()>0) {
            NumEntitiesGlobal[0] += 1;
        }

        NumEntitiesGlobal[1] = shortEdges->getEntityMap()->getMaxAllGlobalIndex();
        if (shortEdges->getEntityMap()->lib()==Xpetra::UseEpetra || shortEdges->getEntityMap()->getGlobalNumElements()>0) {
            NumEntitiesGlobal[1] += 1;
        }

        NumEntitiesGlobal[2] = straightEdges->getEntityMap()->getMaxAllGlobalIndex();
        if (straightEdges->getEntityMap()->lib()==Xpetra::UseEpetra || straightEdges->getEntityMap()->getGlobalNumElements()>0) {
            NumEntitiesGlobal[2] += 1;
        }

        NumEntitiesGlobal[3] = edges->getEntityMap()->getMaxAllGlobalIndex();
        if (edges->getEntityMap()->lib()==Xpetra::UseEpetra || edges->getEntityMap()->getGlobalNumElements()>0) {
            NumEntitiesGlobal[3] += 1;
        }

        NumEntitiesGlobal[4] = faces->getEntityMap()->getMaxAllGlobalIndex();
        if (faces->getEntityMap()->lib()==Xpetra::UseEpetra || faces->getEntityMap()->getGlobalNumElements()>0) {
            NumEntitiesGlobal[4] += 1;
        }

        if (Comm->getRank()==0) {

            cout << "\n\
            --------------------------------------------\n\
            # vertices:       --- " << NumEntitiesGlobal.at(0) << "\n\
            # shortEdges:     --- " << NumEntitiesGlobal.at(1) << "\n\
            # straightEdges:  --- " << NumEntitiesGlobal.at(2) << "\n\
            # edges:          --- " << NumEntitiesGlobal.at(3) << "\n\
            # faces:          --- " << NumEntitiesGlobal.at(4) << "\n\
            --------------------------------------------\n\
            Coarse space:\n\
            --------------------------------------------\n\
            vertices: translations      --- " << true << "\n\
            shortEdges: translations    --- " << true << "\n\
            shortEdges: rotations       --- " << true << "\n\
            straightEdges: translations --- " << true << "\n\
            straightEdges: rotations    --- " << true << "\n\
            edges: translations         --- " << true << "\n\
            edges: rotations            --- " << true << "\n\
            faces: translations         --- " << true << "\n\
            faces: rotations            --- " << true << "\n\
            --------------------------------------------\n";
        }

        Comm->barrier(); if (Comm->getRank()==0) cout << "\n#############\n# Finished! #\n#############" << endl;
    }
    return(EXIT_SUCCESS);
}
