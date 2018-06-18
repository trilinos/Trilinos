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
#include "Epetra_MpiComm.h"

#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

#include "Teuchos_CommandLineProcessor.hpp"

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsMatrix.hpp"

#include <FROSch_DDInterface_def.hpp>

#include <FROSch_Tools_decl.hpp>

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode;
typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	{

		Epetra_MpiComm CommWorld(MPI_COMM_WORLD);

		CommandLineProcessor My_CLP;

		int M = 4;
		My_CLP.setOption("M",&M,"H / h.");
		int Dimension = 3;
		My_CLP.setOption("DIM",&Dimension,"Dimension.");

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

		if (color==0) {

			if (Comm->MyPID()==0) cout << "ASSEMBLY...";

			ParameterList GalerList;
			GalerList.set("nx", N*M);
			GalerList.set("ny", N*M);
			GalerList.set("nz", N*M);
			GalerList.set("mx", N);
			GalerList.set("my", N);
			GalerList.set("mz", N);

			RCP<Epetra_Map> UniqueMap;
			RCP<Epetra_CrsMatrix> K;
			if (Dimension==2) {
				UniqueMap.reset(Galeri::CreateMap("Cartesian2D", *Comm, GalerList));
				K.reset(Galeri::CreateCrsMatrix("Laplace2D", UniqueMap.get(), GalerList));
			} else if (Dimension==3) {
				UniqueMap.reset(Galeri::CreateMap("Cartesian3D", *Comm, GalerList));
				K.reset(Galeri::CreateCrsMatrix("Laplace3D", UniqueMap.get(), GalerList));
			}
			EpetraCrsMatrixT<GO,NO> xK(K);
			RCP<CrsMatrix<SC,LO,GO,NO> > xCrsMat = rcpFromRef(xK);
			RCP<Matrix<SC,LO,GO,NO> > xMat = rcp(new CrsMatrixWrap<SC,LO,GO,NO>(xCrsMat));

			if (Comm->MyPID()==0) cout << "done" << endl << "CONSTRUCTING REPEATEDMAP...";
			RCP<Map<LO,GO,NO> > RepeatedMap = BuildRepeatedMap<SC,LO,GO,NO>(xMat);
			if (Comm->MyPID()==0) cout << "done" << endl;

			RCP<EntitySet<SC,LO,GO,NO> > vertices,shortEdges,straightEdges,edges,faces,interface,interior;
			//RCP<Map<LO,GO,NO> > verticesMap,shortEdgesMap,straightEdgesMap,edgesMap,facesMap;

			DDInterface<SC,LO,GO,NO> dDInterface(Dimension,1,RepeatedMap);
			dDInterface.divideUnconnectedEntities(xMat);
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

			if (Comm->MyPID()==0) {

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

		}

		MPI_Comm_free(&COMM);

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
