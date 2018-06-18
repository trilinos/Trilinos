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

#include <Teuchos_DefaultSerialComm.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsMatrix.hpp"

#include <FROSch_GDSWInterfacePartitionOfUnity_def.hpp>

#include <FROSch_Tools_decl.hpp>

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode;
typedef EpetraNode NO;

int main(int argc, char *argv[])
{
	using namespace std;
	using namespace Teuchos;
	using namespace Xpetra;
	using namespace FROSch;

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

			RCP<Map<LO,GO,NO> > RepeatedNodesMap;
			ArrayRCP<RCP<Map<LO,GO,NO> > > RepeatedDofMaps;
			BuildDofMaps(RepeatedMap,1,NodeWise,RepeatedNodesMap,RepeatedDofMaps);

			RCP<const Teuchos::Comm<int> > SerialComm = createSerialComm<int>();

			RCP<ParameterList> parameterList = getParametersFromXmlFile("ParametersIPOU.xml");;
			RCP<InterfacePartitionOfUnity<SC,LO,GO,NO> > IPOU(new GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>(RepeatedMap->getComm(),SerialComm,Dimension,1,RepeatedNodesMap,RepeatedDofMaps,parameterList));
			IPOU->removeDirichletNodes();
			IPOU->sortInterface(xMat);
			IPOU->computePartitionOfUnity();
		}

		MPI_Comm_free(&COMM);

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
