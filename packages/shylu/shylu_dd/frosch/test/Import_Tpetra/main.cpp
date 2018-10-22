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

#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
//#include <BelosPseudoBlockGmresSolMgr.hpp>

#include "FROSch_Tools_def.hpp"

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace Belos;

using namespace std;

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	{

		Epetra_MpiComm CommWorld(MPI_COMM_WORLD);

		RCP<Epetra_MpiComm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
		RCP<const Teuchos::Comm<int> > TeuchosComm = rcp(new MpiComm<int> (MPI_COMM_WORLD));

		assert(TeuchosComm->getSize()==8);

		ParameterList GalerList;
		GalerList.set("nx", 8);
		GalerList.set("ny", 8);
		GalerList.set("nz", 8);
		GalerList.set("mx", 2);
		GalerList.set("my", 2);
		GalerList.set("mz", 2);

		RCP<Epetra_Map> UniqueMapEpetra;
		RCP<Epetra_CrsMatrix> KEpetra;

		UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian3D", *Comm, GalerList));
		KEpetra.reset(Galeri::CreateCrsMatrix("Laplace3D", UniqueMapEpetra.get(), GalerList));

		ArrayView<GO> uniqueMapArrayView(UniqueMapEpetra->MyGlobalElements(),UniqueMapEpetra->NumMyElements());
		RCP<Map<LO,GO,NO> > UniqueMap = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,uniqueMapArrayView,0,TeuchosComm);
		RCP<Matrix<SC,LO,GO,NO> > K = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
		for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
			LO numEntries;
			GO* indices;
			SC* values;
			KEpetra->ExtractMyRowView(i,numEntries,values,indices);

			Array<GO> indicesArray(numEntries);
			ArrayView<SC> valuesArrayView(values,numEntries);
			for (LO j=0; j<numEntries; j++) {
				indicesArray[j] = KEpetra->ColMap().GID(indices[j]);
			}
			K->insertGlobalValues(UniqueMapEpetra->GID(i),indicesArray(),valuesArrayView);
		}
		K->fillComplete(UniqueMap,UniqueMap);

		Teuchos::RCP<Xpetra::Map<LO,GO,NO> > uniqueMap = Xpetra::MapFactory<LO,GO,NO>::Build(UniqueMap,1);
		Teuchos::RCP<Xpetra::Map<LO,GO,NO> > overlappingMap = uniqueMap;
		FROSch::ExtendOverlapByOneLayer<SC,LO,GO,NO>(K,overlappingMap);

		Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > tmpMatrix = K;
		K = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(overlappingMap,2*tmpMatrix->getGlobalMaxNumRowEntries());
		Teuchos::RCP<Xpetra::Export<LO,GO,NO> > gather = Xpetra::ExportFactory<LO,GO,NO>::Build(overlappingMap,uniqueMap);
		TeuchosComm->barrier(); if (TeuchosComm->getRank()==0) std::cout << "BEFORE IMPORT\n";
		K->doImport(*tmpMatrix,*gather,Xpetra::ADD);
		TeuchosComm->barrier(); if (TeuchosComm->getRank()==0) std::cout << "AFTER IMPORT\n";

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
