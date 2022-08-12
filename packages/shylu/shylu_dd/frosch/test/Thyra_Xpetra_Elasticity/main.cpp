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

#define INIT_KOKKOS_HANDLES
#if defined(INIT_KOKKOS_HANDLES)
#include "Kokkos_Core.hpp"
#include "KokkosKernels_config.h"
#include "KokkosKernels_Controls.hpp"

#include "KokkosBlas_trtri.hpp"
#endif

using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Thyra;

// --------------------------------------------------------------------- //
static GO
min3(GO a, GO b, GO c) {
  return std::min(a, std::min(b, c));
}

static GO
max3(GO a, GO b, GO c) {
  return std::max(a, std::max(b, c));
}

static void
cubic_radical_search(GO n, GO & x, GO & y, GO & z) {
  double best = 0.0;

  for (GO f1 = (GO)(pow(n,1.0/3.0)+0.5); f1 > 0; --f1)
    if (n % f1 == 0) {
      GO n1 = n/f1;
      for (GO f2 = (GO)(pow(n1,0.5)+0.5); f2 > 0; --f2)
        if (n1 % f2 == 0) {
          GO f3 = n1 / f2;
          double current = (double)min3(f1, f2, f3)/max3(f1, f2, f3);
          if (current > best) {
            best = current;
            x = f1;
            y = f2;
            z = f3;
          }
        }
    }
}
// --------------------------------------------------------------------- //

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
    bool useIrregularGrid = false;
    My_CLP.setOption("useIrregularGrid","useRegularGrid",&useIrregularGrid,"Use Irregular Process Grid");
    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

for (int count = 0; count < 2; count++) {
    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Thyra Elasticity Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    int color=1;
    int N = 0;
    GO Nx = 0, Ny = 0, Nz = 0;
    GO Mx = 0, My = 0, Mz = 0;
    if (useIrregularGrid) {
        if (Dimension == 3) {
            color=0;
            cubic_radical_search(CommWorld->getSize(), Nx, Ny, Nz);
            Mx = (M+Nx-1)/Nx;
            My = (M+Ny-1)/Ny;
            Mz = (M+Nz-1)/Nz;
        } else {
            assert(false);
        }
    } else {
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
        if (useIrregularGrid) {
            GaleriList.set("nx", GO(M));
            GaleriList.set("ny", GO(M));
            GaleriList.set("nz", GO(M));
            GaleriList.set("mx", GO(Nx));
            GaleriList.set("my", GO(Ny));
            GaleriList.set("mz", GO(Nz));
            if (Comm->getRank()==0) {
              std::cout << Comm->getSize() << " MPIs with M = " << M << " : " 
                        << " #  (nx * ny * nz) = (" << Nx << "x" << Ny << "x" << Nz << "), (mx * my * mz) = (" << Mx << "x" << My << "x" << Mz << ")"
                        << std::endl << std::endl;
            }
        } else {
            GaleriList.set("nx", GO(N*M));
            GaleriList.set("ny", GO(N*M));
            GaleriList.set("nz", GO(N*M));
            GaleriList.set("mx", GO(N));
            GaleriList.set("my", GO(N));
            GaleriList.set("mz", GO(N));
        }
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

Kokkos::fence();
/*if (Comm->getRank() == 10) {
  printf("\n > main <\n" );
  {
    ArrayView<const LO> indices;
    ArrayView<const SC> values;
    K->getLocalRowView(0,indices,values);
    for(size_t k = 0; k < indices.size(); k++) {
      printf("%d: %d %d (%d %d) : %e\n",k,(int)0,(int)indices[k], K->getRowMap()->getGlobalElement(0), K->getColMap()->getGlobalElement(indices[k]), values[k] );
    }
    printf("\n");
  }
}*/
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

#if defined(INIT_KOKKOS_HANDLES)
  KokkosKernels::Experimental::Controls controls;
  #if defined (KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    if(CommWorld->getRank() == 0) std::cout << " > getCuBlasHandle()" << std::endl;
  controls.getCublasHandle();
  #endif
  #if defined (KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
    if(CommWorld->getRank() == 0) std::cout << " > getCuSparseHandle()" << std::endl;
  controls.getCusparseHandle();
  #endif
  #if defined (KOKKOSKERNELS_ENABLE_TPL_MAGMA)
    if(CommWorld->getRank() == 0) std::cout << " > getMagmaHandle()" << std::endl;
  KokkosBlas::Impl::MagmaSingleton & s = KokkosBlas::Impl::MagmaSingleton::singleton();
  {
    Kokkos::View< SC**, Kokkos::DefaultExecutionSpace > dViewL ("Wup",10,10);
    KokkosBlas::trtri("L", "N", dViewL);
  }
  #endif
  Kokkos::fence();
#endif

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
}

    return(EXIT_SUCCESS);

}
