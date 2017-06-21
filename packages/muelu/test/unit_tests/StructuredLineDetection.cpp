// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_StructuredLineDetectionFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_GeneralGeometricPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredLineDetectionFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    out << "version: " << MueLu::Version() << std::endl;

    RCP<StructuredLineDetectionFactory> lineDetectionFact = rcp(new StructuredLineDetectionFactory);
    TEST_EQUALITY(lineDetectionFact != Teuchos::null, true);

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredLineDetectionFactory, LabelLines, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LO maxLevels = 3;
    LO maxIter   = 10;
    GO nx = 4;
    GO ny = 4;
    GO nz = 4;
    GO numPoints = nx*ny*nz;
    Array<GO> gNodesPerDim(3);
    gNodesPerDim[0] = 4;
    gNodesPerDim[1] = 4;
    gNodesPerDim[2] = 4;
    Array<LO> lNodesPerDim(3);
    if(comm->getSize() == 1) {
      lNodesPerDim[0] = nx;
      lNodesPerDim[1] = ny;
      lNodesPerDim[2] = nz;
    } else if(comm->getSize() == 4) {
      lNodesPerDim[0] = nx;
      lNodesPerDim[1] = ny;
      lNodesPerDim[2] = nz/4;
    }

    const RCP<const Map> map = MapFactory::Build(lib,
                                                 gNodesPerDim[0]*gNodesPerDim[1]*gNodesPerDim[2],
                                                 lNodesPerDim[0]*lNodesPerDim[1]*lNodesPerDim[2],
                                                 0, comm);

    Teuchos::ParameterList problemParamList;
    problemParamList.set("nx",gNodesPerDim[0]);
    problemParamList.set("ny",gNodesPerDim[1]);
    problemParamList.set("nz",gNodesPerDim[2]);
    // problemParamList.set("keepBCs", true);

    // create Poisson problem and matrix
    Galeri::Xpetra::Laplace3DProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector> PoissonOnCube(problemParamList, map);
    RCP<Matrix> Op = PoissonOnCube.BuildMatrix();

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (Scalar) 1.0);
    Teuchos::Array<magnitude_type> norms(1);
    nullSpace->norm1(norms);
    if (comm->getRank() == 0) {
      out << "||NS|| = " << norms[0] << std::endl;
    }

    // build coordinates on rank 0
    global_size_t numGlobalElements = numPoints;
    size_t numLocalElements = (comm->getRank() == 0 ? numPoints : 0);
    RCP<const Map> exportMap = MapFactory::Build(lib, numGlobalElements, numLocalElements, 0, comm);
    RCP<Xpetra::MultiVector<double,LO,GO,NO> > source_coordinates = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(exportMap, 3);
    if (comm->getRank() == 0) {
      for(LO k = 0; k < gNodesPerDim[2]; ++k) {
        for(LO j = 0; j < gNodesPerDim[1]; ++j) {
          for(LO i = 0; i < gNodesPerDim[0]; ++i) {
            source_coordinates->getDataNonConst(0)[k*ny*nx+j*nx+i] = i / Teuchos::as<double>(gNodesPerDim[0]-1);
            source_coordinates->getDataNonConst(1)[k*ny*nx+j*nx+i] = j / Teuchos::as<double>(gNodesPerDim[1]-1);
            source_coordinates->getDataNonConst(2)[k*ny*nx+j*nx+i] = k / Teuchos::as<double>(gNodesPerDim[2]-1);
          }
        }
      }
    }
    RCP<Xpetra::MultiVector<double,LO,GO,NO> > coordinates = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(map,3);
    RCP<Xpetra::Export<LO,GO,Node> > coords_exporter = Xpetra::ExportFactory<LO,GO,Node>::Build(exportMap, map);
    coordinates->doExport(*source_coordinates, *coords_exporter, Xpetra::INSERT);

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    // create the factory manager and the factories
    FactoryManager M;

    RCP<Factory>      Pfact  = rcp( new GeneralGeometricPFactory() );
    RCP<Factory>      Rfact  = rcp( new TransPFactory() );
    RCP<Factory>      LDfact = rcp( new StructuredLineDetectionFactory() );
    RCP<Factory>      Tfact  = rcp( new CoordinatesTransferFactory() );
    RCP<RAPFactory>   Acfact = rcp( new RAPFactory() );
    RCP<Factory>      NSfact = rcp( new NullspaceFactory() );

    // Set paramters needed by the factories
    Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
    Pfact->SetParameter("axisPermutation", Teuchos::ParameterEntry(std::string("{0,1,2}")));
    Pfact->SetParameter("order", Teuchos::ParameterEntry(1));

    Tfact->SetParameter("Geometric", Teuchos::ParameterEntry(true));

    // Set interfactory dependencies
    NSfact->SetFactory("Nullspace", Pfact);

    Acfact->AddTransferFactory(Tfact);

    // Set default factories in the manager
    M.SetFactory("P",           Pfact);
    M.SetFactory("R",           Rfact);
    M.SetFactory("A",           Acfact);
    M.SetFactory("Nullspace",   NSfact);
    M.SetFactory("gNodesPerDim", Tfact);
    M.SetFactory("lNodesPerDim", Tfact);
    M.SetFactory("Coordinates",  Tfact);
    M.SetFactory("coarseCoordinates",  Pfact);
    M.SetFactory("gCoarseNodesPerDim", Pfact);
    M.SetFactory("lCoarseNodesPerDim", Pfact);
    M.SetFactory("CoarseNumZLayers",          LDfact);
    M.SetFactory("LineDetection_Layers",      LDfact);
    M.SetFactory("LineDetection_VertLineIds", LDfact);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Jacobi");
    smootherParamList.set("relaxation: sweeps", (LocalOrdinal) 1);
    smootherParamList.set("relaxation: damping factor", (Scalar) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("LINESMOOTHING_BANDEDRELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

    // Set smoothers and coarse solver in the manager
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarseSolveFact);

    // Populate data on the finest level of the hierarchy
    RCP<Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",Op);                      // set fine level matrix
    Finest->Set("Nullspace",nullSpace);       // set null space information for finest level
    Finest->Set("Coordinates", coordinates);  // set fine level coordinates
    Finest->Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
    Finest->Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info
    Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

    // Setup the hierarchy
    H->SetMaxCoarseSize(10);
    H->Setup(M, 0, maxLevels);

    // Define RHS and solution vector
    RCP<MultiVector> X   = MultiVectorFactory::Build(map,1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

    X->putScalar(1.0);
    X->norm2(norms);
    Op->apply(*X,*RHS,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    {
      X->putScalar( (Scalar) 0.0);
      H->Iterate(*RHS,*X,maxIter);
      X->norm2(norms);
      if (comm->getRank() == 0) {out << "||RHS|| = " << norms << std::endl;}
    }

  }

#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredLineDetectionFactory,Constructor,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredLineDetectionFactory,LabelLines,Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>

}
