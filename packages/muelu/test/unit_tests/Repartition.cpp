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
#include <vector>
#include "Teuchos_UnitTestHarness.hpp"

#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_ExportFactory.hpp"
#include "Xpetra_MatrixFactory.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"

#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraMaps.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_ParameterListInterpreter.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  TEUCHOS_UNIT_TEST(Repartition, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    TEST_EQUALITY(repart != Teuchos::null, true);

  } // Constructor

  TEUCHOS_UNIT_TEST(Repartition, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests build of the permutation matrix for repartitioning." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    if (numProcs != 4) {
      std::cout << "\nThis test must be run on 4 processors!\n" << std::endl;
      return;
    }

    const int nx = 5, ny = 3;

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",      nx);
    matrixList.set("ny",      ny);
    matrixList.set("keepBCs", false);

    // Describes the initial layout of matrix rows across processors.
    const GO     numGlobalElements = nx*ny; // 24
    const GO     indexBase     = 0;
    size_t numMyElements = 0;
    switch (myRank) {
      case 0:  numMyElements = 6; break;
      case 1:  numMyElements = 5; break;
      case 2:  numMyElements = 0; break;
      case 3:  numMyElements = 4; break;
    }

    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map, false);
    RCP<MultiVector>                  coords        = MultiVectorFactory::Build(map, 2);
    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength())
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Assign the partition that each unknown belongs to. In this case: part0 has 6, part1 has 4, part2 has 5
    const int numPartitions = 3;
    switch (myRank)  {
      case 0:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[3] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 2;
        partitionThisDofBelongsTo[4] = 2;
        partitionThisDofBelongsTo[5] = 2;
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[2] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[4] = 1;
        partitionThisDofBelongsTo[3] = 2;
        break;
      case 2:
        break;
      case 3:
        partitionThisDofBelongsTo[1] = 0;
        partitionThisDofBelongsTo[2] = 0;
        partitionThisDofBelongsTo[3] = 1;
        partitionThisDofBelongsTo[0] = 2;
        break;
    }
    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to Repartition. It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside Repartition.
    // Furthermore, that same instance must be supplied to MueLu::Repartition. Coordinates must be
    // provided, even though they are not used, as the Zoltan interface checks for them.
    Level level;
    level.SetLevelID(1);

    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);

    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());

    level.Set<RCP<Matrix> >     ("A",                    A);
    level.Set<RCP<MultiVector> >("Coordinates",          coords);
    level.Set<GO>               ("number of partitions", numPartitions);

    level.Request("Partition", zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition", decomposition, zoltan.get());

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    Teuchos::ParameterList paramList;
    paramList.set("repartition: start level",       1);
    paramList.set("repartition: min rows per proc", 1);
    paramList.set("repartition: max imbalance",     1.2);
    paramList.set("repartition: remap parts",       false);
    repart->SetParameterList(paramList);
    repart->SetFactory("Partition", zoltan);

    // Build
    level.Request("Importer", repart.get());
    repart->Build(level);

    // The local resulting vector should contain the partition number that this pid owns.
    // In this case, pid 0 owns partition 0
    //               pid 1 owns partition 1
    //               pid 2 owns partition 2
    //               pid 3 does not own a partition
    RCP<const Import> importer;
    level.Get("Importer", importer, repart.get());

    RCP<Xpetra::Vector<GO,LO,GO,NO> > result = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(importer->getTargetMap(), false);
    result->doImport(*decomposition, *importer, Xpetra::INSERT);
    Teuchos::ArrayRCP<GO> resultData;
    if (result->getLocalLength() > 0)
      resultData = result->getDataNonConst(0);

    // Local entries in the resulting vector should be equal to the partition number
    int pidFailed = -1;
    for (int i = 0; i < resultData.size(); i++) {
      if (resultData[i] != myRank) {
        pidFailed = myRank;
        break;
      }
    }
    resultData = Teuchos::null;

    // Test passes if all pid's return -1. If a pid hits an error, it sets
    // thisPidFailed equal to comm->getRank().  In this way, you can tell with --details=ALL
    // which pid failed.
    int whichPidFailed;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, pidFailed, Teuchos::outArg(whichPidFailed));
    TEST_EQUALITY(whichPidFailed, -1);

  } // Build

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement1)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << "The matrix is distributed across 4 processors, and there are 4 partitions." << std::endl;
    out << "Process 0 ends up owning a partition for which it originally has no unknowns." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    if (numProcs != 4) {
      std::cout << "\nThis test must be run on 4 processors!\n" << std::endl;
      return;
    }

    const int nx = 4, ny = 6;

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",      nx);
    matrixList.set("ny",      ny);
    matrixList.set("keepBCs", true); // keep Dirichlet rows

    // Describes the initial layout of matrix rows across processors.
    const GO     numGlobalElements = nx*ny; // 24
    const size_t numMyElements = 6;
    const GO     indexBase     = 0;

    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map, false);
    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Assign the partition that each unknown belongs to. In this case: part0 has 6, part1 has 6, part2 has 9, part3 has 6
    const int numPartitions = 4;
    switch (myRank)  {
      case 0:                                       //  nnz by row    nnz by partition
        partitionThisDofBelongsTo[2] = 0;           //      1
        partitionThisDofBelongsTo[3] = 0;           //      1             3
        partitionThisDofBelongsTo[4] = 0;           //      1
        partitionThisDofBelongsTo[5] = 1;           //      5             5
        partitionThisDofBelongsTo[0] = 2;           //      1
        partitionThisDofBelongsTo[1] = 2;           //      1             2
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;           //      5
        partitionThisDofBelongsTo[1] = 0;           //      1             6
        partitionThisDofBelongsTo[2] = 1;           //      1
        partitionThisDofBelongsTo[3] = 1;           //      5             6
        partitionThisDofBelongsTo[4] = 2;           //      5
        partitionThisDofBelongsTo[5] = 2;           //      1             6
        break;
      case 2:
        partitionThisDofBelongsTo[1] = 0;           //      5             5
        partitionThisDofBelongsTo[0] = 2;           //      1
        partitionThisDofBelongsTo[2] = 2;           //      5             8
        partitionThisDofBelongsTo[3] = 2;           //      1
        partitionThisDofBelongsTo[4] = 2;           //      1
        partitionThisDofBelongsTo[5] = 3;           //      5             5
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 1;           //      5
        partitionThisDofBelongsTo[1] = 1;           //      1             7
        partitionThisDofBelongsTo[2] = 1;           //      1
        partitionThisDofBelongsTo[3] = 2;           //      1             1
        partitionThisDofBelongsTo[4] = 3;           //      1
        partitionThisDofBelongsTo[5] = 3;           //      1             2
        break;
    }
    partitionThisDofBelongsTo = Teuchos::null;

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());

    bool keepProc0 = false;
    repart->DeterminePartitionPlacement(*A, *decomposition, numPartitions, keepProc0);

    Teuchos::ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
    switch (myRank) {
      case 0:
        TEST_EQUALITY(decompEntries[2], 1);
        TEST_EQUALITY(decompEntries[3], 1);
        TEST_EQUALITY(decompEntries[4], 1);
        TEST_EQUALITY(decompEntries[5], 3);
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        break;
      case 1:
        TEST_EQUALITY(decompEntries[0], 1);
        TEST_EQUALITY(decompEntries[1], 1);
        TEST_EQUALITY(decompEntries[2], 3);
        TEST_EQUALITY(decompEntries[3], 3);
        TEST_EQUALITY(decompEntries[4], 2);
        TEST_EQUALITY(decompEntries[5], 2);
        break;
      case 2:
        TEST_EQUALITY(decompEntries[1], 1);
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[2], 2);
        TEST_EQUALITY(decompEntries[3], 2);
        TEST_EQUALITY(decompEntries[4], 2);
        TEST_EQUALITY(decompEntries[5], 0);
        break;
      case 3:
        TEST_EQUALITY(decompEntries[0], 3);
        TEST_EQUALITY(decompEntries[1], 3);
        TEST_EQUALITY(decompEntries[2], 3);
        TEST_EQUALITY(decompEntries[3], 2);
        TEST_EQUALITY(decompEntries[4], 0);
        TEST_EQUALITY(decompEntries[5], 0);
        break;
    }
  } // DeterminePartitionPlacement1

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement2)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    if (numProcs != 4) {
      std::cout << "\nThis test must be run on 4 processors!\n" << std::endl;
      return;
    }

    const int nx = 4, ny = 6;

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",      nx);
    matrixList.set("ny",      ny);
    matrixList.set("keepBCs", false);

    // Describes the initial layout of matrix rows across processors.
    const GO     numGlobalElements = nx*ny; // 24
    const size_t numMyElements = 6;
    const GO     indexBase     = 0;

    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map, false);
    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Assign the partition that each unknown belongs to. In this case: part0 has 6, part1 has 6, parth2 has 9, part3 has 6
    const int numPartitions = 4;
    switch (myRank)  {
      case 0:                                       //  nnz by row    nnz by partition
        partitionThisDofBelongsTo[0] = 0;           //       3
        partitionThisDofBelongsTo[1] = 0;           //       4            7
        partitionThisDofBelongsTo[2] = 1;           //       4
        partitionThisDofBelongsTo[3] = 1;           //       3            7
        partitionThisDofBelongsTo[4] = 3;           //       4
        partitionThisDofBelongsTo[5] = 3;           //       5            9
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;           //       5
        partitionThisDofBelongsTo[1] = 0;           //       4            9
        partitionThisDofBelongsTo[2] = 1;           //       4
        partitionThisDofBelongsTo[3] = 1;           //       5            9
        partitionThisDofBelongsTo[4] = 2;           //       5
        partitionThisDofBelongsTo[5] = 2;           //       4            9
        break;
      case 2:
        partitionThisDofBelongsTo[0] = 0;           //       4
        partitionThisDofBelongsTo[1] = 0;           //       5            9
        partitionThisDofBelongsTo[2] = 2;           //       5
        partitionThisDofBelongsTo[3] = 2;           //       4            9
        partitionThisDofBelongsTo[4] = 3;           //       4
        partitionThisDofBelongsTo[5] = 3;           //       5            9
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 1;           //       5
        partitionThisDofBelongsTo[1] = 1;           //       4            9
        partitionThisDofBelongsTo[2] = 2;           //       3
        partitionThisDofBelongsTo[3] = 2;           //       4            7
        partitionThisDofBelongsTo[4] = 3;           //       4
        partitionThisDofBelongsTo[5] = 3;           //       3            7
        break;
    } //switch
    partitionThisDofBelongsTo = Teuchos::null;

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());

    bool keepProc0 = false;
    repart->DeterminePartitionPlacement(*A, *decomposition, numPartitions, keepProc0);

    Teuchos::ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
    switch (myRank)  {
      case 0:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        TEST_EQUALITY(decompEntries[2], 3);
        TEST_EQUALITY(decompEntries[3], 3);
        TEST_EQUALITY(decompEntries[4], 0);
        TEST_EQUALITY(decompEntries[5], 0);
        break;
      case 1:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        TEST_EQUALITY(decompEntries[2], 3);
        TEST_EQUALITY(decompEntries[3], 3);
        TEST_EQUALITY(decompEntries[4], 1);
        TEST_EQUALITY(decompEntries[5], 1);
        break;
      case 2:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        TEST_EQUALITY(decompEntries[2], 1);
        TEST_EQUALITY(decompEntries[3], 1);
        TEST_EQUALITY(decompEntries[4], 0);
        TEST_EQUALITY(decompEntries[5], 0);
        break;
      case 3:
        TEST_EQUALITY(decompEntries[0], 3);
        TEST_EQUALITY(decompEntries[1], 3);
        TEST_EQUALITY(decompEntries[2], 1);
        TEST_EQUALITY(decompEntries[3], 1);
        TEST_EQUALITY(decompEntries[4], 0);
        TEST_EQUALITY(decompEntries[5], 0);
        break;
    }
  } // DeterminePartitionPlacement2

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement3)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << "Matrix is initially placed on pids 0 and 1, but there are 4 partitions." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    if (numProcs != 4) {
      std::cout << "\nThis test must be run on 4 processors!\n" << std::endl;
      return;
    }

    const int nx = 5, ny = 3;

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",      nx);
    matrixList.set("ny",      ny);
    matrixList.set("keepBCs", false);

    // Describes the initial layout of matrix rows across processors.
    const GO     numGlobalElements = nx*ny; // 24
    const GO     indexBase     = 0;
    size_t numMyElements = 0;
    switch (myRank) {
      case 0:  numMyElements = 8; break;
      case 1:  numMyElements = 7; break;
    }

    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map, false);
    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength())
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Assign the partition that each unknown belongs to. In this case: part0 has 3, part1 has 4, parth2 has 3, part3 has 5
    const int numPartitions = 4;
    switch (myRank)  {
      case 0:                                       //  nnz by row    nnz by partition
        partitionThisDofBelongsTo[0] = 0;           //      3
        partitionThisDofBelongsTo[1] = 0;           //      4             7
        partitionThisDofBelongsTo[2] = 1;           //      4
        partitionThisDofBelongsTo[3] = 1;           //      4            11
        partitionThisDofBelongsTo[4] = 1;           //      3
        partitionThisDofBelongsTo[5] = 2;           //      4
        partitionThisDofBelongsTo[6] = 2;           //      5             9
        partitionThisDofBelongsTo[7] = 3;           //      5             5
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;           //      5             5
        partitionThisDofBelongsTo[1] = 1;           //      4             4
        partitionThisDofBelongsTo[2] = 2;           //      3             3
        partitionThisDofBelongsTo[3] = 3;           //      4
        partitionThisDofBelongsTo[4] = 3;           //      4            16
        partitionThisDofBelongsTo[5] = 3;           //      4
        partitionThisDofBelongsTo[6] = 3;           //      3
        break;
    }
    partitionThisDofBelongsTo = Teuchos::null;

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());

    bool keepProc0 = false;
    repart->DeterminePartitionPlacement(*A, *decomposition, numPartitions, keepProc0);

    Teuchos::ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
    switch (myRank)  {
      case 0:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        TEST_EQUALITY(decompEntries[2], 0);
        TEST_EQUALITY(decompEntries[3], 0);
        TEST_EQUALITY(decompEntries[4], 0);
        TEST_EQUALITY(decompEntries[5], 3);
        TEST_EQUALITY(decompEntries[6], 3);
        TEST_EQUALITY(decompEntries[7], 1);
        break;
      case 1:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 0);
        TEST_EQUALITY(decompEntries[2], 3);
        TEST_EQUALITY(decompEntries[3], 1);
        TEST_EQUALITY(decompEntries[4], 1);
        TEST_EQUALITY(decompEntries[5], 1);
        TEST_EQUALITY(decompEntries[6], 1);
        break;
    }
  } // DeterminePartitionPlacement3

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement4)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << "Matrix is distributed across all four processors, but there are only 3 partitions." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    if (numProcs != 4) {
      std::cout << "\nThis test must be run on 4 processors!\n" << std::endl;
      return;
    }

    const int nx = 5, ny = 3;

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",      nx);
    matrixList.set("ny",      ny);
    matrixList.set("keepBCs", false);

    // Describes the initial layout of matrix rows across processors.
    const GO     numGlobalElements = nx*ny; // 24
    const GO     indexBase     = 0;
    size_t numMyElements = 0;
    switch (myRank) {
      case 0:  numMyElements = 3; break;
      case 1:  numMyElements = 4; break;
      case 2:  numMyElements = 3; break;
      case 3:  numMyElements = 5; break;
    }

    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map, false);
    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength())
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Assign the partition that each unknown belongs to. In this case: part0 has 6, part1 has 6, parth2 has 3
    const int numPartitions = 3;
    switch (myRank)  {
      case 0:                                       //  nnz by row    nnz by partition
        partitionThisDofBelongsTo[0] = 0;           //      3             3
        partitionThisDofBelongsTo[1] = 1;           //      4             4
        partitionThisDofBelongsTo[2] = 2;           //      4             4
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;           //      4             4
        partitionThisDofBelongsTo[1] = 1;           //      3             3
        partitionThisDofBelongsTo[2] = 2;           //      4
        partitionThisDofBelongsTo[3] = 2;           //      5             9
        break;
      case 2:
        partitionThisDofBelongsTo[0] = 0;           //      5
        partitionThisDofBelongsTo[1] = 0;           //      5            10
        partitionThisDofBelongsTo[2] = 1;           //      4             4
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 0;           //      3
        partitionThisDofBelongsTo[1] = 0;           //      4             7
        partitionThisDofBelongsTo[2] = 1;           //      4
        partitionThisDofBelongsTo[3] = 1;           //      4
        partitionThisDofBelongsTo[4] = 1;           //      3            11
        break;
    }
    partitionThisDofBelongsTo = Teuchos::null;

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());

    bool keepProc0 = false;
    repart->DeterminePartitionPlacement(*A, *decomposition, numPartitions, keepProc0);

    Teuchos::ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
    switch (myRank)  {
      case 0:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 3);
        TEST_EQUALITY(decompEntries[2], 1);
        break;
      case 1:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 3);
        TEST_EQUALITY(decompEntries[2], 1);
        TEST_EQUALITY(decompEntries[3], 1);
        break;
      case 2:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        TEST_EQUALITY(decompEntries[2], 3);
        break;
      case 3:
        TEST_EQUALITY(decompEntries[0], 2);
        TEST_EQUALITY(decompEntries[1], 2);
        TEST_EQUALITY(decompEntries[2], 3);
        TEST_EQUALITY(decompEntries[3], 3);
        TEST_EQUALITY(decompEntries[4], 3);
        break;
    }
  } // DeterminePartitionPlacement4

  TEUCHOS_UNIT_TEST(Repartition, Correctness)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests application of the permutation matrix to matrix A." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    if (numProcs != 4) {
      std::cout << "\nThis test must be run on 4 processors!\n" << std::endl;
      return;
    }

    const int nx = 3, ny = 5;

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",      nx);
    matrixList.set("ny",      ny);
    matrixList.set("keepBCs", false);

    // Describes the initial layout of matrix rows across processors.
    const GO     numGlobalElements = nx*ny; // 24
    const GO     indexBase     = 0;
    size_t numMyElements = 0;
    switch (myRank) {
      case 0:  numMyElements = 6; break;
      case 1:  numMyElements = 5; break;
      case 2:  numMyElements = 0; break;
      case 3:  numMyElements = 4; break;
    }

    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map, false);
    RCP<MultiVector>                  coords        = MultiVectorFactory::Build(map, 2);
    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength())
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Assign the partition that each unknown belongs to. In this case: part0 has 6, part1 has 4, part2 has 5
    const int numPartitions = 3;
    switch (myRank)  {
      case 0:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 2;
        partitionThisDofBelongsTo[3] = 0;
        partitionThisDofBelongsTo[4] = 2;
        partitionThisDofBelongsTo[5] = 2;
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 0;
        partitionThisDofBelongsTo[3] = 2;
        partitionThisDofBelongsTo[4] = 1;
        break;
      case 2:
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 2;
        partitionThisDofBelongsTo[1] = 0;
        partitionThisDofBelongsTo[2] = 0;
        partitionThisDofBelongsTo[3] = 1;
        break;
    }
    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to Repartition.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // Repartition.  Furthermore, that same instance must be supplied to MueLu::Repartition.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    Level level;
    level.SetLevelID(1);

    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);

    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());

    level.Set<RCP<Matrix> >     ("A",                    A);
    level.Set<RCP<MultiVector> >("Coordinates",          coords);
    level.Set<GO>               ("number of partitions", numPartitions);

    level.Request("Partition", zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition", decomposition, zoltan.get());

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    Teuchos::ParameterList paramList;
    paramList.set("repartition: start level",       1);
    paramList.set("repartition: min rows per proc", 1);
    paramList.set("repartition: max imbalance",     1.2);
    paramList.set("repartition: remap parts",       false);
    repart->SetParameterList(paramList);
    repart->SetFactory("Partition", zoltan);

    // Build
    level.Request("Importer", repart.get());
    repart->Build(level);

    RCP<const Import> importer;
    level.Get("Importer", importer, repart.get());

    RCP<Matrix>        permutedA = MatrixFactory::Build(importer->getTargetMap(), A->getGlobalMaxNumRowEntries());
    permutedA->doImport(*A, *importer, Xpetra::INSERT);
    permutedA->fillComplete(A->getDomainMap(), importer->getTargetMap());

    // Calculate vectors y1 = Perm[A*v] and y2 = Perm[A]*v for random vector v. They should be identical (to machine precision)
    RCP<Vector> randVec = VectorFactory::Build(A->getDomainMap(), false);
    RCP<Vector> workVec = VectorFactory::Build(A->getRangeMap(),  false);
    randVec->randomize();

    RCP<Vector> P_Av = VectorFactory::Build(importer->getTargetMap(), false);
    RCP<Vector> PA_v = VectorFactory::Build(permutedA->getRangeMap(), false);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    permutedA->apply(*randVec, *PA_v, Teuchos::NO_TRANS, one, zero);

    A->apply(*randVec, *workVec, Teuchos::NO_TRANS, one, zero);
    P_Av->doImport(*workVec,*importer,Xpetra::INSERT);

    RCP<MultiVector> diff = VectorFactory::Build(permutedA->getRangeMap());
    diff->update(one, *P_Av, -one, *PA_v, zero);

    Teuchos::Array<STS::magnitudeType> norms(1);
    diff->norm2(norms);
    out << "||diff|| = " << norms[0] << std::endl;
    TEST_EQUALITY(norms[0] < 1e-14, true);

  } // Correctness

  GlobalOrdinal myrandom (GlobalOrdinal i) { return std::rand()%i;}

  TEUCHOS_UNIT_TEST(Repartition, CoordinateMap)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests that repartitioning is invariant to map specified in coordinates." << std::endl;
    out << std::endl;

    /*
       This test checks that MueLu successfully ignores the map of the coordinate MultiVector (MV).
       MueLu treats the coordinate data as if the MV is consistent with the linear system A.
       */

    // Create a matrix and coordinates.
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO nx = 20, ny = 20;

    // Describes the initial layout of matrix rows across processors.
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    galeriList.set("ny", ny);
    RCP<const Map> map = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);

    //build coordinates before expanding map (nodal coordinates, not dof-based)
    RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", map, galeriList);
    map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 2); //expand map for 2 DOFs per node

    galeriList.set("right boundary" , "Neumann");
    galeriList.set("bottom boundary", "Neumann");
    galeriList.set("top boundary"   , "Neumann");
    galeriList.set("front boundary" , "Neumann");
    galeriList.set("back boundary"  , "Neumann");
    galeriList.set("keepBCs",             false);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
    RCP<Matrix> A = Pr->BuildMatrix();
    A->SetFixedBlockSize(2);

    Utils::Write("A.mm", *A);
    comm->barrier();

    RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter("testCoordinates.xml", *comm));
    RCP<Hierarchy> H = mueLuFactory->CreateHierarchy();
    H->GetLevel(0)->Set("A", A);
    H->GetLevel(0)->Set("Coordinates", coordinates);
    mueLuFactory->SetupHierarchy(*H);
    double cplx1 = H->GetOperatorComplexity();

    //build a map that is a "randomly" permuted version of the correct
    //coordinate map.  This map will be used to build the "bad" coordinates.
    RCP<const Map> coordMap = coordinates->getMap();
    std::srand(Teuchos::as<unsigned int>(comm->getRank()*31415));
    Teuchos::ArrayView<const GO> correctLocalElts = coordMap->getNodeElementList();
    std::vector<GO> eltsToShuffle;
    for (size_t i=0; i < Teuchos::as<size_t>(correctLocalElts.size()); ++i)
      eltsToShuffle.push_back(correctLocalElts[i]);
    std::random_shuffle(eltsToShuffle.begin(),eltsToShuffle.end(),myrandom);
    Teuchos::Array<GO> eltList(eltsToShuffle);
    RCP<const Map> badMap = MapFactory::Build(TestHelpers::Parameters::getLib(), coordMap->getGlobalNumElements(), eltList(), coordMap->getIndexBase(), comm);

    Teuchos::Array<Teuchos::ArrayView<const Scalar> > coordVals;
    Teuchos::ArrayRCP<const Scalar> xcoords = coordinates->getData(0);
    Teuchos::ArrayRCP<const Scalar> ycoords = coordinates->getData(1);
    coordVals.push_back(xcoords());
    coordVals.push_back(ycoords());
    RCP<MultiVector> badCoordinates = MultiVectorFactory::Build(badMap, coordVals(), coordinates->getNumVectors());
    xcoords = Teuchos::null;
    ycoords = Teuchos::null;

    mueLuFactory = rcp(new ParameterListInterpreter("testCoordinates.xml", *comm));
    H = mueLuFactory->CreateHierarchy();
    H->GetLevel(0)->Set("A", A);
    H->GetLevel(0)->Set("Coordinates", badCoordinates);
    mueLuFactory->SetupHierarchy(*H);
    double cplx2 = H->GetOperatorComplexity();

    TEST_EQUALITY(cplx1, cplx2);

  } // CoordinateMap

} // namespace MueLuTests
