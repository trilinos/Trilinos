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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
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

#include "Galeri_XpetraUtils.hpp"

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(Repartition, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    TEST_EQUALITY(repart != Teuchos::null, true);

  } //Constructor

  // -----------------------------------------------------------------------

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement1)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << "The matrix is distributed across 4 processors, and there are 4 partitions." << std::endl;
    out << "Process 0 ends up owning a partition for which it originally has no unknowns." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=4;
    int ny=6;
    GO numGlobalElements = nx*ny; //24

    // Describes the initial layout of matrix rows across processors.
    size_t numMyElements=0;
    switch(mypid) {
       case 0:
         numMyElements = 6;
         break;
       case 1:
         numMyElements = 6;
         break;
       case 2:
         numMyElements = 6;
         break;
       case 3:
         numMyElements = 6;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);
    matrixList.set("keepBCs",true); //keeps Dirichlet rows

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D",map,matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();
    level.Set<RCP<Matrix> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions", 4);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 6 unknowns
       partition 1 has 6 unknowns
       partition 2 has 9 unknowns
       partition 3 has 3 unknowns
    */

    switch (mypid)  {
      case 0:                                       //  nnz by row    nnz by partition
        partitionThisDofBelongsTo[2] = 0;           //           1
        partitionThisDofBelongsTo[3] = 0;           //           1          3
        partitionThisDofBelongsTo[4] = 0;           //           1
        partitionThisDofBelongsTo[5] = 1;           //           5          5
        partitionThisDofBelongsTo[0] = 2;           //           1
        partitionThisDofBelongsTo[1] = 2;           //           1          2
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;           //           5
        partitionThisDofBelongsTo[1] = 0;           //           1          6
        partitionThisDofBelongsTo[2] = 1;           //           1
        partitionThisDofBelongsTo[3] = 1;           //           5          6
        partitionThisDofBelongsTo[4] = 2;           //           5
        partitionThisDofBelongsTo[5] = 2;           //           1          6
        break;
      case 2:
        partitionThisDofBelongsTo[1] = 0;           //           5          5
        partitionThisDofBelongsTo[0] = 2;           //           1
        partitionThisDofBelongsTo[2] = 2;           //           5          8
        partitionThisDofBelongsTo[3] = 2;           //           1
        partitionThisDofBelongsTo[4] = 2;           //           1
        partitionThisDofBelongsTo[5] = 3;           //           5          5
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 1;           //           5
        partitionThisDofBelongsTo[1] = 1;           //           1          7
        partitionThisDofBelongsTo[2] = 1;           //           1
        partitionThisDofBelongsTo[3] = 2;           //           1          1
        partitionThisDofBelongsTo[4] = 3;           //           1
        partitionThisDofBelongsTo[5] = 3;           //           1          2
        break;
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to RepartitionFactory.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // RepartitionFactory.  Furthermore, that same instance must be supplied to MueLu::RepartitionFactory.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    RCP<MultiVector> coords = MultiVectorFactory::Build(map,2);
    level.Set<RCP<MultiVector> >("Coordinates",coords);
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
    level.Request("Partition",zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition, zoltan.get());
    level.SetLevelID(2); //partitioning by default won't happen unless level >= 1 // TODO: we should just need to change an option of the factory to do this.

    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    Teuchos::ParameterList paramList;
    paramList.set("startLevel", 1);
    paramList.set("minRowsPerProcessor", 1000);
    paramList.set("nonzeroImbalance", 1.2);
    paramList.set("diffusiveHeuristic", 10);
    paramList.set< RCP<const FactoryBase> >("number of partitions", Teuchos::null); // use user-defined #partitions
    repart->SetParameterList(paramList);
    repart->SetFactory("Partition", zoltan);

    GO myPartitionNumber;
    Array<int> partitionOwners;
    level.Request(*repart);
    repart->DeterminePartitionPlacement(level,myPartitionNumber,partitionOwners);
    level.Release(*repart);

    TEST_EQUALITY(partitionOwners[0],1);
    TEST_EQUALITY(partitionOwners[1],3);
    TEST_EQUALITY(partitionOwners[2],2);
    TEST_EQUALITY(partitionOwners[3],0);
  } //DeterminePartitionPlacement1

  // -----------------------------------------------------------------------

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement2)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=4;
    int ny=6;
    GO numGlobalElements = nx*ny; //24

    // Describes the initial layout of matrix rows across processors.
    size_t numMyElements=0;
    switch(mypid) {
       case 0:
         numMyElements = 6;
         break;
       case 1:
         numMyElements = 6;
         break;
       case 2:
         numMyElements = 6;
         break;
       case 3:
         numMyElements = 6;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D",map,matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();
    level.Set<RCP<Matrix> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions",4);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 6 unknowns
       partition 1 has 6 unknowns
       partition 2 has 8 unknowns
       partition 3 has 4 unknowns
    */

    switch (mypid)  {
      case 0:                                    //nnz in this row       nnz by partition
        partitionThisDofBelongsTo[0] = 0;        //       3
        partitionThisDofBelongsTo[1] = 0;        //       4                     7
        partitionThisDofBelongsTo[2] = 1;        //       4
        partitionThisDofBelongsTo[3] = 1;        //       3                     7
        partitionThisDofBelongsTo[4] = 3;        //       4
        partitionThisDofBelongsTo[5] = 3;        //       5                     9
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;        //       5
        partitionThisDofBelongsTo[1] = 0;        //       4                     9
        partitionThisDofBelongsTo[2] = 1;        //       4
        partitionThisDofBelongsTo[3] = 1;        //       5                     9
        partitionThisDofBelongsTo[4] = 2;        //       5
        partitionThisDofBelongsTo[5] = 2;        //       4                     9
        break;
      case 2:
        partitionThisDofBelongsTo[0] = 0;        //       4
        partitionThisDofBelongsTo[1] = 0;        //       5                     9
        partitionThisDofBelongsTo[2] = 2;        //       5
        partitionThisDofBelongsTo[3] = 2;        //       4                     9
        partitionThisDofBelongsTo[4] = 3;        //       4
        partitionThisDofBelongsTo[5] = 3;        //       5                     9
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 1;        //       5
        partitionThisDofBelongsTo[1] = 1;        //       4                     9
        partitionThisDofBelongsTo[2] = 2;        //       3
        partitionThisDofBelongsTo[3] = 2;        //       4                     7
        partitionThisDofBelongsTo[4] = 3;        //       4
        partitionThisDofBelongsTo[5] = 3;        //       3                     7
        break;
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to Repartition.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // Repartition.  Furthermore, that same instance must be supplied to MueLu::Repartition.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    RCP<MultiVector> coords = MultiVectorFactory::Build(map,2);
    level.Set<RCP<MultiVector> >("Coordinates",coords);
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
    level.Request("Partition",zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition, zoltan.get());
    level.SetLevelID(2); //partitioning by default won't happen unless level >= 1
    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());

    Teuchos::ParameterList paramList;
    paramList.set("startLevel", 1);
    paramList.set("minRowsPerProcessor", 1000);
    paramList.set("nonzeroImbalance", 1.2);
    paramList.set("diffusiveHeuristic", 10);
    paramList.set< RCP<const FactoryBase> >("number of partitions", Teuchos::null); // use user-defined #partitions
    repart->SetParameterList(paramList);

    repart->SetFactory("Partition", zoltan);

    GO myPartitionNumber;
    Array<int> partitionOwners;
    level.Request(*repart);
    repart->DeterminePartitionPlacement(level,myPartitionNumber,partitionOwners);
    level.Release(*repart);

    TEST_EQUALITY(partitionOwners[0],2);
    TEST_EQUALITY(partitionOwners[1],3);
    TEST_EQUALITY(partitionOwners[2],1);
    TEST_EQUALITY(partitionOwners[3],0);
  } //DeterminePartitionPlacement2

  // -----------------------------------------------------------------------

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement3)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << "Matrix is initially placed on pids 0 and 1, but there are 4 partitions." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=5;
    int ny=3;
    GO numGlobalElements = nx*ny; //15

    // Describes the initial layout of matrix rows across processors.
    size_t numMyElements=0;
    switch(mypid) {
       case 0:
         numMyElements = 8;
         break;
       case 1:
         numMyElements = 7;
         break;
       case 2:
         numMyElements = 0;
         break;
       case 3:
         numMyElements = 0;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D",map,matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();
    level.Set<RCP<Matrix> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions",4);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 3 unknowns
       partition 1 has 4 unknowns
       partition 2 has 3 unknowns
       partition 3 has 5 unknowns
    */

    switch (mypid)  {
      case 0:                                //nnz in this row        nnz by partition
        partitionThisDofBelongsTo[0] = 0;    //       3
        partitionThisDofBelongsTo[1] = 0;    //       4                     7
        partitionThisDofBelongsTo[2] = 1;    //       4
        partitionThisDofBelongsTo[3] = 1;    //       4                    11
        partitionThisDofBelongsTo[4] = 1;    //       3
        partitionThisDofBelongsTo[5] = 2;    //       4
        partitionThisDofBelongsTo[6] = 2;    //       5                     9
        partitionThisDofBelongsTo[7] = 3;    //       5                     5
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;    //       5                     5
        partitionThisDofBelongsTo[1] = 1;    //       4                     4
        partitionThisDofBelongsTo[2] = 2;    //       3                     3
        partitionThisDofBelongsTo[3] = 3;    //       4
        partitionThisDofBelongsTo[4] = 3;    //       4                    16
        partitionThisDofBelongsTo[5] = 3;    //       4
        partitionThisDofBelongsTo[6] = 3;    //       3
        break;
      case 2:
        break;
      case 3:
        break;
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to Repartition.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // Repartition.  Furthermore, that same instance must be supplied to MueLu::Repartition.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    RCP<MultiVector> coords = MultiVectorFactory::Build(map,2);
    level.Set<RCP<MultiVector> >("Coordinates",coords);
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
    level.Request("Partition",zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition, zoltan.get());
    level.SetLevelID(2); //partitioning by default won't happen unless level >= 1
    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    Teuchos::ParameterList paramList;
    paramList.set("startLevel", 1);
    paramList.set("minRowsPerProcessor", 1000);
    paramList.set("nonzeroImbalance", 1.2);
    paramList.set("diffusiveHeuristic", 10);
    paramList.set< RCP<const FactoryBase> >("number of partitions", Teuchos::null); // use user-defined #partitions
    repart->SetParameterList(paramList);
    repart->SetFactory("Partition", zoltan);

    GO myPartitionNumber;
    Array<int> partitionOwners;
    level.Request(*repart);
    repart->DeterminePartitionPlacement(level,myPartitionNumber,partitionOwners);
    level.Release(*repart);

    TEST_EQUALITY(partitionOwners[0],2);
    TEST_EQUALITY(partitionOwners[1],0);
    TEST_EQUALITY(partitionOwners[2],3);
    TEST_EQUALITY(partitionOwners[3],1);
  } //DeterminePartitionPlacement3

  // -----------------------------------------------------------------------

  TEUCHOS_UNIT_TEST(Repartition, DeterminePartitionPlacement4)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
    out << "Matrix is distributed across all four processors, but there are only 3 partitions." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=5;
    int ny=3;
    GO numGlobalElements = nx*ny; //15

    // Describes the initial layout of matrix rows across processors.
    size_t numMyElements=0;
    switch(mypid) {
       case 0:
         numMyElements = 3;
         break;
       case 1:
         numMyElements = 4;
         break;
       case 2:
         numMyElements = 3;
         break;
       case 3:
         numMyElements = 5;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D",map,matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();
    level.Set<RCP<Matrix> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions",3);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 3 unknowns
       partition 1 has 4 unknowns
       partition 2 has 3 unknowns
       partition 3 has 5 unknowns
    */

    switch (mypid)  {
      case 0:                              // nnz in this row      nnz in this partition
        partitionThisDofBelongsTo[0] = 0;  //       3                        3
        partitionThisDofBelongsTo[1] = 1;  //       4                        4
        partitionThisDofBelongsTo[2] = 2;  //       4                        4
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;  //       4                        4
        partitionThisDofBelongsTo[1] = 1;  //       3                        3
        partitionThisDofBelongsTo[2] = 2;  //       4
        partitionThisDofBelongsTo[3] = 2;  //       5                        9
        break;
      case 2:
        partitionThisDofBelongsTo[0] = 0;  //       5
        partitionThisDofBelongsTo[1] = 0;  //       5                       10
        partitionThisDofBelongsTo[2] = 1;  //       4                        4
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 0;  //       3
        partitionThisDofBelongsTo[1] = 0;  //       4                        7
        partitionThisDofBelongsTo[2] = 1;  //       4
        partitionThisDofBelongsTo[3] = 1;  //       4
        partitionThisDofBelongsTo[4] = 1;  //       3                       11
        break;
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that given via Level to Repartition.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // Repartition.  Furthermore, that same instance must be supplied to MueLu::Repartition.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    RCP<MultiVector> coords = MultiVectorFactory::Build(map,2);
    level.Set<RCP<MultiVector> >("Coordinates",coords);
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
    level.Set<GO>("number of partitions",3);
    level.Request("Partition",zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition, zoltan.get());
    level.SetLevelID(2); //partitioning by default won't happen unless level >= 1
    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    Teuchos::ParameterList paramList;
    paramList.set("startLevel", 1);
    paramList.set("minRowsPerProcessor", 1000);
    paramList.set("nonzeroImbalance", 1.2);
    paramList.set("diffusiveHeuristic", 10);
    paramList.set< RCP<const FactoryBase> >("number of partitions", Teuchos::null); // use user-defined #partitions
    repart->SetParameterList(paramList);
    repart->SetFactory("Partition", zoltan);

    GO myPartitionNumber;
    Array<int> partitionOwners;
    level.Request(*repart);
    repart->DeterminePartitionPlacement(level,myPartitionNumber,partitionOwners);
    level.Release(*repart);

    TEST_EQUALITY(partitionOwners[0],2);
    TEST_EQUALITY(partitionOwners[1],3);
    TEST_EQUALITY(partitionOwners[2],1);
  } //DeterminePartitionPlacement4

  // -----------------------------------------------------------------------

  TEUCHOS_UNIT_TEST(Repartition, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests build of the permutation matrix for repartitioning." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=3;
    int ny=5;
    GO numGlobalElements = nx*ny; //15

    // Describes the initial layout of unknowns across processors. Note that PID 2 initially has 0 unknowns.
    size_t numMyElements=0;
    switch(mypid) {
       case 0:
         numMyElements = 6;
         break;
       case 1:
         numMyElements = 5;
         break;
       case 2:
         numMyElements = 0;
         break;
       case 3:
         numMyElements = 4;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D",map,matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();
    level.Set<RCP<Matrix> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions",3);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 6 unknowns
       partition 1 has 4 unknowns
       partition 2 has 5 unknowns
       there is no partition 3
    */

    switch (mypid)  {
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
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to Repartition.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // Repartition.  Furthermore, that same instance must be supplied to MueLu::Repartition.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    RCP<MultiVector> coords = MultiVectorFactory::Build(map,2);
    level.Set<RCP<MultiVector> >("Coordinates",coords);
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
    level.Set<GO>("number of partitions",3);
    level.Request("Partition",zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition, zoltan.get());
    level.SetLevelID(2); //partitioning by default won't happen unless level >= 1
    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    Teuchos::ParameterList paramList;
    paramList.set("startLevel", 1);
    paramList.set("minRowsPerProcessor", 1000);
    paramList.set("nonzeroImbalance", 1.2);
    paramList.set("diffusiveHeuristic", 10);
    paramList.set< RCP<const FactoryBase> >("number of partitions", Teuchos::null); // use user-defined #partitions
    repart->SetParameterList(paramList);
    repart->SetFactory("Partition", zoltan);
    level.Request("Importer",repart.get());  // request permutation matrix

    repart->Build(level);

    // Test permutation matrix by multiplying against vector that describes partition numbering.
    // The local resulting vector should contain the partition number that this pid owns.
    // In this case, pid 0 owns partition 2
    //               pid 1 owns partition 1
    //               pid 2 does not own a partition
    //               pid 3 owns partition 0
    RCP<const Import> importer;
    level.Get("Importer",importer,repart.get());
    RCP<Vector> result = VectorFactory::Build(importer->getTargetMap(),false);
    result->doImport(*decompositionAsScalar,*importer,Xpetra::INSERT);
    int thisPidFailed=-1;
    // local entries in the resulting vector should be equal to the partition number
    Teuchos::ArrayRCP<SC> resultData;
    if (result->getLocalLength() > 0)
      resultData = result->getDataNonConst(0);
    for (int i=0; i<resultData.size(); ++i) {
      switch(mypid) {
        case 0:
           if (resultData[i] != 2)
             thisPidFailed = mypid;
           break;
        case 1:
           if (resultData[i] != 1)
             thisPidFailed = mypid;
           break;
        case 2: //pid 2 should have no data after permutation
           thisPidFailed = mypid;
           break;
        case 3:
           if (resultData[i] != 0)
             thisPidFailed = mypid;
           break;
      } //switch
    }
    resultData = Teuchos::null;

    // Test passes if all pid's return -1. If a pid hits an error, it sets
    // thisPidFailed equal to comm->getRank().  In this way, you can tell with --details=ALL
    // which pid failed.
    int whichPidFailed;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,thisPidFailed,Teuchos::outArg(whichPidFailed));
    TEST_EQUALITY(whichPidFailed,-1);

  } //Build

  TEUCHOS_UNIT_TEST(Repartition, Correctness)
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;

    out << "version: " << MueLu::Version() << std::endl;
    out << "Tests application of the permutation matrix to matrix A." << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=3;
    int ny=5;
    GO numGlobalElements = nx*ny; //15

    // Describes the initial layout of unknowns across processors. Note that PID 2 initially has 0 unknowns.
    size_t numMyElements=0;
    switch(mypid) {
       case 0:
         numMyElements = 6;
         break;
       case 1:
         numMyElements = 5;
         break;
       case 2:
         numMyElements = 0;
         break;
       case 3:
         numMyElements = 4;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D",map,matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();
    level.Set<RCP<Matrix> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions",3);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 6 unknowns
       partition 1 has 4 unknowns
       partition 2 has 5 unknowns
       there is no partition 3
    */

    switch (mypid)  {
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
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    // This test uses a made up partitioning  that is given via Level to Repartition.  It must be
    // associated with an instance of the ZoltanInterface so that it can be found inside
    // Repartition.  Furthermore, that same instance must be supplied to MueLu::Repartition.
    // Coordinates must be provided, even though they are not used, as the Zoltan interface checks for them.
    RCP<MultiVector> coords = MultiVectorFactory::Build(map,2);
    level.Set<RCP<MultiVector> >("Coordinates",coords);
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
    level.Set<GO>("number of partitions",3);
    level.Request("Partition",zoltan.get());
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition, zoltan.get());
    level.SetLevelID(2); //partitioning by default won't happen unless level >= 1
    RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
    repart->SetFactory("Partition", zoltan);
    //RCP<RepartitionFactory> repart = rcp(new RepartitionFactory(zoltan, Teuchos::null, 1000, 1.2, 1, 10/*useDiffusiveHeuristic*/, -1));
    //level.Request("Permutation",repart.get());  // request permutation matrix
    level.Request("Importer",repart.get());  // request permutation matrix

    repart->Build(level);

    RCP<const Import> importer;
    level.Get("Importer",importer,repart.get());
    //TODO this next bit needs to be put into Xpetra::Matrix or Xpetra::Utils
    RCP<Matrix> PermTimesA = MatrixFactory::Build(importer->getTargetMap(), Op->getGlobalMaxNumRowEntries());
    RCP<CrsMatrixWrap> crsOp = rcp_dynamic_cast<CrsMatrixWrap>(PermTimesA);
    RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
    RCP<CrsMatrixWrap> origOp = rcp_dynamic_cast<CrsMatrixWrap>(Op);
    RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
    crsMtx->doImport(*origMtx, *importer,Xpetra::INSERT);
    crsMtx = Teuchos::null;
    PermTimesA->fillComplete(Op->getDomainMap(),importer->getTargetMap());
    //TODO end of this next bit

    // Calculate vectors y1=Perm*(A*v) and y2 = (Perm*A)*v for random vector v. They should be identical.
    RCP<Vector> randomVec = VectorFactory::Build(Op->getDomainMap(),false);
    randomVec->randomize();
    RCP<Vector> workVec = VectorFactory::Build(Op->getRangeMap(),false);
    RCP<Vector> P_Av = VectorFactory::Build(importer->getTargetMap(),false);
    RCP<Vector> PA_v = VectorFactory::Build(PermTimesA->getRangeMap(),false);
    PermTimesA->apply(*randomVec,*PA_v,Teuchos::NO_TRANS,1,0);
    Op->apply(*randomVec,*workVec,Teuchos::NO_TRANS,1,0);
    P_Av->doImport(*workVec,*importer,Xpetra::INSERT);

    RCP<MultiVector> diff = VectorFactory::Build(PermTimesA->getRangeMap());
    //diff = P_Av + (-1.0)*(PA_v) + 0*diff
    diff->update(1.0,*P_Av,-1.0,*PA_v,0.0);

    Teuchos::Array<ST::magnitudeType> norms(1);
    diff->norm2(norms);
    out << "||diff|| = " << norms[0] << std::endl;
    TEST_EQUALITY(norms[0]<1e-15, true);

  } //Correctness

}//namespace MueLuTests
