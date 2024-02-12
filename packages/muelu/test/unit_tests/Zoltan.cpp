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
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_ZoltanInterface.hpp>

#include <Galeri_XpetraUtils.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Zoltan, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<ZoltanInterface> zoltan         = rcp(new ZoltanInterface());
  TEST_EQUALITY(zoltan != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Zoltan, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << std::endl;
  out << "This tests that the partitioning produced by Zoltan is \"reasonable\" for a matrix" << std::endl;
  out << "that has a random number of nonzeros per row.  Good results have been precomputed" << std::endl;
  out << "for up to 5 processors.  The results are the number of nonzeros in the local matrix" << std::endl;
  out << "once the Zoltan repartitioning has been applied." << std::endl;
  out << "The results can be viewed in Paraview by enabling code guarded by the macro MUELU_VISUALIZE_REPARTITIONING" << std::endl;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  if (comm->getSize() > 5) {
    out << std::endl;
    out << "This test must be run on 1 to 5 processes." << std::endl;
    TEST_EQUALITY(true, true);
    return;
  }

  Level level;
  RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
  level.SetFactoryManager(factoryHandler);
  GO nx                                   = 7;
  GO ny                                   = nx;
  Xpetra::global_size_t numGlobalElements = nx * ny;
  size_t maxEntriesPerRow                 = 30;

  // Populate CrsMatrix with random number of entries (up to maxEntriesPerRow) per row.
  RCP<const Map> map                                       = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, 0 /*indexBase*/, comm);
  const size_t numMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();
  // memory on the fly.  While not super efficient, this
  // ensures that no zeros are being stored.  Thus, from
  // Zoltan's perspective the matrix is imbalanced.
  // Create a vector with random integer entries in [1,maxEntriesPerRow].
  ST::seedrandom(666 * comm->getRank());
  LO maxNnzPerRow                                    = 0;
  RCP<Xpetra::Vector<LO, LO, GO, NO> > entriesPerRow = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(map, false);
  Teuchos::ArrayRCP<LO> eprData                      = entriesPerRow->getDataNonConst(0);
  for (typename Teuchos::ArrayRCP<LO>::iterator i = eprData.begin(); i != eprData.end(); ++i) {
    *i = static_cast<LO>(std::floor(((Teuchos::ScalarTraits<double>::random() + 1) * 0.5 * maxEntriesPerRow) + 1));
    if (maxNnzPerRow < *i) {
      maxNnzPerRow = *i;
    }
  }
  RCP<Matrix> A = rcp(new CrsMatrixWrap(map, maxNnzPerRow));  // Force underlying linear algebra library to allocate more

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(-1);

  GO zero = Teuchos::OrdinalTraits<GO>::zero();
  GO one  = Teuchos::OrdinalTraits<GO>::one();

  Teuchos::Array<Scalar> vals(maxEntriesPerRow);
  Teuchos::Array<GO> cols(maxEntriesPerRow);
  for (size_t i = 0; i < numMyElements; ++i) {
    Teuchos::ArrayView<SC> av(&vals[0], eprData[i]);
    Teuchos::ArrayView<GO> iv(&cols[0], eprData[i]);
    // stick in ones for values
    for (LO j = 0; j < eprData[i]; ++j) vals[j] = ST::one();
    // figure out valid column indices
    GO start = std::max(myGlobalElements[i] - eprData[i] + one, zero);
    for (LO j = 0; j < eprData[i]; ++j) cols[j] = start + j;
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete();
  level.Set("A", A);

  // build coordinates
  RCP<const Map> rowMap = A->getRowMap();
  Teuchos::ParameterList list;
  list.set("nx", nx);
  list.set("ny", ny);
  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> double_multivector_type;
  RCP<double_multivector_type> XYZ = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, double_multivector_type>("2D", rowMap, list);
  level.Set("Coordinates", XYZ);

  int numPartitions = comm->getSize();
  level.Set("number of partitions", numPartitions);
  RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
  // zoltan->SetNumberOfPartitions(numPartitions);
  // zoltan->SetOutputLevel(0); //options are 0=none, 1=summary, 2=every pid prints
  level.Request("Partition", zoltan.get());
  zoltan->Build(level);

  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = level.Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >("Partition", zoltan.get());
  /* //TODO temporary code to have the trivial decomposition (no change)
     ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
     for (ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
   *i = comm->getRank();
   decompEntries=Teuchos::null;
   */
  // TODO end of temporary code

  // Create vector whose local length is the global number of partitions.
  // This vector will record the local number of nonzeros associated with each partition.
  Teuchos::Array<GO> parts(numPartitions);
  for (int i = 0; i < numPartitions; ++i) parts[i] = i;
  Teuchos::ArrayView<GO> partsView(&parts[0], numPartitions);
  RCP<const Map> partitionMap                        = MapFactory::Build(TestHelpers::Parameters::getLib(),
                                                                         Teuchos::OrdinalTraits<global_size_t>::invalid(), partsView,
                                                                         map->getIndexBase(), comm);
  RCP<Xpetra::Vector<LO, LO, GO, NO> > localPartsVec = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(partitionMap);

  // For the local rows in each partition, tally up the number of nonzeros.  This is what
  // Zoltan should be load-balancing.
  Teuchos::ArrayRCP<LO> lpvData          = localPartsVec->getDataNonConst(0);
  Teuchos::ArrayRCP<const GO> decompData = decomposition->getData(0);
  for (size_t i = 0; i < decomposition->getLocalLength(); ++i) {
    Teuchos::ArrayView<const LO> c;
    Teuchos::ArrayView<const SC> v;
    A->getLocalRowView(i, c, v);
    lpvData[decompData[i]] += v.size();
  }

  lpvData    = Teuchos::null;
  decompData = Teuchos::null;

  // localPartsVec->describe(*fos,Teuchos::VERB_EXTREME);

  // Send the local nnz tallies to pid 0, which can report the global sums.
  size_t mysize = 1;
  if (comm->getRank() == 0) mysize = numPartitions;
  RCP<const Map> globalTallyMap                       = MapFactory::Build(TestHelpers::Parameters::getLib(),
                                                                          Teuchos::OrdinalTraits<global_size_t>::invalid(),
                                                                          mysize,
                                                                          map->getIndexBase(),
                                                                          comm);
  RCP<Xpetra::Vector<LO, LO, GO, NO> > globalTallyVec = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(globalTallyMap);
  RCP<const Export> exporter                          = ExportFactory::Build(partitionMap, globalTallyMap);
  globalTallyVec->doExport(*localPartsVec, *exporter, Xpetra::ADD);

  ArrayRCP<GO> expectedResults(numPartitions);
  switch (comm->getSize()) {
    case 1:
      expectedResults[0] = 807;
      break;

    case 2:
      expectedResults[0] = 364;
      expectedResults[1] = 363;
      break;

    case 3:
      expectedResults[0] = 277;
      expectedResults[1] = 261;
      expectedResults[2] = 269;
      break;

    case 4:
      expectedResults[0] = 195;
      expectedResults[1] = 186;
      expectedResults[2] = 177;
      expectedResults[3] = 168;
      break;

    case 5:
      expectedResults[0] = 161;
      expectedResults[1] = 145;
      expectedResults[2] = 148;
      expectedResults[3] = 159;
      expectedResults[4] = 157;
      break;

    default:
      break;
  };

  // FIXME cool ... this next line causes a hang if locally the globalyTallyVec has no data.
  // FIXME I get around this by making mysize (above) 1 instead of 0. Is this a bug or feature
  // FIXME in getData?
  ArrayRCP<const LO> gtvData = globalTallyVec->getData(0);

#ifdef __linux__
  out << "Checking results..." << std::endl;
  for (int i = 0; i < numPartitions; ++i) {
    if (comm->getRank() == 0) TEST_EQUALITY(expectedResults[i], gtvData[i]);
  }
#endif

#ifdef MUELU_VISUALIZE_REPARTITIONING
  //
  // Now write everything to a comma-separate list that ParaView can grok
  //
  Teuchos::ArrayRCP<const double> X = XYZ->getData(0);
  Teuchos::ArrayRCP<const double> Y = XYZ->getData(1);
  Teuchos::ArrayRCP<const GO> D     = decomposition->getData(0);
  RCP<std::ofstream> outFile;
  std::string fileName = "zoltanResults.csv";

  // write header information
  if (comm->getRank() == 0) {
    outFile = rcp(new std::ofstream(fileName.c_str()));
    *outFile << "x coord, y coord, z coord, scalar" << std::endl;
  }
  comm->barrier();

  // append coordinates
  for (int j = 0; j < comm->getSize(); ++j) {
    int mypid = comm->getRank();
    if (mypid == j) {
      outFile = rcp(new std::ofstream(fileName.c_str(), std::ios::app));
      for (int i = 0; i < D.size(); ++i) {
        *outFile << X[i] << ", " << Y[i] << ", " << ST::zero() << ", " << D[i] << std::endl;
      }
    }
  }  // for (int i=0; i<comm->getSize(); ++i)

  out << std::endl;
  out << "You can view the Zoltan decomposition in ParaView 3.10.1 or later:" << std::endl;
  out << "   1) Load the data file " << fileName << "." << std::endl;
  out << "   2) Run the filter Filters/ Alphabetical/ Table To Points." << std::endl;
  out << "   3) Tell ParaView what columns are the X, Y and Z coordinates." << std::endl;
  out << "   4) Split screen horizontally (Icon, top right)." << std::endl;
  out << "   5) Click on the eyeball in the Pipeline Browser to see the points." << std::endl;
  out << "   6) Under the Display tab, you can color points by scalar value and resize them." << std::endl;
  out << std::endl;
  out << " To display row weights next to each point:" << std::endl;
  out << "   1) Click the \"Select Points Through\" button (2nd row) and select all points." << std::endl;
  out << "   2) Under View pull-down menu, choose the \"Selection Inspector\"." << std::endl;
  out << "   3) Under the Point Label, check the Visible box and set the Label Mode to \"row weight\"." << std::endl;
#endif

}  // Build

#ifdef DISABLED  // JG: FIXME: coordinates format
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Zoltan, Build3PDEs, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef Teuchos::ScalarTraits<Scalar> ST;

#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << std::endl;
  out << "This tests that the partitioning produced by Zoltan is \"reasonable\" for a matrix" << std::endl;
  out << "that has a random number of nonzeros per row and 3 DOFs per mesh point.  Good results have been precomputed" << std::endl;
  out << "for up to 5 processors.  The results are the number of nonzeros in the local matrix" << std::endl;
  out << "once the Zoltan repartitioning has been applied." << std::endl;
  out << "The results can be viewed in Paraview by enabling code guarded by the macro MUELU_VISUALIZE_REPARTITIONING" << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  if (comm->getSize() > 5) {
    out << std::endl;
    out << "This test must be run on 1 to 5 processes." << std::endl;
    TEST_EQUALITY(true, true);
    return;
  }

  Level level;
  RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
  level.SetFactoryManager(factoryHandler);
  int nx                  = 9;
  int ny                  = nx;
  int dofsPerNode         = 3;
  GO numGlobalElements    = nx * ny * dofsPerNode;
  size_t maxEntriesPerRow = 30;

  RCP<const Map> map;
  int numMyNodes = numGlobalElements / dofsPerNode;
  if (comm->getSize() > 1) {
    // In parallel, make sure that the dof's associated with a node all
    // reside on the same processor.
    int numNodes = numGlobalElements / dofsPerNode;
    TEUCHOS_TEST_FOR_EXCEPTION((numGlobalElements - numNodes * dofsPerNode) != 0, MueLu::Exceptions::RuntimeError,
                               "Number of matrix rows is not divisible by #dofs");
    int nproc = comm->getSize();
    if (comm->getRank() < nproc - 1)
      numMyNodes = numNodes / nproc;
    else
      numMyNodes = numNodes - (numNodes / nproc) * (nproc - 1);
    map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyNodes * dofsPerNode, 0 /*indexBase*/, comm);
  } else {
    map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, 0 /*indexBase*/, comm);
  }

  const size_t numMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();
  RCP<Matrix> A                                            = rcp(new CrsMatrixWrap(map, 1));  // Force underlying linear algebra library to allocate more
  // memory on the fly.  While not super efficient, this
  // ensures that no zeros are being stored.  Thus, from
  // Zoltan's perspective the matrix is imbalanced.
  // Populate CrsMatrix with random number of entries (up to maxEntriesPerRow) per row.
  // Create a vector with random integer entries in [1,maxEntriesPerRow].
  ST::seedrandom(666 * comm->getRank());
  RCP<Xpetra::Vector<LO, LO, GO, NO> > entriesPerRow = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(map, false);
  Teuchos::ArrayRCP<LO> eprData                      = entriesPerRow->getDataNonConst(0);
  for (Teuchos::ArrayRCP<LO>::iterator i = eprData.begin(); i != eprData.end(); ++i) {
    *i = (LO)(std::floor(((ST::random() + 1) * 0.5 * maxEntriesPerRow) + 1));
  }

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(-1);

  Teuchos::Array<Scalar> vals(maxEntriesPerRow);
  Teuchos::Array<GO> cols(maxEntriesPerRow);
  for (size_t i = 0; i < numMyElements; ++i) {
    Teuchos::ArrayView<SC> av(&vals[0], eprData[i]);
    Teuchos::ArrayView<GO> iv(&cols[0], eprData[i]);
    // stick in ones for values
    for (LO j = 0; j < eprData[i]; ++j) vals[j] = ST::one();
    // figure out valid column indices
    GO start = std::max(myGlobalElements[i] - eprData[i] + 1, 0);
    for (LO j = 0; j < eprData[i]; ++j) cols[j] = start + j;
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete();

  // Now treat the matrix as if it has 3 DOFs per node.
  A->SetFixedBlockSize(dofsPerNode);
  level.Set("A", A);

  // build coordinates
  Teuchos::ParameterList list;
  list.set("nx", nx);
  list.set("ny", ny);
  RCP<const Map> coalescedMap = MapFactory::Build<NO>(TestHelpers::Parameters::getLib(), numGlobalElements / dofsPerNode, numMyNodes, 0 /*indexBase*/, comm);
  RCP<MultiVector> XYZ        = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", coalescedMap, list);

  // XYZ are the "coalesce" coordinates as it has been generated for 1 DOF/node and we are using them for 3 DOFS/node
  // level.Set("Coordinates",XYZ); "Coordinates" == uncoalesce. "X,Y,ZCoordinates" == coalesce
  {
    RCP<MultiVector> coordinates = XYZ;

    // making a copy because I don't want to keep 'open' the Xpetra_MultiVector
    if (coordinates->getNumVectors() >= 1) {
      Teuchos::ArrayRCP<const SC> coord = coordinates->getData(0);
      Teuchos::ArrayRCP<SC> coordCpy(coord.size());
      for (int i = 0; i < coord.size(); i++) {
        coordCpy[i] = coord[i];
      }
      level.Set("XCoordinates", coordCpy);
      // std::cout << coordCpy << std::endl;
    }

    if (coordinates->getNumVectors() >= 2) {
      Teuchos::ArrayRCP<const SC> coord = coordinates->getData(1);
      Teuchos::ArrayRCP<SC> coordCpy(coord.size());
      for (int i = 0; i < coord.size(); i++) {
        coordCpy[i] = coord[i];
      }
      level.Set("YCoordinates", coordCpy);
    }

    /*if (coordinates->getNumVectors() >= 3) {
      Teuchos::ArrayRCP<const SC> coord = coordinates->getData(2);
      Teuchos::ArrayRCP<SC> coordCpy(coord.size());
      for(int i=0; i<coord.size(); i++) {
      coordCpy[i] = coord[i];
      }
      level.Set("ZCoordinates", coordCpy);
      }*/
  }

  // coalescedMap->describe(*fos,Teuchos::VERB_EXTREME);
  // sleep(1); comm->barrier();
  // XYZ->describe(*fos,Teuchos::VERB_EXTREME);

  int numPartitions = comm->getSize();
  level.Set("number of partitions", numPartitions);
  RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
  // zoltan->SetOutputLevel(0); //options are 0=none, 1=summary, 2=every pid prints
  level.Request("Partition", zoltan.get());
  zoltan->Build(level);

  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = level.Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >("Partition", zoltan.get());
  /* //temporary code to have the trivial decomposition (no change)
     ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
     for (ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
   *i = comm->getRank();
   decompEntries=Teuchos::null;
   */

  // Create vector whose local length is the global number of partitions.
  // This vector will record the local number of nonzeros associated with each partition.
  Teuchos::Array<GO> parts(numPartitions);
  for (int i = 0; i < numPartitions; ++i) parts[i] = i;
  Teuchos::ArrayView<GO> partsView(&parts[0], numPartitions);
  RCP<const Map> partitionMap                        = MapFactory::Build(TestHelpers::Parameters::getLib(),
                                                                         Teuchos::OrdinalTraits<global_size_t>::invalid(), partsView,
                                                                         map->getIndexBase(), comm);
  RCP<Xpetra::Vector<LO, LO, GO, NO> > localPartsVec = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(partitionMap);

  RCP<Xpetra::Vector<LO, LO, GO, NO> > nnzPerRow = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(A->getRowMap());
  Teuchos::ArrayRCP<GO> nnzData                  = nnzPerRow->getDataNonConst(0);
  // For the local rows in each partition, tally up the number of nonzeros.  This is what
  // Zoltan should be load-balancing.
  Teuchos::ArrayRCP<GO> lpvData          = localPartsVec->getDataNonConst(0);
  Teuchos::ArrayRCP<const GO> decompData = decomposition->getData(0);
  for (size_t i = 0; i < decomposition->getLocalLength(); ++i) {
    Teuchos::ArrayView<const LO> c;
    Teuchos::ArrayView<const SC> v;
    A->getLocalRowView(i, c, v);
    lpvData[decompData[i]] += v.size();
    nnzData[i] = v.size();
  }

  lpvData    = Teuchos::null;
  decompData = Teuchos::null;
  nnzData    = Teuchos::null;

  /*
     if (comm->getRank() == 0)
     std::cout << "nnz per row" << std::endl;
     nnzPerRow->describe(*fos,Teuchos::VERB_EXTREME);

     if (comm->getRank() == 0)
     std::cout << "Row-to-partition assignment (from Zoltan)" << std::endl;
     decomposition->describe(*fos,Teuchos::VERB_EXTREME);

     if (comm->getRank() == 0)
     std::cout << "#nonzeros per partition" << std::endl;
     localPartsVec->describe(*fos,Teuchos::VERB_EXTREME);
     */

  // Send the local nnz tallies to pid 0, which can report the global sums.
  size_t mysize = 1;
  if (comm->getRank() == 0) mysize = numPartitions;
  RCP<const Map> globalTallyMap                       = MapFactory::Build(TestHelpers::Parameters::getLib(),
                                                                          Teuchos::OrdinalTraits<global_size_t>::invalid(),
                                                                          mysize,
                                                                          map->getIndexBase(),
                                                                          comm);
  RCP<Xpetra::Vector<LO, LO, GO, NO> > globalTallyVec = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(globalTallyMap);
  RCP<const Export> exporter                          = ExportFactory::Build(partitionMap, globalTallyMap);
  globalTallyVec->doExport(*localPartsVec, *exporter, Xpetra::ADD);

  ArrayRCP<GO> expectedResults(numPartitions);
  switch (comm->getSize()) {
    case 1:
      expectedResults[0] = 3951;
      break;

    case 2:
      expectedResults[0] = 1955;
      expectedResults[1] = 1910;
      break;

    case 3:
      expectedResults[0] = 1326;
      expectedResults[1] = 1340;
      expectedResults[2] = 1321;
      break;

    case 4:
      expectedResults[0] = 950;
      expectedResults[1] = 922;
      expectedResults[2] = 908;
      expectedResults[3] = 936;
      break;

    case 5:
      expectedResults[0] = 774;
      expectedResults[1] = 735;
      expectedResults[2] = 726;
      expectedResults[3] = 771;
      expectedResults[4] = 759;
      break;

    default:
      break;
  };

  ArrayRCP<const LO> gtvData = globalTallyVec->getData(0);

#ifdef __linux__
  out << "Checking results..." << std::endl;
  for (int i = 0; i < numPartitions; ++i) {
    if (comm->getRank() == 0) TEST_EQUALITY(expectedResults[i], gtvData[i]);
  }
#endif

#ifdef MUELU_VISUALIZE_REPARTITIONING
  //
  // Now write everything to a comma-separate list that ParaView can grok
  //
  Teuchos::ArrayRCP<const double> X = XYZ->getData(0);
  Teuchos::ArrayRCP<const double> Y = XYZ->getData(1);
  Teuchos::ArrayRCP<const GO> D     = decomposition->getData(0);
  RCP<std::ofstream> outFile;
  std::string fileName = "zoltanResults.csv";

  // write header information
  if (comm->getRank() == 0) {
    outFile = rcp(new std::ofstream(fileName.c_str()));
    *outFile << "x coord, y coord, z coord, partition, row weight" << std::endl;
  }
  comm->barrier();

  // append coordinates
  nnzData = nnzPerRow->getDataNonConst(0);
  for (int j = 0; j < comm->getSize(); ++j) {
    int mypid = comm->getRank();
    if (mypid == j) {
      outFile       = rcp(new std::ofstream(fileName.c_str(), std::ios::app));
      int blockSize = A->GetFixedBlockSize();
      // Coordinates are for coalesced system, D is for uncoalesced
      for (int i = 0; i < D.size() / blockSize; ++i) {
        int nnz = 0;
        for (int k = 0; k < blockSize; ++k) nnz += nnzData[i * blockSize + k];
        *outFile << X[i] << ", " << Y[i] << ", " << ST::zero() << ", "
                 << D[i * blockSize] << ", " << nnz << std::endl;
      }
    }
  }  // for (int i=0; i<comm->getSize(); ++i)

  out << std::endl;
  out << "You can view the Zoltan decomposition in ParaView 3.10.1 or later:" << std::endl;
  out << "   1) Load the data file " << fileName << "." << std::endl;
  out << "   2) Run the filter Filters/ Alphabetical/ Table To Points." << std::endl;
  out << "   3) Tell ParaView what columns are the X, Y and Z coordinates." << std::endl;
  out << "   4) Split screen horizontally (Icon, top right)." << std::endl;
  out << "   5) Click on the eyeball in the Pipeline Browser to see the points." << std::endl;
  out << "   6) Under the Display tab, you can color points by scalar value and resize them." << std::endl;
  out << std::endl;
  out << " To display row weights next to each point:" << std::endl;
  out << "   1) Click the \"Select Points Through\" button (2nd row) and select all points." << std::endl;
  out << "   2) Under View pull-down menu, choose the \"Selection Inspector\"." << std::endl;
  out << "   3) Under the Point Label, check the Visible box and set the Label Mode to \"row weight\"." << std::endl;
#endif

}  // Build3PDEs
#endif  // TMP

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Zoltan, Constructor, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Zoltan, Build, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
