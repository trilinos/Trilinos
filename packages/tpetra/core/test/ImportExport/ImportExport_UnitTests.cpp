// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include <iterator>
#include <sstream>

#include <Tpetra_Distributor.hpp>
#include "Teuchos_FancyOStream.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Tpetra_Import_Util.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::global_size_t;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::tuple;
  using Teuchos::Range1D;
  using Tpetra::Map;
  using Tpetra::Import;
  using Tpetra::Export;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
  using Tpetra::ADD;
  using std::ostream_iterator;
  using std::endl;

  using Tpetra::createContigMapWithNode;

  // bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, basic, LO, GO, NT ) {
    const Tpetra::global_size_t INVALID =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    // create Maps
    RCP<const Map<LO, GO, NT> > source =
      createContigMapWithNode<LO, GO, NT> (INVALID, 10, comm);
    RCP<const Map<LO, GO, NT> > target =
      createContigMapWithNode<LO, GO, NT> (INVALID,  5, comm);
    // create Import object
    RCP<const Import<LO, GO, NT> > importer =
      Tpetra::createImport<LO, GO, NT> (source, target);

    auto same = importer->getNumSameIDs();
    auto permute = importer->getNumPermuteIDs();
    auto remote = importer->getNumRemoteIDs();
    auto sum = same + permute + remote;
    auto expectedSum = target->getLocalNumElements();
    TEST_EQUALITY( sum, expectedSum );

    bool isvalid=Tpetra::Import_Util::checkImportValidity(*importer);
    TEST_EQUALITY(isvalid,true);
    TEST_EQUALITY(importer->isLocallyFitted(), source->isLocallyFitted(*target));
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ImportExport, GetNeighborsForward, Scalar, LO, GO, Node )
  {
    // import with the importer to duplicate
    // export with the exporter to add and reduce
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const size_t numLocal  = 1;
    const size_t numVectors = 5;
    // my neighbors: myImageID-1, me, myImageID+1
    Array<GO> neighbors;
    if (myImageID != 0) neighbors.push_back(myImageID-1);
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) neighbors.push_back(myImageID+1);
    // two maps: one has one entries per process, the other is the 1-D neighbors
    RCP<const Map<LO,GO,Node> >
      smap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm),
      tmap = rcp(new Map<LO,GO,Node>(INVALID,neighbors(),0,comm) );
    for (size_t tnum=0; tnum < 2; ++tnum) {
      RCP<MV> mvMine, mvWithNeighbors;
      // for tnum=0, these are contiguously allocated multivectors
      // for tnum=1, these are non-contiguous views of multivectors
      if (tnum == 0) {
        mvMine = rcp(new MV(smap,numVectors));
        mvWithNeighbors = rcp(new MV(tmap,numVectors));
      }
      else {
        MV mineParent(smap,2+numVectors);
        MV neigParent(tmap,2+numVectors);
        TEUCHOS_TEST_FOR_EXCEPTION(numVectors != 5, std::logic_error, "Test assumption broken.");
        mvMine = mineParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
        mvWithNeighbors = neigParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
      }
      // mvMine = [myImageID  myImageID+numImages ... myImageID+4*numImages]
      for (size_t j=0; j<numVectors; ++j) {
        mvMine->replaceLocalValue(0,j,static_cast<Scalar>(myImageID + j*numImages));
      }
      // create Import from smap to tmap, Export from tmap to smap, test them
      RCP<const Import<LO,GO,Node> > importer =
        Tpetra::createImport<LO,GO,Node>(smap,tmap);
      RCP<const Export<LO,GO,Node> > exporter =
        Tpetra::createExport<LO,GO,Node>(tmap,smap);
      bool local_success = true;
      // importer testing
      TEST_EQUALITY_CONST( importer->getSourceMap() == smap, true );
      TEST_EQUALITY_CONST( importer->getTargetMap() == tmap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( importer->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( importer->getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->isLocallyFitted(), tmap->isLocallyFitted(*smap));
      // exporter testing
      TEST_EQUALITY_CONST( exporter->getSourceMap() == tmap, true );
      TEST_EQUALITY_CONST( exporter->getTargetMap() == smap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( exporter->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( exporter->isLocallyFitted(), tmap->isLocallyFitted(*smap));
      // import neighbors, test their proper arrival
      //                   [ 0    n     2n    3n    4n ]
      // mvWithNeighbors = [...  ....  ....  ....  ....]
      //                   [n-1  2n-1  3n-1  4n-1  5n-1]
      mvWithNeighbors->doImport(*mvMine,*importer,REPLACE);
      if (myImageID == 0) {
        for (size_t j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)); // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(j*numImages)+ST::one()); // neighbor
        }
      }
      else if (myImageID == numImages-1) {
        for (size_t j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
        }
      }
      else {
        for (size_t j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),2,static_cast<Scalar>(myImageID+j*numImages)+ST::one()); // neighbor
        }
      }

      // export values, test
      mvMine->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      mvMine->doExport(*mvWithNeighbors,*exporter,ADD);
      if (myImageID == 0 || myImageID == numImages-1) {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and one neighbor: double original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(2.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      else {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and two neighbors: triple original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(3.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      success &= local_success;
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ImportExport, GetNeighborsBackward, Scalar, LO, GO, Node )
  {
    // import with the exporter to duplicate
    // export with the importer to add and reduce
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const size_t numLocal = 1;
    const size_t numVectors = 5;
    // my neighbors: myImageID-1, me, myImageID+1
    Array<GO> neighbors;
    if (myImageID != 0) neighbors.push_back(myImageID-1);
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) neighbors.push_back(myImageID+1);
    // two maps: one has one entries per node, the other is the 1-D neighbors
    auto smap = createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm);
    auto tmap = rcp (new Map<LO, GO, Node> (INVALID, neighbors (), 0, comm));
    for (size_t tnum=0; tnum < 2; ++tnum) {
      RCP<MV> mvMine, mvWithNeighbors;
      // for tnum=0, these are contiguously allocated multivectors
      // for tnum=1, these are non-contiguous views of multivectors
      if (tnum == 0) {
        mvMine = rcp(new MV(smap,numVectors));
        mvWithNeighbors = rcp(new MV(tmap,numVectors));
      }
      else {
        MV mineParent(smap,2+numVectors),
           neigParent(tmap,2+numVectors);
        TEUCHOS_TEST_FOR_EXCEPTION(numVectors != 5, std::logic_error, "Test assumption broken.");
        mvMine = mineParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
        mvWithNeighbors = neigParent.subViewNonConst(tuple<size_t>(0,6,3,4,5));
      }
      // mvMine = [myImageID  myImageID+numImages ... myImageID+4*numImages]
      for (size_t j=0; j<numVectors; ++j) {
        mvMine->replaceLocalValue(0,j,static_cast<Scalar>(myImageID + j*numImages));
      }
      // create Import from smap to tmap, Export from tmap to smap, test them
      auto importer = Tpetra::createImport<LO, GO, Node> (smap, tmap);
      auto exporter = Tpetra::createExport<LO, GO, Node> (tmap, smap);
      bool local_success = true;
      // importer testing
      TEST_EQUALITY_CONST( importer->getSourceMap() == smap, true );
      TEST_EQUALITY_CONST( importer->getTargetMap() == tmap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( importer->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( importer->getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->isLocallyFitted(), tmap->isLocallyFitted(*smap));
      // exporter testing
      TEST_EQUALITY_CONST( exporter->getSourceMap() == tmap, true );
      TEST_EQUALITY_CONST( exporter->getTargetMap() == smap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( exporter->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( exporter->isLocallyFitted(), tmap->isLocallyFitted(*smap));
      // import neighbors, test their proper arrival
      //                   [ 0    n     2n    3n    4n ]
      // mvWithNeighbors = [...  ....  ....  ....  ....]
      //                   [n-1  2n-1  3n-1  4n-1  5n-1]
      mvWithNeighbors->doImport(*mvMine,*exporter,REPLACE);
      if (myImageID == 0) {
        for (size_t j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)); // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(j*numImages)+ST::one()); // neighbor
        }
      }
      else if (myImageID == numImages-1) {
        for (size_t j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
        }
      }
      else {
        for (size_t j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),0,static_cast<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),1,static_cast<Scalar>(myImageID+j*numImages));           // me
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors->getData(j),2,static_cast<Scalar>(myImageID+j*numImages)+ST::one()); // neighbor
        }
      }
      // export values, test
      mvMine->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      mvMine->doExport(*mvWithNeighbors,*importer,ADD);
      if (myImageID == 0 || myImageID == numImages-1) {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and one neighbor: double original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(2.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      else {
        for (size_t j=0; j<numVectors; ++j) {
          // contribution from me and two neighbors: triple original value
          TEST_EQUALITY(mvMine->getData(j)[0],static_cast<Scalar>(3.0)*static_cast<Scalar>(myImageID+j*numImages));
        }
      }
      success &= local_success;
    }
    //
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, AbsMax, LO, GO, Node )
  {
    using Tpetra::createContigMapWithNode;
    using Tpetra::createNonContigMapWithNode;

    // test ABSMAX CombineMode
    //
    // The test includes both local and remote entries, to exercise
    // both copying and permuting, and unpacking and combining.
    typedef Tpetra::Vector<>::scalar_type SC;
    typedef Tpetra::Vector<SC,LO,GO,Node> Vec;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 2) return;

    auto smap = createContigMapWithNode<LO, GO, Node> (INVALID, 1, comm);
    const GO myOnlyGID = smap->getGlobalElement (0);
    auto dmap = createNonContigMapWithNode<LO, GO, Node> (tuple<GO> (myOnlyGID, (myOnlyGID+1) % numImages), comm);
    RCP<Vec> srcVec = rcp (new Vec (smap));
    srcVec->putScalar (-1.0);
    RCP<Vec> dstVec = rcp (new Vec (dmap));
    dstVec->putScalar (-3.0);
    // first item of dstVec is local (w.r.t. srcVec), while the second is remote
    // ergo, during the import:
    // - the first will be over-written (by 1.0) from the source, while
    // - the second will be "combined", i.e., abs(max(1.0,3.0)) = 3.0 from the dest
    auto importer = Tpetra::createImport<LO, GO, Node> (smap, dmap);
    TEST_EQUALITY( importer->isLocallyFitted(), dmap->isLocallyFitted(*smap));
    dstVec->doImport (*srcVec,*importer,Tpetra::ABSMAX);
    TEST_COMPARE_ARRAYS( tuple<SC>(-1.0,3.0), dstVec->get1dView() )
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


 TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, ExportReverse, LO, GO, Node )
  {
    // This test reproduces an issue seen in Github Issue #114.
    // As of time of checkin, this test will fail on CUDA but pass on other platforms.
    // This is intentional
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    typedef Tpetra::Map<LO,GO,Node> Tpetra_Map;
    typedef Tpetra::Import<LO,GO,Node> Tpetra_Import;
    typedef Tpetra::Vector<int, LO, GO,Node> IntVector;

    int NumProcs = comm->getSize();
    int MyPID    = comm->getRank();

    // This problem only works on 4 procs
    if(NumProcs!=4) {TEST_EQUALITY(true,true);return;}

    // Problem setup
    int num_per_proc;
    if(MyPID==0) num_per_proc=7;
    else num_per_proc=6;

    GO from_gids_p0[7] = {0,1,2,3,4,5,6};
    GO to_gids_p0[7]   = {0,4,8,12,16,20,24};

    GO from_gids_p1[6] = {7,8,9,10,11,12};
    GO to_gids_p1[6]   = {1,5,9,13,17,21};

    GO from_gids_p2[6] = {13,14,15,16,17,18};
    GO to_gids_p2[6]   = {2,6,10,14,18,22};

    GO from_gids_p3[6] = {19,20,21,22,23,24};
    GO to_gids_p3[6]   = {3,7,11,15,19,23};

    // Correctness check array
    int who_owns[25];
    for(int i=0; i<7; i++)
      who_owns[to_gids_p0[i]] = 0;
    for(int i=0; i<6; i++) {
      who_owns[to_gids_p1[i]] = 1;
      who_owns[to_gids_p2[i]] = 2;
      who_owns[to_gids_p3[i]] = 3;
    }

    GO *from_ptr, *to_ptr;
    if(MyPID==0)      {from_ptr=&from_gids_p0[0]; to_ptr=&to_gids_p0[0];}
    else if(MyPID==1) {from_ptr=&from_gids_p1[0]; to_ptr=&to_gids_p1[0];}
    else if(MyPID==2) {from_ptr=&from_gids_p2[0]; to_ptr=&to_gids_p2[0];}
    else if(MyPID==3) {from_ptr=&from_gids_p3[0]; to_ptr=&to_gids_p3[0];}
    else exit(-1);

    Teuchos::ArrayView<GO> myfromgids(from_ptr,num_per_proc);
    Teuchos::ArrayView<GO> mytogids(to_ptr,num_per_proc);

    // FromMap (from.getRowMap() from Zoltan2)
    RCP<Tpetra_Map> FromMap = rcp(new Tpetra_Map(INVALID,myfromgids,0,comm));

    // ToMap (tmap from Zoltan2)
    RCP<Tpetra_Map> ToMap = rcp(new Tpetra_Map(INVALID,mytogids,0,comm));

    // Importer
    Tpetra_Import Importer(FromMap,ToMap);

    TEST_EQUALITY( Importer.isLocallyFitted(), ToMap->isLocallyFitted(*FromMap));

    // Duplicating what Zoltan2/Tpetra Does
    IntVector FromVector(FromMap);
    IntVector ToVector(ToMap);
    ToVector.putScalar(MyPID);
    FromVector.putScalar(-666);

    FromVector.doExport(ToVector,Importer,Tpetra::REPLACE);

    Teuchos::ArrayRCP<const int> f_rcp = FromVector.getData();
    Teuchos::ArrayView<const int> f_view = f_rcp();
    Teuchos::ArrayRCP<const int> t_rcp = ToVector.getData();
    Teuchos::ArrayView<const int> t_view = t_rcp();

    // Check the "FromAnswer" answer against who_owns
    bool all_is_well=true;
    for(size_t i=0; i<FromMap->getLocalNumElements(); i++) {
      if(f_view[i] != who_owns[FromMap->getGlobalElement(i)]) {
        std::cerr<<"["<<MyPID<<"] ERROR: Ownership of GID"<<FromMap->getGlobalElement(i)<<" is incorrect!"<<std::endl;
        all_is_well=false;
      }
    }
    TEST_EQUALITY(all_is_well,true);


    bool isvalid=Tpetra::Import_Util::checkImportValidity(Importer);
    if(!isvalid) {
      std::ostringstream oss;
      Importer.print(oss);

      std::cout<<oss.str()<<std::endl;
    }

    TEST_EQUALITY(isvalid,true);

  }

  //
  // INSTANTIATIONS
  //


#define UNIT_TEST_3( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, basic, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, AbsMax, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, ExportReverse, LO, GO, NT)

  #define UNIT_TEST_4( SCALAR, LO, GO, NODE )  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ImportExport, GetNeighborsForward,  SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ImportExport, GetNeighborsBackward, SCALAR, LO, GO, NODE )


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_3 )

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_4 )

} // namespace (anonymous)


