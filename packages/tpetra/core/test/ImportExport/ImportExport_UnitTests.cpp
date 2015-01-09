/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>
#include "Teuchos_UnitTestHarness.hpp"
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include <iterator>

namespace {
  using Tpetra::TestingUtilities::getNode;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::DefaultPlatform;
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

  using Tpetra::createContigMap;
  using Tpetra::createContigMapWithNode;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ImportExport, basic, Ordinal, Node ) {
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    RCP<Node> node = getNode<Node>();
    // create Maps
    RCP<const Map<Ordinal,Ordinal,Node> > source = createContigMapWithNode<Ordinal,Ordinal>(INVALID,10,comm,node),
                                          target = createContigMapWithNode<Ordinal,Ordinal>(INVALID, 5,comm,node);
    // create Import object
    RCP<const Import<Ordinal,Ordinal,Node> > importer = Tpetra::createImport<Ordinal,Ordinal,Node>(source, target);

    Ordinal same = importer->getNumSameIDs();
    Ordinal permute = importer->getNumPermuteIDs();
    Ordinal remote = importer->getNumRemoteIDs();
    Ordinal sum = same + permute + remote;
    Ordinal expectedSum = target->getNodeNumElements();
    TEST_EQUALITY( sum, expectedSum );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, GetNeighborsForward, Ordinal, Scalar, Node )
  {
    // import with the importer to duplicate
    // export with the exporter to add and reduce
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    RCP<Node> node = getNode<Node>();
    const int numImages = comm->getSize(),
              myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const size_t numLocal  = 1,
                 numVectors = 5;
    // my neighbors: myImageID-1, me, myImageID+1
    Array<Ordinal> neighbors;
    if (myImageID != 0) neighbors.push_back(myImageID-1);
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) neighbors.push_back(myImageID+1);
    // two maps: one has one entries per node, the other is the 1-D neighbors
    RCP<const Map<Ordinal,Ordinal,Node> >
      smap = createContigMapWithNode<Ordinal,Ordinal>(INVALID,numLocal,comm,node),
      tmap = rcp(new Map<Ordinal,Ordinal,Node>(INVALID,neighbors(),0,comm,node) );
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
      RCP<const Import<Ordinal,Ordinal,Node> > importer =
        Tpetra::createImport<Ordinal,Ordinal,Node>(smap,tmap);
      RCP<const Export<Ordinal,Ordinal,Node> > exporter =
        Tpetra::createExport<Ordinal,Ordinal,Node>(tmap,smap);
      bool local_success = true;
      // importer testing
      TEST_EQUALITY_CONST( importer->getSourceMap() == smap, true );
      TEST_EQUALITY_CONST( importer->getTargetMap() == tmap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( importer->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( importer->getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      // exporter testing
      TEST_EQUALITY_CONST( exporter->getSourceMap() == tmap, true );
      TEST_EQUALITY_CONST( exporter->getTargetMap() == smap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( exporter->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
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


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, GetNeighborsBackward, Ordinal, Scalar, Node )
  {
    // import with the exporter to duplicate
    // export with the importer to add and reduce
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    RCP<Node> node = getNode<Node>();
    const int numImages = comm->getSize(),
              myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const size_t numLocal = 1,
               numVectors = 5;
    // my neighbors: myImageID-1, me, myImageID+1
    Array<Ordinal> neighbors;
    if (myImageID != 0) neighbors.push_back(myImageID-1);
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) neighbors.push_back(myImageID+1);
    // two maps: one has one entries per node, the other is the 1-D neighbors
    RCP<const Map<Ordinal,Ordinal,Node> > smap = createContigMapWithNode<Ordinal,Ordinal,Node>(INVALID,numLocal,comm,node),
                                          tmap = rcp(new Map<Ordinal,Ordinal,Node>(INVALID,neighbors(),0,comm,node) );
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
      RCP<const Import<Ordinal,Ordinal,Node> > importer = Tpetra::createImport<Ordinal,Ordinal,Node>(smap,tmap);
      RCP<const Export<Ordinal,Ordinal,Node> > exporter = Tpetra::createExport<Ordinal,Ordinal,Node>(tmap,smap);
      bool local_success = true;
      // importer testing
      TEST_EQUALITY_CONST( importer->getSourceMap() == smap, true );
      TEST_EQUALITY_CONST( importer->getTargetMap() == tmap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( importer->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
      TEST_EQUALITY( importer->getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      TEST_EQUALITY( importer->getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
      // exporter testing
      TEST_EQUALITY_CONST( exporter->getSourceMap() == tmap, true );
      TEST_EQUALITY_CONST( exporter->getTargetMap() == smap, true );
      TEST_EQUALITY( importer->getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
      TEST_EQUALITY( exporter->getNumPermuteIDs(), static_cast<size_t>(myImageID == 0 ? 0 : 1) );
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


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ImportExport, AbsMax, Node )
  {
    // test ABSMAX CombineMode
    // test with local and remote entries, as copyAndPermute() and unpackAndCombine() both need to be tested
    typedef Tpetra::Vector<double,int,int,Node> Vec;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    RCP<Node> node = getNode<Node>();
    if (numImages < 2) return;
    // create a Map
    RCP<const Map<int,int,Node> > smap =
      Tpetra::createContigMapWithNode<int,int>(INVALID,1,comm,node);
    const int myOnlyGID = smap->getGlobalElement(0);
    RCP<const Map<int,int,Node> > dmap =
      Tpetra::createNonContigMapWithNode<int,int>(tuple<int>(myOnlyGID, (myOnlyGID+1)%numImages), comm, node);
    RCP<Vec> srcVec = Tpetra::createVector<double>(smap);
    srcVec->putScalar(-1.0);
    RCP<Vec> dstVec = Tpetra::createVector<double>(dmap);
    dstVec->putScalar(-3.0);
    // first item of dstVec is local (w.r.t. srcVec), while the second is remote
    // ergo, during the import:
    // - the first will be over-written (by 1.0) from the source, while
    // - the second will be "combined", i.e., abs(max(1.0,3.0)) = 3.0 from the dest
    RCP<const Tpetra::Import<int,int,Node> > importer =
      Tpetra::createImport<int, int, Node> (smap, dmap);
    dstVec->doImport(*srcVec,*importer,Tpetra::ABSMAX);
    TEST_COMPARE_ARRAYS( tuple<double>(-1.0,3.0), dstVec->get1dView() )
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  //
  // INSTANTIATIONS
  //

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT_NODE(ORDINAL, NODE)  \
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE(ORDINAL, ComplexFloat, NODE)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE_NODE(ORDINAL, NODE) \
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE(ORDINAL, ComplexDouble, NODE)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT_NODE(ORDINAL, NODE)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE_NODE(ORDINAL, NODE)
#endif

#define UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE )         \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, GetNeighborsForward,  ORDINAL, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, GetNeighborsBackward, ORDINAL, SCALAR, NODE )

#define UNIT_TEST_GROUP_ORDINAL_NODE( ORDINAL, NODE )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ImportExport, basic, ORDINAL, NODE ) \
  UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( ORDINAL, double, NODE)

#define UNIT_TEST_GROUP_NODE( NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ImportExport, AbsMax, NODE ) \
  UNIT_TEST_GROUP_ORDINAL_NODE( int, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
  TPETRA_INSTANTIATE_N( UNIT_TEST_GROUP_NODE )
#else
  TPETRA_INSTANTIATE_N_NOGPU( UNIT_TEST_GROUP_NODE )
#endif
}


