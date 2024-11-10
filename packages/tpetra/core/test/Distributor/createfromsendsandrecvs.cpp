// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

#include <Tpetra_ConfigDefs.hpp>
#include "Teuchos_UnitTestHarness.hpp"

#include <Tpetra_Distributor.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace { // (anonymous)

  using Teuchos::ArrayView;

  Teuchos::RCP<Tpetra::CrsGraph<> >
  getTpetraGraph (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Tpetra::Map<> map_type;
    typedef Tpetra::CrsGraph<> graph_type;
    typedef Tpetra::Map<>::local_ordinal_type LO;
    typedef Tpetra::Map<>::global_ordinal_type GO;
    typedef Tpetra::global_size_t GST;

    const LO lclNumRows = 101; // prime
    const GST gblNumRows = static_cast<GST> ( 101 * comm->getSize ());
    const GO indexBase = 0;
    const size_t numEntPerRow = 33;

    // A Map describes a distribution of data over MPI processes.
    // This "row Map" will describe the distribution of rows of the
    // sparse graph that we will create.
    RCP<const map_type> rowMap =
      rcp (new map_type (gblNumRows, static_cast<size_t> (lclNumRows),
                         indexBase, comm));
    const GO gblNumCols = static_cast<GO> (rowMap->getGlobalNumElements ());
    // Create the graph structure of the sparse matrix.
    RCP<graph_type> G =
      rcp (new graph_type (rowMap, numEntPerRow));
    // Fill in the sparse graph.
    Teuchos::Array<GO> gblColInds (numEntPerRow);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) { // for each of my rows
      const GO gblInd = rowMap->getGlobalElement (lclRow);
      // Just put some entries in the graph.  The actual column
      // indices don't matter so much, as long as they make the
      // resulting matrix square and don't go out of bounds.
      for (LO k = 0; k < static_cast<LO> (numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
      }
      G->insertGlobalIndices (gblInd, gblColInds ());
    }
    // Make the graph ready for use by (Block)CrsMatrix.
    G->fillComplete ();
    return G;
  }


  // Get a Tpetra::BlockCrsMatrix for use in benchmarks.
  // This method takes the result of getTpetraGraph() (above)
  Teuchos::RCP<Tpetra::BlockCrsMatrix<> >
  getTpetraBlockCrsMatrix (const Teuchos::RCP<const Tpetra::CrsGraph<> >& graph)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Tpetra::BlockCrsMatrix<> matrix_type;
    typedef matrix_type::impl_scalar_type SC;
    typedef Tpetra::Map<>::local_ordinal_type LO;
    //typedef Tpetra::Map<>::global_ordinal_type GO;
    typedef matrix_type::device_type device_type;
    typedef Kokkos::View<SC**, Kokkos::LayoutRight,
                         device_type>::HostMirror block_type;

    const auto& meshRowMap = * (graph->getRowMap ());
    // Contrary to expectations, asking for the graph's number of
    // columns, or asking the column Map for the number of entries,
    // won't give the correct number of columns in the graph.
    // const GO gblNumCols = graph->getDomainMap ()->getGlobalNumElements ();
    const LO lclNumRows = meshRowMap.getLocalNumElements ();
    const LO blkSize = 101;

    RCP<matrix_type> A = rcp (new matrix_type (*graph, blkSize));

    // Create a "prototype block" of values to use when filling the
    // block sparse matrix.  We don't care so much about the values;
    // we just want them not to be Inf or NaN, in case the processor
    // makes the unfortunate choice to handle arithmetic with those
    // via traps.
    block_type curBlk ("curBlk", blkSize, blkSize);
    for (LO j = 0; j < blkSize; ++j) {
      for (LO i = 0; i < blkSize; ++i) {
        curBlk(i,j) = 1.0;
      }
    }

    // Fill in the block sparse matrix.
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) { // for each of my rows
      typename matrix_type::local_inds_host_view_type lclColInds;
      graph->getLocalRowView (lclRow, lclColInds);

      // Put some entries in the matrix.
      for (LO k = 0; k < static_cast<LO> (lclColInds.size ()); ++k) {
        const LO lclColInd = lclColInds[k];
        const LO err =
          A->replaceLocalValues (lclRow, &lclColInd, curBlk.data (), 1);
        TEUCHOS_TEST_FOR_EXCEPTION(err != 1, std::logic_error, "Bug");
      }
    }

    return A;
  }

} // namespace (anonymous)


TEUCHOS_UNIT_TEST( Distributor, createfromsendsandrecvs)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::TimeMonitor;
  using std::endl;
  //typedef Tpetra::Vector<>::scalar_type SC;

  auto comm = Tpetra::getDefaultComm ();
  int my_proc = comm->getRank();
  //int nprocs = comm->getSize(); // unused

  // Set debug = true if you want immediate debug output to stderr.
  const bool debug = true;
  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    debug ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::rcpFromRef (out);
  Teuchos::FancyOStream& myOut = *outPtr;

  myOut << "Distributor createfromsendsandrecvs" << endl;
  Teuchos::OSTab tab1 (myOut);

  myOut << "Create CrsGraph, BlockCrsMatrix, and Vectors" << endl;
  auto G = getTpetraGraph (comm);
  auto A = getTpetraBlockCrsMatrix (G);
  Tpetra::Vector<> X (A->getDomainMap ());
  Tpetra::Vector<> Y (A->getRangeMap ());

  myOut << "Get the CrsGraph's Import object" << endl;
  RCP<const Tpetra::Import<> > importer = G->getImporter ();
  if (importer.is_null ()) {
    TEST_EQUALITY_CONST( comm->getSize (), 1 );
    myOut << "The CrsGraph's Import object is null";
    if (success) {
      myOut << ".  This is to be expected when the communicator only has 1 "
        "process.  We'll say this test succeeded and be done with it." << endl;
    }
    else {
      myOut << ", but the communicator has " << comm->getSize () << " != 1 "
        "processes.  That means we didn't construct the test graph correctly.  "
        "It makes no sense to continue this test beyond this point." << endl;
    }
    return;
  }
  auto dist = importer->getDistributor();

  myOut << "Build up arrays to construct equivalent Distributor" << endl;

  const ArrayView<const int> procF = dist.getProcsFrom();
  const ArrayView<const int> procT = dist.getProcsTo();
  const ArrayView<const size_t> lenF = dist.getLengthsFrom();
  const ArrayView<const size_t> lenT = dist.getLengthsTo();
  // This section takes the consolidated procF and procT with the length and re-builds
  // the un-consolidated lists of processors from and to that
  // This is needed because in Tpetra::constructExpert, the unconsolidated procsFrom and ProcsTo
  // will be used.

  Teuchos::Array<int> nuF;
  Teuchos::Array<int> nuT;
  int sumLenF=0;
  for ( ArrayView<const size_t>::iterator b = lenF.begin(); b!=lenF.end(); ++b)
    sumLenF+=(*b);
  int sumLenT=0;
  for ( ArrayView<const size_t>::iterator b = lenT.begin(); b!=lenT.end(); ++b)
    sumLenT+=(*b);
  nuF.resize(sumLenF);
  nuT.resize(sumLenT);

  size_t p=0;
  for ( size_t j = 0; j<(size_t)procF.size(); ++j) {
    size_t lend = p+lenF[j];
    for (size_t i = p ; i < lend ; ++i)
      nuF[i]=procF[j];
    p+=lenF[j];
  }
  p=0;
  for ( size_t j = 0; j<(size_t) procT.size(); ++j) {
    size_t lend = p+lenT[j];
    for (size_t i = p ; i < lend ; ++i)
      nuT[i]=procT[j];
    p+=lenT[j];
  }

  myOut << "Create a new Distributor using createFromSendsAndRecvs" << endl;

  Tpetra::Distributor newdist(comm);
  TEST_NOTHROW( newdist.createFromSendsAndRecvs(nuT,nuF) );
  {
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess) );
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      myOut << "Test FAILED on some process; giving up early" << endl;
    }
  }

  myOut << "Test resulting Distributor" << endl;

  if(dist.getNumReceives()!=newdist.getNumReceives()) {
    myOut << "ProcID "<<my_proc <<" getNumReceives does not match " << endl;
    success = false;
  }
  if(dist.getNumSends()!=newdist.getNumSends()) {
    myOut << "ProcID "<<my_proc <<" getNumSends does not match " << endl;
    success = false;
  }
  if(dist.hasSelfMessage()!=newdist.hasSelfMessage()) {
    myOut << "ProcID "<<my_proc <<" hasSelfMessage does not match " << endl;
    success = false;
  }
  if(dist.getMaxSendLength()!=newdist.getMaxSendLength()) {
    myOut << "ProcID "<<my_proc <<" getMaxSendLength does not match " << endl;
    success = false;
  }
  if(dist.getTotalReceiveLength()!=newdist.getTotalReceiveLength()) {
    myOut << "ProcID "<<my_proc <<" getTotalReceiveLength does not match " << endl;
    success = false;
  }
  ArrayView<const int> a = dist.getProcsFrom();
  ArrayView<const int> b = newdist.getProcsFrom();
  if(a.size()!=b.size()) {
    myOut << "ProcID "<<my_proc <<" getProcsFrom size does not match " << endl;

    success = false;
  }
  for(unsigned ui = 0;ui<a.size();++ui)
    if(a[ui]!=b[ui]) {
      myOut << "ProcID "<<my_proc <<" old getProcsFrom"<<a<<endl;
      myOut << "ProcID "<<my_proc <<" new getProcsFrom"<<b<<endl;
      success = false;
      break;
    }
  ArrayView<const int> at = dist.getProcsTo();
  ArrayView<const int> bt = newdist.getProcsTo();
  if(at.size()!=bt.size()) {
    myOut << "ProcID "<<my_proc <<" getProcsTo size does not match " << endl;
    success = false;
  }
  for(unsigned ui = 0;ui<at.size();++ui)
    if(at[ui]!=bt[ui]) {
      myOut << "ProcID "<<my_proc <<" getProcsTo old "<<at<<endl;
      myOut << "ProcID "<<my_proc <<" getProcsTo new "<<bt<<endl;
      success = false;
    }
  ArrayView<const long unsigned int> c = dist.getLengthsFrom();
  ArrayView<const long unsigned int> d = newdist.getLengthsFrom();
  if(c.size()!=d.size()) {
    myOut << "ProcID "<<my_proc <<" getLengthsFrom does not match "<<b<<endl;
    success = false;
  }
  for(unsigned ui = 0;ui<c.size();++ui)
    if(c[ui]!=d[ui]) {
      myOut << "ProcID "<<my_proc <<" lengthfrom old "<<c<<endl;
      myOut << "ProcID "<<my_proc <<" lengthsfrom new "<<d<<endl;
      success = false;
      break;
    }
  ArrayView<const long unsigned int> ct = dist.getLengthsTo();
  ArrayView<const long unsigned int> dt = newdist.getLengthsTo();
  if(ct.size()!=dt.size()) {
    myOut << "ProcID "<<my_proc <<" getLengthsTo size does not match " << endl;
    success = false;
  }
  for(unsigned ui = 0;ui<ct.size();++ui)
    if(ct[ui]!=dt[ui]) {
      myOut << "ProcID "<<my_proc <<" lengthTo old "<<ct<<endl;
      myOut << "ProcID "<<my_proc <<" lengthsTo new "<<dt<<endl;
      success = false;
      break;
    }

  if(newdist.howInitialized()!=Tpetra::Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS)
    {
      myOut << "ProcID "<<my_proc <<"howInitialized() from distributor initialized with createFromSendsAndRecvs is incorrect" << endl;
      success = false;
    }

  int globalSuccess_int = -1;
  Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );

  // return EXIT_SUCCESS;
}

