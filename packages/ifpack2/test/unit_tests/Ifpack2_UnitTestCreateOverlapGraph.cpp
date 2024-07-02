// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestCreateOverlapGraph.cpp

\brief Ifpack2 Unit test.

This file unit-tests the createOverlapGraph function.

*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_CreateOverlapGraph.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2CreateOverlapGraph, OverlapGraphTest0, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  global_size_t num_rows_per_proc = 5;

//Create a Tpetra::CrsGraph:

  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_tridiag_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  TEST_EQUALITY( crsgraph->getMap()->getLocalNumElements(), num_rows_per_proc)

  LocalOrdinal overlap_levels = 2;

  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > overlapgraph =
    Ifpack2::createOverlapGraph<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >(crsgraph, overlap_levels);

  const int numProcs = overlapgraph->getMap()->getComm()->getSize();
  const int myProc   = overlapgraph->getMap()->getComm()->getRank();

  //Now test how many local rows there are in the overlapped-graph.
  //For a tri-diagonal input-graph with 5 local rows on each proc:
  //'end' procs (procs 0 and numProcs-1) should have an extra row in the
  //overlapped graph for each level of overlap. Other procs should
  //have an extra 2 rows for each level of overlap.

  //Special case: if numProcs==1 then overlapgraph should have the same
  //number of rows as the input-graph.

  if (numProcs == 1) {
    TEST_EQUALITY(overlapgraph->getMap()->getLocalNumElements(), num_rows_per_proc)
  }
  else {
    if (myProc == 0 || myProc == numProcs-1) {
      TEST_EQUALITY(overlapgraph->getMap()->getLocalNumElements(), num_rows_per_proc+overlap_levels)
    }
    else {
      TEST_EQUALITY(overlapgraph->getMap()->getLocalNumElements(), num_rows_per_proc+overlap_levels*2)
    }
  }
}

#define UNIT_TEST_GROUP_LO_GO(LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2CreateOverlapGraph, OverlapGraphTest0, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of LocalOrdinal (LO)
// and GlobalOrdinal (GO) types.

IFPACK2_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

} // namespace (anonymous)

