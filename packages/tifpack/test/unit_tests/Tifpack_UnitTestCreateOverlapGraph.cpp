// ***********************************************************************
// 
//      Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Tifpack_UnitTestTemplate.cpp

\brief Tifpack Unit testing template.

This file demonstrates how you create a unit test for template code.

*/


#include <Teuchos_ConfigDefs.hpp>
#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>

#include <Tifpack_UnitTestHelpers.hpp>
#include <Tifpack_CreateOverlapGraph.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(TifpackCreateOverlapGraph, OverlapGraphTest0, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  global_size_t num_rows_per_proc = 5;

//Create a Tpetra::CrsGraph:

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_tridiag_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  TEUCHOS_TEST_EQUALITY( crsgraph->getMap()->getNodeNumElements(), num_rows_per_proc, out, success)

  LocalOrdinal overlap_levels = 2;

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > overlapgraph =
    Tifpack::CreateOverlapGraph<LocalOrdinal,GlobalOrdinal,Node>(crsgraph, overlap_levels);

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
    TEUCHOS_TEST_EQUALITY(overlapgraph->getMap()->getNodeNumElements(), num_rows_per_proc, out, success)
  }
  else {
    if (myProc == 0 || myProc == numProcs-1) {
      TEUCHOS_TEST_EQUALITY(overlapgraph->getMap()->getNodeNumElements(), num_rows_per_proc+overlap_levels, out, success)
    }
    else {
      TEUCHOS_TEST_EQUALITY(overlapgraph->getMap()->getNodeNumElements(), num_rows_per_proc+overlap_levels*2, out, success)
    }
  }
}

#define UNIT_TEST_GROUP_ORDINAL(LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( TifpackCreateOverlapGraph, OverlapGraphTest0, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_ORDINAL(int, int)

}//namespace <anonymous>

