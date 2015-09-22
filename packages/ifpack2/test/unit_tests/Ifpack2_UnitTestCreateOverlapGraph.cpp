/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


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

  TEST_EQUALITY( crsgraph->getMap()->getNodeNumElements(), num_rows_per_proc)

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
    TEST_EQUALITY(overlapgraph->getMap()->getNodeNumElements(), num_rows_per_proc)
  }
  else {
    if (myProc == 0 || myProc == numProcs-1) {
      TEST_EQUALITY(overlapgraph->getMap()->getNodeNumElements(), num_rows_per_proc+overlap_levels)
    }
    else {
      TEST_EQUALITY(overlapgraph->getMap()->getNodeNumElements(), num_rows_per_proc+overlap_levels*2)
    }
  }
}

#define UNIT_TEST_GROUP_ORDINAL(LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2CreateOverlapGraph, OverlapGraphTest0, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_ORDINAL(int, int)

}//namespace <anonymous>

