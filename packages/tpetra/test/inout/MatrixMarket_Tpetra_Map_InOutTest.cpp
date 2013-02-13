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

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::ParameterList;
using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

namespace {
  template<class NodeType>
  Teuchos::RCP<NodeType>
  getNode() {
    Teuchos::ParameterList defaultParams;
    return Teuchos::rcp (new NodeType (defaultParams));
  }
} // namespace (anonymous)

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MapOutput, ContiguousUniform, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  out << "Test: output a contiguous uniform Tpetra::Map to a Matrix Market file" << endl;

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  // Forestall compiler warnings about unused variables.
  (void) myRank;
  (void) numProcs;

  const size_t localNumElts = 10;
  const global_size_t globalNumElts = localNumElts * as<global_size_t> (numProcs);
  const GO indexBase = 0;

  map_type map (globalNumElts, indexBase, comm, Tpetra::GloballyDistributed, node);

  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;

  writer_type::writeMap (out, map);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MapOutput, ContiguousNonuniform, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  out << "Test: output a contiguous nonuniform Tpetra::Map to a Matrix Market file" << endl;

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();

  // Proc 0 gets a different number of local elements, so that the Map is nonuniform.
  const size_t localNumElts = (myRank == 0) ? 11 : 10;
  const global_size_t globalNumElts = 10 * as<global_size_t> (numProcs) + 1;
  const GO indexBase = 0;

  map_type map (globalNumElts, localNumElts, indexBase, comm, node);

  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (out, map);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MapOutput, Noncontiguous, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename Array<GO>::size_type size_type;

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();

  out << "Test: output a noncontiguous Tpetra::Map to a Matrix Market file" << endl;

  // Just to make sure that the Map is really stored noncontiguously,
  // we only have it contain even-numbered GIDs.
  const size_t localNumElts = 10;
  const global_size_t globalNumElts = 10 * as<global_size_t> (numProcs);
  Array<GO> localGids (localNumElts);
  const GO indexBase = 0;

  const GO localStartGid = 2 * as<GO> (myRank) * as<GO> (localNumElts);
  for (size_t k = 0; k < localNumElts; ++k) {
    localGids[k] = localStartGid + as<GO> (2*k);
  }
  map_type map (globalNumElts, localGids (), indexBase, comm, node);

  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (out, map);
}


// Unit test macro isn't smart enough to deal with namespace qualifications.
typedef Kokkos::DefaultNode::DefaultNodeType the_node_type;

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MapOutput, ContiguousUniform, int, long, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MapOutput, ContiguousNonuniform, int, long, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MapOutput, Noncontiguous, int, long, the_node_type )
