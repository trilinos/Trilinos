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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include "Tpetra_Map.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Directory_def.hpp"
#endif

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;

////
TEUCHOS_UNIT_TEST( Map, Bug5378_GoodGIDs )
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<int,long> > map = Tpetra::createContigMap<int,long>(10,10,comm);
  Array<long> lookup_gids(  tuple<long>(1,3,5) );
  Array<int> expected_ids(  tuple<int>( 0,0,0) );
  Array<int> expected_lids( tuple<int>( 1,3,5) );
  Array<int> nodeIDs( lookup_gids.size() ),
             nodeLIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs(), nodeLIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::AllIDsPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );
  TEST_COMPARE_ARRAYS( nodeLIDs(), expected_lids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

////
TEUCHOS_UNIT_TEST( Map, Bug5378_BadGIDs )
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<int,long> > map = Tpetra::createContigMap<int,long>(10,10,comm);
  Array<long> lookup_gids(  tuple<long>(1,10,5) );
  Array<int> expected_ids(  tuple<int>( 0,-1,0) );
  Array<int> expected_lids( tuple<int>( 1,-1,5) );
  Array<int> nodeIDs( lookup_gids.size() ),
             nodeLIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs(), nodeLIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );
  TEST_COMPARE_ARRAYS( nodeLIDs(), expected_lids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

////
TEUCHOS_UNIT_TEST( Map, Bug5378_GoodGIDsNoLIDs )
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<int,long> > map = Tpetra::createContigMap<int,long>(10,10,comm);
  Array<long> lookup_gids(  tuple<long>(1,3,5) );
  Array<int> expected_ids(  tuple<int>( 0,0,0) );
  Array<int> nodeIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::AllIDsPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

////
TEUCHOS_UNIT_TEST( Map, Bug5378_BadGIDsNoLIDs )
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<int,long> > map = Tpetra::createContigMap<int,long>(10,10,comm);
  Array<long> lookup_gids(  tuple<long>(1,10,5) );
  Array<int> expected_ids(  tuple<int>( 0,-1,0) );
  Array<int> nodeIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}


