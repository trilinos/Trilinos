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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>

#include <iostream>

using Teuchos::RCP;
using Teuchos::ArrayView;
using Teuchos::tuple;

RCP<Tpetra::Vector<int> >
TestTpetra (const Teuchos::ArrayView<const Tpetra::Map<>::global_ordinal_type>& srcGID,
            const Teuchos::ArrayView<const Tpetra::Map<>::global_ordinal_type>& destGID)
{
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::global_size_t GST;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  const GO indexBase = 0;
  const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const Tpetra::Map<> >  srcMap (new Tpetra::Map<> (INV, srcGID (), indexBase, comm));
  RCP<const Tpetra::Map<> > destMap (new Tpetra::Map<> (INV, destGID (), indexBase, comm));

  RCP<Tpetra::Vector<int> >  srcVector (new Tpetra::Vector<int> (srcMap));
  RCP<Tpetra::Vector<int> > destVector (new Tpetra::Vector<int> (destMap));
  destVector->putScalar (-1);

  Tpetra::Export<> exporter (srcMap, destMap);
  destVector->doExport (*srcVector, exporter, Tpetra::INSERT);

  Teuchos::FancyOStream out (Teuchos::rcpFromRef (std::cout));
  destVector->describe (out, Teuchos::VERB_EXTREME);

  return destVector;
}

////
TEUCHOS_UNIT_TEST( DistObject, SubMapExport1 )
{
  typedef Tpetra::Vector<int>::global_ordinal_type GO;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int MyPid = comm->getRank();
  /**********************************************************************************/
  // Maps of the export operation:
  // -----------------------------
  // SRC Map Processor 0: Global IDs = 2 3 4 5 6 7 8
  //         Processor 1: Global IDs =           7 8 9 10 11 12 13
  //
  // DEST Map  Processor 0: Global IDs = 0 1 2 3 4 5 6
  //           Processor 1: Global IDs =                    9 10 11 12 13 14 15
  //
  //
  // Vectors before export operation:
  // --------------------------------
  // srcVector  = [ 0  0  ...  0 ]
  // destVector = [-1 -1  ... -1 ]
  //
  // Expected result:
  // ----------------
  // destVector Processor 0: Values = [ -1 -1 0 0 0 0 0 ]
  //            Processor 1: Values =           [ 0 0 0 0 0 -1 -1 ]
  RCP<Tpetra::Vector<int> > destVector;
  if (MyPid == 0) {
    TEST_NOTHROW( destVector = TestTpetra (tuple<GO> (2,3,4,5,6,7,8), tuple<GO> (0,1,2,3,4,5,6) ) );
    TEST_COMPARE_ARRAYS( tuple<int>(-1,-1,0,0,0,0,0), destVector->get1dView() )
  }
  else {
    TEST_NOTHROW( destVector = TestTpetra (tuple<GO>(7,8,9,10,11,12,13), tuple<GO>(9,10,11,12,13,14,15) ) );
    TEST_COMPARE_ARRAYS( tuple<int>(0,0,0,0,0,-1,-1), destVector->get1dView() )
  }
}

////
TEUCHOS_UNIT_TEST( DistObject, SubMapExport2 )
{
  typedef Tpetra::Vector<int>::global_ordinal_type GO;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int MyPid = comm->getRank();
  /**********************************************************************************/
  // SRC Map  Processor 0: Global IDs =
  //          Processor 1: Global IDs = 0 1
  //
  // DEST Map Processor 0: Global IDs = 0 1 2
  //          Processor 1: Global IDs =   1 2
  //
  // Vectors before export operation:
  // --------------------------------
  // srcVector  = [] [0  0]
  // destVector = [-1 -1 -1] [-1 -1]
  //
  // Expected result:
  // destVector Processor 0: Values = 0 -1 -1
  //            Processor 1: Values =    0 -1
  RCP<Tpetra::Vector<int> > destVector;
  if (MyPid == 0) {
    TEST_NOTHROW( destVector = TestTpetra (ArrayView<GO>(), tuple<GO>(0,1,2) ) )
    TEST_COMPARE_ARRAYS( tuple<int>(0,-1,-1), destVector->get1dView() )
  }
  else {
    TEST_NOTHROW( destVector = TestTpetra (tuple<GO>(0,1), tuple<GO>(1,2) ) )
    TEST_COMPARE_ARRAYS( tuple<int>(0,-1), destVector->get1dView() )
  }
}

////
TEUCHOS_UNIT_TEST( DistObject, SubMapExport3 )
{
  using GO = Tpetra::Vector<int>::global_ordinal_type;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int MyPid = comm->getRank();
  /**********************************************************************************/
  // SRC Map  Processor 0: Global IDs = 0 1
  //          Processor 1: Global IDs = 0 1
  //
  // DEST Map Processor 0: Global IDs = 2 3
  //          Processor 1: Global IDs = 2 3
  //
  // Vectors before export operation:
  // --------------------------------
  // srcVector  = 0 0
  // destVector = -1 -1
  //
  // Expected result:
  // destVector Processor 0: Values = -1 -1
  //            Processor 1: Values = -1 -1
  RCP<Tpetra::Vector<int> > destVector;
  if (MyPid == 0) {
    TEST_NOTHROW( destVector = TestTpetra (tuple<GO> (0,1), tuple<GO> (2,3) ) )
  }
  else {
    TEST_NOTHROW( destVector = TestTpetra (tuple<GO> (0,1), tuple<GO> (2,3)) )
  }
  TEST_COMPARE_ARRAYS( tuple<int>(-1,-1), destVector->get1dView() )
}


