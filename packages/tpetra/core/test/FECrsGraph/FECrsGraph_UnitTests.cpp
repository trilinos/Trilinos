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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_FECrsGraph.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Assembly_Helpers.hpp"
#include "Tpetra_Details_getNumDiags.hpp"

namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Tpetra::createContigMapWithNode;
using Tpetra::StaticProfile;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::outArg;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using std::endl;
typedef Tpetra::global_size_t GST;


#define GLOBAL_SUCCESS_CHECK(out,comm,success)  \
  { \
    int lclSuccess = success ? 1 : 0; \
    int gblSuccess = 1; \
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MIN, lclSuccess, Teuchos::outArg (gblSuccess)); \
    if (gblSuccess == 1) { \
      out << "Succeeded on all processes!" << endl; \
    } else { \
      out << "FAILED on at least one process!" << endl; \
    } \
    TEST_EQUALITY_CONST(gblSuccess, 1);  \
  }

//
// UNIT TESTS
//



////
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FECrsGraph, Diagonal, LO, GO, Node )
{
    typedef Tpetra::FECrsGraph<LO,GO,Node> FEG;
    typedef Tpetra::CrsGraph<LO,GO,Node> CG;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    // create a Map
    const size_t numLocal = 10;
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);


    // Trivial test that makes sure a diagonal graph can be built
    CG g1(map,1,StaticProfile);
    FEG g2(map,map,1);

    Tpetra::beginFill(g2);
    for(size_t i=0; i<numLocal; i++) {
      GO gid = map->getGlobalElement(i);
      g1.insertGlobalIndices(gid,1,&gid);
      g2.insertGlobalIndices(gid,1,&gid);
    }
    Tpetra::endFill(g2);
    g1.fillComplete();
    g2.fillComplete();

    // FIXME: Use graph comparison here


    GLOBAL_SUCCESS_CHECK(out,comm,true)
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FECrsGraph, Diagonal, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // end namespace (anonymous)


