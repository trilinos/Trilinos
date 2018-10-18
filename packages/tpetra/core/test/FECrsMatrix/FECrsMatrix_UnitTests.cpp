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
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_Details_getNumDiags.hpp"

// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Tpetra::createContigMapWithNode;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::outArg;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::NO_TRANS;
//using Teuchos::TRANS;
using Teuchos::CONJ_TRANS;
using std::endl;
typedef Tpetra::global_size_t GST;

#define STD_TESTS(matrix) \
{ \
    using Teuchos::outArg; \
    RCP<const Comm<int> > STCOMM = matrix.getComm(); \
    ArrayView<const GO> STMYGIDS = matrix.getRowMap()->getNodeElementList(); \
    ArrayView<const LO> loview; \
    ArrayView<const Scalar> sview; \
    size_t STMAX = 0; \
    for (size_t STR=0; STR < matrix.getNodeNumRows(); ++STR) { \
        const size_t numEntries = matrix.getNumEntriesInLocalRow(STR); \
        TEST_EQUALITY( numEntries, matrix.getNumEntriesInGlobalRow( STMYGIDS[STR] ) ); \
        matrix.getLocalRowView(STR,loview,sview); \
        TEST_EQUALITY( static_cast<size_t>(loview.size()), numEntries ); \
        TEST_EQUALITY( static_cast<size_t>( sview.size()), numEntries ); \
        STMAX = std::max( STMAX, numEntries ); \
    } \
    TEST_EQUALITY( matrix.getNodeMaxNumRowEntries(), STMAX ); \
    GST STGMAX; \
    Teuchos::reduceAll<int,GST>( *STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg(STGMAX) ); \
    TEST_EQUALITY( matrix.getGlobalMaxNumRowEntries(), STGMAX ); \
}



//
// UNIT TESTS
//



////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FECrsMatrix, ZeroMatrix, LO, GO, Scalar, Node )
{
    typedef Tpetra::FECrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Teuchos::ScalarTraits<Scalar>          ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType             Mag;
    typedef Teuchos::ScalarTraits<Mag>             MT;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);

    // SCAFFOLDING - just some prints to keep the -Wunused-variable happy (for now)
    std::cout << "numVecs: " << numVecs << std::endl;

    // create the zero matrix
    /*
    MAT zero(map,0);
    zero.fillComplete();

    //
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();
    mvres.putScalar(1);
    zero.apply(mvrand,mvres);
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm1(norms());
    if (ST::isOrdinal)
    {
        TEST_COMPARE_ARRAYS(norms,zeros);
    }
    else
    {
        TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
    */
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FECrsMatrix, ZeroMatrix, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // end namespace (anonymous)


