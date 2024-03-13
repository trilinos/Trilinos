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
#include "Tpetra_BlockVector.hpp"
#include "Tpetra_BlockView.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TypeNameTraits.hpp"


namespace { // anonymous

  using Kokkos::ALL;
  using Kokkos::subview;
  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::TypeNameTraits;
  using std::endl;
  
  
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockVector, Issue11931, Scalar, LO, GO, Node ) {
    // This particular issued exposed the fact that BlockVector did not std::move
    // correctly, due to an issue with the internal pointMap
    
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    using BV = Tpetra::BlockVector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;
    
    const int numGlobalEntries = 15*comm->getSize();
    RCP<const map_t> Map = Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));
    
    
    int block_size = 5;
    
    // To test this, we make a BlockVector, std::move it to a second one
    // and then have the first guy go out of scope.
    RCP<BV> Vec2;
    {
      BV Vec1 (*Map,block_size);     
      Vec2 = Teuchos::rcp(new BV(std::move(Vec1)));
    }
    
    Tpetra::global_size_t block_map_size = Vec2->getMap()->getGlobalNumElements();
    Tpetra::global_size_t point_map_size = Vec2->getVectorView().getMap()->getGlobalNumElements();
    TEST_EQUALITY_CONST( block_map_size * block_size, point_map_size );
    
  }





//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockVector, Issue11931, SCALAR, LO, GO, NODE ) 


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)


