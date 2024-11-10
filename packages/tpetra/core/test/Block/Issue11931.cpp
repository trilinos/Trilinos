// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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


