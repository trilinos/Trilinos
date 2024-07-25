// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Include here just the bare minimum needed for the outer #if.
// Only include more once we know we need the test.
#include "TpetraCore_config.h"

#if defined(HAVE_TPETRA_SYCL)

#include "Tpetra_Test_CrsMatrix_WithGraph.hpp"

namespace Tpetra {
namespace Test {

//
// INSTANTIATIONS
//

TPETRA_ETI_MANGLING_TYPEDEFS()

// Declare a colon- and comma-free typedef, to avoid macro issues.
typedef Tpetra::KokkosCompat::KokkosSYCLWrapperNode sycl_node_type;

#define UNIT_TEST_GROUP_SYCL( SCALAR, LO, GO ) \
  UNIT_TEST_GROUP( SCALAR, LO, GO, sycl_node_type )

TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SYCL )

} // namespace Test
} // namespace Tpetra

#endif // defined(HAVE_TPETRA_SYCL)

