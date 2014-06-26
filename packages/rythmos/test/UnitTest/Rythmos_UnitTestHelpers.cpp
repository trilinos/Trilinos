//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_UnitTestHelpers.hpp"

namespace Rythmos {
  
// Acceptance test for creating default Thyra Vectors
TEUCHOS_UNIT_TEST( Rythmos_UnitTestHelpers, createDefaultVectorSpace ) {
  int N = 5;
  double value = 2.0;
  Teuchos::RCP<Thyra::VectorBase<double> > vec = createDefaultVector<double>(N,value);
  TEST_COMPARE( vec, !=, Teuchos::null );
  TEST_EQUALITY( vec->space()->dim(), N );
  for (int i=0 ; i<N ; ++i) {
    TEST_EQUALITY( Thyra::get_ele(*vec,i), value );
  }
}

#ifdef HAVE_RYTHMOS_DEBUG
TEUCHOS_UNIT_TEST( Rythmos_UnitTestHelpers, debug_mode_is_ON ) {
  TEST_ASSERT(true);
}
#else // HAVE_RYTHMOS_DEBUG
TEUCHOS_UNIT_TEST( Rythmos_UnitTestHelpers, debug_mode_is_OFF ) {
  TEST_ASSERT(true);
}
#endif // HAVE_RYTHMOS_DEBUG

} // namespace Rythmos

