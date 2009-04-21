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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Rythmos {

#ifdef RYTHMOS_BROKEN_TEST
TEUCHOS_UNIT_TEST( Rythmos_Thyra, clone_v_detail ) {
  Teuchos::RCP<Thyra::VectorBase<double> > x0; 
  int dim = 1;
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs = 
    Teuchos::rcp(new Thyra::DefaultSpmdVectorSpace<double>(dim));
  Teuchos::RCP<Thyra::VectorBase<double> > x1 = Thyra::createMember(vs); 
  {
    // x1 owns a false RCP to vs
    TEST_EQUALITY_CONST( x1->space().has_ownership(), false );
    // RCP<> for x1 owns a true RCP to vs
    std::string label = "VectorSpaceBase";
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > extra_data_x1 = 
      Teuchos::get_extra_data<Teuchos::RCP<const Thyra::VectorSpaceBase<double> >, Thyra::VectorBase<double> >(x1, label );
    TEST_EQUALITY_CONST( extra_data_x1.has_ownership(), true );
  }
  x0 = x1->clone_v();
  {
    // x0 owns a false RCP to vs
    TEST_EQUALITY_CONST( x0->space().has_ownership(), false );
    // RCP<> for x0 owns a true RCP to a _DIFFERENT_ VectorSpaceBase
    // object because the one used to clone x1 is a false RCP, so the
    // VectorSpaceBase was cloned and that is the one that was set on the RCP.
    std::string label = "VectorSpaceBase";
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > extra_data_x0 = 
      Teuchos::get_extra_data<Teuchos::RCP<const Thyra::VectorSpaceBase<double> >, Thyra::VectorBase<double> >(x0, label );
    TEST_EQUALITY_CONST( extra_data_x0.has_ownership(), true );
    TEST_COMPARE( extra_data_x0.ptr(), !=, vs.ptr() );
  }
  vs = Teuchos::null; // vs still around because x1's RCP owns it
  x1 = Teuchos::null; // vs deleted 
  {
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs_old = x0->space();
    TEST_EQUALITY_CONST( vs_old->dim(), 1 ); // INVALID READ => VALGRIND ERROR!
  }
  
}

TEUCHOS_UNIT_TEST( Rythmos_Thyra, clone_v ) {
  Teuchos::RCP<Thyra::VectorBase<double> > x0;
  {
    int dim = 1;
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs = 
      Teuchos::rcp(new Thyra::DefaultSpmdVectorSpace<double>(dim));
    Teuchos::RCP<Thyra::VectorBase<double> > x1 = Thyra::createMember(vs);
    V_S(Teuchos::outArg(*x1),2.0);
    x0 = x1->clone_v();
  }
  TEST_ASSERT(!is_null(x0->space()));
}
#endif // RYTHMOS_BROKEN_TEST

} // namespace Rythmos


