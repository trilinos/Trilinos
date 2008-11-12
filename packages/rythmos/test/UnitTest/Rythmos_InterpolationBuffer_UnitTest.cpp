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

#include "Rythmos_UnitTestHelpers.hpp"

#include "Rythmos_InterpolationBuffer.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, newBuffer ) {
  InterpolationBuffer<double> ib;
  TEST_EQUALITY_CONST( ib.getStorage(), 2 );
  TEST_EQUALITY_CONST( ib.getTimeRange().isValid(), false ); 
  TEST_EQUALITY_CONST( ib.getOrder(), 1 ); // linear interpolator by default
  TEST_EQUALITY( ib.getParameterList(), Teuchos::null );
}

// Verify we can append entries to an IB correctly
//   * try passing a range of points that includes the last point in the IB
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, invalidAddVectors ) {
  Array<double> time_vec;
  time_vec.push_back(2.0);
  time_vec.push_back(3.0);

  Array<RCP<const VectorBase<double> > > v_vec;
  RCP<VectorBase<double> > v = createDefaultVector(2,2.0);
  v_vec.push_back(v);
  v = createDefaultVector(2,3.0);
  v_vec.push_back(v);

  Array<RCP<const VectorBase<double> > > vdot_vec;
  RCP<VectorBase<double> > vdot = createDefaultVector(2,4.0);
  vdot_vec.push_back(vdot);
  vdot = createDefaultVector(2,5.0);
  vdot_vec.push_back(vdot);

  InterpolationBuffer<double> ib;
  ib.addPoints(time_vec,v_vec,vdot_vec);

  Array<double> newTime_vec;
  newTime_vec.push_back(3.0);
  newTime_vec.push_back(4.0);

  Array<RCP<const VectorBase<double> > > newV_vec;
  RCP<VectorBase<double> > newV = createDefaultVector(2,6.0);
  newV_vec.push_back(newV);
  newV = createDefaultVector(2,7.0);
  newV_vec.push_back(newV);

  Array<RCP<const VectorBase<double> > > newVdot_vec;
  RCP<VectorBase<double> > newVdot = createDefaultVector(2,8.0);
  newVdot_vec.push_back(newVdot);
  newVdot = createDefaultVector(2,9.0);
  newVdot_vec.push_back(newVdot);

  TEST_THROW( ib.addPoints(newTime_vec, newV_vec, newVdot_vec), std::logic_error);
}

// Verify that the IB is copying the vectors rather than just storing the pointers
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, copyVectors ) {
  double t = 2.5;
  Array<double> time_vec;
  time_vec.push_back(t);

  RCP<VectorBase<double> > v = createDefaultVector(2,2.0);
  Array<RCP<const VectorBase<double> > > v_vec;
  v_vec.push_back(v);

  RCP<VectorBase<double> > v_dot = createDefaultVector(2,3.0);
  Array<RCP<const VectorBase<double> > > v_dot_vec;
  v_dot_vec.push_back(v_dot);

  InterpolationBuffer<double> ib;
  ib.addPoints(time_vec, v_vec, v_dot_vec);

  Thyra::V_S(&*v, 4.0);
  Thyra::V_S(&*v_dot, 5.0);

  Array<RCP<const VectorBase<double> > > v_vec_out;
  Array<RCP<const VectorBase<double> > > v_dot_vec_out;
  Array<double> accuracy_vec_out;
  ib.getPoints(time_vec, &v_vec_out, &v_dot_vec_out, &accuracy_vec_out);

  TEST_EQUALITY_CONST( get_ele(*(v_vec_out[0]),0), 2.0 );
  TEST_EQUALITY_CONST( get_ele(*(v_dot_vec_out[0]),0), 3.0 );
}

// TODO Test the storage limit and the buffer policies

} // namespace Rythmos



