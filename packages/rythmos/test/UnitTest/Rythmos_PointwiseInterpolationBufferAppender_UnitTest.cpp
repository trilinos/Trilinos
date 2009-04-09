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

#include "Rythmos_Types.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {

using Teuchos::tuple;
using Teuchos::as;

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, create ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  TEST_ASSERT(!is_null(piba));
}

/*
TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_all_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();

  int N = 2;
  Teuchos::Array<double> time_vec;
  RCP<InterpolationBuffer<double> > ibSource = interpolationBuffer<double>();
  ibSource->setStorage(N);
  {
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    int dim = 2;
    for (int i=0 ; i<N ; ++i) {
      time_vec.push_back(0.0+1.0*i);
      x_vec.push_back(createDefaultVector(dim,1.0+1.0*i)); 
      xdot_vec.push_back(createDefaultVector(dim,3.0+1.0*i)); 
    }
    ibSource->addPoints(time_vec,x_vec,xdot_vec);
  }
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  //TimeRange<double> sourceRange(0.0,1.0); 
  {
    //Array<double> time_vec_in;
    //ibSource->getNodes(&time_vec_in);

    //Array<double> time_vec_out;
    //selectPointsInTimeRange(&time_vec_out,time_vec_in,sourceRange);

    Array<RCP<const Thyra::VectorBase<double> > > x_vec;
    Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    ibSource->getPoints(time_vec, &x_vec, &xdot_vec, &accuracy_vec);
    ibSink->addPoints(time_vec, x_vec, xdot_vec);
  }
//  piba->append(*ibSource,sourceRange,Teuchos::outArg(*ibSink));
//  // Verify both points got added:
//  Array<double> time_vec_out(N);
//  ibSink->getNodes(&time_vec_out);
//  TEST_COMPARE_ARRAYS( time_vec_out, time_vec );
}
*/

/*
RCP<InterpolationBuffer<double> > createIBSource(
    const Array<double>& t_vals, 
    const Array<double>& x_vals, 
    const Array<double>& xdot_vals)
{
  int N = t_vals.size();
  TEUCHOS_ASSERT( N == x_vals.size() );
  TEUCHOS_ASSERT( N == xdot_vals.size() );
  RCP<InterpolationBuffer<double> > ibSource = interpolationBuffer<double>(Teuchos::null,N);
  Teuchos::Array<double> time_vec(N);
  Array<RCP<const VectorBase<double> > > x_vec(N);
  Array<RCP<const VectorBase<double> > > xdot_vec(N);
  int dim = 2;
  for (int i=0 ; i<N ; ++i) {
    time_vec[i] = t_vals[i];
    x_vec[i] = createDefaultVector(dim,x_vals[i]); 
    xdot_vec[i] = createDefaultVector(dim,xdot_vals[i]); 
  }
  ibSource->addPoints(time_vec,x_vec,xdot_vec);
  return ibSource;
}

RCP<InterpolationBuffer<double> > createDefaultIBSource()
{
  int N = 2;
  Array<double> time_vec(N);
  time_vec[0] = 0.0;
  time_vec[1] = 1.0;
  Array<double> x_vec(N);
  x_vec[0] = 1.0;
  x_vec[1] = 2.0;
  Array<double> xdot_vec(N);
  xdot_vec[0] = 3.0;
  xdot_vec[1] = 4.0;
  RCP<InterpolationBuffer<double> > ib = createIBSource( time_vec, x_vec, xdot_vec );
  return ib;
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_all_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();

  RCP<InterpolationBuffer<double> > ibSource = createDefaultIBSource();
  Array<double> time_vec;
  ibSource->getNodes(&time_vec);
  int N = time_vec.size();

  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> sourceRange(0.0,1.0); 
  piba->append(*ibSource,sourceRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out(N);
  ibSink->getNodes(&time_vec_out);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_first_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int N = 2;
  Array<double> time_vec = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      time_vec,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> sourceRange(0.0,0.5); 
  piba->append(*ibSource,sourceRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out(N);
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid(1);
  time_vec_valid.push_back(time_vec[0]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_last_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int N = 2;
  Array<double> time_vec = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      time_vec,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> sourceRange(0.5,1.0); 
  piba->append(*ibSource,sourceRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out(N);
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid(1);
  time_vec_valid.push_back(time_vec[1]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_none_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int N = 2;
  Array<double> time_vec = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      time_vec,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> sourceRange(0.25,0.75); 
  piba->append(*ibSource,sourceRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out(N);
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid(0);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}
*/

// ADD MORE TESTS!!!!!!!!!!!!!!!


} // namespace Rythmos


