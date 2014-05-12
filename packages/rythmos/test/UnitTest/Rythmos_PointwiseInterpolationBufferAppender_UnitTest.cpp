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

RCP<InterpolationBuffer<double> > createIBSource(
    const RCP<const Thyra::VectorSpaceBase<double> >& vs, // FOR BUG 4388
    const Array<double>& t_vals, 
    const Array<double>& x_vals, 
    const Array<double>& xdot_vals)
{
  int N = t_vals.size();
  TEUCHOS_ASSERT( N == as<int>(x_vals.size()) );
  TEUCHOS_ASSERT( N == as<int>(xdot_vals.size()) );
  RCP<InterpolationBuffer<double> > ibSource = interpolationBuffer<double>(Teuchos::null,N);
  Teuchos::Array<double> time_vec(N);
  Array<RCP<const VectorBase<double> > > x_vec(N);
  Array<RCP<const VectorBase<double> > > xdot_vec(N);
  for (int i=0 ; i<N ; ++i) {
    time_vec[i] = t_vals[i];
    x_vec[i] = createDefaultVector(vs,x_vals[i]); 
    xdot_vec[i] = createDefaultVector(vs,xdot_vals[i]); 
  }
  ibSource->addPoints(time_vec,x_vec,xdot_vec);
  return ibSource;
}

RCP<InterpolationBuffer<double> > createDefaultIBSource(
    const RCP<const Thyra::VectorSpaceBase<double> >& vs // FOR BUG 4388
    )
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
  RCP<InterpolationBuffer<double> > ib = createIBSource( vs, time_vec, x_vec, xdot_vec );
  return ib;
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_all_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();

  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388
  RCP<InterpolationBuffer<double> > ibSource = createDefaultIBSource(vs);
  Array<double> time_vec;
  ibSource->getNodes(&time_vec);
  int N = time_vec.size();

  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> appendRange(0.0,1.0); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out(N);
  ibSink->getNodes(&time_vec_out);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_first_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> appendRange(0.0,0.5); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify only the first point got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  time_vec_valid.push_back(time_vec[0]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_last_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> appendRange(0.5,1.0); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify only the last point got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  time_vec_valid.push_back(time_vec[1]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_none_to_empty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  RCP<InterpolationBuffer<double> > ibSink = interpolationBuffer<double>();
  TimeRange<double> appendRange(0.25,0.75); 
#ifdef HAVE_RYTHMOS_DEBUG
  // no time points to add, will throw (04/21/09 tscoffe:  Is that the right behavior?)
  TEST_THROW(piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink)), std::logic_error);
#else // HAVE_RYTHMOS_DEBUG
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify no points got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
#endif // HAVE_RYTHMOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_last_to_nonempty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(3);
  Array<double> time_vec_source = tuple<double>(1.0,2.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec_source,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  TimeRange<double> appendRange(1.0,2.0); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  time_vec_valid.push_back(time_vec_sink[0]);
  time_vec_valid.push_back(time_vec_sink[1]);
  time_vec_valid.push_back(time_vec_source[1]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_first_to_nonempty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(1.0,2.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(3);
  Array<double> time_vec_source = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec_source,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  TimeRange<double> appendRange(0.0,1.0); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  time_vec_valid.push_back(time_vec_source[0]);
  time_vec_valid.push_back(time_vec_sink[0]);
  time_vec_valid.push_back(time_vec_sink[1]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, append_all_to_nonempty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(0.0,1.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(4);
  Array<double> time_vec_source = tuple<double>(2.0,3.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec_source,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  TimeRange<double> appendRange(2.0,3.0); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  time_vec_valid.push_back(time_vec_sink[0]);
  time_vec_valid.push_back(time_vec_sink[1]);
  time_vec_valid.push_back(time_vec_source[0]);
  time_vec_valid.push_back(time_vec_source[1]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, prepend_all_to_nonempty ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(3.0,4.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(4);
  Array<double> time_vec_source = tuple<double>(1.0,2.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec_source,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  TimeRange<double> appendRange(1.0,2.0); 
  piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
  // Verify both points got added:
  Array<double> time_vec_out;
  ibSink->getNodes(&time_vec_out);
  Array<double> time_vec_valid;
  time_vec_valid.push_back(time_vec_source[0]);
  time_vec_valid.push_back(time_vec_source[1]);
  time_vec_valid.push_back(time_vec_sink[0]);
  time_vec_valid.push_back(time_vec_sink[1]);
  TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
}

// Test appendRange overlaps sinkRange (both prepend and append) (throw, throw)
TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, invalid_appendRange ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(3.0,4.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(3);
  Array<double> time_vec_source = tuple<double>(1.0,2.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec_source,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  {
    TimeRange<double> appendRange(2.0,3.5); 
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW(
      piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink)),
      std::logic_error
      );
#else // HAVE_RYTHMOS_DEBUG
    piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
    // Verify one points got added:
    Array<double> time_vec_out;
    ibSink->getNodes(&time_vec_out);
    Array<double> time_vec_valid = tuple<double>(2.0,3.0,4.0);
    TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
#endif // HAVE_RYTHMOS_DEBUG
  }
  {
    TimeRange<double> appendRange(3.5,5.0); 
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW(
      piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink)),
      std::logic_error
      );
#else // HAVE_RYTHMOS_DEBUG
    piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
    // Verify no points got added:
    Array<double> time_vec_out;
    ibSink->getNodes(&time_vec_out);
    Array<double> time_vec_valid = tuple<double>(2.0,3.0,4.0);
    TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
#endif // HAVE_RYTHMOS_DEBUG
  }
}
// Test appendRange does not sit inside sourceRange (throw)
TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, invalid_appendRange2 ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(3.0,4.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(6);
  {
    Array<double> time_vec_source = tuple<double>(1.0,2.0);
    RCP<InterpolationBuffer<double> > ibSource = createIBSource(
        vs,
        time_vec_source,
        tuple<double>(1.0,2.0),
        tuple<double>(3.0,4.0)
        );
    TimeRange<double> appendRange(1.0,2.5); 
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW(
      piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink)),
      std::logic_error
      );
#else // HAVE_RYTHMOS_DEBUG
    piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
    // Verify both points got added:
    Array<double> time_vec_out;
    ibSink->getNodes(&time_vec_out);
    Array<double> time_vec_valid = tuple<double>(1.0,2.0,3.0,4.0);
    TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
#endif // HAVE_RYTHMOS_DEBUG
  }
  {
    Array<double> time_vec_source = tuple<double>(5.0,6.0);
    RCP<InterpolationBuffer<double> > ibSource = createIBSource(
        vs,
        time_vec_source,
        tuple<double>(1.0,2.0),
        tuple<double>(3.0,4.0)
        );
    TimeRange<double> appendRange(4.5,6.0); 
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW(
      piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink)),
      std::logic_error
      );
#else // HAVE_RYTHMOS_DEBUG
    piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
    // Verify both points got added:
    Array<double> time_vec_out;
    ibSink->getNodes(&time_vec_out);
    Array<double> time_vec_valid = tuple<double>(1.0,2.0,3.0,4.0,5.0,6.0);
    TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
#endif // HAVE_RYTHMOS_DEBUG
  }
}

// Test that sourceRange can sit inside sinkRange as long as appendRange does not
TEUCHOS_UNIT_TEST( Rythmos_PointwiseInterpolationBufferAppender, valid_append ) {
  RCP<PointwiseInterpolationBufferAppender<double> > piba = 
    pointwiseInterpolationBufferAppender<double>();
  int dim = 2;
  RCP<const Thyra::VectorSpaceBase<double> > vs = createDefaultVectorSpace<double>(dim); // FOR BUG 4388

  Array<double> time_vec_sink = tuple<double>(3.0,4.0);
  RCP<InterpolationBuffer<double> > ibSink = createIBSource(
      vs,
      time_vec_sink,
      tuple<double>(1.0,2.0),
      tuple<double>(3.0,4.0)
      );
  ibSink->setStorage(6);
  Array<double> time_vec_source = tuple<double>(1.0,2.0,3.0,4.0,5.0,6.0);
  RCP<InterpolationBuffer<double> > ibSource = createIBSource(
      vs,
      time_vec_source,
      tuple<double>(1.0,2.0,3.0,4.0,5.0,6.0),
      tuple<double>(7.0,8.0,9.0,10.0,11.0,12.0)
      );
  {
    TimeRange<double> appendRange(1.0,2.5); 
    piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
    // Verify two points got added:
    Array<double> time_vec_out;
    ibSink->getNodes(&time_vec_out);
    Array<double> time_vec_valid;
    time_vec_valid.push_back(time_vec_source[0]);
    time_vec_valid.push_back(time_vec_source[1]);
    time_vec_valid.push_back(time_vec_sink[0]);
    time_vec_valid.push_back(time_vec_sink[1]);
    TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
  }
  {
    TimeRange<double> appendRange(4.5,6.0); 
    piba->append(*ibSource,appendRange,Teuchos::outArg(*ibSink));
    // Verify two points got added:
    Array<double> time_vec_out;
    ibSink->getNodes(&time_vec_out);
    Array<double> time_vec_valid;
    time_vec_valid.push_back(1.0);
    time_vec_valid.push_back(2.0);
    time_vec_valid.push_back(3.0);
    time_vec_valid.push_back(4.0);
    time_vec_valid.push_back(5.0);
    time_vec_valid.push_back(6.0);
    TEST_COMPARE_ARRAYS( time_vec_out, time_vec_valid );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBufferAppenderBase, selectPointsInTimeRange ) {
  Array<double> points_in = tuple<double>(1.0,2.0,3.0,4.0,5.0);
  {
    Array<double> points_out;
    TimeRange<double> tr(1.0,2.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(1.0,2.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange<double> tr(0.0,3.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(1.0,2.0,3.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange<double> tr(1.5,3.3);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(2.0,3.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange<double> tr(0.0,0.5);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points;
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange<double> tr(6.0,7.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points;
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange<double> tr(4.0,7.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(4.0,5.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
}
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBufferAppenderBase, selectPointsInTimeRange_oc ) {
  Array<double> points_in = tuple<double>(1.0,2.0,3.0,4.0,5.0);
  {
    Array<double> points_out;
    TimeRange_oo<double> tr(2.0,5.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(3.0,4.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange_oc<double> tr(1.0,3.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(2.0,3.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange_co<double> tr(3.0,5.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(3.0,4.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
  {
    Array<double> points_out;
    TimeRange_cc<double> tr(2.0,3.0);
    selectPointsInTimeRange(points_in,tr,Teuchos::outArg(points_out));
    Array<double> valid_points = tuple<double>(2.0,3.0);
    TEST_COMPARE_ARRAYS(points_out,valid_points);
  }
}

} // namespace Rythmos


