// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Transformation.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>

#include <cmath>

namespace krino {

// Unit coverage for the Quaternion rotation helper and the rigid-body
// Transformation (translational + rotational velocity integrated over time).
// These are pure math with no mesh dependence.

namespace {
const double halfPi = 0.5*M_PI;
}

//--------------------------------------------------------------------------
// Quaternion
//--------------------------------------------------------------------------

TEST(Quaternion, defaultIsIdentityRotation)
{
  Quaternion q;
  const stk::math::Vector3d v(1.,2.,3.);
  expect_eq(v, q.rotate_3d_vector(v));
  expect_eq(v, q.reverse_rotate_3d_vector(v));
}

TEST(Quaternion, zeroRotationVectorIsIdentity)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d::ZERO);
  const stk::math::Vector3d v(1.,2.,3.);
  expect_eq(v, q.rotate_3d_vector(v));
}

TEST(Quaternion, ninetyDegreesAboutZ_rotatesInPlane)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(0.,0.,halfPi));
  expect_eq(stk::math::Vector3d(0.,1.,0.), q.rotate_3d_vector(stk::math::Vector3d(1.,0.,0.)));
  expect_eq(stk::math::Vector3d(-1.,0.,0.), q.rotate_3d_vector(stk::math::Vector3d(0.,1.,0.)));
  expect_eq(stk::math::Vector3d(0.,0.,1.), q.rotate_3d_vector(stk::math::Vector3d(0.,0.,1.))); // axis is invariant
}

TEST(Quaternion, oneEightyDegreesAboutX_negatesPerpendicularComponents)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(M_PI,0.,0.));
  expect_eq(stk::math::Vector3d(0.,-1.,0.), q.rotate_3d_vector(stk::math::Vector3d(0.,1.,0.)));
  expect_eq(stk::math::Vector3d(0.,0.,-1.), q.rotate_3d_vector(stk::math::Vector3d(0.,0.,1.)));
  expect_eq(stk::math::Vector3d(1.,0.,0.), q.rotate_3d_vector(stk::math::Vector3d(1.,0.,0.))); // axis is invariant
}

TEST(Quaternion, reverseRotateIsInverseOfRotate)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(0.3,-0.7,1.1));
  const stk::math::Vector3d v(1.,-2.,0.5);
  expect_eq(v, q.reverse_rotate_3d_vector(q.rotate_3d_vector(v)));
  expect_eq(v, q.rotate_3d_vector(q.reverse_rotate_3d_vector(v)));
}

TEST(Quaternion, reverseRotateNinetyAboutZ_rotatesOppositeDirection)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(0.,0.,halfPi));
  expect_eq(stk::math::Vector3d(0.,-1.,0.), q.reverse_rotate_3d_vector(stk::math::Vector3d(1.,0.,0.)));
}

TEST(Quaternion, rotationPreservesLength)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(0.5,1.2,-0.8));
  const stk::math::Vector3d v(3.,-1.,2.);
  EXPECT_NEAR(v.length(), q.rotate_3d_vector(v).length(), 1.e-9);
}

//--------------------------------------------------------------------------
// Transformation
//--------------------------------------------------------------------------

TEST(Transformation, defaultApplyIsNoOp)
{
  Transformation transformation;
  stk::math::Vector3d x(1.,2.,3.);
  transformation.apply(x);
  expect_eq(stk::math::Vector3d(1.,2.,3.), x);
}

TEST(Transformation, pureTranslationVelocity_translatesByVelocityTimesTime)
{
  Transformation transformation;
  const stk::math::Vector3d velocity(1.,-2.,3.);
  transformation.set_translational_velocity(velocity);
  transformation.initialize();

  const double time = 2.;
  transformation.update(time);

  stk::math::Vector3d x(10.,10.,10.);
  transformation.apply(x);
  expect_eq(stk::math::Vector3d(12.,6.,16.), x); // (10,10,10) + 2*(1,-2,3)
}

TEST(Transformation, pureRotationVelocityAboutOrigin_rotatesPoint)
{
  Transformation transformation;
  transformation.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi)); // pi/2 per unit time about z
  transformation.set_reference_point(stk::math::Vector3d::ZERO);
  transformation.initialize();

  transformation.update(1.); // rotate by pi/2

  stk::math::Vector3d x(1.,0.,0.);
  transformation.apply(x);
  expect_eq(stk::math::Vector3d(0.,1.,0.), x);
}

TEST(Transformation, rotationAboutReferencePoint_leavesReferencePointFixed)
{
  const stk::math::Vector3d refPt(5.,-3.,2.);
  Transformation transformation;
  transformation.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi));
  transformation.set_reference_point(refPt);
  transformation.initialize();

  transformation.update(1.);

  stk::math::Vector3d x = refPt;
  transformation.apply(x);
  expect_eq(refPt, x);
}

TEST(Transformation, updateIsCachedForRepeatedTime)
{
  Transformation transformation;
  transformation.set_translational_velocity(stk::math::Vector3d(1.,0.,0.));
  transformation.initialize();

  transformation.update(3.);
  stk::math::Vector3d x1(0.,0.,0.);
  transformation.apply(x1);

  transformation.update(3.); // same time: must not advance the transform again
  stk::math::Vector3d x2(0.,0.,0.);
  transformation.apply(x2);

  expect_eq(stk::math::Vector3d(3.,0.,0.), x1);
  expect_eq(x1, x2);
}

TEST(Transformation, initializeWithDisplacement_appliesOffset)
{
  Transformation transformation;
  const stk::math::Vector3d displacement(1.,2.,3.);
  transformation.initialize(displacement, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO);

  stk::math::Vector3d x(4.,5.,6.);
  transformation.apply(x);
  expect_eq(stk::math::Vector3d(5.,7.,9.), x);
}

}
