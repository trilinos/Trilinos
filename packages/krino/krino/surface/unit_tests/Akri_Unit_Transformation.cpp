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
// Transformation (where a motion with constant velocities puts a body at a time).
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

TEST(Quaternion, rotationVector_isTheVectorItWasSetFrom)
{
  Quaternion q;
  const stk::math::Vector3d rotationVector(0.3,-0.7,1.1);
  q.set_from_rotation_vector(rotationVector);
  expect_eq(rotationVector, q.rotation_vector());
}

TEST(Quaternion, rotationVectorOfIdentity_isZero)
{
  expect_eq(stk::math::Vector3d::ZERO, Quaternion().rotation_vector());
}

TEST(Quaternion, rotationVectorOfMoreThanHalfATurn_isTheSameRotationTheShortWay)
{
  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(0.,0.,1.5*M_PI));
  expect_eq(stk::math::Vector3d(0.,0.,-halfPi), q.rotation_vector());
}

TEST(Quaternion, rotationVectorOfAComposition_givesTheSameRotation)
{
  Quaternion aboutX, aboutY;
  aboutX.set_from_rotation_vector(stk::math::Vector3d(halfPi,0.,0.));
  aboutY.set_from_rotation_vector(stk::math::Vector3d(0.,halfPi,0.));
  const Quaternion xThenY = aboutY.after(aboutX);
  Quaternion fromRotationVector;
  fromRotationVector.set_from_rotation_vector(xThenY.rotation_vector());
  const stk::math::Vector3d v(1.,2.,3.);
  expect_eq(xThenY.rotate_3d_vector(v), fromRotationVector.rotate_3d_vector(v));
}

//--------------------------------------------------------------------------
// Transformation
//--------------------------------------------------------------------------

TEST(Transformation, default_isIdentityAtAnyTime)
{
  const Transformation motion;
  const stk::math::Vector3d x(1.,2.,3.);
  expect_eq(x, motion.rigid_transformation_at(0.).transform_point(x));
  expect_eq(x, motion.rigid_transformation_at(2.).transform_point(x));
}

TEST(Transformation, pureTranslationVelocity_translatesByVelocityTimesTime)
{
  Transformation motion;
  motion.set_translational_velocity(stk::math::Vector3d(1.,-2.,3.));
  expect_eq(stk::math::Vector3d(12.,6.,16.), motion.rigid_transformation_at(2.).transform_point(stk::math::Vector3d(10.,10.,10.))); // (10,10,10) + 2*(1,-2,3)
}

TEST(Transformation, pureRotationVelocityAboutOrigin_rotatesPoint)
{
  Transformation motion;
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi)); // pi/2 per unit time about z
  expect_eq(stk::math::Vector3d(0.,1.,0.), motion.rigid_transformation_at(1.).transform_point(stk::math::Vector3d(1.,0.,0.)));
}

TEST(Transformation, rotationAboutReferencePoint_leavesReferencePointFixed)
{
  const stk::math::Vector3d refPt(5.,-3.,2.);
  Transformation motion;
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi));
  motion.set_reference_point(refPt);
  expect_eq(refPt, motion.rigid_transformation_at(1.).transform_point(refPt));
}

TEST(Quaternion, after_rotatesByFirstThenByThis)
{
  Quaternion first, second;
  first.set_from_rotation_vector(stk::math::Vector3d(0.3,-0.7,1.1));
  second.set_from_rotation_vector(stk::math::Vector3d(-0.5,0.2,0.4));
  const stk::math::Vector3d v(1.,-2.,0.5);
  expect_eq(second.rotate_3d_vector(first.rotate_3d_vector(v)), second.after(first).rotate_3d_vector(v));
}

//--------------------------------------------------------------------------
// Transformation::rigid_transformation_at: where the motion puts the body
//--------------------------------------------------------------------------

TEST(Transformation, rigidTransformationAt_referencePointMovesWithDisplacementAndVelocity)
{
  const stk::math::Vector3d refPt(5.,-3.,2.);
  const stk::math::Vector3d displacement(1.,2.,3.);
  const stk::math::Vector3d velocity(-1.,0.5,2.);
  Transformation motion;
  motion.set_reference_point(refPt);
  motion.set_initial_displacement(displacement);
  motion.set_initial_rotation(stk::math::Vector3d(0.3,-0.7,1.1));
  motion.set_translational_velocity(velocity);
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi));
  for (const double time : {0., 1., 2.5})
  {
    expect_eq(refPt+displacement+time*velocity, motion.rigid_transformation_at(time).transform_point(refPt));
  }
}

TEST(Transformation, rigidTransformationAt_initialRotationAboutReferencePointThenDisplacement)
{
  const stk::math::Vector3d refPt(1.,1.,0.);
  const stk::math::Vector3d displacement(0.,0.,3.);
  Transformation motion;
  motion.set_reference_point(refPt);
  motion.set_initial_displacement(displacement);
  motion.set_initial_rotation(stk::math::Vector3d(0.,0.,halfPi));
  expect_eq(stk::math::Vector3d(1.,2.,3.), motion.rigid_transformation_at(0.).transform_point(stk::math::Vector3d(2.,1.,0.)));
}

TEST(Transformation, rigidTransformationAt_turnsAboutTheMovingReferencePoint)
{
  const stk::math::Vector3d refPt(1.,0.,0.);
  Transformation motion;
  motion.set_reference_point(refPt);
  motion.set_translational_velocity(stk::math::Vector3d(0.,1.,0.));
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi));
  // After a unit time the reference point is at (1,1,0), and the body has made a quarter turn about it.
  expect_eq(stk::math::Vector3d(1.,2.,0.), motion.rigid_transformation_at(1.).transform_point(stk::math::Vector3d(2.,0.,0.)));
  // The same at time 2, after a half turn, however the motion got there.
  expect_eq(stk::math::Vector3d(0.,2.,0.), motion.rigid_transformation_at(2.).transform_point(stk::math::Vector3d(2.,0.,0.)));
}

TEST(Transformation, rigidTransformationAt_initialRotationThenRotationSinceTimeZero)
{
  const stk::math::Vector3d initialRotation(0.3,-0.7,1.1);
  const stk::math::Vector3d rotationalVelocity(-0.5,0.2,0.4);
  Transformation motion;
  motion.set_initial_rotation(initialRotation);
  motion.set_rotational_velocity(rotationalVelocity);
  const double time = 1.5;
  Quaternion first, second;
  first.set_from_rotation_vector(initialRotation);
  second.set_from_rotation_vector(time*rotationalVelocity);
  const stk::math::Vector3d x(1.,-2.,0.5);
  expect_eq(second.rotate_3d_vector(first.rotate_3d_vector(x)), motion.rigid_transformation_at(time).transform_point(x));
}

TEST(Transformation, rigidTransformation_inverseUndoesTransformAndVectorsAreOnlyRotated)
{
  Transformation motion;
  motion.set_reference_point(stk::math::Vector3d(5.,-3.,2.));
  motion.set_initial_displacement(stk::math::Vector3d(1.,2.,3.));
  motion.set_initial_rotation(stk::math::Vector3d(0.3,-0.7,1.1));
  motion.set_translational_velocity(stk::math::Vector3d(-1.,0.5,2.));
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,halfPi));
  const RigidTransformation transformation = motion.rigid_transformation_at(0.7);
  const stk::math::Vector3d x(1.,-2.,0.5);
  expect_eq(x, transformation.inverse_transform_point(transformation.transform_point(x)));
  expect_eq(x, transformation.inverse_transform_vector(transformation.transform_vector(x)));
  const stk::math::Vector3d y(-4.,1.,3.);
  expect_eq(transformation.transform_point(y)-transformation.transform_point(x), transformation.transform_vector(y-x));
}

}

