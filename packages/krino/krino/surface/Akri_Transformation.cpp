// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Transformation.hpp>

#include <cmath>
#include <limits>

namespace krino{

void
Quaternion::set_from_rotation_vector(const stk::math::Vector3d & v)
{
  // the angle theta is the length of the rotation vector omega
  const double theta = v.length();
  const double real_min = 10.0*std::numeric_limits<double>::min();

  double coef;
  if ( theta > real_min ) {
    coef = std::sin(0.5*theta)/theta;
  } else {
    coef = 0.5;
  }

  q[0] = std::cos( 0.5*theta );
  q[1] = coef*v[0];
  q[2] = coef*v[1];
  q[3] = coef*v[2];

  const double inv_length = 1.0/std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  for (auto && val : q) val *= inv_length;
}

stk::math::Vector3d use_quaternion_to_rotate_3d_vector(const std::array<double,4> & q, const stk::math::Vector3d & v)
{
  return stk::math::Vector3d(
    ((2. * ( q[0] * q[0] + q[1] * q[1] ) - 1.) * v[0] +
     (2. * ( q[1] * q[2] - q[0] * q[3] )     ) * v[1] +
     (2. * ( q[1] * q[3] + q[0] * q[2] )     ) * v[2]),
    ((2. * ( q[1] * q[2] + q[0] * q[3] )     ) * v[0] +
     (2. * ( q[0] * q[0] + q[2] * q[2] ) - 1.) * v[1] +
     (2. * ( q[2] * q[3] - q[0] * q[1] )     ) * v[2]),
    ((2. * ( q[1] * q[3] - q[0] * q[2] )     ) * v[0] +
     (2. * ( q[2] * q[3] + q[0] * q[1] )     ) * v[1] +
     (2. * ( q[0] * q[0] + q[3] * q[3] ) - 1.) * v[2]));
}


stk::math::Vector3d
Quaternion::rotate_3d_vector(const stk::math::Vector3d & v) const
{
  return use_quaternion_to_rotate_3d_vector(q, v);
}

stk::math::Vector3d
Quaternion::reverse_rotate_3d_vector(const stk::math::Vector3d & v) const
{
  return use_quaternion_to_rotate_3d_vector(std::array<double,4>{q[0],-q[1],-q[2],-q[3]}, v);
}

Quaternion
Quaternion::after(const Quaternion & first) const
{
  // The Hamilton product q*first.
  const std::array<double,4> & p = first.q;
  Quaternion product;
  product.q = {{q[0]*p[0] - q[1]*p[1] - q[2]*p[2] - q[3]*p[3],
      q[0]*p[1] + q[1]*p[0] + q[2]*p[3] - q[3]*p[2],
      q[0]*p[2] - q[1]*p[3] + q[2]*p[0] + q[3]*p[1],
      q[0]*p[3] + q[1]*p[2] - q[2]*p[1] + q[3]*p[0]}};
  return product;
}

stk::math::Vector3d
Quaternion::rotation_vector() const
{
  // q and -q are the same rotation: take the one with q[0] >= 0, whose angle is at most pi.
  const double sign = (q[0] < 0.) ? -1. : 1.;
  const stk::math::Vector3d axisTimesSinHalfAngle(sign*q[1], sign*q[2], sign*q[3]);
  const double sinHalfAngle = axisTimesSinHalfAngle.length();
  if (sinHalfAngle == 0.)
    return stk::math::Vector3d::ZERO;
  const double angle = 2.*std::atan2(sinHalfAngle, sign*q[0]);
  return (angle/sinHalfAngle) * axisTimesSinHalfAngle;
}

RigidTransformation
Transformation::rigid_transformation_at(const double time) const
{
  Quaternion initialRotation;
  initialRotation.set_from_rotation_vector(myInitialRotation);
  Quaternion rotationSinceTimeZero;
  rotationSinceTimeZero.set_from_rotation_vector(time*my_rotational_velocity);
  const Quaternion rotation = rotationSinceTimeZero.after(initialRotation);
  const stk::math::Vector3d translation = my_reference_point + myInitialDisplacement + time*my_translational_velocity - rotation.rotate_3d_vector(my_reference_point);
  return RigidTransformation(rotation, translation);
}

} // namespace krino
