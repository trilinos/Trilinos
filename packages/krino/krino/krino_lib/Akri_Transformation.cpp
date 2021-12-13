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
Quaternion::set_from_rotation_vector(const Vector3d & v)
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

Vector3d use_quaternion_to_rotate_3d_vector(const std::array<double,4> & q, const Vector3d & v)
{
  return Vector3d(
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


Vector3d
Quaternion::rotate_3d_vector(const Vector3d & v) const
{
  return use_quaternion_to_rotate_3d_vector(q, v);
}

Vector3d
Quaternion::reverse_rotate_3d_vector(const Vector3d & v) const
{
  return use_quaternion_to_rotate_3d_vector(std::array<double,4>{q[0],-q[1],-q[2],-q[3]}, v);
}

void
Transformation::initialize( const Vector3d & initial_displacement, const Vector3d & initial_rotation, const Vector3d & reference_point )
{
  my_reference_point = reference_point;
  my_update_orientation.set_from_rotation_vector(initial_rotation);
  my_update_offset = initial_displacement - my_update_orientation.rotate_3d_vector(my_reference_point);
  my_reference_point += initial_displacement;
}

void
Transformation::initialize()
{
  // this form of initialization assumes that set_initial_displacement() and set_initial_rotation() have been called or are zero
  const Vector3d initial_displacement = my_update_offset;
  my_update_offset = initial_displacement - my_update_orientation.rotate_3d_vector(my_reference_point);
  my_reference_point += initial_displacement;
}

void
Transformation::update( const double time ) const
{
  if (time == my_last_update)
  {
    return;
  }
  if (my_last_update < 0.0)
  {
    my_last_update = time;
    return;
  }

  const double dt = time - my_last_update;
  const Vector3d update_rotation_angle = dt*my_rotational_velocity;
  my_update_orientation.set_from_rotation_vector(update_rotation_angle);
  my_update_offset = my_reference_point - my_update_orientation.rotate_3d_vector(my_reference_point) + my_translational_velocity * dt;
  my_reference_point += my_translational_velocity * dt;

  my_last_update = time;
}

void
Transformation::apply( Vector3d & x0 ) const
{
  x0 = my_update_orientation.rotate_3d_vector(x0) + my_update_offset;
}

} // namespace krino
