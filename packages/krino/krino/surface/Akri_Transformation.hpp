// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Transformation_h
#define Akri_Transformation_h

#include <array>
#include <stk_math/StkVector.hpp>

namespace krino {

class Quaternion {
  public:
    Quaternion() : q{{ 1.0, 0.0, 0.0, 0.0 }} {}

    // The rotation by the angle |v|, in radians, about the axis v, by the right hand rule: a rotation
    // vector of (0,0,pi/2) turns x to y.
    void set_from_rotation_vector(const stk::math::Vector3d & v);
    stk::math::Vector3d rotate_3d_vector(const stk::math::Vector3d & v) const;
    // The inverse rotation.
    stk::math::Vector3d reverse_rotate_3d_vector(const stk::math::Vector3d & v) const;
    // The rotation by first, then by this one.
    Quaternion after(const Quaternion & first) const;
    // The rotation vector of this rotation, with an angle from 0 to pi (set_from_rotation_vector
    // gives the same rotation from it).
    stk::math::Vector3d rotation_vector() const;

  private:
    std::array<double,4> q; ///< q[0] stores angle, (q[1],q[2],a[3]) stores axis.
};

// The rigid transformation x -> R x + t.
class RigidTransformation {
public:
  RigidTransformation() = default;
  RigidTransformation(const Quaternion & rotation, const stk::math::Vector3d & translation)
  : myRotation(rotation), myTranslation(translation) {}

  stk::math::Vector3d transform_point(const stk::math::Vector3d & x) const { return myRotation.rotate_3d_vector(x) + myTranslation; }
  stk::math::Vector3d inverse_transform_point(const stk::math::Vector3d & x) const { return myRotation.reverse_rotate_3d_vector(x - myTranslation); }
  stk::math::Vector3d transform_vector(const stk::math::Vector3d & v) const { return myRotation.rotate_3d_vector(v); }
  stk::math::Vector3d inverse_transform_vector(const stk::math::Vector3d & v) const { return myRotation.reverse_rotate_3d_vector(v); }

private:
  Quaternion myRotation;
  stk::math::Vector3d myTranslation{stk::math::Vector3d::ZERO};
};

// A motion of a body, given at time 0, with constant velocities: at time t its point x is at
//   R(rotationalVelocity t) R(initialRotation) (x - referencePoint) + referencePoint + initialDisplacement + translationalVelocity t,
// so it is turned by its initial rotation about the reference point and displaced, then turns about
// the reference point as it moves with it.  The rotations are rotation vectors (see Quaternion).
class Transformation {
public:
  void set_translational_velocity(const stk::math::Vector3d & v) { my_translational_velocity = v; }
  void set_rotational_velocity(const stk::math::Vector3d & v) { my_rotational_velocity = v; }
  void set_reference_point(const stk::math::Vector3d & v) { my_reference_point = v; }

  void set_initial_displacement(const stk::math::Vector3d & v) { myInitialDisplacement = v; }
  void set_initial_rotation(const stk::math::Vector3d & v) { myInitialRotation = v; }

  const stk::math::Vector3d & get_translational_velocity() const { return my_translational_velocity; }
  const stk::math::Vector3d & get_rotational_velocity() const { return my_rotational_velocity; }
  const stk::math::Vector3d & get_reference_point() const { return my_reference_point; }
  const stk::math::Vector3d & get_initial_displacement() const { return myInitialDisplacement; }
  const stk::math::Vector3d & get_initial_rotation() const { return myInitialRotation; }

  // Where the motion puts the body at time.
  RigidTransformation rigid_transformation_at(const double time) const;

private:
  stk::math::Vector3d my_translational_velocity{stk::math::Vector3d::ZERO};
  stk::math::Vector3d my_rotational_velocity{stk::math::Vector3d::ZERO};
  stk::math::Vector3d my_reference_point{stk::math::Vector3d::ZERO};
  stk::math::Vector3d myInitialDisplacement{stk::math::Vector3d::ZERO};
  stk::math::Vector3d myInitialRotation{stk::math::Vector3d::ZERO};
};


} // namespace krino

#endif // Akri_Transformation_h
