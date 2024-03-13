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

    void set_from_rotation_vector(const stk::math::Vector3d & v);
    stk::math::Vector3d rotate_3d_vector(const stk::math::Vector3d & v) const;
    stk::math::Vector3d reverse_rotate_3d_vector(const stk::math::Vector3d & v) const;

  private:
    std::array<double,4> q; ///< q[0] stores angle, (q[1],q[2],a[3]) stores axis.
};

class Transformation {
public:
  Transformation()
    : my_translational_velocity(stk::math::Vector3d::ZERO), my_rotational_velocity(stk::math::Vector3d::ZERO),
      my_reference_point(stk::math::Vector3d::ZERO), my_last_update(0.0), my_update_orientation(), my_update_offset(stk::math::Vector3d::ZERO) {}
  virtual ~Transformation() {}

  void set_translational_velocity(const stk::math::Vector3d & v) { my_translational_velocity = v; }
  void set_rotational_velocity(const stk::math::Vector3d & v) { my_rotational_velocity = v; }
  void set_reference_point(const stk::math::Vector3d & v) { my_reference_point = v; }

  // temporary storage until initialize()
  void set_initial_displacement(const stk::math::Vector3d & v) { my_update_offset = v; }
  void set_initial_rotation(const stk::math::Vector3d & v) { my_update_orientation.set_from_rotation_vector(v); }

  const stk::math::Vector3d & get_translational_velocity() const { return my_translational_velocity; }
  const stk::math::Vector3d & get_rotational_velocity() const { return my_rotational_velocity; }
  const stk::math::Vector3d & get_reference_point() const { return my_reference_point; }

  void initialize( const stk::math::Vector3d & initial_displacement, const stk::math::Vector3d & initial_rotation, const stk::math::Vector3d & reference_point );
  void initialize(); // this form assumes that set_initial_displacement() and set_initial_rotation() have been called or are zero
  void update( const double time ) const;
  void apply( stk::math::Vector3d & x0 ) const;

protected:
  stk::math::Vector3d my_translational_velocity;
  stk::math::Vector3d my_rotational_velocity;
  mutable stk::math::Vector3d my_reference_point;
  mutable double my_last_update;
  mutable Quaternion my_update_orientation;
  mutable stk::math::Vector3d my_update_offset;
};


} // namespace krino

#endif // Akri_Transformation_h
