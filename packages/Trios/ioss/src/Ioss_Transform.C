/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Transform.h>
#include <Ioss_Field.h>
#include <vector>
#include <string>

namespace Ioss {

  Transform::Transform() {}
  Transform::~Transform() {}

  bool Transform::execute(const Ioss::Field &field, void *data)
  {
    return internal_execute(field, data);
  }

  void Transform::set_property(const std::string&, int) {}
  void Transform::set_property(const std::string&, double) {}
  void Transform::set_properties(const std::string&, const std::vector<int>&) {}
  void Transform::set_properties(const std::string&, const std::vector<double>&) {}
}
