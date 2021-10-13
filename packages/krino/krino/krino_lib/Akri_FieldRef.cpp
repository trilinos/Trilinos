// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_FieldRef.hpp>
#include <sstream>                      // for operator<<, stringstream, etc

namespace krino{

bool FieldRef::operator < ( const FieldRef & rhs ) const
{
  if ( my_field == rhs.my_field || rhs.my_field == NULL )
  {
    return false; // identical
  }
  else if (my_field == NULL)
  {
    return true;
  }
  else if ( my_field->field_state(stk::mesh::StateNone) == rhs.field().field_state(stk::mesh::StateNone) )
  {
    return my_field->state() < rhs.state();
  }
  return name() < rhs.name();
}

std::string
state_string(stk::mesh::FieldState state)
{
  static const char * const local_state_name[] = {
    "STATE_NONE/NEW/NP1",
    "STATE_OLD/N",
    "STATE_NM1",
    "STATE_NM2",
    "STATE_NM3",
    "STATE_NM4"
  };

  if ((unsigned) state < (unsigned)sizeof(local_state_name)/(unsigned)sizeof(local_state_name[0]))
    return local_state_name[state];
  else {
    std::stringstream strout;
    strout << "<state invalid>(" << (unsigned) state << ")";
    return strout.str().c_str();
  }
}

//----------------------------------------------------------------------

} // namespace krino
