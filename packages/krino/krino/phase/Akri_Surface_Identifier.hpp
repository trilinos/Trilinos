// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Surface_Identifier_h
#define Akri_Surface_Identifier_h

#include <ostream>

namespace krino {

class Surface_Identifier {
public:
  explicit Surface_Identifier(const unsigned id) : my_id(id) {}

  bool operator < ( const Surface_Identifier & RHS ) const { return my_id < RHS.my_id; }
  bool operator == ( const Surface_Identifier & RHS ) const { return my_id == RHS.my_id; }
  bool operator != ( const Surface_Identifier & RHS ) const { return my_id != RHS.my_id; }

  unsigned get() const { return my_id; }

  friend std::ostream& operator<<(std::ostream & os, const Surface_Identifier & id)
  {
    os << id.get();
    return os;
  }
private:
  unsigned my_id;
};

} // namespace krino

#endif // Akri_Surface_Identifier_h
