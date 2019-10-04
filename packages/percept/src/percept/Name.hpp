// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Name_hpp
#define percept_Name_hpp

#include <string>

  namespace percept
  {

    /// this is to avoid a common bug where the name of the String Function is given instead of function_string,
    ///    ie., first two args are accidentally reversed - this is essentially a model of a "named argument",
    ///    as opposed to a positional one; of course, it is still only a hint (though a strong one) to the user
    
    /// Useful in other places where two strings are passed into a function or constructor
    class Name 
    {
      const std::string m_name;
    public:
      explicit 
      Name(const std::string name) : m_name(name) {}
      const std::string& getName() const { return m_name; }
    };

  }

#endif
