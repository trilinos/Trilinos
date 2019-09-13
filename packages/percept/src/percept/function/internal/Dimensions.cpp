// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/function/internal/Dimensions.hpp>

  namespace percept
  {

    std::ostream &operator<<(std::ostream& out, const Dimensions& dim)
    {
      //out << "Dimensions: ";
      out << "[ ";
      for (unsigned int i = 0; i < dim.size(); i++)
        {
          out << " " << dim[i];
        }
      out << "] ";
      return out;
    }

  }
