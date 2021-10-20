// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SIERRA_Akri_TypeDefs_h
#define SIERRA_Akri_TypeDefs_h

#include <stk_util/util/Array.hpp>
#include <Akri_Vec.hpp>
#include <vector>

namespace krino {

  class NINT{};
  class NPE_VAR{};
  class NPE_COORD{};
  class DIM{};
  
  typedef std::vector< Vector3d > PointVec;

}

#endif // SIERRA_Akri_TypeDefs_h
