//@HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
//
//                Shards : Shared Discretization Tools
//            The Array Dims are a copy of shards array dims.
//
// Copyright 2008 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include <stdlib.h>
#include <stdexcept>
#include <sstream>

#include "Phalanx_DimTag.hpp"

namespace PHX {
  
  DimTag::~DimTag() {}
  
  std::string
  DimTag::to_string(DimTag::size_type /* n */, DimTag::size_type i) const
  {
    //array_traits::check_range( i , n );
    std::ostringstream s ;
    s << i ;
    return s.str();
  }

  DimTag::size_type
  DimTag::to_index( DimTag::size_type /* n */, const std::string & s ) const
  {
    const int i = atoi( s.c_str() );
    //array_traits::check_range( i , n );
    return i ;
  }
  
}
