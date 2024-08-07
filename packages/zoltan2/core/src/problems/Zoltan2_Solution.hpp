// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_Solution.hpp
    \brief Defines the Solution base class.
    \todo The Solution classes are part of the API.  This is problematic
                right now - they use many internal classes and are used
                by Problems and algorithms.  Maybe there
                should be a SolutionAPI object within the Solution which
                is the visible part of the Solution.
           The source for this class could be
                 in the src/input directory.  This is where the input
                 adapters which use that part of the solution reside.
*/

#ifndef _ZOLTAN2_SOLUTION_HPP_
#define _ZOLTAN2_SOLUTION_HPP_

#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! \brief Just a placeholder for now.
*/

class Solution
{
public:
  virtual ~Solution() {}
};

}

#endif
