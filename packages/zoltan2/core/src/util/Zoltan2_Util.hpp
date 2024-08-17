// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_Util.hpp
 *  \brief A gathering of useful namespace methods.
 *  \todo Should each class of utility functions be in a separate source file
 *         instead of having a source file with the unhelpful name of Util?
 */

#ifndef ZOLTAN2_UTIL_HPP
#define ZOLTAN2_UTIL_HPP

#include <Zoltan2_Standards.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <stdexcept>

namespace Zoltan2{

long getProcessKilobytes();

template <typename scalar_t>
  inline bool outsideRegion(scalar_t val, scalar_t mark, double epsilon){
    return ((val < mark-epsilon) || (val > mark+epsilon));
}

static inline void
AssertCondition(bool condition, const std::string &message,
                const char *file = __FILE__, int line = __LINE__) {
    if (!condition) {
      std::ostringstream eMsg;
      eMsg << "Error: " << file << ", " << line << ": " << message
                << std::endl;
      throw std::runtime_error(eMsg.str());
    }
}

} // namespace Zoltan2

#endif
