// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_UPDATE_TYPE_HPP
#define ROL_UPDATE_TYPE_HPP

#include <cstdint>
#include <string>

namespace ROL {

enum class UpdateType : std::uint8_t {
  Initial = 0, // Update has not been called before
  Accept,      // This is the new iterate, trial must be called before
  Revert,      // Revert to the previous iterate, trial must be called before
  Trial,       // This is a candidate for the next iterate
  Temp         // For temporary uses including finite difference computations
};


inline std::string UpdateTypeToString(const UpdateType& type) {
  std::string retString;
  switch(type) {
    case UpdateType::Initial:  retString = "Initial";  break;
    case UpdateType::Accept:   retString = "Accept";   break;
    case UpdateType::Revert:   retString = "Revert";   break;
    case UpdateType::Trial:    retString = "Trial";    break;
    case UpdateType::Temp:     retString = "Temp";     break;
  }
  return retString;
}

} // namespace ROL

#endif
