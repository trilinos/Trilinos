// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_ApplyType.hpp"
#include <stdexcept>

namespace TSQR {
  ApplyType::ApplyType (const std::string& op) :
    type_ (decide_apply_type (op)),
    lapackString_ (enumToLapackString (type_))
  {}

  ApplyType::ApplyType (const ApplyType& rhs) :
    type_ (rhs.type_),
    lapackString_ (rhs.lapackString_)
  {}

  ApplyType& ApplyType::operator= (const ApplyType& rhs) {
    type_ = rhs.type_;
    lapackString_ = rhs.lapackString_;
    return *this;
  }

  const ApplyType ApplyType::NoTranspose = ApplyType ("N");
  const ApplyType ApplyType::Transpose = ApplyType ("T");
  const ApplyType ApplyType::ConjugateTranspose = ApplyType ("C");

  std::string
  ApplyType::enumToLapackString (const ApplyType::ApplyType_ theType)
  {
    if (theType == NoTranspose_)
      return std::string("N");
    else if (theType == Transpose_)
      return std::string("T");
    else if (theType == ConjugateTranspose_)
      return std::string("C");
    else
      throw std::logic_error("Invalid ApplyType: should never get here");
  }

  bool
  ApplyType::decide_transposed (const std::string& op) const
  {
    if (op[0] == 'N' || op[0] == 'n')
      return false;
    else
      {
        const char validTransposeLetters[] = "TtCcHh";
        const int numValidTransposeLetters = 6;

        for (int k = 0; k < numValidTransposeLetters; ++k)
          if (op[0] == validTransposeLetters[k])
            return true;

        throw std::invalid_argument ("Invalid \"op\" argument \"" + op + "\"");
      }
  }

  ApplyType::ApplyType_
  ApplyType::decide_apply_type (const std::string& op) const
  {
    if (op[0] == 'T' || op[0] == 't')
      return Transpose_;
    else if (op[0] == 'N' || op[0] == 'n')
      return NoTranspose_;
    else if (op[0] == 'C' || op[0] == 'c' || op[0] == 'H' || op[0] == 'h')
      return ConjugateTranspose_;
    else
      throw std::invalid_argument ("Invalid \"op\" argument \"" + op + "\"");
  }
}

