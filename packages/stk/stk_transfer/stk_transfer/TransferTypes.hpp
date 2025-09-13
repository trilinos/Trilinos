// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_STK_TRANSFER_STK_TRANSFER_TRANSFERTYPES_HPP_
#define STK_STK_TRANSFER_STK_TRANSFER_TRANSFERTYPES_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <functional>                    // for function
#include <limits>                        // for numeric_limits
#include <string>                        // for string, operator==, basic_st...
#include <utility>                       // for pair
#include <vector>                        // for vector
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

namespace impl {
inline double default_transform(double value) { return value; }
}

using FieldTransform = double (*)(double);

struct InterpolationData {
  unsigned nFields;
  std::vector<double*> fieldPtr;
  std::vector<int> fieldSize;
  std::vector<unsigned> fieldKey;
  std::vector<unsigned> fieldDataIndex;
  bool debug{false};

  InterpolationData(unsigned maxSize)
    : nFields(0)
    , fieldPtr(maxSize)
    , fieldSize(maxSize)
    , fieldKey(maxSize)
    , fieldDataIndex(maxSize)
    , debug(false)
  {
  }
};

} // namespace transfer
} // namespace stk

#endif /* STK_STK_TRANSFER_STK_TRANSFER_TRANSFERTYPES_HPP_ */
