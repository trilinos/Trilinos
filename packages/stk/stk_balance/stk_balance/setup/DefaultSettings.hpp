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

#ifndef STK_BALANCE_DEFAULT_SETTINGS_HPP
#define STK_BALANCE_DEFAULT_SETTINGS_HPP

#include <string>

namespace stk {
namespace balance {

enum class VertexWeightMethod
{
  CONSTANT = 0,
  TOPOLOGY,
  CONNECTIVITY,
  FIELD
};

std::string vertex_weight_method_name(VertexWeightMethod method);

struct DefaultSettings {
  static constexpr const char * logFile {"stk_balance.log"};
  static constexpr const char * outputDirectory {"."};
  static constexpr const char * decompMethod {"parmetis"};

  static constexpr bool useContactSearch {true};
  static constexpr bool fixMechanisms {false};

  static constexpr double faceSearchRelTol {0.15};
  static constexpr double faceSearchAbsTol {0.0001};

  static constexpr double particleSearchTol {3.0};

  static constexpr VertexWeightMethod vertexWeightMethod {VertexWeightMethod::CONSTANT};
  static constexpr double graphEdgeWeightMultiplier {1.0};
  static constexpr double faceSearchVertexMultiplier {1.0};
  static constexpr double faceSearchEdgeWeight {1.0};
  static constexpr bool fixSpiders {false};

  static constexpr VertexWeightMethod sdVertexWeightMethod {VertexWeightMethod::CONNECTIVITY};
  static constexpr double sdGraphEdgeWeightMultiplier {10.0};
  static constexpr double sdFaceSearchVertexMultiplier {2.0};
  static constexpr double sdFaceSearchEdgeWeight {1.0};
  static constexpr bool sdFixSpiders {true};

  static constexpr VertexWeightMethod smVertexWeightMethod {VertexWeightMethod::CONSTANT};
  static constexpr double smGraphEdgeWeightMultiplier {1.0};
  static constexpr double smFaceSearchVertexMultiplier {3.0};
  static constexpr double smFaceSearchEdgeWeight {1.0};
  static constexpr bool smFixSpiders {false};

  static constexpr const char * vertexWeightBlockMultiplier {""};
  static constexpr const char * cohesiveElements {""};

};

} }

#endif
