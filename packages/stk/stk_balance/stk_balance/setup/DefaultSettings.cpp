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

#include "DefaultSettings.hpp"

namespace stk {
namespace balance {

std::string vertex_weight_method_name(VertexWeightMethod method) {
  switch (method) {
  case (VertexWeightMethod::CONSTANT) : {
    return "constant";
  }
  case (VertexWeightMethod::TOPOLOGY) : {
    return "topology";
  }
  case (VertexWeightMethod::CONNECTIVITY) : {
    return "connectivity";
  }
  case (VertexWeightMethod::FIELD) : {
    return "field";
  }
  default: {
    return "unknown";
  }
  }
}

// Declaration of static class members does not count as a definition,
// so they must be initialized outside the class.  However, if they are also
// constexpr, then they must be initialized at declaration time.  Defining
// them as empty here works around this conflict, preserving their initialization
// at declaration time.
//
// In C++17, these external definitions are no longer necessary, meaning that
// this entire file can be deleted.

constexpr const char * DefaultSettings::logFile;
constexpr const char * DefaultSettings::outputDirectory;
constexpr const char * DefaultSettings::decompMethod;

constexpr bool DefaultSettings::useContactSearch;
constexpr bool DefaultSettings::fixMechanisms;

constexpr double DefaultSettings::faceSearchRelTol;
constexpr double DefaultSettings::faceSearchAbsTol;

constexpr double DefaultSettings::particleSearchTol;

constexpr VertexWeightMethod DefaultSettings::vertexWeightMethod;
constexpr double DefaultSettings::graphEdgeWeightMultiplier;
constexpr double DefaultSettings::faceSearchVertexMultiplier;
constexpr double DefaultSettings::faceSearchEdgeWeight;
constexpr bool DefaultSettings::fixSpiders;

constexpr VertexWeightMethod DefaultSettings::sdVertexWeightMethod;
constexpr double DefaultSettings::sdGraphEdgeWeightMultiplier;
constexpr double DefaultSettings::sdFaceSearchVertexMultiplier;
constexpr double DefaultSettings::sdFaceSearchEdgeWeight;
constexpr bool DefaultSettings::sdFixSpiders;

constexpr VertexWeightMethod DefaultSettings::smVertexWeightMethod;
constexpr double DefaultSettings::smGraphEdgeWeightMultiplier;
constexpr double DefaultSettings::smFaceSearchVertexMultiplier;
constexpr double DefaultSettings::smFaceSearchEdgeWeight;
constexpr bool DefaultSettings::smFixSpiders;

constexpr const char * DefaultSettings::vertexWeightBlockMultiplier;
constexpr const char * DefaultSettings::cohesiveElements;
}
}
