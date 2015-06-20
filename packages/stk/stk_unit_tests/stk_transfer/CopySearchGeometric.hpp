// Copyright (c) 2015, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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


#ifndef  STK_COPYSEARCHGEOMETRIC_HPP
#define  STK_COPYSEARCHGEOMETRIC_HPP

#include "CopySearchBase.hpp"
#include <stk_search/IdentProc.hpp>
#include "CopyTransferMeshBase.hpp"

namespace stk {
namespace transfer {

class CopySearchGeometric : public CopySearchBase {
public:
  CopySearchGeometric() : m_radius(1.0e-6) {}
  virtual ~CopySearchGeometric() {}
  virtual void intialize(const CopyTransferMeshBase & mesha, const CopyTransferMeshBase & meshb) {}
  virtual void do_search(const CopyTransferMeshBase & mesha, const CopyTransferMeshBase & meshb, KeyToTargetProcessor & key_to_target_processor);
  virtual const MeshIDSet & get_remote_keys() const { return m_remote_keys; }
  void set_bounding_box_radius(float radius_in) { m_radius = radius_in; }
  float get_bounding_box_radius() const { return m_radius; }

private:
  MeshIDSet m_remote_keys;
  float m_radius;
};

template< typename BoundingBox>
struct BoundingBoxCompare {
    bool operator()(const BoundingBox & a, const BoundingBox & b) const
    {
      return a.second.id() < b.second.id();
    }
};

typedef stk::search::IdentProc<CopyTransferMeshBase::Mesh_ID> MeshIDProc;

}  } // namespace transfer stk

#endif //  STK_COPYSEARCHGEOMETRIC_HPP
