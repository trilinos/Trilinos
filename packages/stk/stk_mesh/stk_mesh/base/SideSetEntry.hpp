// Copyright (c) 2013, Sandia Corporation.
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


#ifndef SIDESETENTRY_HPP_
#define SIDESETENTRY_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk
{
namespace mesh
{

struct SideSetEntry
{
  SideSetEntry() : element(stk::mesh::Entity()), side(stk::mesh::INVALID_CONNECTIVITY_ORDINAL){};
  SideSetEntry(stk::mesh::Entity in_element, stk::mesh::ConnectivityOrdinal in_side)
    : element(in_element),
      side(in_side)
  {  }
  SideSetEntry(stk::mesh::Entity in_element, int in_side)
    : SideSetEntry(in_element, static_cast<stk::mesh::ConnectivityOrdinal>(in_side))
  {  }

  bool operator==(const SideSetEntry &rhs) const
  {
      return ((element == rhs.element) && (side == rhs.side));
  }

  bool operator<(const SideSetEntry &rhs) const
  {
      if(element < rhs.element)
          return true;
      else if (element == rhs.element && side < rhs.side)
          return true;
      else return false;
  }

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal side;
};

typedef std::vector<SideSetEntry> SideSet;
typedef std::vector<SideSet*> SideSetVector;

}
}

#endif /* SIDESETENTRY_HPP_ */
