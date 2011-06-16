// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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

#ifndef IOSS_Ioss_EntityBlock_h
#define IOSS_Ioss_EntityBlock_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_GroupingEntity.h>
#include <string>

namespace Ioss {
  class ElementTopology;
  class ElementBlock;
  class DatabaseIO;

  class EntityBlock : public GroupingEntity {
  public:
    virtual Property
      get_implicit_property(const std::string& name) const = 0;

    // Describes the contained entities topology
    const ElementTopology *topology() const {return topology_;}

    // Describes the contained entities element block topology
    const ElementTopology *parent_element_topology() const
      { return parentTopology_;}

    // For faceblock, edgeblock, if they are split by element block, then this
    // will be non-NULL and is a pointer to the parent element block for this
    // faceblock or edgeblock. Has no meaning for other EntityBlock types or split
    // types.
    const ElementBlock *parent_element_block() const
      { return parentElementBlock_;}
    void set_parent_element_block(ElementBlock *element_block)
    { parentElementBlock_ = element_block; }
    
    // For faceblock, edgeblock, return whether the surface is applied
    // to the same face/edge for all elements in the surface. If not,
    // return 0; otherwise return the consistent face number.
    int  get_consistent_side_number() const;
    void set_consistent_side_number(int side) {consistentSideNumber = side;}

  protected:
    EntityBlock(const DatabaseIO *io_database,
		const std::string& name, const std::string& entity_type,
		const std::string& parent_topology_type,
		size_t entity_count);

  private:
    EntityBlock(const EntityBlock&); // do not implement
    EntityBlock& operator=(const EntityBlock&); // do not implement
    ElementTopology *topology_;
    ElementTopology *parentTopology_; // Topology of parent element (if any)
    ElementBlock    *parentElementBlock_;
  protected:
    mutable int consistentSideNumber;
  };
}
#endif
