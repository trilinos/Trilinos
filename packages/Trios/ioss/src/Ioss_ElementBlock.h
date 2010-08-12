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

#ifndef IOSS_Ioss_ElementBlock_h
#define IOSS_Ioss_ElementBlock_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_EntityBlock.h>
#include <string>
#include <assert.h>

namespace Ioss {
  class DatabaseIO;

  class ElementBlock : public EntityBlock {
  public:
    ElementBlock(const DatabaseIO *io_database,
		 const std::string& name, const std::string& element_type,
		 size_t number_elements, size_t number_attributes);

    ~ElementBlock();

    std::string type_string() const {return "ElementBlock";}
    EntityType type() const {return ELEMENTBLOCK;}

    /// Handle implicit properties -- These are calcuated from data stored
    /// in the grouping entity instead of having an explicit value assigned.
    /// An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string& name) const;

    bool contains(size_t local_id) const
      { return idOffset < local_id && local_id <= idOffset + elementCount; }

    void  set_offset(int offset) {idOffset = offset;}
    int get_offset() const {return idOffset;}

    void get_block_adjacencies(std::vector<std::string> &block_adjacency_list) const;

  protected:
    int internal_get_field_data(const Field& field,
				void *data, size_t data_size) const;

    int internal_put_field_data(const Field& field,
				void *data, size_t data_size) const;
  private:
    void count_attributes() const;

    /**
     * The 'offset' is used to map an element location within an
     * element block to the element 'file descriptor'.
     * For example, the file descriptor of the 37th element in the 4th
     * block is calculated by:
     *
     * file_descriptor = offset of block 4 + 37
     *
     * This can also be used to determine which element block
     * an element with a file_descriptor maps into. An particular
     * element block contains all elements in the range:
     *
     * offset < file_descriptor <= offset+number_elements_per_block
     */
    size_t idOffset;
    
    size_t elementCount; ///< stored locally to avoid looking up property
    mutable size_t attributeCount;
  };
}
#endif
