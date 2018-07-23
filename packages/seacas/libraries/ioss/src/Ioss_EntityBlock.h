// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//     * Neither the name of NTESS nor the names of its
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

#include <Ioss_GroupingEntity.h> // for GroupingEntity
#include <Ioss_Property.h>       // for Property
#include <cstddef>               // for size_t
#include <string>                // for string
namespace Ioss {
  class DatabaseIO;
} // namespace Ioss
namespace Ioss {
  class ElementTopology;
} // namespace Ioss

namespace Ioss {
  class ElementBlock;

  /** \brief Base class for all 'block'-type grouping entities, which means all
   *         members of the block are similar or have the same topology.
   *
   *   The following derived classes are typical:
   *
   *   -- NodeBlock -- grouping of 'similar' nodes (same degree of freedom, ...)
   *
   *   -- ElementBlock -- grouping of 'similar' elements (same element topology,
   *                      attributes, ...)
   *      0d, 1d, 2d, 3d topology possible -- e.g., sphere, bar, quad, hex
   */
  class EntityBlock : public GroupingEntity
  {
  public:
    Property get_implicit_property(const std::string &my_name) const override = 0;

    /** \brief Get the topology of the entities in the block.
     *
     *  \returns The topology.
     */
    const ElementTopology *topology() const { return topology_; }

    /** \brief Determine whether the block contains the entity with a given id.
     *
     *  \param[in] local_id The id to check.
     *  \returns True if the block contains the entity.
     */
    bool contains(size_t local_id) const
    {
      return idOffset < local_id && local_id <= idOffset + entityCount;
    }
    /** \brief Set the 'offset' for the block.
     *
     *  The 'offset' is used to map an element location within an
     *  element block to the element 'file descriptor'.
     *  For example, the file descriptor of the 37th element in the 4th
     *  block is calculated by:
     *
     *  file_descriptor = offset of block 4 + 37
     *
     *  This can also be used to determine which element block
     *  an element with a file_descriptor maps into. An particular
     *  element block contains all elements in the range:
     *
     *  offset < file_descriptor <= offset+number_elements_per_block
     */
    void set_offset(size_t offset) { idOffset = offset; }

    /** \brief Get the 'offset' for the block.
     *
     *  The 'offset' is used to map an element location within an
     *  element block to the element 'file descriptor'.
     *  For example, the file descriptor of the 37th element in the 4th
     *  block is calculated by:
     *
     *  file_descriptor = offset of block 4 + 37
     *
     *  This can also be used to determine which element block
     *  an element with a file_descriptor maps into. An particular
     *  element block contains all elements in the range:
     *
     *  offset < file_descriptor <= offset+number_elements_per_block
     */
    size_t get_offset() const { return idOffset; }

  protected:
    EntityBlock(DatabaseIO *io_database, const std::string &my_name, const std::string &entity_type,
                size_t entity_cnt);

    EntityBlock(const EntityBlock &) = delete;
    EntityBlock &operator=(const EntityBlock &) = delete;

    ElementTopology *topology_;

  protected:
    size_t idOffset;
  };
} // namespace Ioss
#endif
