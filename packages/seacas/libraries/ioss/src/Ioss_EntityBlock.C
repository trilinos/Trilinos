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

#include <Ioss_EntityBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_VariableType.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>

#include <string>
#include <assert.h>

Ioss::EntityBlock::EntityBlock(const Ioss::DatabaseIO *io_database,
			       const std::string& my_name,
			       const std::string& entity_type,
			       size_t entity_count)
  : Ioss::GroupingEntity(io_database, my_name, entity_count), 
    idOffset(0)
  
{
  topology_ = ElementTopology::factory(entity_type);
  assert(topology_ != NULL);

  if (topology()->master_element_name() != entity_type &&
      topology()->name() != entity_type) {
    // Maintain original element type on output database if possible.
    properties.add(Ioss::Property("original_topology_type", entity_type));
  }

  properties.add(Ioss::Property(this, "topology_node_count",
			       Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this, "topology_type",
			       Ioss::Property::STRING));
  fields.add(Ioss::Field("connectivity", Ioss::Field::INTEGER,
			topology_->name(),
			Ioss::Field::MESH, entity_count));

  // Returns connectivity in local id space
  fields.add(Ioss::Field("connectivity_raw", Ioss::Field::INTEGER,
			 topology()->name(),
			 Ioss::Field::MESH, entity_count));
}

Ioss::Property Ioss::EntityBlock::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "topology_node_count")
    return Ioss::Property(my_name, topology()->number_nodes());
  else if (my_name == "topology_type")
    return Ioss::Property(my_name, topology()->name());
  else
    return Ioss::GroupingEntity::get_implicit_property(my_name);
}

