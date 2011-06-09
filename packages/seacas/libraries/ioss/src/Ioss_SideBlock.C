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

#include <Ioss_SideBlock.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_Utils.h>
#include <Ioss_ElementTopology.h>
#include <string>

Ioss::SideBlock::SideBlock(const Ioss::DatabaseIO *io_database,
			   const std::string& my_name,
			   const std::string& side_type,
			   const std::string& element_type,
			   size_t side_count)
  : Ioss::EntityBlock(io_database, my_name, side_type, element_type, side_count), owner_(NULL)
{
  properties.add(Ioss::Property(this, "distribution_factor_count",
				Ioss::Property::INTEGER));
  
  fields.add(Ioss::Field("element_side",
			 Ioss::Field::INTEGER, "pair",
			 Ioss::Field::MESH, side_count));

  // Same as element_side except that the element id are the local
  // element position (1-based) and not the global element id.
  fields.add(Ioss::Field("element_side_raw",
			 Ioss::Field::INTEGER, "pair",
			 Ioss::Field::MESH, side_count));
  
  // Distribution factors are optional...
}

int Ioss::SideBlock::internal_get_field_data(const Ioss::Field& field,
					     void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::SideBlock::internal_put_field_data(const Ioss::Field& field,
					     void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::SideBlock::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "distribution_factor_count") {
    if (field_exists("distribution_factors")) {
      int nnodes = topology()->number_nodes();
      int nside  = get_property("entity_count").get_int();
      return Ioss::Property(my_name, nnodes*nside);
    } else {
      return Ioss::Property(my_name, 0);
    }
  }
  else
    return Ioss::EntityBlock::get_implicit_property(my_name);
}

void Ioss::SideBlock::block_membership(std::vector<std::string> &block_members)
{
  // Simplest case.  If the surfaces are split by element block, then this will return non-null
  // and we are done.
  const Ioss::ElementBlock *eb = parent_element_block();
  if (eb != NULL) {
    block_members.push_back(eb->name());
    return;
  }

  if (blockMembership.empty()) {
    get_database()->compute_block_membership(this, blockMembership);
  } 
  block_members = blockMembership;
}
