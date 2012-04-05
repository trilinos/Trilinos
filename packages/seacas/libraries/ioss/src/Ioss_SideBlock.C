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

#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_EntityBlock.h>
#include <Ioss_Field.h>
#include <Ioss_Property.h>
#include <Ioss_SideBlock.h>
#include <assert.h>
#include <stddef.h>
#include <string>
#include <vector>

#include "Ioss_FieldManager.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_PropertyManager.h"

Ioss::SideBlock::SideBlock(Ioss::DatabaseIO *io_database,
			   const std::string& my_name,
			   const std::string& side_type,
			   const std::string& element_type,
			   size_t side_count)
  : Ioss::EntityBlock(io_database, my_name, side_type, side_count), owner_(NULL),
    parentTopology_(NULL), parentElementBlock_(NULL), consistentSideNumber(-1)
{
  parentTopology_ = ElementTopology::factory(element_type);
  assert(parentTopology_ != NULL);

  properties.add(Ioss::Property(this, "parent_topology_type",
			       Ioss::Property::STRING));

  properties.add(Ioss::Property(this, "distribution_factor_count",
				Ioss::Property::INTEGER));
  
  fields.add(Ioss::Field("element_side",
			 field_int_type(), "pair",
			 Ioss::Field::MESH, side_count));

  // Same as element_side except that the element id are the local
  // element position (1-based) and not the global element id.
  fields.add(Ioss::Field("element_side_raw",
			 field_int_type(), "pair",
			 Ioss::Field::MESH, side_count));
  
  // Distribution factors are optional...
}

int64_t Ioss::SideBlock::internal_get_field_data(const Ioss::Field& field,
					     void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int64_t Ioss::SideBlock::internal_put_field_data(const Ioss::Field& field,
					     void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::SideBlock::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "distribution_factor_count") {
    if (field_exists("distribution_factors")) {
      int64_t nnodes = topology()->number_nodes();
      int64_t nside  = get_property("entity_count").get_int();
      return Ioss::Property(my_name, nnodes*nside);
    } else {
      return Ioss::Property(my_name, 0);
    }
  }
  else if (my_name == "parent_topology_type")
    return Ioss::Property(my_name, parent_element_topology()->name());
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

namespace {
  template <typename INT>
  int internal_consistent_side_number(std::vector<INT> &element_side)
  {
    size_t ecount = element_side.size();
    int side = ecount > 0 ? element_side[1] : 0;
    for (size_t i=3; i < ecount; i+=2) {
      int this_side = element_side[i];
      if (this_side != side) {
	side = 999; // Indicates the sides are not consistent ;
	break;
      }
    }
    return side;
  }
}

int  Ioss::SideBlock::get_consistent_side_number() const
{
  if (consistentSideNumber == -1) {
    // It wasn't calculated during the metadata reading of the surfaces.
    // Determine it now...
    if (field_exists("element_side")) {
      int side = 0;
      if (get_database()->int_byte_size_api() == 8) {
	std::vector<int64_t> element_side;
	get_field_data("element_side", element_side);
	side = internal_consistent_side_number(element_side);
      } else {
	std::vector<int> element_side;
	get_field_data("element_side", element_side);
	side = internal_consistent_side_number(element_side);
      }

      int side_max = get_database()->util().global_minmax(side, Ioss::ParallelUtils::DO_MAX);
      if (side_max != 999)
	consistentSideNumber = side_max;
      else
	consistentSideNumber = 0;
    } else {
      consistentSideNumber = 0;
    }
  }
  return consistentSideNumber;
}
