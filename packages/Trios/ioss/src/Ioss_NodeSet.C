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

#include <Ioss_NodeSet.h>

#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_Utils.h>
#include <string>

namespace {
  const std::string SCALAR() { return std::string("scalar");}
}

Ioss::NodeSet::NodeSet(const Ioss::DatabaseIO *io_database, const std::string& my_name,
		       size_t number_nodes)
  : Ioss::GroupingEntity(io_database, my_name)
{
  properties.add(Ioss::Property("entity_count", static_cast<int>(number_nodes)));
  properties.add(Ioss::Property("distribution_factor_count",
				static_cast<int>(number_nodes)));
      // Add the standard fields...
  fields.add(Ioss::Field("ids", Ioss::Field::INTEGER, SCALAR(),
			 Ioss::Field::MESH, number_nodes));

  fields.add(Ioss::Field("distribution_factors",
			 Ioss::Field::REAL, SCALAR(),
			 Ioss::Field::MESH, number_nodes));

}

int Ioss::NodeSet::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::NodeSet::internal_put_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property Ioss::NodeSet::get_implicit_property(const std::string& my_name) const
{
  return Ioss::GroupingEntity::get_implicit_property(my_name);
}

void Ioss::NodeSet::block_membership(std::vector<std::string> &block_members)
{
  block_members.push_back("nodeblock_1");
}
