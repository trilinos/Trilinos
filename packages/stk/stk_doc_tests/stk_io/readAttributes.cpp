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
//

#include <gtest/gtest.h>                // for AssertHelper, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <iomanip>                      // for operator<<
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_io/AttributeFields.hpp"

namespace {

class ExodusFileWithAttributes : public stk::unit_test_util::MeshFixture { };

stk::mesh::FieldVector get_attribute_fields_for_part(const stk::mesh::MetaData &meta, const stk::mesh::Part *ioPart)
{
  stk::mesh::FieldVector attributes;

  for(stk::mesh::FieldBase *field : meta.get_fields())
  {
    const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
    if(fieldRole != nullptr && *fieldRole == Ioss::Field::ATTRIBUTE)
    {
      for(const stk::mesh::FieldBase::Restriction &restriction : field->restrictions())
      {
        const stk::mesh::Selector &selector = restriction.selector();
        if(selector(ioPart))
          attributes.push_back(field);
      }
    }
  }
  return attributes;
}

//-BEGIN
std::vector<double> get_attributes_of_first_element(const stk::mesh::BulkData &bulk, const stk::mesh::Part *ioPart)
{
  stk::mesh::FieldVector attributeFields = get_attribute_fields_for_part(bulk.mesh_meta_data(), ioPart);

  stk::mesh::EntityVector elements;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, *ioPart, elements);

  std::vector<double> attributes;
  if(!elements.empty()) {
    for(const stk::mesh::FieldBase *field : attributeFields) {
      unsigned numAttribute = stk::mesh::field_scalars_per_entity(*field, elements[0]);
      double *dataForElement = static_cast<double*> (stk::mesh::field_data(*field, elements[0]));
      for(unsigned i=0; i<numAttribute; ++i)
        attributes.push_back(dataForElement[i]);
    }
  }
  return attributes;
}

TEST_F(ExodusFileWithAttributes, readAttributes_haveFieldsWithAttributes)
{
  setup_mesh("hex_spider.exo", stk::mesh::BulkData::AUTO_AURA);

  const stk::mesh::Part *partBlock2  = get_meta().get_part("block_2");
  const stk::mesh::Part *partBlock10 = get_meta().get_part("block_10");

  EXPECT_EQ(1u, get_attributes_of_first_element(get_bulk(), partBlock2).size());
  EXPECT_EQ(7u, get_attributes_of_first_element(get_bulk(), partBlock10).size());
}

void mark_field_as_attribute(stk::mesh::FieldBase &field)
{
  stk::io::set_field_role(field, Ioss::Field::ATTRIBUTE);
}

TEST_F(ExodusFileWithAttributes, addAttribute_haveFieldsWithAttribute)
{
  allocate_bulk(stk::mesh::BulkData::AUTO_AURA);

  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(get_bulk());
  stkIo.add_mesh_database("hex_spider.exo", stk::io::READ_MESH);
  stkIo.create_input_mesh();

  double initialValue = 0.0;
  auto &newAttrField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "newAttr");
  mark_field_as_attribute(newAttrField);

  const stk::mesh::Part *partBlock10 = get_meta().get_part("block_10");
  stk::mesh::put_field_on_mesh(newAttrField, *partBlock10, &initialValue);

  stkIo.populate_bulk_data();

  EXPECT_EQ(8u, get_attributes_of_first_element(get_bulk(), partBlock10).size());
}
//-END
}
