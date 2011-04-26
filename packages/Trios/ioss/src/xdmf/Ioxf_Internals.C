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

#include <xdmf/Ioxf_Internals.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_SideBlock.h>
#include <algorithm>

#include <assert.h>

#include <string>

using namespace Ioxf;

Block::Block(const Ioss::ElementBlock &other)
{
  id = other.get_property("id").get_int();
  elementCount = other.get_property("entity_count").get_int();
  nodesPerElement = other.get_property("topology_node_count").get_int();
  attributeCount = other.get_property("attribute_count").get_int();
  offset_ = other.get_offset();
  std::string el_type = other.get_property("topology_type").get_string();
  if (other.property_exists("original_element_type")) {
    el_type = other.get_property("original_element_type").get_string();
  }

  std::strncpy(elType, el_type.c_str(), MAX_STR_LENGTH);
  elType[MAX_STR_LENGTH] = 0;

  // Fixup an exodusII kludge.  For triangular elements, the same
  // name is used for 2D elements and 3D shell elements.  Convert
  // to unambiguous names for the IO Subsystem.  The 2D name
  // stays the same, the 3D name becomes 'trishell#'
  // Here, we need to map back to the 'triangle' name...
  if (std::strncmp(elType, "trishell", 8) == 0)
    std::strncpy(elType, "triangle", 8);

  std::string io_name = other.name();
  std::strncpy(name, io_name.c_str(), MAX_STR_LENGTH);
  name[MAX_STR_LENGTH] = 0;
}

Block& Block::operator=(const Block& other)
{
  id = other.id;
  elementCount = other.elementCount;
  nodesPerElement = other.nodesPerElement;
  attributeCount = other.attributeCount;
  offset_ = other.offset_;
  std::strcpy(elType, other.elType);
  std::strcpy(name,   other.name);
  return *this;
}

bool Block::operator==(const Block& other) const
{
  return id == other.id &&
    elementCount == other.elementCount &&
    nodesPerElement == other.nodesPerElement &&
    attributeCount == other.attributeCount &&
    std::strcmp(elType, other.elType) == 0 &&
    std::strcmp(name, other.name) == 0;
}

NodeSet::NodeSet(const Ioss::NodeSet &other)
{
  id = other.get_property("id").get_int();
  nodeCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();

  std::string io_name = other.name();
  std::strncpy(name, io_name.c_str(), MAX_STR_LENGTH);
  name[MAX_STR_LENGTH] = 0;
}

bool NodeSet::operator==(const NodeSet& other) const
{
  return id == other.id &&
    nodeCount == other.nodeCount &&
    dfCount == other.dfCount &&
    std::strcmp(name, other.name) == 0;
}

SideSet::SideSet(const Ioss::SideBlock &other)
{
  id = other.get_property("id").get_int();
  sideCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  std::string io_name = other.name();
  std::strncpy(name, io_name.c_str(), MAX_STR_LENGTH);
  name[MAX_STR_LENGTH] = 0;
}

bool SideSet::operator==(const SideSet& other) const
{
  return id == other.id &&
    sideCount == other.sideCount &&
    dfCount == other.dfCount &&
    std::strcmp(name, other.name) == 0;
}

bool CommunicationMap::operator==(const CommunicationMap& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    type == other.type;
}
