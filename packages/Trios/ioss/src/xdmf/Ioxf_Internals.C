// Copyright 2001 Sandia Corporation, Albuquerque, NM.

#include <xdmf/Ioxf_Internals.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_EdgeBlock.h>
#include <Ioss_FaceBlock.h>
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

SideSet::SideSet(const Ioss::FaceBlock &other)
{
  id = other.get_property("id").get_int();
  sideCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  std::string io_name = other.name();
  std::strncpy(name, io_name.c_str(), MAX_STR_LENGTH);
  name[MAX_STR_LENGTH] = 0;
}

SideSet::SideSet(const Ioss::EdgeBlock &other)
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
