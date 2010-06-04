/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_EdgeSet.h>
#include <Ioss_EdgeBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>

#include <string>
#include <utility>
#include <algorithm>

static const std::string SCALAR("scalar");

Ioss::EdgeSet::EdgeSet(const Ioss::DatabaseIO *io_database, const std::string& my_name)
  : Ioss::GroupingEntity(io_database, my_name)
{
  properties.add(Ioss::Property(this,
				"edge_block_count", Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this,
				"block_count", Ioss::Property::INTEGER));
}

Ioss::EdgeSet::~EdgeSet()
{
  try {
    Ioss::EdgeBlockContainer::const_iterator i = edgeBlocks.begin();
    while (i != edgeBlocks.end()) {
      delete *i++;
    }
  } catch (...) {
  }
}

const Ioss::EdgeBlockContainer& Ioss::EdgeSet::get_edge_blocks() const
{ return edgeBlocks; }

Ioss::EntityBlock* Ioss::EdgeSet::get_block(size_t which) const
{
  if (which < edgeBlocks.size())
    return edgeBlocks[which];
  else
    return NULL;
}

Ioss::EdgeBlock* Ioss::EdgeSet::get_edge_block(const std::string& my_name) const
{
  Ioss::EdgeBlock *ge = NULL;
  Ioss::EdgeBlockContainer::const_iterator i = edgeBlocks.begin();
  while (i != edgeBlocks.end()) {
    if ((*i)->name() == my_name) {
      ge = *i;
      break;
    }
    ++i;
  }
  return ge;
}

bool Ioss::EdgeSet::add(Ioss::EdgeBlock *edge_block)
{
    edgeBlocks.push_back(edge_block);
    edge_block->owner_ = this;
    return true;
}

int Ioss::EdgeSet::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::EdgeSet::internal_put_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::EdgeSet::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "edge_block_count")
    return Ioss::Property(my_name, (int)edgeBlocks.size());
  else if (my_name == "block_count")
    return Ioss::Property(my_name, (int)edgeBlocks.size());
  else
    return Ioss::GroupingEntity::get_implicit_property(my_name);
}

void Ioss::EdgeSet::block_membership(std::vector<std::string> &block_members)
{
  if (blockMembership.empty()) {
    Ioss::EdgeBlockContainer::iterator i = edgeBlocks.begin();
    while (i != edgeBlocks.end()) {
      std::vector<std::string> blocks;
      (*i)->block_membership(blocks);
      blockMembership.insert(blockMembership.end(), blocks.begin(), blocks.end());
      ++i;
    }
    std::sort(blockMembership.begin(), blockMembership.end());
    blockMembership.erase(std::unique(blockMembership.begin(), blockMembership.end()), blockMembership.end());
  }
  block_members = blockMembership;
}
