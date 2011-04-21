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

#include <Ioss_SideSet.h>

#include <Ioss_SideBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Region.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_ElementTopology.h>

#include <string>
#include <algorithm>

static const std::string SCALAR("scalar");

Ioss::SideSet::SideSet(const Ioss::DatabaseIO *io_database, const std::string& my_name)
  : Ioss::GroupingEntity(io_database, my_name)
{
  properties.add(Ioss::Property(this,
			       "side_block_count", Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this,
			       "block_count", Ioss::Property::INTEGER));
}

Ioss::SideSet::~SideSet()
{
  try {
    Ioss::SideBlockContainer::const_iterator i = sideBlocks.begin();
    while (i != sideBlocks.end()) {
      delete *i++;
    }
  } catch (...) {
  }
}

const Ioss::SideBlockContainer& Ioss::SideSet::get_side_blocks() const
{ return sideBlocks; }

Ioss::SideBlock* Ioss::SideSet::get_block(size_t which) const
{
  if (which < sideBlocks.size())
    return sideBlocks[which];
  else
    return NULL;
}

Ioss::SideBlock* Ioss::SideSet::get_side_block(const std::string& my_name) const
{
  Ioss::SideBlock *ge = NULL;
  Ioss::SideBlockContainer::const_iterator i = sideBlocks.begin();
  while (i != sideBlocks.end()) {
    if ((*i)->name() == my_name) {
      ge = *i;
      break;
    }
    ++i;
  }
  return ge;
}

bool Ioss::SideSet::add(Ioss::SideBlock *side_block)
{
    sideBlocks.push_back(side_block);
    side_block->owner_ = this;
    return true;
}

int Ioss::SideSet::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::SideSet::internal_put_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::SideSet::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "side_block_count")
    return Ioss::Property(my_name, (int)sideBlocks.size());
  else if (my_name == "block_count")
    return Ioss::Property(my_name, (int)sideBlocks.size());
  else
    return Ioss::GroupingEntity::get_implicit_property(my_name);
}

int Ioss::SideSet::max_parametric_dimension() const
{
  int max_par_dim = 0;

  Ioss::SideBlockContainer::const_iterator i = sideBlocks.begin();
  while (i != sideBlocks.end()) {
    Ioss::SideBlock * const sideblock = *i ;
    int parametric_dim = sideblock->topology()->parametric_dimension();
    if (parametric_dim > max_par_dim)
      max_par_dim = parametric_dim;
    ++i;
  }
  if (max_par_dim == 0) {
    // If the sideset is empty, return the maximum that the parametric dimension could be...
    // Faces for 3D model; Edges for 2D model
    const Ioss::Region *reg = get_database()->get_region();
    max_par_dim = reg->get_property("spatial_dimension").get_int()-1;
  }
  return max_par_dim;
}

void Ioss::SideSet::block_membership(std::vector<std::string> &block_members)
{
  if (blockMembership.empty()) {
    Ioss::SideBlockContainer::iterator i = sideBlocks.begin();
    while (i != sideBlocks.end()) {
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
