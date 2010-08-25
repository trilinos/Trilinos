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

#include <Ioss_FaceSet.h>

#include <Ioss_FaceBlock.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_Property.h>
#include <Ioss_Field.h>

#include <string>
#include <algorithm>

static const std::string SCALAR("scalar");

Ioss::FaceSet::FaceSet(const Ioss::DatabaseIO *io_database, const std::string& my_name)
  : Ioss::GroupingEntity(io_database, my_name)
{
  properties.add(Ioss::Property(this,
			       "face_block_count", Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this,
			       "block_count", Ioss::Property::INTEGER));
}

Ioss::FaceSet::~FaceSet()
{
  try {
    Ioss::FaceBlockContainer::const_iterator i = faceBlocks.begin();
    while (i != faceBlocks.end()) {
      delete *i++;
    }
  } catch (...) {
  }
}

const Ioss::FaceBlockContainer& Ioss::FaceSet::get_face_blocks() const
{ return faceBlocks; }

Ioss::EntityBlock* Ioss::FaceSet::get_block(size_t which) const
{
  if (which < faceBlocks.size())
    return faceBlocks[which];
  else
    return NULL;
}

Ioss::FaceBlock* Ioss::FaceSet::get_face_block(const std::string& my_name) const
{
  Ioss::FaceBlock *ge = NULL;
  Ioss::FaceBlockContainer::const_iterator i = faceBlocks.begin();
  while (i != faceBlocks.end()) {
    if ((*i)->name() == my_name) {
      ge = *i;
      break;
    }
    ++i;
  }
  return ge;
}

bool Ioss::FaceSet::add(Ioss::FaceBlock *face_block)
{
    faceBlocks.push_back(face_block);
    face_block->owner_ = this;
    return true;
}

int Ioss::FaceSet::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->get_field(this, field, data, data_size);
}

int Ioss::FaceSet::internal_put_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  return get_database()->put_field(this, field, data, data_size);
}

Ioss::Property
Ioss::FaceSet::get_implicit_property(const std::string& my_name) const
{
  if (my_name == "face_block_count")
    return Ioss::Property(my_name, (int)faceBlocks.size());
  else if (my_name == "block_count")
    return Ioss::Property(my_name, (int)faceBlocks.size());
  else
    return Ioss::GroupingEntity::get_implicit_property(my_name);
}

void Ioss::FaceSet::block_membership(std::vector<std::string> &block_members)
{
  if (blockMembership.empty()) {
    Ioss::FaceBlockContainer::iterator i = faceBlocks.begin();
    while (i != faceBlocks.end()) {
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
