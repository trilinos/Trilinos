/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
