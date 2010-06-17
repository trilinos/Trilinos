/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_FaceSet_h
#define IOSS_Ioss_FaceSet_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_GroupingEntity.h>
#include <string>

#include <vector>

namespace Ioss {
  class DatabaseIO;
  class FaceBlock;
  class EntityBlock;

  typedef std::vector<FaceBlock*> FaceBlockContainer;

  class FaceSet : public GroupingEntity {
  public:
    FaceSet(const DatabaseIO *io_database, const std::string& name);
    ~FaceSet();

    std::string type_string() const {return "FaceSet";}
    EntityType type() const {return FACESET;}

    bool add(FaceBlock    *face_block);
    const FaceBlockContainer&    get_face_blocks() const;
    FaceBlock*    get_face_block(const std::string& name) const;
    size_t face_block_count() const {return faceBlocks.size();}
    
    size_t block_count() const {return faceBlocks.size();}
    EntityBlock* get_block(size_t which) const;

    void block_membership(std::vector<std::string> &block_members);

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string& name) const;

  protected:
    int internal_get_field_data(const Field& field,
				void *data, size_t data_size) const;

    int internal_put_field_data(const Field& field,
				void *data, size_t data_size) const;

  private:
    FaceBlockContainer      faceBlocks;
    std::vector<std::string> blockMembership; // What element blocks do the
                                             // elements in this faceset belong to.
  };
}
#endif
