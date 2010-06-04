/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_ElementBlock_h
#define SIERRA_Ioss_ElementBlock_h

#include <Ioss_CodeTypes.h>
#include <Ioss_Property.h>
#include <Ioss_EntityBlock.h>
#include <string>
#include <assert.h>

namespace Ioss {
  class DatabaseIO;

  class ElementBlock : public EntityBlock {
  public:
    ElementBlock(const DatabaseIO *io_database,
		 const std::string& name, const std::string& element_type,
		 size_t number_elements, size_t number_attributes);

    ~ElementBlock();

    std::string type_string() const {return "ElementBlock";}
    EntityType type() const {return ELEMENTBLOCK;}

    /// Handle implicit properties -- These are calcuated from data stored
    /// in the grouping entity instead of having an explicit value assigned.
    /// An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string& name) const;

    bool contains(size_t local_id) const
      { return idOffset < local_id && local_id <= idOffset + elementCount; }

    void  set_offset(int offset) {idOffset = offset;}
    int get_offset() const {assert(idOffset >= 0); return idOffset;}

    void get_block_adjacencies(std::vector<std::string> &block_adjacency_list) const;

  protected:
    int internal_get_field_data(const Field& field,
				void *data, size_t data_size) const;

    int internal_put_field_data(const Field& field,
				void *data, size_t data_size) const;
  private:
    void count_attributes() const;

    /**
     * The 'offset' is used to map an element location within an
     * element block to the element 'file descriptor'.
     * For example, the file descriptor of the 37th element in the 4th
     * block is calculated by:
     *
     * file_descriptor = offset of block 4 + 37
     *
     * This can also be used to determine which element block
     * an element with a file_descriptor maps into. An particular
     * element block contains all elements in the range:
     *
     * offset < file_descriptor <= offset+number_elements_per_block
     */
    size_t idOffset;
    
    size_t elementCount; ///< stored locally to avoid looking up property
    mutable size_t attributeCount;
  };
}
#endif
