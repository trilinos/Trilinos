/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// -*- Mode: c++ -*-
#ifndef SIERRA_Ioss_NodeBlock_h
#define SIERRA_Ioss_NodeBlock_h

#include <Ioss_CodeTypes.h>
#include <Ioss_EntityBlock.h>
#include <string>

namespace Ioss {
  class DatabaseIO;

  class NodeBlock : public EntityBlock {
  public:
    NodeBlock(const DatabaseIO *io_database,
	      const std::string& name,
	      size_t node_count,
	      size_t degrees_of_freedom);

    ~NodeBlock();

    std::string type_string() const {return "NodeBlock";}
    EntityType type() const {return NODEBLOCK;}

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string& name) const;

  protected:
    int internal_get_field_data(const Field& field,
				void *data, size_t data_size) const;

    int internal_put_field_data(const Field& field,
				void *data, size_t data_size) const;

  };
}
#endif
