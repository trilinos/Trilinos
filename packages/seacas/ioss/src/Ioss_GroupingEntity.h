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

#ifndef IOSS_Ioss_GroupingEntity_h
#define IOSS_Ioss_GroupingEntity_h

#include <Ioss_CodeTypes.h>
#include <string>
#include <Ioss_State.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_FieldManager.h>
#include <Ioss_EntityType.h>

namespace Ioss {

  // Base class for all 'grouping' entities. The following derived classes
  // are typical:
  // -- NodeList  -- grouping of nodes (0d topology)
  // -- EdgeList  -- grouping of edges (1d topology)
  // -- FaceList  -- grouping of faces (2d topology) [Surface]
  //
  // Similarly, there is:
  // -- NodeBlock -- grouping of 'similar' nodes (same degree of freedom, ...)
  // -- ElementBlock -- grouping of 'similar' elements (same element topology,
  //                    attributes, ...)
  //    0d, 1d, 2d, 3d topology possible -- e.g., sphere, bar, quad, hex
  //
  // A Region is also a grouping entity, except that its list of subentites
  // are other GroupingEntities. That is, it maintains a list of NodeBlocks,
  // ElementBlocks, NodeLists, CommLists and Surfaces. [Similar to the
  // "Composite Patter" in Design Patterns]  All interface to GroupingEntities
  // is through the Region class; clients of the IO subsystem have no direct
  // access to the underlying GroupingEntities (other than the Region).
  //
  // Each GroupingEntity contains:
  // -- name
  // -- MeshEntities of the specified topological dimension
  // -- Optional attributes, either global (applied to the groupingentity), or
  //    unique value(s) to be applied to each subentity.
  // -- Data items
  //
  //========================================================================
  class EntityBlock;
  
  class GroupingEntity {
  public:
    friend class Property;

    explicit GroupingEntity(const DatabaseIO *io_database = NULL,
			    const std::string& name="");
    virtual ~GroupingEntity();

    State get_state() const;

    const DatabaseIO* get_database() const;
    void  set_database(const DatabaseIO *io_database);
    virtual void delete_database();

    //: Return name of entity.  This short-circuits the process of
    //  getting the name via the property.  Returns the same information
    //  as: entity->get_property("name").get_string()
    const std::string & name() const {return entityName;}
    void set_name(const std::string &new_name) {entityName = new_name;}
    
    //: Returns true if 'name' is an alias for this entity.
    bool is_alias(const std::string &name) const;

    //: Return list of blocks that the entities in this GroupingEntity "touch"
    //: For a FaceSet or EdgeSet, returns a list of the element blocks that the
    //: elements in the set belong to.  For a nodeset, returns "nodeblock_1".
    //: For others, it returns an empty vector.
    //: Entries are pushed onto the "block_members" vector, so it will be
    //: appended to if it is not empty at entry to the function.
    virtual void block_membership(std::vector<std::string> &block_members) {}
      
    std::string get_filename() const;
    virtual std::string type_string() const = 0;
    virtual EntityType type() const = 0;
      
    // ========================================================================
    //                                PROPERTIES
    // ========================================================================
    // Property-related information....
    // Just forward it through to the property manager...
    void property_add(const Property& new_prop);
    void property_erase(const std::string& property_name);
    bool property_exists(const std::string& property_name) const;
    Property get_property(const std::string& property_name) const;
    int property_describe(NameList* names) const;
    size_t property_count() const;

    // ========================================================================
    //                                FIELDS
    // ========================================================================
    // Just forward these through to the field manager...
    void field_add(const Field& new_field);
    void field_erase(const std::string& field_name);
    bool field_exists(const std::string& field_name) const;
    Field get_field(const std::string& field_name) const;
    const Field &get_fieldref(const std::string& field_name) const;
    int field_describe(NameList* names) const;
    int field_describe(Field::RoleType role, NameList* names) const;
    size_t field_count() const;

    // Put this fields data into 'data'.
    // Returns number of entities for which the field was read.
    // Assumes 'data' is large enough to hold all values.
    int get_field_data(const std::string& field_name, void *data,
		       size_t data_size) const;

    int put_field_data(const std::string& field_name, void *data,
		       size_t data_size) const;

    // Put this fields data into the specified std::vector space.
    // Returns number of entities for which the field was read.
    // Resizes 'data' to size needed to hold all values.
    int get_field_data(const std::string& field_name, std::vector<char>    &data) const;
    int get_field_data(const std::string& field_name, std::vector<double>  &data) const;
    int get_field_data(const std::string& field_name, std::vector<int>     &data) const;
    int get_field_data(const std::string& field_name, std::vector<Complex> &data) const;

    int put_field_data(const std::string& field_name, std::vector<char>    &data) const;
    int put_field_data(const std::string& field_name, std::vector<double>  &data) const;
    int put_field_data(const std::string& field_name, std::vector<int>     &data) const;
    int put_field_data(const std::string& field_name, std::vector<Complex> &data) const;


  protected:
    bool set_state(State new_state) {
      entityState = new_state;
      return true;
    }

    // Protected to give access to Region which is the only
    // class that should delete the database. May have to make
    // private and provide friend...
    void really_delete_database();

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    // Note that even though this is a pure virtual function, an implementation
    // is provided to return properties that are common to all grouping entities.
    // Derived classes should call 'GroupingEntity::get_implicit_property'
    // if the requested property is not specific to their type.
    virtual Property
      get_implicit_property(const std::string& name) const = 0;

    PropertyManager properties;
    FieldManager    fields;

    virtual int internal_get_field_data(const Field& field,
					void *data, size_t data_size=0) const = 0;
    virtual int internal_put_field_data(const Field& field,
					void *data, size_t data_size=0) const = 0;

  private:
    GroupingEntity(const GroupingEntity&); // do not implement
    GroupingEntity& operator=(const GroupingEntity&); // do not implement

    std::string entityName;
    const DatabaseIO* database_;

    State entityState;
  };
}
inline void
Ioss::GroupingEntity::property_add(const Ioss::Property& new_prop)
{properties.add(new_prop);}

inline void
Ioss::GroupingEntity::property_erase(const std::string& property_name)
{properties.erase(property_name);}

inline bool
Ioss::GroupingEntity::property_exists(const std::string& property_name) const
{return properties.exists(property_name);}

inline Ioss::Property
Ioss::GroupingEntity::get_property(const std::string& property_name) const
{return properties.get(property_name);}

inline int
Ioss::GroupingEntity::property_describe(NameList* names) const
{return properties.describe(names);}

inline size_t
Ioss::GroupingEntity::property_count() const
{return properties.count();}

// ------------------------------------------------------------------------

inline void
Ioss::GroupingEntity::field_erase(const std::string& field_name)
{fields.erase(field_name);}

inline bool
Ioss::GroupingEntity::field_exists(const std::string& field_name) const
{return fields.exists(field_name);}

inline Ioss::Field
Ioss::GroupingEntity::get_field(const std::string& field_name) const
{return fields.get(field_name);}

inline const Ioss::Field&
Ioss::GroupingEntity::get_fieldref(const std::string& field_name) const
{return fields.getref(field_name);}

inline int
Ioss::GroupingEntity::field_describe(NameList* names) const
{return fields.describe(names);}

inline int
Ioss::GroupingEntity::field_describe(Ioss::Field::RoleType role, NameList* names) const
{return fields.describe(role, names);}

inline size_t
Ioss::GroupingEntity::field_count() const
{return fields.count();}

#endif
