// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
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
// 

#ifndef stk_mesh_MetaData_hpp
#define stk_mesh_MetaData_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for NULL, size_t
#include <string.h>                     // for strlen
#include <sys/types.h>                  // for int64_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <map>                          // for map, map<>::value_compare
#include <stk_mesh/base/CellTopology.hpp>  // for CellTopology
#include <stk_mesh/base/EntityKey.hpp>  // for EntityKey
#include <stk_mesh/base/PartField.hpp>    // for PartField
#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/PropertyBase.hpp>  // for Property
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/Types.hpp>      // for EntityRank, PropertyBase, etc
#include <stk_mesh/baseImpl/FieldRepository.hpp>  // for FieldRepository, etc
#include <stk_mesh/baseImpl/PartRepository.hpp>  // for PartRepository
#include <stk_topology/topology.hpp>    // for topology, topology::rank_t, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/util/string_case_compare.hpp>  // for equal_case
#include <string>                       // for string, char_traits
#include <typeinfo>                     // for type_info
#include <utility>                      // for pair
#include <vector>                       // for vector, vector<>::size_type
#include "Shards_CellTopology.hpp"      // for operator<, CellTopology
#include "Shards_CellTopologyTraits.hpp"  // for getCellTopologyData
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits (ptr only), etc
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/FieldState.hpp"  // for ::MaximumFieldStates, etc
#include "stk_mesh/baseImpl/PartImpl.hpp"  // for PartImpl
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsgIf, etc
#include "stk_util/util/CSet.hpp"       // for CSet

namespace shards { class ArrayDimTag; }
namespace shards { class CellTopologyManagedData; }
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Ghosting; } }
namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Print an entity key for this meta data */
std::ostream &
print_entity_key( std::ostream & os, const MetaData & meta_data, const EntityKey & key);

std::string
print_entity_key( const MetaData & meta_data, const EntityKey & key );

bool is_topology_root_part(const Part & part);

/** set a cell_topology on a part */
void set_cell_topology( Part &part, const CellTopology cell_topology);

/** set a cell_topology on a part */
template<class Topology>
inline void set_cell_topology(Part & part)
{
  stk::mesh::set_cell_topology(part, CellTopology(shards::getCellTopologyData<Topology>()));
}

stk::topology get_topology(const MetaData& meta_data, EntityRank entity_rank, const std::pair<const unsigned*, const unsigned*>& supersets);


/** set a stk::topology on a part */
void set_topology(Part &part, stk::topology topology);

/** get the stk::topology given a Shards Cell Topology */
stk::topology get_topology(CellTopology shards_topology, int spatial_dimension = 3);

/** Get the Shards Cell Topology given a stk::topology  */
CellTopology get_cell_topology(stk::topology topo);

//----------------------------------------------------------------------
/** \brief  The manager of an integrated collection of
 *          \ref stk::mesh::Part "parts" and
 *          \ref stk::mesh::Field "fields".
 *
 *  Mesh meta data must be identical on all processors.
 *
 * The FEM features include the concept of spatial dimension with
 * entity ranks tied to the given spatial dimension, topology
 * mapping to parts along with induced topology membership
 * through part subsetting, and many additional invariants that are
 * enforced.
 *
 * Invariants for MetaData:
 * 1.  Each topology has one and only one root topology part.  The
 *     root topology part for a topology is a unique part that can be
 *     subsetted to induce topology on the subset part.
 *     -> Enforced by register_cell_topology is the only function that modifies
 *     the PartCellTopologyVector private data on MetaData.
 * 2.  Root cell topology parts cannot be subsets of parts with cell topologies
 *     -> Enforced by declare_part_subset
 * 3.  Incompatible cell topologies are prohibited.  I.e. parts with
 *     different cell topologies of the same rank (as the part) cannot be subsets
 *     of each other.
 *     -> Enforced by declare_part_subset
 */
class MetaData {
public:

  /** \} */
  //------------------------------------
  /** \name  Meta data manager construction and destruction
   *  \{
   */

  /// CellTopologyPartEntityRankMap maps each Cell Topology to its root cell topology part and its associated rank
  typedef std::map<CellTopology, std::pair<Part *, EntityRank> > CellTopologyPartEntityRankMap;
  /// PartCellTopologyVector is a fast-lookup vector of size equal to the number of parts
  typedef std::vector<CellTopology> PartCellTopologyVector;


  inline static MetaData & get( const Part & part ) { return part.meta_data(); }
  inline static MetaData & get( const FieldBase & field ) { return field.meta_data(); }
  inline static MetaData & get( const PropertyBase & property ) { return property.meta_data(); }

  static MetaData & get( const BulkData & bulk_data );
  static MetaData & get( const Bucket & bucket );
  static MetaData & get( const Ghosting & ghost );

  /** \brief  Construct a meta data manager to own parts and fields.  */
  explicit MetaData(size_t spatial_dimension, const std::vector<std::string>& rank_names = std::vector<std::string>());

  /** \brief  Construct a meta data manager to own parts and fields.  */
  MetaData();

  /** \brief  Destroy the meta data manager and
   *          all of the parts and fields that it owns.
   */
  ~MetaData();

  /** Standard usage associates just one BulkData with a MetaData.
   * An error is thrown if this method is called with a non-NULL BulkData while another
   * non-NULL BulkData has already been set through a previous call.
   * If you wish to replace the BulkData with a different one, you must clear the
   * first one by setting it to NULL, then set the new one.
   */
  void set_mesh_bulk_data(BulkData* bulk)
  {
      ThrowRequireMsg(m_bulk_data == NULL || m_bulk_data == bulk || bulk == NULL, "MetaData::set_mesh_bulk_data ERROR, trying to set mesh when it's already set.");
      m_bulk_data = bulk;
  }

  BulkData& mesh_bulk_data() {
      ThrowRequireMsg(m_bulk_data != NULL, "MetaData::mesh_bulk_data() ERROR, mesh not set yet.");
    return *m_bulk_data;
  }

  const BulkData& mesh_bulk_data() const {
      ThrowRequireMsg(m_bulk_data != NULL, "MetaData::mesh_bulk_data() ERROR, mesh not set yet.");
    return *m_bulk_data;
  }

  bool has_mesh() const {
      return m_bulk_data != NULL;
  }

  //------------------------------------
  /** \name Predefined Parts
   *  \{
   */

  /** \brief  Universal subset for the problem domain.
   *          All other parts are a subset of the universal part.
   */
  Part & universal_part() const { return *m_universal_part; }

  /** \brief  Subset for the problem domain that is owned by the
   *          local process.  Ghost entities are not members of this part.
   */
  Part & locally_owned_part()  const { return *m_owns_part ; }

  /** \brief  Subset for the problem domain that is shared with another
   *          process.  Ghost entities are not members of this part.
   */
  Part & globally_shared_part() const { return *m_shares_part ; }

  /** \brief  Subset for the problem domain that is ghosted from another
   *          process.
   */
  Part & aura_part() const { return *m_aura_part ; }

  /** \} */
  //------------------------------------
  /** \name  Declare and query parts
   *  \{
   */

  /** \brief  Get an existing part by its application-defined text name.
   *
   *  Return NULL if not present.
   *  If not present throws an exception
   */
  /// \todo REFACTOR remove required_by argument
  Part * get_part( const std::string & p_name,
                   const char * required_by = NULL ) const ;

  /** \brief  Get an existing part by its ordinal */
  Part & get_part( unsigned ord ) const ;

  /** \brief  Query all parts of the mesh ordered by the parts' ordinal. */
  const PartVector & get_parts() const { return m_part_repo.get_all_parts(); }

  /** \brief  Query non-internal parts of the mesh ordered by the parts' ordinal. */
  const PartVector get_mesh_parts() const { return m_part_repo.get_mesh_parts(); }

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part will have member entities that are of equal or lesser rank.
   *  When an entity of equal rank becomes a member
   *  then all related entities of lesser rank also become members.
   *
   *  Normally, specifying a valid rank for a part will turn on inducability, but
   *  we want to allow the client to decide whether they want induction or not.
   *  For now, provide a flag to allow turn inducability off since it can be
   *  expensive (bucket fragmentation). Eventually, we should decouple induction
   *  from part rank.
   */
  Part & declare_part( const std::string & p_name, EntityRank rank, bool arg_force_no_induce = false );

  /** \brief  Declare a part of the given name
   *          Redeclaration returns the previously declared part.
   *
   *  This part does not have an entity rank.
   */
  Part & declare_part( const std::string & p_name);

  /** \brief  Declare a part with a given stk topology
   */
  Part &declare_part_with_topology( const std::string &name, stk::topology::topology_t topology, bool arg_force_no_induce = false )
  {
    ThrowRequireMsg(is_initialized(),"MetaData::declare_part: initialize() must be called before this function");
    Part &root_part = get_cell_topology_root_part(stk::mesh::get_cell_topology(topology));
    EntityRank primary_entity_rank = root_part.primary_entity_rank();
    Part & part = declare_part(name, primary_entity_rank, arg_force_no_induce);
    declare_part_subset(root_part, part);
    return part;
  }

  void force_no_induce(Part& part)
  {
    if (part.primary_entity_rank() != InvalidEntityRank) {
      declare_part( part.name(), part.primary_entity_rank(), true /*force no induce*/);
    }
  }

  void set_part_id(Part& part, int64_t lid)
  { part.m_partImpl.set_id(lid); }

  /** \brief  Declare a superset-subset relationship between parts */
  void declare_part_subset( Part & superset , Part & subset );

    /** \brief  Declare an attribute on a part.
   *          Return the attribute of that type,
   *          which may be an already existing value.
   * \todo REFACTOR  Should be using a shared pointer in this interface.
   *       declare_attribute( Part & , shared_ptr<const T> & );
   */
  template<class T>
  const T * declare_attribute_with_delete( Part & part, const T * attribute);
  template<class T>
  const T * declare_attribute_no_delete( Part & part, const T * attribute);
  template<class T>
  bool remove_attribute( Part & part, const T * attribute);

  /** \} */
  //------------------------------------
  /** \name  Entity-ranks
   *  \{
   */

  /** \brief initialize
   *
   */
  void initialize(size_t spatial_dimension, const std::vector<std::string> &rank_names = std::vector<std::string>());

  bool is_initialized() const
  { return !m_entity_rank_names.empty(); }

  EntityRank entity_rank( const std::string &name ) const;

  const std::vector<std::string> & entity_rank_names() const
    { return m_entity_rank_names ; }

  std::vector<std::string>::size_type entity_rank_count() const
    { return m_entity_rank_names.size(); }

  const std::string & entity_rank_name( EntityRank entity_rank ) const ;

  /** \brief Returns the side rank which changes depending on spatial dimension
   */
  EntityRank side_rank() const
  {
    switch (m_spatial_dimension)
    {
    case 1 : return stk::topology::NODE_RANK;
    case 2 : return stk::topology::EDGE_RANK;
    case 3 : return stk::topology::FACE_RANK;
    default: break;
    }
    return stk::topology::INVALID_RANK;
  }

  /**
   * Return true if rank is valid.
   */
  bool check_rank(EntityRank rank) const;

  /** \} */
  //------------------------------------
  /** \name  Declare and query fields
   *  \{
   */

  /** \brief  Get a field by name, return NULL if it does not exist.
   *
   * Note that this is a case-insensitive name search.
   * E.g., 'TEMPERATURE' is the same as 'temperature'.
   *
   *  \exception std::runtime_error
   *    If the field exits and the
   *    \ref stk::mesh::Field "field_type" does not match or
   *    if a field of that name is not found.
   */
  template< class field_type >
  field_type * get_field( stk::mesh::EntityRank entity_rank, const std::string & name ) const ;

  /**
   * \brief Get a field by name with unknown type, NULL if does not exist
   *
   * Note that this is a case-insensitive name search.
   * E.g., 'TEMPERATURE' is the same as 'temperature'.
   */
  FieldBase* get_field( stk::mesh::EntityRank entity_rank, const std::string& name ) const;

  /** \brief  Get/Set the coordinate field */
  FieldBase const* coordinate_field() const;
  void set_coordinate_field(FieldBase* coord_field) { m_coord_field = coord_field; }

  /** \brief  Get all defined fields */
  const FieldVector & get_fields() const {
    return m_field_repo.get_fields() ;
  }

  const FieldVector & get_fields(stk::topology::rank_t rank) const {
    return m_field_repo.get_fields(rank) ;
  }

  /** \brief  Declare a field of the given
   *          \ref stk::mesh::Field "field_type", test name,
   *          and number of states.
   *
   *  A compatible redeclaration returns the previously declared field.
   *  \exception std::runtime_error  If a redeclaration is incompatible
   *
   *  See Field.hpp for a full discussion of Fields.
   */

  template< class field_type >
  field_type & declare_field( stk::topology::rank_t arg_entity_rank,
                              const std::string & name ,
                              unsigned number_of_states = 1 );

  /** \brief  Declare an attribute on a field.
   *          Return the attribute of that type,
   *          which may be an already existing value.
   * \todo REFACTOR  Should be using a shared pointer in this interface.
   */
  template<class T>
  const T * declare_attribute_with_delete( FieldBase & field, const T * attribute);
  template<class T>
  const T * declare_attribute_no_delete( FieldBase & field, const T * attribute);
  template<class T>
  bool remove_attribute( FieldBase & field, const T * attribute);

  template< class data_type >
  PartField<data_type>& declare_part_field(unsigned itemsPerPart = 1)
  {
    unsigned part_field_index = m_part_fields.size();
    PartField<data_type>* new_part_field = new PartField<data_type>(this, part_field_index, itemsPerPart);
    m_part_fields.push_back(new_part_field);
    return *new_part_field;
  }

  /** \} */
  //------------------------------------

  template<class T>
  const T * get_attribute() const ;

  /** \brief  Declare an attribute on the meta data.
   *          Return the attribute of that type,
   *          which may be an already existing value.
   */
  template<class T>
  const T * declare_attribute_with_delete( const T * attribute);

  template<class T>
  const T * declare_attribute_no_delete( const T * attribute);

  template<class T>
  bool remove_attribute( const T * );

  //------------------------------------
  /** \} */
  /** \name  Declare and query properties associated with parts
   *  \{
   */

  /** \brief  Get a property, return NULL if it does not exist.
   *
   *  \exception std::runtime_error
   *    If the property exits and the
   *    \ref stk::mesh::Property "type" does not match or
   */
  template< typename DataType >
  Property<DataType> * get_property( const std::string & name ) const ;

  /** \brief  Get all defined properties */
  const std::vector< PropertyBase * > & get_properties() const
    { return m_properties ; }

  /** \brief  Declare a property of the given
   *          \ref stk::mesh::Property "type", name, and dimensions.
   *
   *  A compatible redeclaration returns the previously declared property.
   *  \exception std::runtime_error  If a redeclaration is incompatible
   */
  template< typename DataType >
  Property<DataType> & declare_property( const std::string & name ,
                                         unsigned size = 1 );

  /** \brief  Put a property on the given part */
  void put_property( PropertyBase & property, Part & part);

  /** get the spatial-dimension. */
  unsigned spatial_dimension() const { return m_spatial_dimension; }

  /** \brief  Commit the part and field declarations so that the
   *          meta data manager can be used to create
   *          \ref stk::mesh::BulkData "mesh bulk data".
   *
   *  Verifies consistency of the meta data and clean out redundant
   *  field data allocation rules.
   *  Once committed no further part or field declarations can be made.
   */
  void commit();

  /** \brief  Query if the meta data manager is committed */
  bool is_commit() const { return m_commit ; }

  /** \} */
  //------------------------------------

  /** \name  Field declaration with weak type information;
   *         direct use in application code is strongly discouraged.
   *  \{
   */

  /** \brief  Declare a field via runtime type information */
  FieldBase * declare_field_base(
    const std::string & arg_name,
    stk::topology::rank_t arg_entity_rank,
    const DataTraits  & arg_traits ,
    unsigned            arg_rank ,
    const shards::ArrayDimTag * const * arg_dim_tags ,
    unsigned arg_num_states )
  {
    require_not_committed();

    return m_field_repo.declare_field(
                  arg_name, arg_entity_rank, arg_traits, arg_rank, arg_dim_tags,
                  arg_num_states, this
                 );
  }



  /** \brief  Declare a field restriction via runtime type information.
   */
  void declare_field_restriction( FieldBase      & arg_field ,
                                  const Part     & arg_part ,
                                  const unsigned   arg_num_scalars_per_entity ,
                                  const unsigned   arg_first_dimension ,
                                  const void*      arg_init_value = NULL );

  /** \brief  Declare a field restriction via runtime type information.
   */
  void declare_field_restriction( FieldBase      & arg_field ,
                                  const Selector & arg_selector ,
                                  const unsigned   arg_num_scalars_per_entity ,
                                  const unsigned   arg_first_dimension ,
                                  const void*      arg_init_value = NULL );

  /** \brief This function is used to register new cell topologies and their associated ranks with MetaData.
   * Currently, several shards Cell Topologies are registered with appropriate ranks at initialization time.
   * See:  internal_declare_known_cell_topology_parts for the whole list.
   *
   * Note:  This function also creates the root cell topology part which is accessible from get_cell_topology_root_part
   */
  void register_cell_topology(const CellTopology cell_topology, EntityRank in_entity_rank);

  shards::CellTopology register_superelement_cell_topology(stk::topology t);

  /** \brief Return the root cell topology part associated with the given cell topology.
   * This Part is created in register_cell_topology
   */

  Part &get_cell_topology_root_part(const CellTopology cell_topology) const;

  /** \brief Return the topology part given a stk::topology.
   */
  Part &get_topology_root_part(stk::topology topology) const
  { return get_cell_topology_root_part(stk::mesh::get_cell_topology(topology)); }

  /** \brief Return the cell topology associated with the given part.
   * The cell topology is set on a part through part subsetting with the root
   * cell topology part.
   */
  CellTopology get_cell_topology( const Part & part) const;
  stk::topology get_topology(const Part & part) const;

  CellTopology get_cell_topology( const std::string & topology_name) const;

  void dump_all_meta_info(std::ostream& out = std::cout) const;

  void set_mesh_on_fields(BulkData* bulk);

  /** \} */
private:
  // Functions

  MetaData( const MetaData & );                ///< \brief  Not allowed
  MetaData & operator = ( const MetaData & );  ///< \brief  Not allowed

  void add_new_part_in_part_fields();
  void synchronize_part_fields_with_parts();

  Part & declare_internal_part( const std::string & p_name);

  Part & declare_internal_part( const std::string & p_name, EntityRank rank);

  void internal_declare_known_cell_topology_parts();

  void internal_declare_part_subset( Part & superset , Part & subset );

  void assign_cell_topology( Part & part, CellTopology topo);

  // Members

  BulkData* m_bulk_data;
  bool   m_commit ;
  impl::PartRepository m_part_repo ;
  CSet   m_attributes ;

  Part * m_universal_part ;
  Part * m_owns_part ;
  Part * m_shares_part ;
  Part * m_aura_part ;

  impl::FieldRepository        m_field_repo ;
  mutable FieldBase* m_coord_field;

  std::vector< PropertyBase* > m_properties ;
  std::vector< std::string >   m_entity_rank_names ;
  std::vector<shards::CellTopologyManagedData*> m_created_topologies;

  unsigned m_spatial_dimension;
  EntityRank m_side_rank;

  std::vector<PartFieldBase*> m_part_fields;

  /// Used to store mapping between Cell Topologies and their associated root parts and specified ranks:
  CellTopologyPartEntityRankMap m_cellTopologyPartEntityRankMap;
  /// Fast-lookup vector that maps part ordinals to Cell Topologies.
  PartCellTopologyVector        m_partCellTopologyVector;

  /** \name  Invariants/preconditions for MetaData.
   * \{
   */
  void require_committed() const ;

  void require_not_committed() const ;

  void require_same_mesh_meta_data( const MetaData & rhs ) const ;

  void require_valid_entity_rank( EntityRank rank) const ;

  void require_not_relation_target( const Part * const part ) const ;
  /** \} */
  //------------------------------------

  Property<void> * get_property_base( const std::string & ,
                                      const std::type_info & ,
                                      unsigned = 0 ) const ;

  void clean_field_restrictions();
};

template< typename part_field_type >
typename part_field_type::PartFieldDataType* part_field_data(part_field_type& part_field, const Part& part)
{
  unsigned part_ordinal = part.mesh_meta_data_ordinal();
  return part_field.data(part_ordinal);
}

/** \brief  Verify that the meta data is identical on all processors */
void verify_parallel_consistency( const MetaData & , ParallelMachine );

/** Get default entity rank names */
const std::vector<std::string>& entity_rank_names();

/** \name  Declare field data allocation rules
 *  \{
 */

/** \brief  Declare a field to exist for a given entity type and Part.
 *
 * See Field.hpp for a full discussion of field restrictions.
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Selector & selector ,
                        const void* init_value = NULL);

/** \brief Declare a field to exist for a given entity type and Part. The
 *         extra unsigned arguments specify the size of a dimension. So,
 *         put_field( field, rank, part, 3, 3 ) would create a 3x3 2D field.
 *         Fields of up to seven dimensions are supported.
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Selector & selector ,
                        unsigned     n1 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Selector & selector ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Selector & selector ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        unsigned     n6 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        unsigned     n6 ,
                        unsigned     n7 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field_on_entire_mesh_with_initial_value(field_type & field, const typename FieldTraits<field_type>::data_type *initial_value)
{
    return put_field(field, field.mesh_meta_data().universal_part(), initial_value);
}

template< class field_type >
field_type & put_field_on_entire_mesh(field_type & field)
{
    typename FieldTraits<field_type>::data_type* init_value = NULL;
    return put_field_on_entire_mesh_with_initial_value(field, init_value);
}

/** \} */
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

namespace stk {
namespace mesh {

inline
Part & MetaData::get_part( unsigned ord ) const
{ return * m_part_repo.get_all_parts()[ord] ; }

template< class field_type >
inline
field_type * MetaData::get_field( stk::mesh::EntityRank arg_entity_rank, const std::string & name ) const
{
  typedef FieldTraits< field_type > Traits ;

  const DataTraits & dt = data_traits< typename Traits::data_type >();

  const shards::ArrayDimTag * tags[8] ;

  Traits::assign_tags( tags );

  FieldBase * const field =
    m_field_repo.get_field( static_cast<stk::topology::rank_t>(arg_entity_rank), name , dt , Traits::Rank , tags , 0 );
  if (field == NULL) {
    return static_cast<field_type*>(NULL);
  }
  else {
    return dynamic_cast< field_type * >( field );
  }
}


template< class field_type >
inline
field_type & MetaData::declare_field( stk::topology::rank_t arg_entity_rank,
                                      const std::string & name ,
                                      unsigned number_of_states )
{

  typedef FieldTraits< field_type > Traits ;

  const DataTraits & traits = data_traits< typename Traits::data_type >();

  const shards::ArrayDimTag * dim_tags[8] ;

  Traits::assign_tags( dim_tags );

  static const char* reserved_state_suffix[6] = {
    "_STKFS_OLD",
    "_STKFS_N",
    "_STKFS_NM1",
    "_STKFS_NM2",
    "_STKFS_NM3",
    "_STKFS_NM4"
  };

  // Check that the name does not have a reserved suffix

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    const int len_name   = name.size();
    const int len_suffix = std::strlen( reserved_state_suffix[i] );
    const int offset     = len_name - len_suffix ;
    if ( 0 <= offset ) {
      const char * const name_suffix = name.c_str() + offset ;
      ThrowErrorMsgIf( equal_case( name_suffix , reserved_state_suffix[i] ),
          "For name = \"" << name_suffix <<
          "\" CANNOT HAVE THE RESERVED STATE SUFFIX \"" <<
          reserved_state_suffix[i] << "\"" );
    }
  }

  // Check that the field of this name has not already been declared

  field_type * f[ MaximumFieldStates ] ;

  f[0] = dynamic_cast<field_type*>(m_field_repo.get_field(
      arg_entity_rank , name ,
      traits , Traits::Rank , dim_tags , number_of_states
      ));

  if ( NULL != f[0] ) {
    for ( unsigned i = 1 ; i < number_of_states ; ++i ) {
      f[i] = &f[0]->field_of_state(static_cast<FieldState>(i));
    }
  }
  else {
    // Field does not exist then create it

    std::string field_names[ MaximumFieldStates ];

    field_names[0] = name ;

    if ( 2 == number_of_states ) {
      field_names[1] = name ;
      field_names[1].append( reserved_state_suffix[0] );
    }
    else {
      for ( unsigned i = 1 ; i < number_of_states ; ++i ) {
        field_names[i] = name ;
        field_names[i].append( reserved_state_suffix[i] );
      }
    }

    for ( unsigned i = 0 ; i < number_of_states ; ++i ) {

      f[i] = new field_type(
          this,
          arg_entity_rank,
          m_field_repo.get_fields().size() ,
          field_names[i] ,
          traits ,
          Traits::Rank,
          dim_tags,
          number_of_states ,
          static_cast<FieldState>(i)
          );

      m_field_repo.add_field( f[i] );
    }

    for ( unsigned i = 0 ; i < number_of_states ; ++i ) {
      f[i]->m_impl.set_field_states( f );
    }
  }

  f[0]->set_mesh(m_bulk_data);

  return *f[0] ;
}


template< class field_type >
inline
field_type & put_field(
  field_type & field ,
  const Part & part ,
  const void* init_value)
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;
  unsigned stride[8] = {0,0,0,0,0,0,0,0};
  Helper::assign(stride);

  unsigned numScalarsPerEntity = 1;
  if(field.field_array_rank() > 0)
  {
      numScalarsPerEntity = stride[0];
  }
  unsigned firstDimension = numScalarsPerEntity;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field(
  field_type & field ,
  const Selector & selector ,
  const void* init_value)
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;
  unsigned stride[8] = {0,0,0,0,0,0,0,0};
  Helper::assign(stride);

  unsigned numScalarsPerEntity = 1;
  if(field.field_array_rank() > 0)
  {
      numScalarsPerEntity = stride[0];
  }
  unsigned firstDimension = numScalarsPerEntity;
  MetaData::get(field).declare_field_restriction( field, selector, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Selector &selector ,
                        unsigned    n1 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, selector, numScalarsPerEntity, firstDimension, init_value);

return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Selector &selector ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, selector, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2*n3;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Selector &selector ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2*n3;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, selector, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2*n3*n4;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        unsigned    n5 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2*n3*n4*n5;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        unsigned    n5 ,
                        unsigned    n6 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2*n3*n4*n5*n6;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        unsigned    n5 ,
                        unsigned    n6 ,
                        unsigned    n7 ,
                        const void* init_value )
{
  unsigned numScalarsPerEntity = n1*n2*n3*n4*n5*n6*n7;
  unsigned firstDimension = n1;
  MetaData::get(field).declare_field_restriction( field, part, numScalarsPerEntity, firstDimension, init_value);

  return field ;
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( const T * a )
{
  return m_attributes.insert_with_delete( a );
}

template<class T>
inline
const T *
MetaData::get_attribute() const
{ return m_attributes.get<T>(); }

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( const T * attribute )
{
  return m_attributes.insert_no_delete( attribute );
}

template<class T>
inline
bool
MetaData::remove_attribute( const T * a )
{
  return m_attributes.remove( a );
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( Part & part , const T * attribute )
{
  return m_part_repo.declare_attribute_with_delete( part, attribute );
}

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( Part & part , const T * attribute )
{
  return m_part_repo.declare_attribute_no_delete( part, attribute );
}

template<class T>
inline
bool
MetaData::remove_attribute( Part & part , const T * attribute )
{
  return m_part_repo.remove_attribute(part, attribute);
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( FieldBase & field , const T * attribute )
{
  return m_field_repo.declare_attribute_with_delete(field, attribute);
}

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( FieldBase & field , const T * attribute )
{
  return m_field_repo.declare_attribute_no_delete(field, attribute);
}

template<class T>
inline
bool
MetaData::remove_attribute( FieldBase & field , const T * attribute )
{
  return m_field_repo.remove_attribute(field, attribute);
}

//----------------------------------------------------------------------


template< typename DataType >
inline
Property<DataType> *
MetaData::get_property( const std::string & name ) const
{
  Property<void> * const pv = get_property_base( name, typeid(DataType) );
  return pv ? pv->property<DataType>() : static_cast<Property<DataType>*>(NULL) ;
}

template< typename DataType >
inline
Property<DataType> &
MetaData::declare_property( const std::string & name , unsigned size )
{
  Property<void> * pv = get_property_base(name,typeid(DataType),size);
  Property<DataType> * prop = NULL ;

  if ( pv != NULL ) {
    prop = pv->property<DataType>();
  }
  else {
    if ( 1 == size ) {
      pv = prop = new Property<DataType>( *this , m_properties.size() , name );
    }
    else {
      pv = prop = new Property< std::vector<DataType> >(
                    *this , m_properties.size() , name , size );
    }
    m_properties.push_back( pv );
  }
  return *prop ;
}

inline
void MetaData::put_property( PropertyBase & property , Part & part )
{
  property.add_property( part.mesh_meta_data_ordinal() );
}

inline
bool MetaData::check_rank(EntityRank rank) const
{
  return rank < m_entity_rank_names.size();
}

inline
bool
is_auto_declared_part(const Part &part)
{
  return stk::mesh::impl::is_internal_part(part);
}

template< class field_type >
field_type * get_field_by_name( const std::string & name, const MetaData & metaData )
{
  field_type* field = NULL;
  unsigned num_nonnull_fields = 0;
  for(stk::topology::rank_t i=stk::topology::NODE_RANK; i<=stk::topology::CONSTRAINT_RANK; ++i) {
    field_type* thisfield = metaData.get_field<field_type>(i, name);
    if (thisfield != NULL) {
      if (field == NULL) {
        field = thisfield;
      }
      ++num_nonnull_fields;
    }
  }

  if (num_nonnull_fields > 1) {
    std::cerr << "get_field_by_name WARNING, found "<<num_nonnull_fields<<" fields with name="<<name
      <<". Returning the first one."<<std::endl;
  }

  return field;
}

FieldBase* get_field_by_name( const std::string& name, const MetaData & metaData );

} // namespace mesh
} // namespace stk

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_MetaData_hpp */
