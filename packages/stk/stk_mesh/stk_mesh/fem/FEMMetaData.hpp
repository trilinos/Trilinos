#ifndef stk_mesh_FEMMetaData_hpp
#define stk_mesh_FEMMetaData_hpp

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/util/string_case_compare.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp> // TODO:  Remove!
#include <stk_mesh/fem/CellTopology.hpp>

#include <vector>
#include <string>

namespace stk {
namespace mesh {
namespace fem {


/** \brief FEMMetaData is a class that implements a Finite Element Method skin
 * on top of the Sierra Tool Kit Meta Data class.  The FEM features include the
 * concept of spatial dimension with entity ranks tied to the given spatial
 * dimension, cell topology mapping to parts along with induced cell topology
 * membership through part subsetting, and many additional invariants that are
 * enforced.
 *
 * Users who are interested in a FEM based mesh database should use FEMMetaData
 * rather than MetaData.
 *
 * Invariants for FEM MetaData:
 * 1.  Each cell topology has one and only one root cell topology part.  The
 *     root cell topology part for a cell topology is a unique part that can be
 *     subsetted to induce cell topology on the subset part.
 *     -> Enforced by register_cell_topology is the only function that modifies
 *     the PartCellTopologyVector private data on FEMMetaData.
 * 2.  Root cell topology parts cannot be subsets of parts with cell topologies
 *     -> Enforced by declare_part_subset
 * 3.  Incompatible cell topologies are prohibited.  I.e. parts with
 *     different cell topologies of the same rank (as the part) cannot be subsets
 *     of each other.
 *     -> Enforced by declare_part_subset
 *
 */

// 02/10/11 FEMMetaData Todo:
// * Implement get_cell_topology for Part.
// * Implement declare_part with cell topology
// Non-critical:
// * Implement stk::mesh::fem::get namespace to include getters for MetaData,
//   BulkData, FEMMetaData, FEMBulkData, CellTopology from things like Part,
//   Bucket, Entity, etc.
// * Create impl class inside the handle classes to hold their parent pointer
//   and a friend to the getter above.

class FEMMetaData {
  public:

  /// CellTopologyPartEntityRankMap maps each Cell Topology to its root cell topology part and its associated rank
  typedef std::map<fem::CellTopology, std::pair<Part *, EntityRank> > CellTopologyPartEntityRankMap;
  /// PartCellTopologyVector is a fast-lookup vector of size equal to the number of parts
  typedef std::vector<fem::CellTopology> PartCellTopologyVector;

  
#ifdef SWIG   //SRK from NLM, this is to avoid pulling in a bunch more headers just to define EntityRank
	enum 
	{
    INVALID_RANK = stk::mesh::InvalidEntityRank,
    NODE_RANK = 0u,
    EDGE_RANK = 1u,
    FACE_RANK = 2u,
    VOLUME_RANK = 3u
	};
#else
  static const EntityRank INVALID_RANK = stk::mesh::InvalidEntityRank;
	static const EntityRank NODE_RANK = 0u;
  static const EntityRank EDGE_RANK = 1u;
  static const EntityRank FACE_RANK = 2u;
  static const EntityRank VOLUME_RANK = 3u;
#endif

  FEMMetaData();
  ~FEMMetaData() {}

  /**
   * \brief Construct and initialize a FEMMetaData
   */
  FEMMetaData(size_t spatial_dimension,
              const std::vector<std::string>& in_entity_rank_names = std::vector<std::string>());


  /// --------------------------------------------------------------------------------
  /// FEMMetaData Specific functions begin:
  /// --------------------------------------------------------------------------------
  /** \brief Initialize the spatial dimension and an optional list of entity rank names associated with each rank
   *
   * This function can only be called once.
   * To determine if a FEMMetaData class has been initialized, call the is_FEM_initialized function.
   */
  void FEM_initialize(
      size_t spatial_dimension,
      const std::vector<std::string>& in_entity_rank_names = std::vector<std::string>()
      );

  /** \brief This function returns whether this class has been initialized or not.
   */
  bool is_FEM_initialized() const
  {
    return m_fem_initialized;
  }

  // NOTE: This is a temporary function that will be removed once a FEMBulkData exists.
  /** \brief Getter for MetaData off of a FEMMetaData object.
   */
  inline static MetaData & get_meta_data( FEMMetaData & fem_meta )
    { return fem_meta.m_meta_data; }

  /** \brief Returns the spatial dimension that was passed in through FEM_initialize.
   */
  size_t spatial_dimension() const
  {
    return m_spatial_dimension;
  }

  /** \brief Returns the node rank, which is always zero.
   */
  EntityRank node_rank() const
  {
    return NODE_RANK;
  }

  /** \brief Returns the edge rank which changes depending on spatial dimension
   */
  EntityRank edge_rank() const
  {
    return EDGE_RANK;
  }

  /** \brief Returns the face rank which changes depending on spatial dimension
   */
  EntityRank face_rank() const
  {
    return FACE_RANK;
  }

  /** \brief Returns the volume rank which changes depending on spatial dimension
   */
  EntityRank volume_rank() const
  {
    return VOLUME_RANK;
  }

  /** \brief Returns the side rank which changes depending on spatial dimension
   */
  EntityRank side_rank() const
  {
    return m_side_rank;
  }

  /** \brief Returns the element rank which is always equal to spatial dimension
   */
  EntityRank element_rank() const
  {
    return m_element_rank;
  }
  //  void check_topo_db();


  /** \brief This function is used to register new cell topologies and their associated ranks with FEMMetaData.
   * Currently, several shards Cell Topologies are registered with appropriate ranks at initialization time.
   * See:  internal_declare_known_cell_topology_parts for the whole list.
   *
   * Note:  This function also creates the root cell topology part which is accessible from get_cell_topology_root_part
   */
  void register_cell_topology(const fem::CellTopology cell_topology, EntityRank in_entity_rank);

  /** \brief Return the root cell topology part associated with the given cell topology.
   * This Part is created in register_cell_topology
   */
  Part &get_cell_topology_root_part(const fem::CellTopology cell_topology) const;

  /** \brief Return the cell topology associated with the given part.
   * The cell topology is set on a part through part subsetting with the root
   * cell topology part.
   */
  fem::CellTopology get_cell_topology( const Part & part) const;

  fem::CellTopology get_cell_topology( const std::string & topology_name) const;

  /** \brief Return the EntityRank that is associated with the given cell
   * topology.  In several cases, this rank is dependent on spatial
   * dimension.
   */
  EntityRank get_entity_rank(const fem::CellTopology cell_topology) const;

  /** \brief Getter for FEMMetaData off of a MetaData object.
   */
  inline static FEMMetaData & get ( const MetaData & meta )
    { return *const_cast<FEMMetaData * >(meta.get_attribute<FEMMetaData>()); }

  /** \brief Getter for FEMMetaData off of a Part object.
   */
  inline static FEMMetaData & get( const Part & part )
    { return FEMMetaData::get(MetaData::get(part)); }

  /** \brief Getter for FEMMetaData off of a FieldBase object.
   */
  inline static FEMMetaData & get( const FieldBase & field )
    { return FEMMetaData::get(MetaData::get(field)); }

  /** \brief Getter for FEMMetaData off of a PropertyBase object.
   */
  inline static FEMMetaData & get( const PropertyBase & property )
    { return FEMMetaData::get(MetaData::get(property)); }

  /** \brief Getter for FEMMetaData off of a BulkData object.
   */
  inline static FEMMetaData & get( const BulkData & bulk_data )
    { return FEMMetaData::get(MetaData::get(bulk_data)); }

  /** \brief Getter for FEMMetaData off of a Bucket object.
   */
  inline static FEMMetaData & get( const Bucket & bucket )
    { return FEMMetaData::get(MetaData::get(bucket)); }

  /** \brief Getter for FEMMetaData off of a Entity object.
   */
  inline static FEMMetaData & get( const Entity & entity )
    { return FEMMetaData::get(MetaData::get(entity)); }

  /** \brief Getter for FEMMetaData off of a Ghosting object.
   */
  inline static FEMMetaData & get( const Ghosting & ghost )
    { return FEMMetaData::get(MetaData::get(ghost)); }

  /** \brief  Declare a part with a given cell topology
   */
  Part &declare_part( const std::string &name, fem::CellTopology cell_topology)
  {
    ThrowRequireMsg(is_FEM_initialized(),"FEMMetaData::declare_part: FEM_initialize() must be called before this function");
    Part &root_part = get_cell_topology_root_part(cell_topology);
    EntityRank primary_entity_rank = root_part.primary_entity_rank();
    Part & part = m_meta_data.declare_part(name, primary_entity_rank);
    declare_part_subset(root_part, part);
    return part;
  }

  /** \brief  Declare a part with a given cell topology
   */
  template< class Top >
  Part &declare_part(const std::string &name) {
    return declare_part(name, shards::getCellTopologyData<Top>());
  }

  /// --------------------------------------------------------------------------------
  /// FEMMetaData Specific functions end
  /// --------------------------------------------------------------------------------

  /// --------------------------------------------------------------------------------
  /// The following functions are call-throughs to the underlying MetaData class:
  /// --------------------------------------------------------------------------------

  //------------------------------------
  /** \name MetaData predefined parts
   */

  /** \brief  Universal subset for the problem domain.
   *          All other parts are a subset of the universal part.
   */
  Part & universal_part() const { return m_meta_data.universal_part(); }

  /** \brief  Subset for the problem domain that is owned by the
   *          local process.  Ghost entities are not members of this part.
   */
  Part & locally_owned_part()  const { return m_meta_data.locally_owned_part(); }

  /** \brief  Subset for the problem domain that is shared with another
   *          process.  Ghost entities are not members of this part.
   */
  Part & globally_shared_part() const { return m_meta_data.globally_shared_part(); }

  //------------------------------------
  /** \name  Declare and query parts
   */

  /** \brief  Get an existing part by its application-defined text name.
   *
   *  Return NULL if not present and required_by == NULL.
   *  If required and not present then throws an exception
   *  with the 'required_by' text.
   */
  /// \todo REFACTOR remove required_by argument
  Part * get_part( const std::string & p_name,
                   const char * required_by = NULL ) const
    { return m_meta_data.get_part(p_name,required_by); }

  /** \brief  Get an existing part by its ordinal */
  Part & get_part( unsigned ord ) const
    { return m_meta_data.get_part(ord); }

  /** \brief  Query all parts of the mesh ordered by the parts' ordinal. */
  const PartVector & get_parts() const
    { return m_meta_data.get_parts(); }

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part will have member entities that are of equal or lesser rank.
   *  When an entity of equal rank becomes a member
   *  then all related entities of lesser rank also become members.
   */
  Part & declare_part( const std::string & p_name, EntityRank rank )
    { return m_meta_data.declare_part(p_name,rank); }

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part does not have an entity type rank.
   */
  Part & declare_part( const std::string & p_name)
    { return m_meta_data.declare_part(p_name); }

  /** \brief  Declare a superset-subset relationship between parts
   *  Note:  Cell Topologies are induced through part subsets.
   *  See the invariants that are enforced by this function in the documentation for FEMMetaData.
   * */
  void declare_part_subset( Part & superset , Part & subset );

  /** \brief  Declare an entity-relationship between parts.
   *
   *  If \ref stk::mesh::Entity "entity" <b> e1 </b> is a member
   *  of <em> root_part </em> and there exists an
   *  \ref stk::mesh::Relation "entity relation"
   *  from <b> e1 </b> to <b> e2 </b> that satisfies the
   *  \ref stk_mesh_relations "relation stencil"
   *  then <b> e2 </b> must be a member of the <em> target_part </em>.
   */
  void declare_part_relation( Part & root_part ,
                              relation_stencil_ptr stencil ,
                              Part & target_part )
    { m_meta_data.declare_part_relation(root_part, stencil, target_part); }

  /** \brief Set the entity rank names in a vector.
   * This also currently sets the maximum entity rank.
   */
  void set_entity_rank_names(const std::vector<std::string> &in_entity_rank_names)
  {
    m_entity_rank_names = in_entity_rank_names;
    m_meta_data.set_entity_rank_names(in_entity_rank_names);
  }

  /** \brief Return the rank for the given name that was provided in set_entity_rank_names
   */
  EntityRank entity_rank( const std::string &name ) const
  {
    EntityRank my_entity_rank = InvalidEntityRank;

    for (size_t i = 0; i < m_entity_rank_names.size(); ++i)
      if (equal_case(name, m_entity_rank_names[i])) {
        my_entity_rank = i;
      break;
      }
    return my_entity_rank;
  }

  /** \brief Return the set of entity rank names specified in set_entity_rank_names
   */
  const std::vector<std::string> & entity_rank_names() const
  {
    return m_entity_rank_names;
  }

  /** \brief Return the maximum entity rank
   */
  std::vector<std::string>::size_type entity_rank_count() const
  {
    return m_entity_rank_names.size();
  }

  /** \brief Return the name for a given entity rank as was specified in set_entity_rank_names
   */
  const std::string & entity_rank_name( EntityRank in_entity_rank ) const
  {
    ThrowErrorMsgIf( in_entity_rank >= m_entity_rank_names.size(),
        "entity-rank " << in_entity_rank <<
        " out of range. Must be in range 0.." << m_entity_rank_names.size());
    return m_entity_rank_names[in_entity_rank];
  }

  /** \brief Return true if the given entity rank is valid.
   */
  bool is_valid_entity_rank(EntityRank rank) const
  {
    return rank < m_entity_rank_names.size();
  }

  //------------------------------------
  /** \name  Declare and query fields
   */

  /** \brief  Get a field, return NULL if it does not exist.
   *
   *  \exception std::runtime_error
   *    If the field exits and the
   *    \ref stk::mesh::Field "field_type" does not match or
   *    if required_by != NULL and a field of that name is not found.
   */
  template< class field_type >
  field_type * get_field( const std::string & name ) const
    { return m_meta_data.get_field<field_type>(name); }

  /** \brief  Get all defined fields */
  const FieldVector & get_fields() const {
    return m_meta_data.get_fields();
  }

  /** \brief  Declare a field of the given
   *          \ref stk::mesh::Field "field_type", test name,
   *          and number of states.
   *
   *  A compatible redeclaration returns the previously declared field.
   *  \exception std::runtime_error  If a redeclaration is incompatible
   */
  template< class field_type >
  field_type & declare_field( const std::string & name ,
                              unsigned number_of_states = 1 )
    { return m_meta_data.declare_field<field_type>( name, number_of_states ); }

  /** \brief Declare a field relation.
   *
   *  The pointer_field's scalar type must be a pointer to the
   *  scalar type of the reference_field.  The following
   *  derived field data relationship maintained.
   *
   *  Let   e_root -> Relation( e_target , ord , kind )
   *  Let   i = stencil( e_root.entity_rank() ,
   *                     e_target.entity_rank() , ord , kind )
   *  Let   Scalar ** ptr = field_data( pointer_field , e_root )
   *  then  ptr[i] = field_data( referenced_field , e_target )
   *
   *  This derived field data relationship is typically used
   *  to support fast access to field data on entities
   *  related to the root entity; e.g. field data associated with
   *  the nodes of an element.
   */
  template< class PointerFieldType , class ReferencedFieldType >
  void declare_field_relation( PointerFieldType & pointer_field ,
                               relation_stencil_ptr stencil ,
                               ReferencedFieldType & referenced_field )
    { return m_meta_data.declare_field_relation( pointer_field, stencil, referenced_field ); }

  /** \brief  Get field relations */
  const std::vector<FieldRelation> & get_field_relations() const
    { return m_meta_data.get_field_relations(); }

  /** \brief  Commit the part and field declarations so that the
   *          meta data manager can be used to create
   *          \ref stk::mesh::BulkData "mesh bulk data".
   *
   *  Verifies consistency of the meta data and clean out redundant
   *  field data allocation rules.
   *  Once committed no further part or field declarations can be made.
   */
  void commit()
    { m_meta_data.commit(); }

  /** \brief  Query if the meta data manager is committed */
  bool is_commit() const
    { return m_meta_data.is_commit(); }

  //------------------------------------

  /** \name  Field declaration with weak type information;
   *         direct use in application code is strongly discouraged.
   */

  /** \brief  Declare a field via runtime type information */
  FieldBase * declare_field_base(
    const std::string & arg_name,
    const DataTraits  & arg_traits ,
    unsigned            arg_rank ,
    const shards::ArrayDimTag * const * arg_dim_tags ,
    unsigned arg_num_states )
    { return m_meta_data.declare_field_base( arg_name, arg_traits, arg_rank, arg_dim_tags, arg_num_states); }

  /** \brief  Declare a field restriction via runtime type information.
   */
  void declare_field_restriction( FieldBase      & arg_field ,
                                  EntityRank       arg_entity_rank ,
                                  const Part     & arg_part ,
                                  const unsigned * arg_stride )
    { m_meta_data.declare_field_restriction(arg_field, arg_entity_rank, arg_part, arg_stride); }

  private: // functions

    void internal_set_spatial_dimension_and_ranks(size_t spatial_dimension);

    void internal_declare_known_cell_topology_parts();

  private: // data
    MetaData                      m_meta_data;
    bool                          m_fem_initialized;
    size_t                        m_spatial_dimension;
    EntityRank                    m_side_rank;
    EntityRank                    m_element_rank;
    std::vector< std::string >    m_entity_rank_names;
    /// Used to store mapping between Cell Topologies and their associated root parts and specified ranks:
    CellTopologyPartEntityRankMap m_cellTopologyPartEntityRankMap;
    /// Fast-lookup vector that maps part ordinals to Cell Topologies.
    PartCellTopologyVector        m_partCellTopologyVector;
};

/** \brief Determine if the given part is a root part for a cell topology.
 */
bool is_cell_topology_root_part(const Part & part);

/** set a cell_topology on a part */
void set_cell_topology( Part &part, const fem::CellTopology cell_topology);

/** set a cell_topology on a part */
template<class Topology>
inline void set_cell_topology(Part & part)
{
  stk::mesh::fem::set_cell_topology(part, fem::CellTopology(shards::getCellTopologyData<Topology>()));
}


/** Get the cell_topology off a bucket */
CellTopology get_cell_topology(const Bucket &bucket);


/** Get the cell_topology off an entity */
inline CellTopology get_cell_topology(const Entity &entity) {
  return get_cell_topology(entity.bucket());
}



std::vector<std::string> entity_rank_names(size_t spatial_dimension);

} // namespace fem
} // namespace mesh
} // namespace stk

#endif //  stk_mesh_FEMMetaData_hpp
