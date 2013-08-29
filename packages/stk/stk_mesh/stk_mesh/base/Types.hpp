/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Types_hpp
#define stk_mesh_Types_hpp

//----------------------------------------------------------------------

#include <stdint.h>
#include <limits>
#include <utility>
#include <vector>

#include <stk_topology/topology.hpp>

#include <stk_util/util/PairIter.hpp>
#include <stk_util/util/NamedPair.hpp>
#include <stk_util/util/TrackingAllocator.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_topology/topology.hpp>

// #ifndef STK_MESH_DISALLOW_DEPRECATED_ENTITY_FNS
// #define STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
// #endif

namespace stk {
namespace mesh {

// Tags used by tracking allocator
struct FieldDataTag {};
struct PartitionTag {};
struct BucketTag {};
struct EntityCommTag {};
struct BucketRelationTag {};
struct DynamicBucketRelationTag {};
struct DynamicBucketNodeRelationTag {};
struct DynamicBucketEdgeRelationTag {};
struct DynamicBucketFaceRelationTag {};
struct DynamicBucketElementRelationTag {};
struct DynamicBucketOtherRelationTag {};
struct AuxRelationTag {};
struct DeletedEntityTag {};

void print_dynamic_connectivity_profile( ParallelMachine parallel, int parallel_rank, std::ostream & out);

void print_max_stk_memory_usage( ParallelMachine parallel, int parallel_rank, std::ostream & out);

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

class MetaData ;  // Meta-data description of a mesh
class Part ;      // Defined subset of the mesh

/** \brief  Collections of \ref stk::mesh::Part "parts" are frequently
 *          maintained as a vector of Part pointers.
 */
typedef std::vector< Part * > PartVector;
typedef std::vector< unsigned > OrdinalVector;
typedef std::vector< unsigned > PermutationIndexVector;

class FieldBase;

template< typename Scalar = void ,
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void >
  class Field ;

/** \brief Maximum
 *  \ref shards::Array "multi-dimenaional array" dimension of a
 *  \ref stk::mesh::Field "field"
 */
enum { MaximumFieldDimension = 7 };

template< typename DataType = void > class Property ;

typedef Property< void > PropertyBase ;

/** \} */

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

class BulkData ; // Bulk-data of a mesh
class Bucket ;   // Homogeneous collection of mesh entitities and their field data
union Entity ;   // Individual entity within the mesh
struct EntityKey ;
class Relation ; // Relation pair of local mesh entities
class Ghosting ;

typedef std::vector<Bucket *> BucketVector;
typedef std::vector<Entity> EntityVector;

/** Change log to reflect change from before 'modification_begin'
  *  to the current status.
  */
enum EntityState { Unchanged = 0 ,
                   Created  = 1 ,
                   Modified = 2 ,
                   Deleted  = 3 };

template< class FieldType > struct EntityArray ;
template< class FieldType > struct BucketArray ;
template< class FieldType > struct FieldTraits ;

//MeshIndex describes an Entity's location in the mesh, specifying which bucket,
//and the offset (ordinal) into that bucket.
//Ultimately we want this struct to contain two ints rather than a pointer and an int...
struct MeshIndex {
  Bucket* bucket;
  size_t bucket_ordinal;
};

typedef unsigned Ordinal;
static const Ordinal InvalidOrdinal = static_cast<Ordinal>(-1); // std::numeric_limits<PartOrdinal>::max();
typedef Ordinal EntityRank ;
typedef Ordinal PartOrdinal;
typedef Ordinal FieldOrdinal;
typedef Ordinal RelationIdentifier;
typedef Ordinal FieldArrayRank;

typedef uint64_t EntityId ;

// Base Entity Rank
// Note:  This BaseEntityRank can be considered the leaf of a tree and it
// represents the furthest out you can go in downward relations.
static const EntityRank BaseEntityRank = 0;
static const EntityRank InvalidEntityRank = stk::topology::INVALID_RANK;
static const PartOrdinal InvalidPartOrdinal = InvalidOrdinal;
static const RelationIdentifier InvalidRelationIdentifier = InvalidOrdinal;
static const int InvalidProcessRank = -1;

  inline unsigned GetInvalidLocalId() {
    static unsigned InvalidLocalId = std::numeric_limits<unsigned int>::max();
    return InvalidLocalId;
  }

/**
* Predefined identifiers for mesh object relationship types.
*/
struct RelationType
{
  enum relation_type_t
  {
    USES      = 0 ,
    USED_BY   = 1 ,
    CHILD     = 2 ,
    PARENT    = 3 ,
    EMBEDDED  = 0x00ff , // 4
    CONTACT   = 0x00ff , // 5
    AUXILIARY = 0x00ff ,
    INVALID   = 10
  };

  RelationType(relation_type_t value = INVALID) : m_value(value) {}

  operator relation_type_t() const { return m_value; }

  relation_type_t m_value;
};

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_bulk_data_parallel
 *  \{
 */

/** \brief  Pairing of an entity with a processor rank */
typedef std::pair<Entity , int> EntityProc ;
typedef std::vector<EntityProc>     EntityProcVec ;

typedef std::pair<EntityKey, int> EntityKeyProc;

/** \brief  Spans of a vector of entity-processor pairs are common.
 *
 */
typedef PairIter< std::vector< EntityProc >::const_iterator >
  PairIterEntityProc ;
#ifndef SWIG
	//NLM SWIG cannot handle this macro

NAMED_PAIR( EntityCommInfo , unsigned , ghost_id , int , proc )

/** \brief  Span of ( communication-subset-ordinal , process-rank ) pairs
 *          for the communication of an entity.
 */
typedef std::vector<EntityCommInfo, tracking_allocator<EntityCommInfo,EntityCommTag > > EntityCommInfoVector;
typedef PairIter<  EntityCommInfoVector::const_iterator >  PairIterEntityComm ;

#endif
/** \} */

//----------------------------------------------------------------------
/** \ingroup stk_mesh_relations
 *  \brief  A relation stencil maps entity relationships to ordinals.
 *
 *  A relation stencil function is the inverse mapping of a contiguous
 *  span of non-negative integers to a template of entity relations.
 *  For example, a triangle-to-vertex relation stencil would map:
 *  -  0 = relation_stencil( Element , Node , 0 )
 *  -  1 = relation_stencil( Element , Node , 1 )
 *  -  2 = relation_stencil( Element , Node , 2 )
 *
 *  If the input entity relationship is within the stencil then
 *  a stencil function returns a non-negative integer;
 *  otherwise a stencil function returns a negative value.
 */
typedef int ( * relation_stencil_ptr )( unsigned  from_type ,
                                        unsigned  to_type ,
                                        unsigned  identifier );

//----------------------------------------------------------------------
/** \brief  Span of a sorted relations for a given domain entity.
 *
 *  The span is sorted by
 *  -# range entity rank,
 *  -# relation identifier, and
 *  -# range entity global identifier.
 */
typedef std::vector<Relation, tracking_allocator<Relation,AuxRelationTag> > RelationVector;

typedef PairIter< RelationVector::const_iterator > PairIterRelation ;

enum ConnectivityType
{
  FIXED_CONNECTIVITY,
  DYNAMIC_CONNECTIVITY,
  INVALID_CONNECTIVITY_TYPE
};

#define EXTRACT_BUCKET_ID(idx) ((idx) >> NUM_BUCKET_ORDINAL_BITS)

#define EXTRACT_BUCKET_ORDINAL(idx) ((idx) & BUCKET_ORDINAL_MASK)

enum ConnectivityOrdinal
{
  INVALID_CONNECTIVITY_ORDINAL = ~0U
};

inline
ConnectivityOrdinal& operator++(ConnectivityOrdinal& ord)
{
  ord = static_cast<ConnectivityOrdinal>(ord + 1);
  return ord;
}

enum Permutation
{
  INVALID_PERMUTATION = ~0U
};

enum ConnectivityId
{
  INVALID_CONNECTIVITY_ID = ~0U
};

//////////////////////////////////////////////////////////////////////////////

template <EntityRank TargetRank>
struct DynamicConnectivityTagSelector
{
  typedef DynamicBucketOtherRelationTag type;
};

template <>
struct DynamicConnectivityTagSelector<stk::topology::NODE_RANK>
{
  typedef DynamicBucketNodeRelationTag type;
};

template <>
struct DynamicConnectivityTagSelector<stk::topology::EDGE_RANK>
{
  typedef DynamicBucketEdgeRelationTag type;
};

template <>
struct DynamicConnectivityTagSelector<stk::topology::FACE_RANK>
{
  typedef DynamicBucketFaceRelationTag type;
};

template <>
struct DynamicConnectivityTagSelector<stk::topology::ELEMENT_RANK>
{
  typedef DynamicBucketElementRelationTag type;
};

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

#endif
