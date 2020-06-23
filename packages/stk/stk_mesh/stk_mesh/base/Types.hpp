// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef stk_mesh_Types_hpp
#define stk_mesh_Types_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for uint64_t
#include <limits>                       // for numeric_limits
#include <stk_util/stk_config.h>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/util/NamedPair.hpp>  // for NAMED_PAIR
#include <stk_util/util/PairIter.hpp>   // for PairIter
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include <set>
#include <map>


namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Relation; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { namespace impl { class EntityRepository; } } }
namespace stk { namespace mesh { struct EntityKey; } }


namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

class BulkData;
class MetaData;
class FieldBase;

/** \brief  Collections of \ref stk::mesh::Part "parts" are frequently
 *          maintained as a vector of Part pointers.
 */
typedef std::vector< Part * >       PartVector;
typedef std::vector< Bucket * >     BucketVector;
typedef std::vector< const Part * > ConstPartVector;
typedef std::vector< FieldBase * >  FieldVector;
typedef std::vector< unsigned >     OrdinalVector;
typedef std::vector< unsigned >     PermutationIndexVector;
typedef std::vector<Entity>         EntityVector;


template< typename Scalar = void ,
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void >
  class Field ;

/** \brief Maximum
 *  \ref "multi-dimensional array" dimension of a
 *  \ref stk::mesh::Field "field"
 */
enum { MaximumFieldDimension = 7 };

enum class Operation
{
  SUM,
  MIN,
  MAX
};


/** \} */

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

/** Change log to reflect change from before 'modification_begin'
  *  to the current status.
  */
enum EntityState : char { Unchanged = 0 ,
                   Created  = 1 ,
                   Modified = 2 ,
                   Deleted  = 3 };

template< class FieldType > struct FieldTraits ;

//MeshIndex describes an Entity's location in the mesh, specifying which bucket,
//and the offset (ordinal) into that bucket.
//Ultimately we want this struct to contain two ints rather than a pointer and an int...
struct MeshIndex
{
  Bucket* bucket;
  unsigned bucket_ordinal;

  MeshIndex(Bucket *bucketIn, size_t ordinal) : bucket(bucketIn), bucket_ordinal(ordinal) {}
};

// Smaller than MeshIndex and replaces bucket pointer with bucket_id to
// remove hop.
struct FastMeshIndex
{
  unsigned bucket_id;
  unsigned bucket_ord;
};

NAMED_PAIR(BucketInfo, unsigned, bucket_id, unsigned, num_entities_this_bucket)

struct BucketIndices
{
  std::vector<BucketInfo> bucket_info;
  std::vector<unsigned> ords;
};

typedef std::vector<BucketIndices> VolatileFastSharedCommMapOneRank;
typedef stk::topology::rank_t EntityRank ;

typedef std::map<std::pair<EntityRank, Selector>, std::pair<size_t, size_t> > SelectorCountMap;
typedef std::map<std::pair<EntityRank, Selector>, BucketVector> SelectorBucketMap;
typedef std::vector<VolatileFastSharedCommMapOneRank> VolatileFastSharedCommMap;

typedef std::map<EntityKey,std::set<int> > EntityToDependentProcessorsMap;
typedef std::vector<std::pair<EntityKey,Entity> >::const_iterator const_entity_iterator;
typedef std::vector<std::pair<EntityKey,Entity> >::iterator entity_iterator;

typedef unsigned Ordinal;
static const Ordinal InvalidOrdinal = static_cast<Ordinal>(-1); // std::numeric_limits<PartOrdinal>::max();

typedef Ordinal PartOrdinal;
typedef Ordinal RelationIdentifier;
typedef Ordinal FieldArrayRank;

typedef uint64_t EntityId ;
static const EntityId InvalidEntityId = std::numeric_limits<stk::mesh::EntityId>::max();

typedef std::vector<EntityId> EntityIdVector;

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
using EntityProc    = std::pair<Entity, int>;
using EntityProcVec = std::vector<EntityProc>;
using EntityProcMap = std::map<Entity, int>;

using EntityIdProc    = std::pair<EntityId, int>;
using EntityIdProcVec = std::vector<EntityIdProc>;
using EntityIdProcMap = std::map<EntityId, int>;

using EntityKeyProc    = std::pair<EntityKey, int>;
using EntityKeyProcVec = std::vector<EntityKeyProc>;
using EntityKeyProcMap = std::map<EntityKey, int>;

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
typedef std::vector<EntityCommInfo> EntityCommInfoVector;
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
typedef int ( * relation_stencil_ptr )( EntityRank  from_type ,
                                        EntityRank  to_type ,
                                        unsigned  identifier );

//----------------------------------------------------------------------
/** \brief  Span of a sorted relations for a given domain entity.
 *
 *  The span is sorted by
 *  -# range entity rank,
 *  -# relation identifier, and
 *  -# range entity global identifier.
 */
typedef std::vector<Relation> RelationVector;

typedef PairIter< RelationVector::const_iterator > PairIterRelation ;

#ifdef SIERRA_MIGRATION
typedef RelationVector::const_iterator   RelationIterator;
#endif // SIERRA_MIGRATION

enum ConnectivityType
{
  FIXED_CONNECTIVITY,
  DYNAMIC_CONNECTIVITY,
  INVALID_CONNECTIVITY_TYPE
};

constexpr unsigned INVALID_BUCKET_ID = std::numeric_limits<unsigned>::max();

#define STK_16BIT_CONNECTIVITY_ORDINAL
#ifdef STK_16BIT_CONNECTIVITY_ORDINAL
using ConnectivityOrdinal = uint16_t;
constexpr ConnectivityOrdinal INVALID_CONNECTIVITY_ORDINAL = 65535;
#else
using ConnectivityOrdinal = uint32_t;
constexpr ConnectivityOrdinal INVALID_CONNECTIVITY_ORDINAL = ~0U;
#endif

enum Permutation
{
  DEFAULT_PERMUTATION = 0,
  INVALID_PERMUTATION = 128
};

} // namespace mesh
} // namespace stk

#endif
