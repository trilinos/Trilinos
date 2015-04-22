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

#ifndef stk_mesh_Types_hpp
#define stk_mesh_Types_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for uint64_t
#include <ostream>
#include <iosfwd>                       // for ostream
#include <limits>                       // for numeric_limits
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/util/NamedPair.hpp>  // for NAMED_PAIR
#include <stk_util/util/PairIter.hpp>   // for PairIter
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include <set>
#include "boost/range/iterator_range_core.hpp"  // for iterator_range
#include "boost/unordered/unordered_map.hpp"  // for unordered_map

namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Relation; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { namespace impl { class EntityRepository; } } }
namespace stk { namespace mesh { struct EntityKey; } }
namespace stk { namespace mesh { template <typename DataType = void> class Property; } }


namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

class MetaData ;  // Meta-data description of a mesh
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


typedef Property< void > PropertyBase ;

/** \} */

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

class BulkData ; // Bulk-data of a mesh
class Ghosting ;

/** Change log to reflect change from before 'modification_begin'
  *  to the current status.
  */
enum EntityState { Unchanged = 0 ,
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
  size_t bucket_ordinal;
};

// Smaller than MeshIndex and replaces bucket pointer with bucket_id to
// remove hop.
struct FastMeshIndex
{
  unsigned bucket_id;
  unsigned bucket_ord;
};

typedef std::vector<std::vector<FastMeshIndex> > VolatileFastSharedCommMapOneRank;
typedef stk::topology::rank_t EntityRank ;

typedef boost::unordered_map<EntityKey, size_t> GhostReuseMap;
typedef std::map<std::pair<EntityRank, Selector>, std::pair<size_t, size_t> > SelectorCountMap;
typedef std::map<std::pair<EntityRank, Selector>, BucketVector> SelectorBucketMap;
typedef std::vector<VolatileFastSharedCommMapOneRank> VolatileFastSharedCommMap;

typedef std::map<EntityKey,std::set<int> > EntityToDependentProcessorsMap;
typedef std::map<EntityKey,int> NewOwnerMap;

typedef unsigned Ordinal;
static const Ordinal InvalidOrdinal = static_cast<Ordinal>(-1); // std::numeric_limits<PartOrdinal>::max();

//typedef Ordinal EntityRank ;
typedef Ordinal PartOrdinal;
typedef Ordinal FieldOrdinal;
typedef Ordinal RelationIdentifier;
typedef Ordinal FieldArrayRank;

typedef uint64_t EntityId ;

// Base Entity Rank
// Note:  This BaseEntityRank can be considered the leaf of a tree and it
// represents the furthest out you can go in downward relations.
static const EntityRank BaseEntityRank = stk::topology::BEGIN_RANK;
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
typedef boost::iterator_range<RelationIterator> RelationRange;
#endif // SIERRA_MIGRATION

enum ConnectivityType
{
  FIXED_CONNECTIVITY,
  DYNAMIC_CONNECTIVITY,
  INVALID_CONNECTIVITY_TYPE
};

inline std::ostream & operator<<(std::ostream &out, ConnectivityType type)
{
  switch (type)
  {
      case ConnectivityType::FIXED_CONNECTIVITY:
          out << "Fixed  ";
          break;
      case ConnectivityType::DYNAMIC_CONNECTIVITY:
          out << "Dynamic";
          break;
      case ConnectivityType::INVALID_CONNECTIVITY_TYPE:
          out << "Invalid";
          break;
      default:
          out << "Unknown";
  }
  return out;
}

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
  DEFAULT_PERMUTATION = 0,
  INVALID_PERMUTATION = 128
};

enum ConnectivityId
{
  INVALID_CONNECTIVITY_ID = ~0U
};

//----------------------------------------------------------------------

// Use macro below to deprecate functions
#ifdef __GNUC__
#define STK_DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define STK_DEPRECATED(func) func
#endif


} // namespace mesh
} // namespace stk

#endif
