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

#include <stk_util/stk_config.h>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/NamedPair.hpp>  // for NAMED_PAIR
#include <stk_util/util/PairIter.hpp>   // for PairIter
#include <cstddef>
#include <cstdint>
#include <limits>                       // for numeric_limits
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include <set>
#include <map>


namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Relation; } }
namespace stk { namespace mesh { struct Entity; } }
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
typedef std::vector<EntityKey>      EntityKeyVector;

enum class Layout : uint8_t
{
  Left,    // Adjacent Entities in memory
  Right,   // Adjacent components in memory
  Auto     // Run-time access to Field data with the correct layout; Not for Field registration
};

inline std::ostream& operator<<(std::ostream& os, Layout layout) {
    switch (layout) {
        case Layout::Left:
            os << "Layout::Left";
            break;
        case Layout::Right:
            os << "Layout::Right";
            break;
        case Layout::Auto:
            os << "Layout::Auto";
            break;
        default:
            os << "Unknown Layout";
            break;
    }
    return os;
}

#ifdef STK_UNIFIED_MEMORY
constexpr Layout DefaultHostLayout = Layout::Left;
#else
constexpr Layout DefaultHostLayout = Layout::Right;
#endif

constexpr Layout DefaultDeviceLayout = Layout::Left;

template <typename Scalar = void, Layout HostLayout = DefaultHostLayout>
class Field;

enum class Operation
{
  SUM,
  MIN,
  MAX
};

enum FieldAccessTag : uint8_t
{
  ReadOnly     = 0, // Sync values to memory space and do not mark as modified; Disallow modification
  ReadWrite    = 1, // Sync values to memory space and mark as modified; Allow modification
  OverwriteAll = 2, // Do not sync values to memory space and mark as modified; Allow modification

  Unsynchronized,      // Do not sync values to memory space and do not mark as modified; Allow modification
  ConstUnsynchronized, // Do not sync values to memory space and do not mark as modified; Disallow modification

  InvalidAccess     // For internal use only.  Not valid for accessing data.
};

constexpr int NumTrackedFieldAccessTags = 3;

inline std::ostream& operator<<(std::ostream& os, FieldAccessTag accessTag) {
    switch (accessTag) {
        case ReadOnly:
            os << "ReadOnly";
            break;
        case ReadWrite:
            os << "ReadWrite";
            break;
        case OverwriteAll:
            os << "OverwriteAll";
            break;
        case Unsynchronized:
            os << "Unsynchronized";
            break;
        case ConstUnsynchronized:
            os << "ConstUnsynchronized";
            break;
        case InvalidAccess:
            os << "InvalidAccess";
            break;
        default:
            os << "Unknown Access Tag";
            break;
    }
    return os;
}

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
//
//MeshIndex describes an Entity's location in the mesh, specifying which bucket,
//and the offset (ordinal) into that bucket.
struct MeshIndex
{
  Bucket* bucket;
  unsigned bucket_ordinal;

  constexpr MeshIndex(Bucket* bucketIn, size_t ordinal) : bucket(bucketIn), bucket_ordinal(ordinal) {}
};

// Smaller than MeshIndex and replaces bucket pointer with bucket_id to
// remove hop.
struct FastMeshIndex
{
  unsigned bucket_id;
  unsigned bucket_ord;
};

KOKKOS_INLINE_FUNCTION
constexpr bool operator<(const FastMeshIndex& lhs, const FastMeshIndex& rhs)
{
  return lhs.bucket_id == rhs.bucket_id ? lhs.bucket_ord < rhs.bucket_ord : lhs.bucket_id < rhs.bucket_id;
}

KOKKOS_INLINE_FUNCTION
constexpr bool operator==(const FastMeshIndex& lhs, const FastMeshIndex& rhs)
{
  return lhs.bucket_id == rhs.bucket_id && lhs.bucket_ord == rhs.bucket_ord;
}

KOKKOS_INLINE_FUNCTION
constexpr bool operator!=(const FastMeshIndex& lhs, const FastMeshIndex& rhs)
{
  return lhs.bucket_id != rhs.bucket_id || lhs.bucket_ord != rhs.bucket_ord;
}

typedef stk::topology::rank_t EntityRank ;

typedef std::map<std::pair<EntityRank, Selector>, std::pair<size_t, size_t> > SelectorCountMap;
typedef std::map<Selector, BucketVector> SelectorBucketMap;

typedef std::map<EntityKey,std::set<int> > EntityToDependentProcessorsMap;

typedef unsigned Ordinal;
static const Ordinal InvalidOrdinal = static_cast<Ordinal>(-1); // std::numeric_limits<PartOrdinal>::max();

typedef Ordinal PartOrdinal;
typedef Ordinal RelationIdentifier;
typedef Ordinal FieldArrayRank;

typedef uint64_t EntityId ;
static constexpr EntityId InvalidEntityId = std::numeric_limits<stk::mesh::EntityId>::max();

typedef std::vector<EntityId> EntityIdVector;

static constexpr EntityRank InvalidEntityRank = stk::topology::INVALID_RANK;
static constexpr PartOrdinal InvalidPartOrdinal = InvalidOrdinal;
static constexpr RelationIdentifier InvalidRelationIdentifier = InvalidOrdinal;
static constexpr int InvalidProcessRank = -1;

constexpr unsigned GetInvalidLocalId()
{
  unsigned InvalidLocalId = std::numeric_limits<unsigned int>::max();
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
    INVALID   = 10
  };

  constexpr RelationType(relation_type_t value = INVALID) : m_value(value) {}

  constexpr operator relation_type_t() const { return m_value; }

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
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Aug 2025
using EntityKeyProcMap STK_DEPRECATED = std::map<EntityKey, int>;
#endif

/** \brief  Spans of a vector of entity-processor pairs are common.
 *
 */
typedef PairIter< std::vector< EntityProc >::const_iterator >
  PairIterEntityProc ;
#ifndef SWIG
	//NLM SWIG cannot handle this macro

NAMED_PAIR( EntityCommInfo , unsigned , ghost_id , int , proc )

inline bool operator>=(const EntityCommInfo& lhs, const EntityCommInfo& rhs)
{
  return lhs.ghost_id >= rhs.ghost_id && lhs.proc >= rhs.proc;
}

/** \brief  Span of ( communication-subset-ordinal , process-rank ) pairs
 *          for the communication of an entity.
 */
typedef std::vector<EntityCommInfo> EntityCommInfoVector;
typedef PairIter<const EntityCommInfo*>  PairIterEntityComm ;

#endif
/** \} */

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
constexpr unsigned INVALID_PARTITION_ID = std::numeric_limits<unsigned>::max();

#define STK_16BIT_CONNECTIVITY_ORDINAL
#ifdef STK_16BIT_CONNECTIVITY_ORDINAL
using ConnectivityOrdinal = uint16_t;
constexpr ConnectivityOrdinal INVALID_CONNECTIVITY_ORDINAL = 65535;
#else
using ConnectivityOrdinal = uint32_t;
constexpr ConnectivityOrdinal INVALID_CONNECTIVITY_ORDINAL = std::numeric_limits<ConnectivityOrdinal>::max();
#endif

#ifdef STK_16BIT_UPWARDCONN_INDEX_TYPE
using UpwardConnIndexType = uint16_t;
constexpr UpwardConnIndexType INVALID_UPWARDCONN_INDEX = 65535;
#else
using UpwardConnIndexType = uint32_t;
constexpr UpwardConnIndexType INVALID_UPWARDCONN_INDEX = std::numeric_limits<UpwardConnIndexType>::max();
#endif

enum Permutation : unsigned char
{
  DEFAULT_PERMUTATION = 0,
  INVALID_PERMUTATION = 128
};

struct FieldMetaData
{
  std::byte* m_data {};
  int m_bytesPerEntity {};
  int m_numComponentsPerEntity {};
  int m_numCopiesPerEntity {};
  int m_bucketSize {};
  int m_bucketCapacity {};
};

struct DeviceFieldMetaData
{
  std::byte* m_data {};
  std::byte* m_hostData {};
  int m_numComponentsPerEntity {};
  int m_numCopiesPerEntity {};
  int m_bucketSize {};
  int m_bucketCapacity {};
};

} // namespace mesh
} // namespace stk

#endif
