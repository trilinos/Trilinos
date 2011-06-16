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

#include <stk_util/util/PairIter.hpp>
#include <stk_util/util/NamedPair.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

class MetaData ;  // Meta-data description of a mesh
class Part ;      // Defined subset of the mesh

/** \brief  Collections of \ref stk::mesh::Part "parts" are frequently
 *          maintained as a vector of Part pointers.
 */
typedef std::vector< Part * > PartVector ;

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
class Bucket ;   // Homogeneous collection of mesh entitities their field data
class Entity ;   // Individual entity within the mesh
class Relation ; // Relation pair of local mesh entities
class Ghosting ;

typedef std::vector<Bucket *> BucketVector;
typedef std::vector<Entity *> EntityVector;

/** Change log to reflect change from before 'modification_begin'
  *  to the current status.
  */
enum EntityModificationLog { EntityLogNoChange = 0 ,
                             EntityLogCreated  = 1 ,
                             EntityLogModified = 2 ,
                             EntityLogDeleted  = 3 };

template< class FieldType > struct EntityArray ;
template< class FieldType > struct BucketArray ;
template< class FieldType > struct FieldTraits ;


typedef unsigned Ordinal;
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
static const EntityRank InvalidEntityRank = static_cast<EntityRank>(-1); // std::numeric_limits<EntityRank>::max();
static const PartOrdinal InvalidPartOrdinal = static_cast<PartOrdinal>(-1); // std::numeric_limits<PartOrdinal>::max();

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_bulk_data_parallel
 *  \{
 */

/** \brief  Pairing of an entity with a processor rank */
typedef std::pair<Entity*,unsigned> EntityProc ;
typedef std::vector<EntityProc>     EntityProcVec ;

/** \brief  Spans of a vector of entity-processor pairs are common.
 *
 */
typedef PairIter< std::vector< EntityProc >::const_iterator >
  PairIterEntityProc ;

NAMED_PAIR( EntityCommInfo , unsigned , ghost_id , unsigned , proc )

/** \brief  Span of ( communication-subset-ordinal , process-rank ) pairs
 *          for the communication of an entity.
 */
typedef PairIter< std::vector< EntityCommInfo >::const_iterator >
  PairIterEntityComm ;

typedef std::vector<EntityCommInfo> EntityCommInfoVector;

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
typedef PairIter< std::vector<Relation>::const_iterator > PairIterRelation ;
typedef std::vector<Relation> RelationVector;

//----------------------------------------------------------------------


} // namespace mesh
} // namespace stk

#endif
