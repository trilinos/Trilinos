#ifndef stk_mesh_Types_hpp
#define stk_mesh_Types_hpp

//----------------------------------------------------------------------

#include <stdint.h>
#include <limits>
#include <utility>
#include <vector>

#include <stk_util/util/PairIter.hpp>
#include <stk_mesh/base/EntityKey.hpp>

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

template< typename Scalar = void ,
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void >
class Field ;

typedef Field< void, void, void, void, void, void, void, void > FieldBase ;

template< class Field, class Relation > 
class GatherField ;


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

template< class FieldType > struct EntityArray ;
template< class FieldType > struct BucketArray ;
template< class FieldType > struct FieldTraits ;

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_bulk_data_parallel
 *  \{
 */

/** \brief  Pairing of an entity with a processor rank */
typedef std::pair<Entity*,unsigned> EntityProc ;

/** \brief  Spans of a vector of entity-processor pairs are common.
 * 
 */
typedef PairIter< std::vector< EntityProc >::const_iterator >
  PairIterEntityProc ;

/** \} */

//----------------------------------------------------------------------
/** \ingroup stk_mesh_relations
 *  \brief  A relation stencil maps entity relationships to ordinals.
 *
 *  A relation stencil function is the inverse mapping of a contiguous
 *  span of non-negative integers to a template of entity relations.
 *  For example, a triangle-to-vertex relation stencil would map:
 *  -  0 = relation_stencil( Element , Node , 0 , 0 )
 *  -  1 = relation_stencil( Element , Node , 1 , 0 )
 *  -  2 = relation_stencil( Element , Node , 2 , 0 )
 *  
 *  If the input entity relationship is within the stencil then
 *  a stencil function returns a non-negative integer;
 *  otherwise a stencil function returns a negative value.
 */
typedef int ( * relation_stencil_ptr )( unsigned  from_type ,
                                        unsigned  to_type ,
                                        unsigned  identifier ,
                                        unsigned  kind );

//----------------------------------------------------------------------


} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

class UnitTestMetaData ;
class UnitTestBulkData ;

} // namespace mesh
} // namespace stk

#endif

