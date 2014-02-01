/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_base_Combo_hpp
#define stk_mesh_base_Combo_hpp

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/BucketConnectivity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/ConnectivityMap.hpp>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/util/PairIter.hpp>

#include <boost/static_assert.hpp>
#include <boost/range.hpp>
#include <boost/type_traits/is_pod.hpp>

#include <stk_topology/topology.hpp>

#include <utility>
#include <vector>
#include <iosfwd>
#include <string>
#include <algorithm>

#include <boost/functional/hash.hpp>


#include <stk_mesh/base/Bucket.tcc> //only place where this file should be included

///////////////////////////////////////////////////////////////////////////////
// Put methods below that could not otherwise be inlined due to cyclic dependencies between
// Relation/Entity/Bucket
///////////////////////////////////////////////////////////////////////////////


#ifdef SIERRA_MIGRATION

namespace stk { namespace mesh { class BulkData;} }

namespace sierra {
namespace Fmwk {

class MeshObjRoster;
class MeshObjSharedAttr;
class MeshBulkData;

namespace detail {
bool set_attributes( MeshBulkData& meshbulk, stk::mesh::Entity , const int , const MeshObjSharedAttr*, const int);
bool set_attributes( MeshBulkData& meshbulk, stk::mesh::Entity , const MeshObjSharedAttr*, const int);
bool update_relation( stk::mesh::Entity, const stk::mesh::RelationIterator ir, const bool back_rel_flag, MeshBulkData& bulk);
}

namespace roster_only {
void destroy_meshobj(stk::mesh::Entity, MeshBulkData& meshbulk );
}

bool insert_relation( stk::mesh::Entity , const stk::mesh::RelationType,  stk::mesh::Entity , const unsigned, const unsigned, const bool, MeshBulkData &);
bool remove_relation(stk::mesh::Entity , const stk::mesh::RelationIterator, MeshBulkData &);
}
}

namespace sierra {
  namespace Fmwk {
extern const stk::mesh::RelationIterator INVALID_RELATION_ITR;
  }
}


#endif


#endif /* stk_mesh_Combo_hpp */
