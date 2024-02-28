#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODALBOUNDINGBOX_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODALBOUNDINGBOX_HPP_

#include <Akri_BoundingBox.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Entity; } }

namespace krino {

class FieldRef;

BoundingBox compute_nodal_bbox( const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector, const FieldRef coordsField );
BoundingBox compute_nodal_bbox( const stk::mesh::BulkData & mesh, const FieldRef coordsField, const std::vector<stk::mesh::Entity> & nodes );

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODALBOUNDINGBOX_HPP_ */
