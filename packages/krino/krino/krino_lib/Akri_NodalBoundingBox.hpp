#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODALBOUNDINGBOX_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODALBOUNDINGBOX_HPP_

#include <Akri_BoundingBox.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }

namespace krino {

class FieldRef;

BoundingBox compute_nodal_bbox( const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector, const FieldRef coordsField );

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODALBOUNDINGBOX_HPP_ */
