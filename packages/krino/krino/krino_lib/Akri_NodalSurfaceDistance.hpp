#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODALSURFACEDISTANCE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODALSURFACEDISTANCE_HPP_

namespace krino { class Composite_Surface; }
namespace krino { class FieldRef; }
namespace stk { namespace mesh { class BulkData; } }

namespace krino {

void compute_nodal_surface_distance(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distanceField, Composite_Surface & surfaces, const double time=0, const double narrowBandSize=0);

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODALSURFACEDISTANCE_HPP_ */
