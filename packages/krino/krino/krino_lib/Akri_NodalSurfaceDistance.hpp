#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODALSURFACEDISTANCE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODALSURFACEDISTANCE_HPP_

namespace krino { class Composite_Surface; }
namespace krino { class FieldRef; }
namespace stk { namespace mesh { class BulkData; } }

namespace krino {

void compute_nodal_surface_distance(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    Composite_Surface & surfaces,
    const double time=0,
    const double narrowBandSize=0);

void compute_nodal_distance_from_spheres(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    const std::vector<std::pair<stk::math::Vector3d,double>> & spheres,
    const int sign = 1);

void compute_nodal_distance_from_plane(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    const stk::math::Vector3d & normal,
    const double offset);

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODALSURFACEDISTANCE_HPP_ */
