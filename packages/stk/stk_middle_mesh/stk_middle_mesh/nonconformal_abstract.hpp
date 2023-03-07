#ifndef NONCONFORMAL_ABSTRACT_H
#define NONCONFORMAL_ABSTRACT_H

#include "field.hpp"
#include "mesh.hpp"
#include "variable_size_field.hpp"
#include "application_interface.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal {
namespace impl {




// TODO: error checking: call create() before other members
class NonconformalAbstract
{
  public:
    NonconformalAbstract() {}

    virtual ~NonconformalAbstract(){};

    // creates and returns the middle grid
    virtual std::shared_ptr<mesh::Mesh> create() = 0;

    // returns the middle grid created by the create() method
    virtual std::shared_ptr<mesh::Mesh> get_mesh_in() = 0;

    // returns a field that maps the elements of mesh_in to mesh1
    virtual mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh1_classification() = 0;

    // returns a field that maps the elements of mesh_in to mesh2
    virtual mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh2_classification() = 0;

    virtual mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_inverse_classification() = 0;

    virtual mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_inverse_classification() = 0;

    // given a set of xi coordinates for middle mesh elements, compute the xi coordinates of the corresponding
    // points on mesh1
    virtual mesh::FieldPtr<utils::Point> compute_points_on_mesh1(std::shared_ptr<XiCoordinates> pts) = 0;

    // similar to compute_points_on_mesh1(), but for mesh2
    virtual mesh::FieldPtr<utils::Point> compute_points_on_mesh2(std::shared_ptr<XiCoordinates> pts) = 0;    
};

} // namespace impl

} // namespace nonconformal
} // namespace middle_mesh
} // namespace stk
#endif
