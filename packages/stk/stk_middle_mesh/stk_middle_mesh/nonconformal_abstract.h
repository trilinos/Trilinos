#ifndef NONCONFORMAL_ABSTRACT_H
#define NONCONFORMAL_ABSTRACT_H

#include "field.h"
#include "mesh.h"
#include "variable_size_field.h"

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
};

} // namespace impl

} // namespace nonconformal
} // namespace middle_mesh
} // namespace stk
#endif
