#ifndef NONCONFORMAL_FOUR_H
#define NONCONFORMAL_FOUR_H

#include "edge_tracer.hpp"
#include "geometry_improver_factory.hpp"
#include "mesh_relational_data.hpp"
#include "nonconformal4_opts.hpp"
#include "nonconformal_abstract.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"
#include "variable_size_field.hpp"
#include "middle_mesh_point_projection.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

// inverts a the field field_ptr, defined on mesh_in, creating a new field on mesh_other.
// field_ptr must give the elements on mesh_other
mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
invert_classification_field(std::shared_ptr<mesh::Mesh> meshIn, std::shared_ptr<mesh::Mesh> meshOther,
                            mesh::FieldPtr<mesh::MeshEntityPtr> fieldPtr);

class Nonconformal4 : public nonconformal::impl::NonconformalAbstract
{
  public:
    Nonconformal4(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                  const NormalProjectionOpts& opts = NormalProjectionOpts())
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_meshIn(mesh::make_empty_mesh(mesh1->get_comm()))
      , m_relationalData(std::make_shared<MeshRelationalData>(mesh1, mesh2, m_meshIn))
      , m_classifier(std::make_shared<predicates::impl::PointClassifierNormalWrapper>(mesh2, opts.classifierTolerances))
      , m_edgeTracerTolerances(opts.edgeTracerTolerances)
      , m_geometryImprovers(opts.geometryImprovers)
    {}

    std::shared_ptr<mesh::Mesh> create() override;

    std::shared_ptr<mesh::Mesh> get_mesh_in() override { return m_meshIn; }

    // returns a field that maps the elements of mesh_in to mesh1
    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh1_classification() override
    {
      return m_relationalData->meshInElementsToMesh1Elements;
    }

    // returns a field that maps the elements of mesh_in to mesh2
    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh2_classification() override
    {
      return m_relationalData->meshInElementsToMesh2Elements;
    }

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_inverse_classification() override
    {
      return impl::invert_classification_field(m_meshIn, m_mesh1, m_relationalData->meshInElementsToMesh1Elements);
    }

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_inverse_classification() override
    {
      return impl::invert_classification_field(m_meshIn, m_mesh2, m_relationalData->meshInElementsToMesh2Elements);
    }

    mesh::FieldPtr<utils::Point> compute_points_on_mesh1(std::shared_ptr<XiCoordinates> pts) override
    {
      MiddleMeshPointProjection projector(m_mesh1, m_mesh2, m_meshIn, m_relationalData, m_classifier);
      return projector.projection_onto_mesh1(pts);
    }

    mesh::FieldPtr<utils::Point> compute_points_on_mesh2(std::shared_ptr<XiCoordinates> pts) override
    {
      MiddleMeshPointProjection projector(m_mesh1, m_mesh2, m_meshIn, m_relationalData, m_classifier);
      return projector.projection_onto_mesh2(pts);
    }

  private:
    void apply_geometry_improvers();

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<MeshRelationalData> m_relationalData;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_classifier;
    middle_mesh::impl::EdgeTracerTolerances m_edgeTracerTolerances;
    std::vector<GeometryImprovers> m_geometryImprovers;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
