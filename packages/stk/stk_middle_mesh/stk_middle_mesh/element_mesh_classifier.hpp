#ifndef ELEMENT_MESH_CLASSIFIER_H
#define ELEMENT_MESH_CLASSIFIER_H

#include "element_mesh_extractor.hpp"
#include "mesh_agglomerator.hpp"
#include "mesh_entity.hpp"
#include "mesh_relational_data.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"

#include <map>
#include <set>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

// given an element mesh and its classification on mesh1, compute the
// classification on mesh2
class ElementMeshClassifier
{
  public:
    ElementMeshClassifier(std::shared_ptr<MeshRelationalData> relationalData,
                          std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> classifier)
      : m_relationalData(relationalData)
      , m_classifier(classifier)
    {}

    void classify(ElementMeshData& elementMeshData, int numConstraintEdges, const std::vector<mesh::MeshEntityPtr>& mesh2Els);

  private:
    template <typename T>
    using SetType = std::set<T, mesh::MeshEntityCompare>;

    template <typename Key, typename Value>
    using MapType = std::map<Key, Value, mesh::MeshEntityCompare>;

    void get_mesh2_element_classification(int numConstraintEdges);

    mesh::impl::MeshAgglomerator setup_agglomerations(int numConstraintEdges);

    void get_mesh2_edges(std::vector<mesh::MeshEntityPtr>& edges2);

    bool is_classified_on_element(const predicates::impl::PointRecord& record);

    bool need_numerical_classification();

    mesh::MeshEntityPtr classification_pass1(const std::vector<mesh::MeshEntityPtr>& mesh2Els,
                                             mesh::impl::MeshAgglomerator& agg,
                                             std::vector<mesh::MeshEntityPtr>& els2List);

    // TODO: DEBUGGING
    mesh::MeshEntityPtr get_element_mesh_entity(mesh::MeshEntityPtr meshInEntity);

    void classification_pass2(const std::vector<mesh::MeshEntityPtr>& els2List, mesh::impl::MeshAgglomerator& agg);

    void classification_pass3(const std::vector<mesh::MeshEntityPtr>& els2List, mesh::impl::MeshAgglomerator& agg);

    void classify_group(mesh::impl::MeshAgglomerator& agg, int group, mesh::MeshEntityPtr el2);

    // returns true if the third vert is on el2
    bool classify_third_vert(mesh::MeshEntityPtr el2, SetType<mesh::MeshEntityPtr>& vertsIn, mesh::MeshEntityPtr elIn);

    // get all verts on mesh_in that are associated with edges
    // of the given element on mesh2
    void get_verts_in(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2, SetType<mesh::MeshEntityPtr>& vertsIn);

    std::vector<int> remove_already_classified_groups(const std::vector<int>& groups);

    // gets elements that have nthres or more vertices from verts_in
    void get_common_elements(const SetType<mesh::MeshEntityPtr>& vertsIn, const unsigned int nthres,
                             std::vector<mesh::MeshEntityPtr>& commonEls);

    mesh::MeshEntityPtr get_most_inside_element(mesh::MeshEntityPtr elIn, const std::vector<mesh::MeshEntityPtr>& els2);

    // classifies all remaining mesh_in elements on el2
    void classify_all(mesh::MeshEntityPtr el2);

    std::shared_ptr<MeshRelationalData> m_relationalData;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_classifier;
    ElementMeshData m_elementMeshData;
    std::vector<mesh::MeshEntityPtr> m_mesh2Els;
    std::vector<bool> m_agglomerationsClassified;
    static constexpr bool M_OUTPUT = false;
};

bool are_all_verts_on_same_edge(std::set<mesh::MeshEntityPtr, mesh::MeshEntityCompare>& vertsIn,
                                mesh::FieldPtr<predicates::impl::PointRecord> vertsInClassOnMesh1Ptr);

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
