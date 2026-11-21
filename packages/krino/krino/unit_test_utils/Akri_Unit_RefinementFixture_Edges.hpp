#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_EDGES_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_EDGES_HPP_

#include <Akri_AllReduce.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace krino
{

// Adds explicit edge entities (not sides) to the 3D mesh of a refinement fixture, and checks them
// after each adaptation.  The fixture builds its mesh first; then create_edges() is called.
template <typename FIXTURE> class WithExplicitEdges : public FIXTURE
{
public:
  using FIXTURE::mMesh;
  using FIXTURE::myRefinement;

  // Creates all edges of the mesh, in edgePart if given, and records them as the original edges.
  void create_edges(stk::mesh::Part * edgePart = nullptr)
  {
    stk::mesh::create_edges(mMesh, mMesh.mesh_meta_data().universal_part(), edgePart);
    myAllOriginalEdgeSegments =
        gather_segments(get_owned_edge_segments(mMesh.mesh_meta_data().universal_part()));
    myHaveEdges = true;
  }

  void check_mesh_after_adaptation() override
  {
    if (myHaveEdges) check_edges_of_active_elements();
  }

  // Edges are only refined, never added: an element edge of an active element has an edge exactly
  // when it lies along an edge of the original mesh.  Every edge is attached to all the elements
  // that have it (check_face_and_edge_relations).
  void check_edges_of_active_elements()
  {
    const stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    const stk::mesh::Selector activeOwnedLeaves = meta.locally_owned_part() &
        this->get_aux_meta().active_part() & !myRefinement.parent_part();
    size_t numMissing = 0;
    size_t numUnexpected = 0;
    size_t numWrongMidNodes = 0;
    std::array<stk::mesh::Entity, 3> elemEdgeNodes;
    for (auto * bucketPtr : mMesh.get_buckets(stk::topology::ELEMENT_RANK, activeOwnedLeaves))
    {
      const stk::topology topo = bucketPtr->topology();
      for (auto elem : *bucketPtr)
      {
        for (unsigned iEdge = 0; iEdge < topo.num_edges(); ++iEdge)
        {
          topo.edge_nodes(mMesh.begin_nodes(elem), iEdge, elemEdgeNodes.data());
          const bool shouldHaveEdge =
              is_along_segment(elemEdgeNodes[0], elemEdgeNodes[1], myAllOriginalEdgeSegments);
          const stk::mesh::Entity edge =
              find_entity_by_ordinal(mMesh, elem, stk::topology::EDGE_RANK, iEdge);
          const bool hasEdge = mMesh.is_valid(edge);
          if (shouldHaveEdge && !hasEdge) ++numMissing;
          if (!shouldHaveEdge && hasEdge) ++numUnexpected;
          // The mid-edge node of a LINE_3 edge is the element's node on that edge.
          if (hasEdge && topo.edge_topology(iEdge).num_nodes() == 3 &&
              mMesh.begin_nodes(edge)[2] != elemEdgeNodes[2])
            ++numWrongMidNodes;
        }
      }
    }
    all_reduce_sum(mMesh.parallel(), numMissing, numUnexpected, numWrongMidNodes);
    EXPECT_EQ(0u, numMissing);
    EXPECT_EQ(0u, numUnexpected);
    EXPECT_EQ(0u, numWrongMidNodes);
    EXPECT_TRUE(check_face_and_edge_relations(mMesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mMesh));
  }

  size_t num_edges_in(const stk::mesh::Selector & selector) const
  {
    std::vector<size_t> counts;
    stk::mesh::count_entities(
        selector & mMesh.mesh_meta_data().locally_owned_part(), mMesh, counts);
    size_t num = counts[stk::topology::EDGE_RANK];
    all_reduce_sum(mMesh.parallel(), num);
    return num;
  }

  size_t num_original_edges() const { return myAllOriginalEdgeSegments.size(); }

protected:
  bool is_along_segment(const stk::mesh::Entity node0,
      const stk::mesh::Entity node1,
      const std::vector<std::array<stk::math::Vector3d, 2>> & segments) const
  {
    const FieldRef coords(mMesh.mesh_meta_data().coordinate_field());
    const stk::math::Vector3d x0 = get_vector_field(mMesh, coords, node0, 3);
    const stk::math::Vector3d x1 = get_vector_field(mMesh, coords, node1, 3);
    for (const auto & [a, b] : segments)
    {
      const stk::math::Vector3d ab = b - a;
      const double len2 = ab.length_squared();
      auto onSegment = [&](const stk::math::Vector3d & x)
      {
        const double t = Dot(x - a, ab) / len2;
        return t > -1.e-12 && t < 1. + 1.e-12 && (a + t * ab - x).length_squared() < 1.e-24 * len2;
      };
      if (onSegment(x0) && onSegment(x1)) return true;
    }
    return false;
  }

  std::vector<std::array<stk::math::Vector3d, 2>>
  gather_segments(const std::vector<std::array<stk::math::Vector3d, 2>> & localSegments) const
  {
    std::vector<double> local;
    for (const auto & [a, b] : localSegments)
      for (const auto & x : {a, b})
        for (int d = 0; d < 3; ++d)
          local.push_back(x[d]);
    const int numProcs = stk::parallel_machine_size(mMesh.parallel());
    std::vector<int> counts(numProcs), displs(numProcs);
    const int localCount = local.size();
    MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, mMesh.parallel());
    int total = 0;
    for (int p = 0; p < numProcs; ++p)
    {
      displs[p] = total;
      total += counts[p];
    }
    std::vector<double> all(total);
    MPI_Allgatherv(local.data(),
        localCount,
        MPI_DOUBLE,
        all.data(),
        counts.data(),
        displs.data(),
        MPI_DOUBLE,
        mMesh.parallel());
    std::vector<std::array<stk::math::Vector3d, 2>> segments;
    for (int i = 0; i + 6 <= total; i += 6)
      segments.push_back({stk::math::Vector3d(all[i], all[i + 1], all[i + 2]),
          stk::math::Vector3d(all[i + 3], all[i + 4], all[i + 5])});
    return segments;
  }

  std::vector<std::array<stk::math::Vector3d, 2>>
  get_owned_edge_segments(const stk::mesh::Selector & selector) const
  {
    const FieldRef coords(mMesh.mesh_meta_data().coordinate_field());
    std::vector<std::array<stk::math::Vector3d, 2>> segments;
    for (auto * bucketPtr : mMesh.get_buckets(
             stk::topology::EDGE_RANK, selector & mMesh.mesh_meta_data().locally_owned_part()))
      for (auto edge : *bucketPtr)
      {
        const stk::mesh::Entity * nodes = mMesh.begin_nodes(edge);
        segments.push_back({get_vector_field(mMesh, coords, nodes[0], 3),
            get_vector_field(mMesh, coords, nodes[1], 3)});
      }
    return segments;
  }

  std::vector<std::array<stk::math::Vector3d, 2>> myAllOriginalEdgeSegments;
  bool myHaveEdges{false};
};

} // namespace krino

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_EDGES_HPP_ */
