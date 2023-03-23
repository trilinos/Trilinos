#ifndef ADJACENCY_SEARCH_H
#define ADJACENCY_SEARCH_H

#include <memory>
#include <queue>
#include <set>
#include <utility>
#include <vector>

#include "field.hpp"
#include "mesh.hpp"
#include "predicates/adjacency_predicates.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// This class provides the ability to iterate over two meshes simultaneously.
// Specifically, each call to getNext() gives an element on mesh1 and
// a vector of elements of mesh2 that are guaranteed to cover the element
// on mesh1 ("cover" mean that every point in the element on mesh1 lies within
// some element in the vector, or on the boundary of the elements in the vector).
// This class does *not* guarantee that vector of elements is the minimal
// set elements required to cover the element on mesh 1, but does guarantee
// that this vector can be generated a time that is independent of the size
// of the mesh.  Thus iterating over all elements on mesh1 and getting
// the covering elements on mesh2 can be done in O(n) time, where n is the
// number of elements in mesh.
//
// It is required that mesh2 cover mesh1 (but the reverse is not required)
class AdjacencySearch
{
  public:
    AdjacencySearch(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2)
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_preds(mesh2)
    {
      m_field                   = mesh::create_field<int_least8_t>(mesh1, FieldShape(0, 0, 1), 1, false);
      auto p                    = get_first_element();
      (*m_field)(p.first, 0, 0) = true;
      m_mesh1Queue.push(p);
    }

    // gets the next element on mesh1, and overwrites els2 with elements
    // on mesh2 that cover the returned element
    MeshEntityPtr get_next(std::vector<MeshEntityPtr>& els2);

  private:
    template <typename T>
    using SetType = std::set<T, MeshEntityCompare>;

    std::pair<MeshEntityPtr, MeshEntityPtr> get_first_element();

    // Inputs:
    //   el1: element on mesh1,
    //   el2: a nearby element on mesh2
    // Outputs:
    //   els2: vector to be overwritten with elements on mesh2 that
    //         cover el1
    void get_next(MeshEntityPtr el1, MeshEntityPtr el2, std::vector<MeshEntityPtr>& els2);

    // Gets an element on mesh2 that has some overlap with e1.
    // Currently checks for an element on mesh2 that contains el1s centroid
    // or has a vertex contained within el1
    // Inputs:
    //   el1: the element on mesh1
    //   el2: a starting guess for el2
    // Output: the element on mesh2
    MeshEntityPtr get_mesh2_element(MeshEntityPtr el1, MeshEntityPtr el2);

    // gets the set of all vertices contained within el1
    // Inputs:
    //   el1: element on mesh1
    //   v: a vertex on mesh2 within el1
    // Inputs/Outputs:
    //   verts: set of all vertices within el1.  Set will be cleared before
    //          adding vertices
    void get_contained_vertices(MeshEntityPtr el1, MeshEntityPtr v, SetType<MeshEntityPtr>& verts);

    // Gets the (unique) set of elements that are upward adjacent of verts
    // Inputs:
    //   verts: set of vertices
    // Outputs:
    //   els: vector to be overwritten with elements
    void get_elements_from_vertices(MeshEntityPtr el1, const SetType<MeshEntityPtr>& verts,
                                    std::vector<MeshEntityPtr>& els);

    // Does an edge-based element adjacency search
    // Overwrites els with the (unique) set of elements that cover
    // el1
    // Inputs:
    //   el1: element on mesh1
    //   edge2: an edge on mesh2 that intersects an edge of el1
    void do_edge_search(MeshEntityPtr el1, MeshEntityPtr edge2, std::vector<MeshEntityPtr>& els);

    void update_mesh1_queue(MeshEntityPtr el1, const std::vector<MeshEntityPtr>& els);

    // Given a set of candidate elements, returns the one with the
    // centroid closest to the given element
    // Inputs:
    //   el1: the element
    //   els2: the candidate elements
    // Outputs:
    //   one of els2
    MeshEntityPtr get_nearest_element(MeshEntityPtr el1, const std::vector<MeshEntityPtr>& els2);

    std::shared_ptr<Mesh> m_mesh1;
    std::shared_ptr<Mesh> m_mesh2;
    std::queue<std::pair<MeshEntityPtr, MeshEntityPtr>> m_mesh1Queue;
    // note: Field<bool> doesn't work because std::vector<bool> has
    // abnormal behavior
    FieldPtr<int_least8_t> m_field; // records whether elements on
                                    // mesh1 have been put into the queue
    predicates::impl::AdjacencyPredicates m_preds;
    const bool m_output = false;
};

// compute distance squared between two points
double compute_d2(const utils::Point& p1, const utils::Point& p2);

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
