#ifndef MESH_QUALITY_IMPROVER_H
#define MESH_QUALITY_IMPROVER_H

#include "mesh_layers.hpp"
#include "optimization_step.hpp"

// #include ".hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

struct RunOpts
{
    double maxDeltaX;
    int itermax;       // itermax
    bool trimValid;    // remove valid elements from set as they are
                       // discovered
    bool requireValid; // require all elements to be valid before returning,
                       // even if max_delta_x is satisfied
                       // (the itermax criteria is always applied)
};

class MeshQualityImprover
{
  public:
    // Tfunc is a function f(MeshEntityPtr) -> bool, returns true if
    // the given MeshEntity should be moved by the quality improver,
    // false otherwise
    // nlayers: determines which vertices will be moved by the quality
    //          improver.  Vertices within nlayers number of layers
    //          from an excluded vert (as determined by the function f)
    //          will be moved.  Set to -1 to move all vertices
    //          Note that f and nlayers work together to determine which
    //          vertices to move: a vertex must have both f(vert) == true
    //          *and* be within nlayers of an excluded vert to be moved
    // tol: if the maximum distance a node is moved in a single iteration
    //      is less than this, the problem is considered converged
    //      Note: the algorithm will continue iterating even if this
    //      criteria is satisifed when there are invalid elements present
    //      The itermax criteria will still terminate the algorithm even
    //      if invalid elements are detected
    // itermax: maximum number of iterations
    template <typename Tfunc>
    MeshQualityImprover(std::shared_ptr<Mesh> mesh, Tfunc filter, const int nlayers,
                        std::shared_ptr<opt::impl::OptimizationStep> qualityOpt, const double tol = 1e-6,
                        const int itermax = 10)
      : m_mesh(mesh)
      , m_quality(qualityOpt)
      , m_distTol(tol)
      , m_itermax(itermax)
    {
      get_active_verts(mesh, filter, nlayers);
      std::cout << "number of active verts = " << m_activeVerts.size() << std::endl;
    }

    void run();

    // returns the number of vertices that have at least one invalid
    // virtual element attached to them
    int count_invalid_points();

  private:
    template <typename Tfunc>
    void get_active_verts(std::shared_ptr<Mesh> mesh, Tfunc filter, const int nlayers)
    {
      std::vector<MeshEntityPtr> roots, verts;
      get_roots(mesh, filter, roots);

      MeshLayers layers;
      if (nlayers == -1)
        layers.get_all_layers(mesh, filter, roots, verts);
      else
        layers.get_layers(mesh, filter, roots, nlayers, verts);

      for (auto& v : verts)
        if (v && filter(v))
          m_activeVerts.emplace_back(v);
    }

    template <typename Tfunc>
    void get_roots(std::shared_ptr<Mesh> mesh, Tfunc filter, std::vector<MeshEntityPtr>& roots)
    {
      roots.clear();
      for (auto v : mesh->get_vertices())
      {
        if (v)
        {
          if (!filter(v))
          { // use the excluded verts as the roots
            roots.push_back(v);
          } else
          {
            opt::impl::ActiveVertData active(v);
            if (m_quality->has_invalid(active))
              roots.push_back(v);
          }
        }
      }

      // if no excluded verts, pick an arbitrary vertex
      if (roots.size() == 0)
        for (auto v : mesh->get_vertices())
          if (v)
          {
            roots.push_back(mesh->get_vertices()[0]);
            break;
          }
    }

    void run_single(std::shared_ptr<opt::impl::OptimizationStep> step, std::vector<opt::impl::ActiveVertData*>& verts,
                    const RunOpts& opts);

    // overwrites invalid_verts with pointers to the invalid verts
    void get_invalid_verts(std::shared_ptr<opt::impl::OptimizationStep> opt,
                           std::vector<opt::impl::ActiveVertData*>& invalidVerts);

    void get_all_verts(std::vector<opt::impl::ActiveVertData*>& verts);

    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<opt::impl::OptimizationStep> m_quality;
    std::vector<opt::impl::ActiveVertData> m_activeVerts;
    double m_distTol = 1e-6; // TODO: needs to be relataive to the size of the element
    int m_itermax    = 10;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
