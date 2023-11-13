#ifndef MESH_QUALITY_IMPROVER_H
#define MESH_QUALITY_IMPROVER_H

#include "optimization_step.hpp"
#include "active_vert_container.hpp"

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

namespace {

class IsVertValidPredicate
{
  public:
    explicit IsVertValidPredicate(std::shared_ptr<Mesh> mesh,
                                  std::shared_ptr<opt::impl::OptimizationStep> quality) :
      m_mesh(mesh),
      m_quality(quality)
    {}

    bool operator()(MeshEntityPtr v)
    {
      assert(v->get_type() == MeshEntityType::Vertex);
      opt::impl::ActiveVertData active(m_mesh, v);
      return !(m_quality->has_invalid(active));
    }

  private:
    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<opt::impl::OptimizationStep> m_quality;
};

}

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
                        const int itermax = 10, bool verboseOutput = false)
      : m_mesh(mesh)
      , m_quality(qualityOpt)
      , m_activeVertContainer(mesh, filter, IsVertValidPredicate(mesh, qualityOpt), nlayers)
      , m_distTol(tol)
      , m_itermax(itermax)
      , m_verboseOutput(verboseOutput)
    {
    }

    void run();

    // returns the number of vertices that have at least one invalid
    // virtual element attached to them
    int count_invalid_points();

    bool verbose_output() const
    {
      return m_verboseOutput;
    }

  private:

    void run_single(std::shared_ptr<opt::impl::OptimizationStep> step, std::vector<opt::impl::ActiveVertData*>& verts,
                    const RunOpts& opts);

    // overwrites invalid_verts with pointers to the invalid verts
    void get_invalid_verts(std::shared_ptr<opt::impl::OptimizationStep> opt,
                           std::vector<opt::impl::ActiveVertData*>& invalidVerts);

    void get_all_verts(std::vector<opt::impl::ActiveVertData*>& verts);

    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<opt::impl::OptimizationStep> m_quality;
    ActiveVertContainer m_activeVertContainer;
    double m_distTol = 1e-6; // TODO: needs to be relataive to the size of the element
    int m_itermax    = 10;
    bool m_verboseOutput = false;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
