#ifndef MESH_QUALITY_IMPROVER_OPTS
#define MESH_QUALITY_IMPROVER_OPTS

namespace stk {
namespace middle_mesh {

struct MeshQualityImproverOpts
{
    int nlayers = -1; // tolerance for when to stop quality improvement.  When the
                      // maximum distance a vertex moves is less than this quantity,
                      // the optimization will terminate, unless invalid elements are
                      // still present.

    double maxDeltaX = 0.01; // Parameter to limit the section of the mesh the quality
                             // improver is applied to.  If the value is positive, vertices
                             // that can be reached from the boundary within nlayers edges
                             // will be moved by the quality improver

    int itermax = 100; // maximum number of steps the quality optimization will take, if
                       // all elements are valid.  The quality optimization can exceed
                       // this number of steps if invalid elements still exist

    double delta = 0.1; // parameter for the mesh quality objective.  Smaller values prioritize
                         // fixing invalid elements over improving the average quality

    bool verboseOutput = false;
};

} // namespace middle_mesh
} // namespace stk
#endif
