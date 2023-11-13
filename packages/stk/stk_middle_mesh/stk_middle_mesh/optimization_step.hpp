#ifndef OPTIMIZATION_STEP_H
#define OPTIMIZATION_STEP_H

#include "active_vert_data.hpp"
#include "patch_energy_objective.hpp"
#include "patch_objective.hpp"

#include "mesh_io.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

// class that does a single pass of the optimization algorithm over
// the provided vertices
class OptimizationStep
{
  public:
    virtual ~OptimizationStep() {}

    // verts is an iterator that yields an ActiveVertData*
    // verts_end is the past-the-end iterator
    // trim_invalid: if true, if a vert is in a valid position after
    // the update, it will be overwritten with the nullptr
    template <typename Titer>
    double improve_quality(Titer vert, Titer vertsEnd, const bool trimInvalid)
    {
      double maxDeltaX = 0;
      while (vert != vertsEnd)
      {
        if (*vert)
        {
          ActiveVertData* activePtr = *vert;

          auto deltaX               = improve_quality_single(*activePtr);
          maxDeltaX                 = std::max(maxDeltaX, deltaX);

          if (trimInvalid && !has_invalid(*activePtr))
            *vert = nullptr;
        }

        vert++;
      }

      return maxDeltaX;
    }

    // returns true if the vert is in an invalid position
    // this method exists as a performance optimization: it avoids
    // the work of setActivePatch
    bool has_invalid(ActiveVertData& active)
    {
      PatchObjective& obj = get_objective();
      obj.set_active_patch(active);
      return obj.has_invalid();
    }

  private:
    virtual double improve_quality_single(ActiveVertData& active);

    // function to compute an updated vertex location
    virtual utils::Point compute_updated_point(ActiveVertData& active) = 0;

    virtual PatchObjective& get_objective() = 0;
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
