#ifndef PATCH_OBJECTIVE_H
#define PATCH_OBJECTIVE_H

#include <vector>

#include "active_vert_data.hpp"
#include "change_of_basis.hpp"
#include "matrix.hpp"
#include "taylor_patch.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

class PatchObjective
{
  public:
    explicit PatchObjective() {}

    virtual ~PatchObjective() {}

    void set_active_patch(ActiveVertData& data)
    {
      // TODO: in many cases this gets calls 2 or 3 times with the same ActiveVertData,
      //       so there is no need to recompute everything.  Unfortunately it sometimes
      //       gets called with an ActiveVertData that lives on the stack (MeshQualityImprover::getRoots()),
      //       which has the same address every time even though the data is different.
      m_data = &data;
      set_parameterization(data);
      set_patch(data);
      set_active_patch_d(data);
    }

    // computes quality when the active vert is at the given (parametric)
    // point
    virtual double compute_quality(const utils::Point& ptIn) = 0;

    // computes derivative of computeQuality wrt the active vert only
    // Computes derivative wrt the parameterized coords
    // pt is the current position of the active vert
    virtual utils::Point compute_quality_rev(const utils::Point& ptIn, double q3Bar = 1) = 0;

    // computes the Hessian of the quality metric wrt the parameterized
    // coordinates of the active vertex
    // pt_in is the current position of the active vert
    virtual utils::impl::Mat2x2<double> compute_hessian(const utils::Point& ptIn) = 0;

    // given a point in xyz space, computes the parameterized coordinates
    utils::Point compute_parameterization(const utils::Point& pt);

    // given a parameterized pt, computes the point in the original
    // xyz coordinate system
    utils::Point compute_inverse_parameterization(const utils::Point& pt);

    bool has_invalid();

  protected:
    // for derived class additional work
    virtual void set_active_patch_d(ActiveVertData& data) {}

    // get the data for the patch under the given parameterization,
    // setting the current vert position to the given coordinates
    void compute_parameterization(const utils::Point& pt, utils::impl::Matrix<utils::Point>& triPts,
                                  std::vector<utils::impl::Mat2x2<double>>& wMats);

    ActiveVertData* m_data = nullptr;

  private:
    // set which plane to project onto
    void set_parameterization(ActiveVertData& data);

    void set_patch(ActiveVertData& data);

    // utils::impl::PlaneProjection m_proj;
    utils::impl::ChangeOfBasis m_proj;
    utils::impl::TaylorPatch m_patch;
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
