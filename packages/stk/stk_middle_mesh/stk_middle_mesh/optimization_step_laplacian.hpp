#ifndef OPTIMIZATION_STEP_LAPLACIAN
#define OPTIMIZATION_STEP_LAPLACIAN

#include "optimization_step.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

class OptimizationStepLaplacian : public OptimizationStep
{
  public:
    virtual ~OptimizationStepLaplacian() {}

  private:
    utils::Point compute_updated_point(ActiveVertData& active) override;

    PatchObjective& get_objective() override { return m_obj; }

    PatchEnergyObjective m_obj;
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
