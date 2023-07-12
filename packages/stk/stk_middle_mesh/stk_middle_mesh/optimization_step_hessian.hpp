#ifndef OPTIMIZATION_STEP_HESSIAN
#define OPTIMIZATION_STEP_HESSIAN

#include "backtracking_line_search.hpp"
#include "optimization_step.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

class OptimizationStepHessian : public OptimizationStep
{
  public:
    OptimizationStepHessian(std::shared_ptr<PatchObjective> obj)
      : m_obj(obj)
      , m_linesearch(0.5, 0.1)
    {}

    virtual ~OptimizationStepHessian() {}

  private:
    // returns updated position of the active vert in parametric space
    utils::Point compute_updated_point(ActiveVertData& active) override;

    PatchObjective& get_objective() override { return *m_obj; }

    // modifies A such that all eigenvalues >= delta
    void modify_eigenvalues(utils::impl::Mat2x2<double>& a, const double delta);

    std::shared_ptr<PatchObjective> m_obj;
    BacktrackingLineSearch m_linesearch;
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
