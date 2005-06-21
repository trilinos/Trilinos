
#ifndef RYTHMOS_STEPPER_FORWARDEULER
#define RYTHMOS_STEPPER_FORWARDEULER

#include "Stepper.hpp"
#include "ModelEvaluator.hpp"

namespace Rythmos {
class ForwardEuler : public Stepper
{
  public:
    ForwardEuler();
    ForwardEuler(ModelEvaluator *model);
    ~ForwardEuler();
    double TakeStep(double dt);
    double get_solution();
  protected:
    double t_;
    double x_;
    double f_;
    ModelEvaluator *model_;
};
} // namespace Rythmos

#endif // RYTHMOS_STEPPER_FORWARDEULER
