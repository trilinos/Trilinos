
#ifndef RYTHMOS_STEPPER_FORWARDEULER
#define RYTHMOS_STEPPER_FORWARDEULER

#include "Stepper.hpp"
#include "ModelEvaluator.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Rythmos {
template<class Scalar>
class ForwardEuler : public Stepper<Scalar>
{
  public:
    ForwardEuler();
    ForwardEuler(ModelEvaluator<Scalar> *model);
    ~ForwardEuler();
    Scalar TakeStep(Scalar dt);
    Scalar get_solution();
  protected:
    Scalar t_;
    Scalar x_;
    Scalar f_;
    ModelEvaluator<Scalar> *model_;
};
} // namespace Rythmos

#endif // RYTHMOS_STEPPER_FORWARDEULER
