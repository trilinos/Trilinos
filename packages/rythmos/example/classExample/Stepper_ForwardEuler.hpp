
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


template<class Scalar>
ForwardEuler<Scalar>::ForwardEuler() 
{
}
template<class Scalar>
ForwardEuler<Scalar>::ForwardEuler(ModelEvaluator<Scalar> *model)
{
  model_ = model;
  t_ = Teuchos::ScalarTraits<Scalar>::zero();
  x_ = model_->get_vector();
}
template<class Scalar>
ForwardEuler<Scalar>::~ForwardEuler() 
{
}
template<class Scalar>
Scalar ForwardEuler<Scalar>::TakeStep(Scalar dt)
{
  f_ = model_->evalModel(x_,t_);
  x_ = x_ + dt*f_;
  t_ = t_ + dt;
  return(dt);
}
template<class Scalar>
Scalar ForwardEuler<Scalar>::get_solution()
{
  return(x_);
}

} // namespace Rythmos


#endif // RYTHMOS_STEPPER_FORWARDEULER
