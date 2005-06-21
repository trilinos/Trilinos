

#include "Stepper_ForwardEuler.hpp"

namespace Rythmos {
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

