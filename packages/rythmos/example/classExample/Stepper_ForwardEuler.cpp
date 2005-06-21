

#include "Stepper_ForwardEuler.hpp"

namespace Rythmos {
ForwardEuler::ForwardEuler() 
{
}
ForwardEuler::ForwardEuler(ModelEvaluator *model)
{
  model_ = model;
  t_ = 0.0;
  x_ = model_->get_vector();
}
ForwardEuler::~ForwardEuler() 
{
}
double ForwardEuler::TakeStep(double dt)
{
  f_ = model_->evalModel(x_,t_);
  x_ = x_ + dt*f_;
  t_ = t_ + dt;
  return(dt);
}
double ForwardEuler::get_solution()
{
  return(x_);
}
} // namespace Rythmos

