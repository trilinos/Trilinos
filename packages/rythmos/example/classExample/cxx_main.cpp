
// This example is to debug why my code isn't compiling.
// 1.  Set up basic class structure here all in one file
//     virtual class Stepper with concrete derived class ForwardEuler
//     virtual class NonlinearModel with concrete derived class ExampleApplicationRythmosInterface
//     in main:  create object of type ExampleApplicationRythmosInterface and pass to constructer of ForwardEuler
//     Done.  My mistake was not including default constructors & destructors for virtual base classes.
// 2.  Add in RefCountPtr
//     Done.  If I pass by const reference, with const object, then the object
//     needs to have const versions of its member functions in order to call
//     them.
// 3.  Separate into namespaces
//     Done.  No errors.
// 4.  Add in templating
//     Done.  No errors.
// 5.  Separate into files

namespace Rythmos {
template<class Scalar>
class Stepper
{
  public:
    Stepper() {};
    virtual ~Stepper() {};
    virtual Scalar TakeStep(Scalar dt) = 0;
    virtual Scalar get_solution() = 0;
};
} // namespace Rythmos

namespace Rythmos {
template<class Scalar>
class ModelEvaluator
{
  public:
    ModelEvaluator() {};
    virtual ~ModelEvaluator() {};
    virtual Scalar evalModel(Scalar x, Scalar t) const = 0;
    virtual Scalar get_vector() const = 0;
};
} // namespace Rythmos

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Rythmos {
template<class Scalar>
class ForwardEuler : public Stepper<Scalar>
{
  public:
    ForwardEuler();
    ForwardEuler(const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model);
    ~ForwardEuler();
    Scalar TakeStep(Scalar dt);
    Scalar get_solution();
  protected:
    Scalar t_;
    Scalar x_;
    Scalar f_;
    Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > model_;
};
template<class Scalar>
ForwardEuler<Scalar>::ForwardEuler() 
{
}
template<class Scalar>
ForwardEuler<Scalar>::ForwardEuler(const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model)
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

class LinearProblem : public Rythmos::ModelEvaluator<double>
{
  public:
    LinearProblem();
    ~LinearProblem();
    double evalModel(double x, double t) const;
    double get_vector() const;
  protected:
    double lambda_;
};
LinearProblem::LinearProblem() 
{ 
  lambda_ = -0.5;
}
LinearProblem::~LinearProblem() 
{ 
}
double LinearProblem::evalModel(double x, double t) const
{
  return(lambda_*x);
}
double LinearProblem::get_vector() const
{
  return(1.0);
}

#include <iostream>
#include <cmath>

int main(int argc, char *argv[])
{
  Teuchos::RefCountPtr<LinearProblem> problem = Teuchos::rcp(new LinearProblem());
  Teuchos::RefCountPtr<Rythmos::ForwardEuler<double> > stepper = Teuchos::rcp(new Rythmos::ForwardEuler<double>(problem));

  double t0 = 0.0;
  double t1 = 1.0;
  int N = 10;
  double dt = (t1-t0)/N;
  for (int i=1 ; i<=N ; ++i)
  {
    stepper->TakeStep(dt);
  }
  std::cout << "Computed x = " << stepper->get_solution() << std::endl;
  std::cout << "Exact    x = " << 1.0*std::exp(-0.5*1.0) << std::endl;
  return(0);
}

