
// This example is to debug why my code isn't compiling.
// 1.  Set up basic class structure here all in one file
//     virtual class Stepper with concrete derived class ForwardEuler
//     virtual class NonlinearModel with concrete derived class ExampleApplicationRythmosInterface
//     in main:  create object of type ExampleApplicationRythmosInterface and pass to constructer of ForwardEuler
// 2.  Add in RefCountPtr
// 3.  Separate into namespaces
// 4.  Add in templating
// 5.  Separate into files


class Stepper
{
  public:
    Stepper() {};
    virtual ~Stepper() {};
    virtual double TakeStep(double dt) = 0;
    virtual double get_solution() = 0;
};

class ModelEvaluator
{
  public:
    ModelEvaluator() {};
    virtual ~ModelEvaluator() {};
    virtual double evalModel(double x, double t) = 0;
    virtual double get_vector() = 0;
};

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

class LinearProblem : public ModelEvaluator
{
  public:
    LinearProblem();
    ~LinearProblem();
    double evalModel(double x, double t);
    double get_vector();
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
double LinearProblem::evalModel(double x, double t) 
{
  return(lambda_*x);
}
double LinearProblem::get_vector()
{
  return(1.0);
}

#include <iostream>
#include <cmath>

int main(int argc, char *argv[])
{
  LinearProblem *problem = new LinearProblem();
  ForwardEuler *stepper = new ForwardEuler(problem);

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
  delete(stepper);
  delete(problem);
  return(0);
}

