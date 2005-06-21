
// This example is to debug why my code isn't compiling.
// 1.  Set up basic class structure here all in one file
//     virtual class Stepper with concrete derived class ForwardEuler
//     virtual class NonlinearModel with concrete derived class LinearProblem
//     in main:  create object of type LinearProblem and pass to constructer of ForwardEuler
//     Done.  My mistake was not including default constructors & destructors for virtual base classes.
// 2.  Add in RefCountPtr
//     Done.  Note:  if I pass by const reference, with const object, then the object
//     needs to have const versions of its member functions in order to call
//     them.
// 3.  Separate into namespaces
//     Done.  No errors.
// 4.  Add in templating
//     Done.  No errors.
// 5.  Separate into files
//     Problems!  
// 6.  Get rid of RefCountPtr
//     Still have problems!
// 7.  Get rid of templates
//     This was the problem!  Now it works.
// 8.  Put templates back and isolate the problem.


#include "Stepper_ForwardEuler.hpp"

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
  LinearProblem *problem = new LinearProblem();
  Rythmos::ForwardEuler<double> *stepper = new Rythmos::ForwardEuler<double>(problem);

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
  delete stepper;
  delete problem;
  return(0);
}

