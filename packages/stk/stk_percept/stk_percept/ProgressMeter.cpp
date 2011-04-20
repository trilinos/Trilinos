#include <stk_percept/ProgressMeter.hpp>

namespace stk
{
namespace percept
{

ProgressMeter::ProgressMeter(Observable<double>& observable) :  Observer<double>(observable) { }

void ProgressMeter::notify(double *data) 
{
  std::cout << "ProgressMeter:: notify" << std::endl;
}

}
}
