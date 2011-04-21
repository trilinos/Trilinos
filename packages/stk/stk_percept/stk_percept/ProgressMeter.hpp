#include <iostream>

#include <stk_percept/Observable.hpp>

namespace stk
{
namespace percept
{

class ProgressMeter : public Observer<double>
{
public:
  ProgressMeter(Observable<double>& observable);

  virtual void notify(double *data) ;

};


}
}
