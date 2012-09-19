#ifndef PROGRESS_METER_HPP
#define PROGRESS_METER_HPP

#include <iostream>
#include <string>

#include <stk_percept/Observable.hpp>

namespace stk
{
namespace percept
{

struct ProgressMeterData
{
  enum STATE
  {
    INIT,
    RUNNING,
    FINI
  };

  ProgressMeterData(STATE state, double data, std::string stage="") : m_state(state), m_data(data), m_stage(stage) {}

  //static std::string *m_stateStrings;
  static const char* m_stateStrings[3];
  STATE m_state;
  double m_data;
  std::string m_stage;
};



class ProgressMeter : public Observer<ProgressMeterData>
{
public:
  ProgressMeter(Observable<ProgressMeterData>& observable);

  virtual void notify(ProgressMeterData *data) ;

};


}
}

#endif
