#include <stk_percept/ProgressMeter.hpp>
#include <unistd.h>

namespace stk
{
namespace percept
{

const char* ProgressMeterData::m_stateStrings[3] = {"INIT", "RUNNING", "FINI"};

ProgressMeter::ProgressMeter(Observable<ProgressMeterData>& observable) :  Observer<ProgressMeterData>(observable) { }

void ProgressMeter::notify(ProgressMeterData *data) 
{
  static std::string meter =  "[===================================================================================================]";  // length 100
  static std::string blanks = "[                                                                                                   ]";

  int percent_done = (int)data->m_data;

  if (data->m_state == ProgressMeterData::INIT)
  {
    percent_done = 0;
  }
  else if (data->m_state == ProgressMeterData::RUNNING)
  {
  }
  else if (data->m_state == ProgressMeterData::FINI)
  {
    percent_done = 100;
  }
  std::cout << "\r" << meter.substr(0,1+percent_done) << blanks.substr(1+percent_done,100-percent_done) << " " << percent_done << " [%] " << data->m_stage;
  if (data->m_state == ProgressMeterData::FINI)
  {
    std::cout << std::endl;
  }

  std::cout << std::flush;
  //sleep(1);

}

}
}
