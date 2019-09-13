// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/ProgressMeter.hpp>
#include <unistd.h>

namespace percept
{

const char* ProgressMeterData::m_stateStrings[3] = {"INIT", "RUNNING", "FINI"};

ProgressMeter::ProgressMeter(Observable<ProgressMeterData>& observable) :  Observer<ProgressMeterData>(observable) { }

void ProgressMeter::notify(ProgressMeterData *data) 
{
  //                           123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 1
  static std::string meter =  "[===================================================================================================]";  // length 101
  static std::string blanks = "[                                                                                                   ]";

  int percent_done = (int)data->m_data;

  if (data->m_state == ProgressMeterData::INIT)
  {
    percent_done = 0;
    std::cout << "Stage: " << data->m_stage << " % done: 0";
  }
  else if (data->m_state == ProgressMeterData::RUNNING)
  {
  }
  else if (data->m_state == ProgressMeterData::FINI)
  {
    percent_done = 100;
    std::cout << " 100";
  }
  if (percent_done < 0) { std::cout << "percent_done= " << percent_done << " data= " << data->m_data << std::endl; percent_done = 0; }
  if (percent_done > 100) { std::cout << "percent_done= " << percent_done << " data= " << data->m_data << std::endl; percent_done = 100; }
  //std::cout << "\r" << meter.substr(0,1+percent_done) << blanks.substr(1+percent_done,100-percent_done) << " " << percent_done << " [%] " << data->m_stage;
  //   size_t meter_end = std::min(1+percent_done, meter.size());
  //   size_t blanks_start = std::min(1+percent_done, blanks.size()-1);
  //   size_t blanks_end = std::min(100-percent_done, blanks.size()-blansk_start);

  //std::cout << "\r" << meter.substr(0,1+percent_done) << blanks.substr(1+percent_done,100-percent_done) << " " << percent_done << " [%] " << data->m_stage;
  if (percent_done > 0 && percent_done < 100) std::cout << " " << percent_done;
  if (data->m_state == ProgressMeterData::FINI)
  {
    std::cout << std::endl;
  }

  std::cout << std::flush;
  //sleep(1);

}

}
