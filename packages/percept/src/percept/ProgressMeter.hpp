// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef PROGRESS_METER_HPP
#define PROGRESS_METER_HPP

#include <iostream>
#include <string>

#include <percept/Observable.hpp>

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

  virtual void notify(ProgressMeterData *data) override;

};


}

#endif
