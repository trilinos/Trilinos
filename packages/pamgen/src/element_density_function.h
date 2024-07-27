// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef element_density_functionH
#define element_density_functionH

#include "RTC_FunctionRTC.hh"
namespace PAMGEN_NEVADA {

class Element_Density_Function{
public:
  Element_Density_Function(const std::string & funcBody,
			   const std::string & direction,
			   std::stringstream & error_stream);
  ~Element_Density_Function();
  double Interpolate(double, std::stringstream & error_stream);
  void Integrate(double start_var, double end_var,std::stringstream & error_stream);
  void Display_Class(std::ostream&, const std::string &indent);
  void deleteRunningSum(){
    if(running_sum){
      delete [] running_sum;
      running_sum = NULL;
    }};
  
private:
  std::string                     _funcBody;
  std::string                     application_direction;
  PG_RuntimeCompiler::Function  _function;
  double * running_sum;
  static const long long running_sum_length = 10000;
  bool integrated; // has this been integrated
  double integral_total;//Total integration sum
  double max_eval;//maximum finction evaluates to during integration
  double min_eval;//minimum ...
  double min_eval_range;//minimum range of evaluation of the function
  double max_eval_range;//maximum ...

};

}// end namespace PAMGEN_NEVADA

#endif
