// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef geometry_transformH
#define geometry_transformH

#include "RTC_FunctionRTC.hh"
#include "pamgen_kokkos_utils.h"

namespace PAMGEN_NEVADA {

class Geometry_Transform{
public:
  Geometry_Transform(const std::string & funcBody,
		   std::stringstream & error_stream);
  ~Geometry_Transform();
  void Display_Class(std::ostream&, const std::string &indent);
  void Operate(double* coords, long long num_nodes,long long dim);
  KOKKOS_INLINE_FUNCTION
  void Operate_Device(View2D<double> coords, long long num_nodes, long long dim);
  
private:
  void parseTransformationMultipliers();
  
  
private:
  std::string                     _funcBody;
  std::string                     application_direction;
  PG_RuntimeCompiler::Function  _function;
  
  // Transformation multipliers for device execution
  double _x_multiplier;
  double _y_multiplier;
  double _z_multiplier;

};
}
#endif
