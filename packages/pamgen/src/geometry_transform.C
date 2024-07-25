// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "geometry_transform.h"
#include <sstream>

namespace PAMGEN_NEVADA {

/*****************************************************************************/
void Geometry_Transform::Display_Class(std::ostream& s, const std::string &/* indent */)
/*****************************************************************************/
{
  s << std::endl 
    << "Evaluation information for user defined run-time compiled " << std::endl 
    << "Geometry Transformation" << std::endl;

}


/*****************************************************************************/
Geometry_Transform::Geometry_Transform(const std::string & funcBody,
				       std::stringstream & error_string)
/*****************************************************************************/
:  _funcBody(funcBody),
  _function(6)//passing in by reference coord(input value) and field(output value)
{
  _function.addVar("double", "inxcoord");
  _function.addVar("double", "inycoord");
  _function.addVar("double", "inzcoord");
  _function.addVar("double", "outxcoord");
  _function.addVar("double", "outycoord");
  _function.addVar("double", "outzcoord");
  bool success = _function.addBody(_funcBody);
  if (!success) {
    error_string << "User_Defined_Geometry_Transform::User_Defined_Geometry_Transform: "
      << "function body error: " << _function.getErrors();
  }
}


/*****************************************************************************/
Geometry_Transform::~Geometry_Transform()
/*****************************************************************************/
{
 
}


/*****************************************************************************/
void Geometry_Transform::Operate(double * coords, long long num_nodes,long long dim)
/*****************************************************************************/
{

  double ivar1,ivar2,ivar3,ovar1,ovar2,ovar3;
  double * pivar1 = &ivar1;
  double * pivar2 = &ivar2;
  double * pivar3 = &ivar3;
  double * povar1 = &ovar1;
  double * povar2 = &ovar2;
  double * povar3 = &ovar3;
  
  for(long long i = 0; i < num_nodes; i ++){
    *pivar1 = coords[i];
    *povar1 = coords[i];
    *pivar2 = coords[i+num_nodes];
    *povar2 = coords[i+num_nodes];
    if(dim==3)   *pivar3 = coords[i+2*num_nodes];
    if(dim==3)   *povar3 = coords[i+2*num_nodes];

    _function.varAddrFill(0,pivar1);
    _function.varAddrFill(1,pivar2);
    _function.varAddrFill(2,pivar3);

    _function.varAddrFill(3,povar1);
    _function.varAddrFill(4,povar2);
    _function.varAddrFill(5,povar3);

      _function.execute();

    coords[i]             = *povar1;
    coords[i+num_nodes]   = *povar2;
    if(dim==3)coords[i+2*num_nodes] = *povar3;
  }
}


}//end namespace
