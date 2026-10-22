// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "geometry_transform.h"
#include "pamgen_kokkos_utils.h"
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
void Geometry_Transform::parseTransformationMultipliers()
/*****************************************************************************/
{
  // Simple parser to extract multipliers from transformation strings
  // Looks for patterns like "outxcoord = inxcoord * 2.0;"
  
  // Default multipliers (identity transformation)
  _x_multiplier = 1.0;
  _y_multiplier = 1.0;
  _z_multiplier = 1.0;
  
  // Parse x multiplier
  size_t x_pos = _funcBody.find("outxcoord = inxcoord * ");
  if (x_pos != std::string::npos) {
    size_t start = x_pos + 23; // length of "outxcoord = inxcoord * "
    size_t end = _funcBody.find_first_of(";\n", start);
    if (end != std::string::npos) {
      std::string multiplier_str = _funcBody.substr(start, end - start);
      try {
        _x_multiplier = std::stod(multiplier_str);
      } catch (...) {
        // Keep default if parsing fails
      }
    }
  }
  
  // Parse y multiplier  
  size_t y_pos = _funcBody.find("outycoord = inycoord * ");
  if (y_pos != std::string::npos) {
    size_t start = y_pos + 23; // length of "outycoord = inycoord * "
    size_t end = _funcBody.find_first_of(";\n", start);
    if (end != std::string::npos) {
      std::string multiplier_str = _funcBody.substr(start, end - start);
      try {
        _y_multiplier = std::stod(multiplier_str);
      } catch (...) {
        // Keep default if parsing fails
      }
    }
  }
  
  // Parse z multiplier
  size_t z_pos = _funcBody.find("outzcoord = inzcoord * ");
  if (z_pos != std::string::npos) {
    size_t start = z_pos + 23; // length of "outzcoord = inzcoord * "
    size_t end = _funcBody.find_first_of(";\n", start);
    if (end != std::string::npos) {
      std::string multiplier_str = _funcBody.substr(start, end - start);
      try {
        _z_multiplier = std::stod(multiplier_str);
      } catch (...) {
        // Keep default if parsing fails
      }
    }
  }
}


/*****************************************************************************/
Geometry_Transform::Geometry_Transform(const std::string & funcBody,
				       std::stringstream & error_string)
/*****************************************************************************/
:  _funcBody(funcBody),
  _x_multiplier(1.0),
  _y_multiplier(1.0),
  _z_multiplier(1.0),
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
  
  // Parse the transformation to extract multipliers for device execution
  parseTransformationMultipliers();
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
  // Create host view that wraps the raw pointer
  // Note: The input coords array uses structure-of-arrays layout (all x, then all y, then all z)
  // We need to create a 2D view with array-of-structures layout for device processing
  HostView2D<double> coords_host("coords_host", num_nodes, dim);
  
  // Copy data from flattened structure-of-arrays to array-of-structures layout
  for (long long i = 0; i < num_nodes; i++) {
    for (long long axis = 0; axis < dim; axis++) {
      coords_host(i, axis) = coords[i + axis * num_nodes];
    }
  }
  
  // Create device view and copy data
  auto coords_device = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), coords_host);
  Kokkos::deep_copy(coords_device, coords_host);
  
  // Execute device computation
  Operate_Device(coords_device, num_nodes, dim);
  
  // Copy results back to host memory
  Kokkos::deep_copy(coords_host, coords_device);
  
  // Copy back from array-of-structures layout to structure-of-arrays layout
  for (long long i = 0; i < num_nodes; i++) {
    for (long long axis = 0; axis < dim; axis++) {
      coords[i + axis * num_nodes] = coords_host(i, axis);
    }
  }
}

/*****************************************************************************/
KOKKOS_INLINE_FUNCTION
void Geometry_Transform::Operate_Device(View2D<double> coords, long long num_nodes, long long dim)
/*****************************************************************************/
{
  // Device version of the geometry transformation
  // Apply the parsed multipliers to each coordinate
  
  // Capture the multipliers for use in the lambda
  double x_mult = _x_multiplier;
  double y_mult = _y_multiplier;
  double z_mult = _z_multiplier;
  
  Kokkos::parallel_for("GeometryTransform", num_nodes, 
    KOKKOS_LAMBDA(const size_t i) {
      // Apply transformation based on parsed multipliers
      if (dim > 0) {
        coords(i, 0) = coords(i, 0) * x_mult; // x coordinate
      }
      if (dim > 1) {
        coords(i, 1) = coords(i, 1) * y_mult; // y coordinate
      }
      if (dim > 2) {
        coords(i, 2) = coords(i, 2) * z_mult; // z coordinate
      }
    });
  
  Kokkos::fence();
}


}//end namespace
