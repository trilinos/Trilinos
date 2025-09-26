// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_BoundingRegion_hpp
#define percept_BoundingRegion_hpp

#include <vector>

#include <percept/FieldTypes.hpp>

namespace stk { namespace mesh {
    struct Entity;  
}}

namespace percept {

class BoundingRegion {
public:
  BoundingRegion(CoordinatesFieldType * coords) : 
    coordinatesField_(coords) {}

  virtual ~BoundingRegion() {}

  virtual bool withinGeometricLimits(stk::mesh::Entity elem) = 0;  

protected:
  CoordinatesFieldType * coordinatesField_;
};

class CylinderBoundingRegion : public BoundingRegion
{
public:
  CylinderBoundingRegion(CoordinatesFieldType * coords, double radius, 
                         std::vector<double> &start, 
                         std::vector<double> &end);
  
  bool withinGeometricLimits(stk::mesh::Entity elem) override;
  
private:
  double radius_, length_;
  std::vector<double> start_, normal_;
};
  
class SphereBoundingRegion : public BoundingRegion
{
public:
  SphereBoundingRegion(CoordinatesFieldType * coords, double radius, 
                       std::vector<double> &center); 
  
  bool withinGeometricLimits(stk::mesh::Entity elem) override;
  
private:
  double radius_;
  std::vector<double> center_;
};
  
class BoxBoundingRegion : public BoundingRegion
{
public:
  BoxBoundingRegion(CoordinatesFieldType * coords, 
                    std::vector<double> &start, 
                    std::vector<double> &end);
  
  bool withinGeometricLimits(stk::mesh::Entity elem) override;
  
private:
  std::vector<double> start_, end_;
};
  
}

#endif
