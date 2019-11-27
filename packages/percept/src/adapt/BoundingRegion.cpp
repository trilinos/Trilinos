// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include<cmath>

#include <adapt/BoundingRegion.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/PerceptUtils.hpp>

namespace percept {

CylinderBoundingRegion::CylinderBoundingRegion(
  CoordinatesFieldType * coords, double radius, 
  std::vector<double> &start, 
  std::vector<double> &end)
  : 
  BoundingRegion(coords),
  radius_(radius),
  length_(0),
  start_(start),
  normal_(3)
{ 
  if (radius<=0.) 
    throw std::runtime_error("CylinderBoundingRegion: cannot use negative radius");

  if (start.size() !=3 || end.size() != 3) 
    throw std::runtime_error("CylinderBoundingRegion: start and end must have length 3");

  for (int i=0; i<3; i++) {
    normal_[i] = end[i]-start_[i];
    length_ += pow(normal_[i],2);
  }
  length_ = std::sqrt(length_);
  for (int i=0; i<3; i++) {
    normal_[i] /= length_;
  }
}

bool 
CylinderBoundingRegion::withinGeometricLimits(stk::mesh::Entity elem)
{
  double centroid[3];
  computeCentroid(elem, centroid, *coordinatesField_);
  
  double proj_to_normal_ = 0;
  for (int i=0; i<3; i++) {
    proj_to_normal_ += (centroid[i]-start_[i])*normal_[i];
  }
  
  double radial_distance = 0;
  for (int i=0; i<3; i++) {
    radial_distance += std::pow(centroid[i] - start_[i] - proj_to_normal_*normal_[i], 2);
  }
  radial_distance = std::sqrt(radial_distance);
  
  // point outside axis of cylinder
  return (proj_to_normal_ <= length_ &&
          proj_to_normal_ >= 0 &&
          radial_distance <= radius_);
}
  
SphereBoundingRegion::SphereBoundingRegion(
  CoordinatesFieldType * coords, double radius, 
  std::vector<double> &center)
  : 
  BoundingRegion(coords),
  radius_(radius),
  center_(center)
{
  if (radius<=0.) 
    throw std::runtime_error("SphereBoundingRegion: cannot use negative radius");
}

bool 
SphereBoundingRegion::withinGeometricLimits(stk::mesh::Entity elem)
{
  double centroid[3];
  computeCentroid(elem, centroid, *coordinatesField_);
  
  double distance = 0;
  for (int i=0; i<3; i++) {
    distance += pow(centroid[i]-center_[i],2);
  }
  distance = std::sqrt(distance);

  return (distance <= radius_);
}

BoxBoundingRegion::BoxBoundingRegion(
  CoordinatesFieldType * coords, 
  std::vector<double> &start,
  std::vector<double> &end)
  : 
  BoundingRegion(coords),
  start_(start),
  end_(end)
{
  if (start.size() !=3 || end.size() != 3) 
    throw std::runtime_error("BoxBoundingRegion: start and end must have length 3");

  if (start[0]>=end[0] || start[1]>=end[1] || start[2]>=end[2]) 
    throw std::runtime_error("BoxBoundingRegion: must have start < end");
}

bool 
BoxBoundingRegion::withinGeometricLimits(stk::mesh::Entity elem)
{
  double centroid[3];
  computeCentroid(elem, centroid, *coordinatesField_);
  
  return (centroid[0]>=start_[0] &&
          centroid[1]>=start_[1] &&
          centroid[2]>=start_[2] &&
          centroid[0]<=end_[0] &&
          centroid[1]<=end_[1] &&
          centroid[2]<=end_[2] 
          );
}

} // namespace percept
