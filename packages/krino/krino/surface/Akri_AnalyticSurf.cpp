// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AnalyticSurf.hpp>
#include <Akri_Transformation.hpp>
#include <Akri_IntersectionUtils.hpp>
#include <Akri_ML_Models.hpp>

#include <cmath>
#include <stdexcept>

namespace krino{

Cylinder::Cylinder(const stk::math::Vector3d endPt1,  // first endpoint of axis
                   const stk::math::Vector3d endPt2,  // second endpoint of axis
                   const double radius,
                   const int sign)
  : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
    myEndPt1(endPt1),
    myEndPt2(endPt2),
    myRadius(radius),
    myDistSign(sign)
{
  set_axis_and_bounding_box();
}

void
Cylinder::prepare_to_compute(const double time, const BoundingBox & /*point_bbox*/, const double /*truncation_length*/)
{
  if (NULL == my_transformation)
  {
    return;
  }

  my_transformation->update(time);
  my_transformation->apply(myEndPt1);
  my_transformation->apply(myEndPt2);

  set_axis_and_bounding_box();
}

void Cylinder::set_axis_and_bounding_box()
{
  myAxis = myEndPt2 - myEndPt1;

  myLength = myAxis.unitize();

  myBoundingBox.clear();
  myBoundingBox.accommodate(myEndPt1);
  myBoundingBox.accommodate(myEndPt2);
  myBoundingBox.pad(myRadius);
}

double
Cylinder::point_signed_distance(const stk::math::Vector3d &x) const
{
  double D = 0;

  const stk::math::Vector3d v = x - myEndPt1;
  const double distAlongAxis = Dot(v, myAxis);
  const double distFromAxis = (v - distAlongAxis*myAxis).length();

  if( distAlongAxis >= 0 && distAlongAxis <= myLength)
  {
    // the projection of x onto the axis of the cylinder is between the endpoints of the cylinder.
    // the signed distance is given by Rx - myRadius unless this point is closer
    // to one of the circular faces

    if ( distFromAxis >= myRadius )
      D = distFromAxis - myRadius; // outside
    else
      D = -std::min(myRadius-distFromAxis, std::min(distAlongAxis, myLength-distAlongAxis)); // inside
  }
  else
  {
    // the projection of x onto the axis of the cylinder is outside of the endpoints of the cylinder.
    const double distAlongAxisFromFace =
        (distAlongAxis <= 0) ?
        (-distAlongAxis) :        // x is behind the cylinder, so the closest point on the cylinder lies on the rear circular face
        (distAlongAxis-myLength); // x is in front of the cylinder, so the closest point on the cylinder lies on the front circular face.

    if ( distFromAxis <= myRadius )
      D = distAlongAxisFromFace;
    else
      D = std::sqrt(distAlongAxisFromFace*distAlongAxisFromFace + (distFromAxis-myRadius)*(distFromAxis-myRadius));
  }

  return myDistSign*D;
}

stk::math::Vector3d
Cylinder::closest_point_normal(const stk::math::Vector3d &x) const
{
  // If closest point is on a surface (end cap or cylindrical surface),
  // then closest_point_normal is the surface normal at the closest point
  // Otherwise (closest point is on the sharp curves of the can),
  // the closest_point_normal is (x - closest_point).unit_vector().
  const stk::math::Vector3d v = x - myEndPt1;
  const double distAlongAxis = Dot(v, myAxis);
  const double distFromAxis = (v - distAlongAxis*myAxis).length();

  const stk::math::Vector3d closestPtOnAxis = myEndPt1 + distAlongAxis*myAxis;

  if (distFromAxis < myRadius) // within infinite cylinder
  {
    if (distAlongAxis <= 0.5*myLength && distAlongAxis <= myRadius-distFromAxis) // closest point on bottom of can
      return -myAxis;
    if (distAlongAxis >= 0.5*myLength && distAlongAxis >= myLength-(myRadius-distFromAxis)) // closest point on bottom of can
      return myAxis;
    return (x-closestPtOnAxis).unit_vector(); // closest point on cylindrical surface
  }

  if( distAlongAxis >= 0 && distAlongAxis <= myLength)
    return (x-closestPtOnAxis).unit_vector(); // closest point on cylindrical surface

  const stk::math::Vector3d closestPointOnInfiniteCylinder = closestPtOnAxis + distFromAxis * (x-closestPtOnAxis).unit_vector();

  if (distAlongAxis < 0)
    return (x - (closestPointOnInfiniteCylinder - distAlongAxis*myAxis)).unit_vector(); // closest point on bottom rim
  return (x - (closestPointOnInfiniteCylinder - (distAlongAxis-myLength)*myAxis)).unit_vector();  // closest point on top rim
}

Point::Point(const stk::math::Vector3d & coords)
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
      my_coords(coords)
{
}

double
Point::point_signed_distance(const stk::math::Vector3d &x) const
{
  return (x-my_coords).length_squared();
}

void Point::insert_into(BoundingBox & bbox) const
{
  bbox.accommodate(my_coords);
}

bool Point::does_intersect(const BoundingBox & bbox) const
{
  return bbox.contains(my_coords);
}

Sphere::Sphere(const stk::math::Vector3d & center,
               const double radius,
               const int sign)
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
      myDistSign(sign),
      myCenter(center),
      myRadius(radius)
{
}

void
Sphere::prepare_to_compute(const double time, const BoundingBox & /*point_bbox*/, const double /*truncation_length*/)
{
  if (nullptr != my_transformation)
  {
    my_transformation->update(time);
    my_transformation->apply(myCenter);
  }
}

void Sphere::insert_into(BoundingBox & bbox) const
{
  const stk::math::Vector3d min(myCenter[0]-myRadius, myCenter[1]-myRadius, myCenter[2]-myRadius);
  const stk::math::Vector3d max(myCenter[0]+myRadius, myCenter[1]+myRadius, myCenter[2]+myRadius);

  bbox.accommodate(min);
  bbox.accommodate(max);
}

bool Sphere::does_intersect(const BoundingBox & bbox) const
{
  const stk::math::Vector3d min(myCenter[0]-myRadius, myCenter[1]-myRadius, myCenter[2]-myRadius);
  const stk::math::Vector3d max(myCenter[0]+myRadius, myCenter[1]+myRadius, myCenter[2]+myRadius);

  BoundingBox sphereBbox(min, max);
  return bbox.intersects(sphereBbox);
}

stk::math::Vector3d Sphere::closest_point_normal(const stk::math::Vector3d &x) const
{
  stk::math::Vector3d normal = x-myCenter;
  normal.unitize();
  return normal;
}

double
Sphere::point_signed_distance(const stk::math::Vector3d &x) const
{
  return (myDistSign*((x-myCenter).length() - myRadius));
}

Ellipsoid::Ellipsoid(
      const std::vector<double> & center,
      const std::vector<double> & semiAxes,
      const std::vector<double> & rotationVec,
      const int sign)
: SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
  mySign(sign)
{
  STK_ThrowAssert(center.size() == 3 && semiAxes.size() == 3);
  mySemiAxesNorm = 0.0;
  for(unsigned i = 0; i < 3; ++i) {
    myCenter[i] = center[i];
    mySemiAxes[i] = semiAxes[i];
    mySemiAxesNorm += mySemiAxes[i]*mySemiAxes[i];
  }
  mySemiAxesNorm = std::sqrt(mySemiAxesNorm);

  if (!rotationVec.empty())
  {
    STK_ThrowAssert(rotationVec.size() == 3);
    stk::math::Vector3d pointRotation(-rotationVec[0], -rotationVec[1], -rotationVec[2]); // myRotation is the rotation used for the query point locations, which is the opposite of the rotation of the ellipsoid
    myRotation = std::make_unique<Quaternion>();
    myRotation->set_from_rotation_vector(pointRotation);
  }

  set_bounding_box();
}

void
Ellipsoid::set_bounding_box()
{
  const stk::math::Vector3d min = myCenter-mySemiAxes;
  const stk::math::Vector3d max = myCenter+mySemiAxes;

  myBoundingBox.clear();
  if (!myRotation)
  {
    myBoundingBox.accommodate(min);
    myBoundingBox.accommodate(max);
  }
  else
  {
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(min[0],min[1],min[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(min[0],min[1],max[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(min[0],max[1],min[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(min[0],max[1],max[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(max[0],min[1],min[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(max[0],min[1],max[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(max[0],max[1],min[2])));
    myBoundingBox.accommodate(myRotation->reverse_rotate_3d_vector(stk::math::Vector3d(max[0],max[1],max[2])));
  }
}

double
Ellipsoid::point_signed_distance(const stk::math::Vector3d &x) const
{
  stk::math::Vector3d delta = x-myCenter;
  if (myRotation)
    delta = myRotation->rotate_3d_vector(delta);
  // Not an exact distance except in the case of a sphere
  const double dx = delta[0]/mySemiAxes[0];
  const double dy = delta[1]/mySemiAxes[1];
  const double dz = delta[2]/mySemiAxes[2];
  return mySign*(mySemiAxesNorm*std::sqrt(dx*dx+dy*dy+dz*dz)-mySemiAxesNorm);
}

Cuboid::Cuboid(const stk::math::Vector3d & center,
    const stk::math::Vector3d & dimensions,
    const stk::math::Vector3d & rotationVec,
    const int sign)
: SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
  mySign(sign),
  myCenter(center),
  myHalfDimensions(0.5*dimensions)
{
  if (!rotationVec.zero_length())
  {
    // myRotation is the rotation used for the cuboid, which is the opposite of the rotation of the query points
    myRotation = std::make_unique<Quaternion>();
    myRotation->set_from_rotation_vector(rotationVec);
  }

  set_bounding_box();
}

void
Cuboid::set_bounding_box()
{
  myBoundingBox.clear();
  if (!myRotation)
  {
    myBoundingBox.accommodate(myCenter-myHalfDimensions);
    myBoundingBox.accommodate(myCenter+myHalfDimensions);
  }
  else
  {
    for (unsigned i=0; i<8; ++i)
      myBoundingBox.accommodate(vertex_location(i));
  }
}

stk::math::Vector3d Cuboid::vertex_location(const unsigned i) const
{
  const std::array<std::array<int,3>,8> vertexDeltas
  {{
    {{-1,-1,-1}}, {{+1,-1,-1}}, {{+1,+1,-1}}, {{-1,+1,-1}},
    {{-1,-1,+1}}, {{+1,-1,+1}}, {{+1,+1,+1}}, {{-1,+1,+1}}
  }};
  const std::array<int,3> & delta = vertexDeltas[i];
  const stk::math::Vector3d locationInCuboidCoords(myCenter[0]+delta[0]*myHalfDimensions[0], myCenter[1]+delta[1]*myHalfDimensions[1], myCenter[2]+delta[2]*myHalfDimensions[2]);
  if (!myRotation)
    return locationInCuboidCoords;
  return myRotation->rotate_3d_vector(locationInCuboidCoords);
}

double
Cuboid::point_signed_distance(const stk::math::Vector3d &x) const
{
  stk::math::Vector3d ptInCuboidCoords = x-myCenter;
  if (myRotation)
    ptInCuboidCoords = myRotation->reverse_rotate_3d_vector(ptInCuboidCoords);
  // This calculation is far from obvious.  It is derived in an YouTube video.
  const stk::math::Vector3d delta(std::abs(ptInCuboidCoords[0])-myHalfDimensions[0], std::abs(ptInCuboidCoords[1])-myHalfDimensions[1], std::abs(ptInCuboidCoords[2])-myHalfDimensions[2]);
  const stk::math::Vector3d clip(std::max(delta[0],0.), std::max(delta[1],0.), std::max(delta[2],0.));
  return mySign*(clip.length() + std::min(0., std::max(delta[0], std::max(delta[1], delta[2]))));
}

void Cuboid::fill_triangle_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,3> & faceNodes, std::vector<stk::math::Vector3d> & intParamCoords) const
{
  intParamCoords.clear();

  if (does_bounding_box_intersect_triangle_3d(myBoundingBox, faceNodes))
  {
    const auto [triArea, triNormal] = CalcTriangle3<double>::area_and_normal(faceNodes);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(0), vertex_location(1), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(1), vertex_location(2), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(2), vertex_location(3), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(3), vertex_location(0), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(4), vertex_location(5), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(5), vertex_location(6), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(6), vertex_location(7), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(7), vertex_location(4), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(0), vertex_location(4), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(1), vertex_location(5), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(2), vertex_location(6), intParamCoords);
    append_triangle_edge_intersection_parametric_coordinates(faceNodes, triArea, triNormal, vertex_location(3), vertex_location(7), intParamCoords);
  }
}

void Cuboid::fill_tetrahedon_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,4> & tetNodes, std::vector<stk::math::Vector3d> & intParamCoords) const
{
  intParamCoords.clear();
  if (does_bounding_box_intersect_tetrahedron(myBoundingBox, tetNodes))
  {
    for (unsigned i=0; i<8; ++i)
      append_tetrahedron_intersection_parametric_coordinates(tetNodes, vertex_location(i), intParamCoords);
  }
}


Plane::Plane(const stk::math::Vector3d & normal,
    const double offset,
    const double multiplier)
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
      myMultiplier(multiplier),
      myNormal(normal),
      myOffset(offset)
{
  myNormal.unitize();
}

double
Plane::point_signed_distance(const stk::math::Vector3d &x) const
{
  return myMultiplier*(Dot(myNormal, x) + myOffset);
}

void Plane::insert_into(BoundingBox & bbox) const
{
  bbox.accommodate(BoundingBox::ENTIRE_DOMAIN);
}

bool Plane::does_intersect(const BoundingBox & /*bbox*/) const
{
  return true; // This could actually test the sidedness of the min and max points of the bbox if desired
}

Random::Random(const unsigned long seed)
  : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
    my_amplitude(1.0)
{
  my_srand(seed);
}

double
Random::point_signed_distance(const stk::math::Vector3d &/*x*/) const
{
  // generate random number between -my_amplitude and my_amplitude
  return my_amplitude * (-1.0 + 2.0 * my_rand());
}

void Random::insert_into(BoundingBox & bbox) const
{
  bbox.accommodate(BoundingBox::ENTIRE_DOMAIN);
}

bool Random::does_intersect(const BoundingBox & /*bbox*/) const
{
  return true;
}

LevelSet_String_Function::LevelSet_String_Function(const std::string & expression)
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
      myExpression(expression),
      myBoundingBox(BoundingBox::ENTIRE_DOMAIN)
{
}

double
LevelSet_String_Function::point_signed_distance(const stk::math::Vector3d &coord) const
{
  return myExpression.evaluate(coord);
}

Implicit_Neural_Representation::Implicit_Neural_Representation(const std::string & modelBaseName)
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign()
{
  myModel.create(modelBaseName);
  if (myModelWithGradient.is_enabled())
    myModelWithGradient.create(modelBaseName);
}

double
Implicit_Neural_Representation::point_signed_distance(const stk::math::Vector3d &coord) const
{
  double dist = 0.;
  myModel.evaluate(3, coord.data(), 1, &dist);
  return dist;
}

stk::math::Vector3d
Implicit_Neural_Representation::closest_point_normal(const stk::math::Vector3d &coord) const
{
  double dist = 0.;
  stk::math::Vector3d grad;
  myModelWithGradient.evaluate(3, coord.data(), 1, &dist, grad.data());
  return grad.unit_vector();
}

} // namespace krino
