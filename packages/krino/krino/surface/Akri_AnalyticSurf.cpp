// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AnalyticSurf.hpp>
#include <Akri_Transformation.hpp>

#include <cmath>

namespace krino{

Cylinder::Cylinder(const double e1[3],  // first endpoint of axis
                   const double e2[3],  // second endpoint of axis
                   const double r,      // radius of cylinder
                   const int sign)
  : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
    dist_sign(sign)
{
  // initialize internal description of cylinder
  radius = r;

  p1 = stk::math::Vector3d(e1);
  p2 = stk::math::Vector3d(e2);

  set_axis_and_bounding_box();
}

void
Cylinder::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  if (NULL == my_transformation)
  {
    return;
  }

  my_transformation->update(time);
  my_transformation->apply(p1);
  my_transformation->apply(p2);

  set_axis_and_bounding_box();
}

void Cylinder::set_axis_and_bounding_box()
{
  xi = p2 - p1;

  length = xi.unitize();

  myBoundingBox.clear();
  myBoundingBox.accommodate(p1);
  myBoundingBox.accommodate(p2);
  myBoundingBox.pad(radius);
}

double
Cylinder::point_signed_distance(const stk::math::Vector3d &x) const
{
  double D = std::numeric_limits<double>::max();

  // convert x to cylindrical coordinates

  // v is the vector from p1 to x.
  stk::math::Vector3d v = x - p1;

  // Xix is the projection of v along the axis of the cylinder
  const double Xix = Dot(v,xi);

  // u is the associated vector
  stk::math::Vector3d u = Xix*xi;

  // Rx is the distance from the axis of the cylinder to x

  const double Rx = (v - u).length();

  if ( Xix <= 0 )
  {
    // x is behind the cylinder.
    // The closest point on the cylinder lies on the rear circular face.

    if ( Rx <= radius )
    {
      D = -Xix;
    }
    else
    {
      D = std::sqrt(Xix*Xix + (Rx-radius)*(Rx-radius));
    }
  }
  else if( Xix > 0  && Xix < length)
  {
    // x is not in front of or behind the cylinder: it is alongside the cylinder.
    // the signed distance is given by Rx - radius unless this point is closer
    // to one of the circular faces

    if ( Rx >= radius )
    {
      D = Rx - radius;
    }
    else
    {
      D = std::max( std::max(Rx - radius, -Xix), Xix-length );
    }
  }
  else
  {
    // x is in front of the cylinder.
    // The closest point on the cylinder lies on the front circular face.

    if ( Rx <= radius )
    {
      D = Xix-length;
    }
    else
    {
      D = std::sqrt((Xix-length)*(Xix-length) + (Rx-radius)*(Rx-radius));
    }
  }

  return dist_sign*D;
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
Sphere::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
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

bool Plane::does_intersect(const BoundingBox & bbox) const
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
Random::point_signed_distance(const stk::math::Vector3d &x) const
{
  // generate random number between -my_amplitude and my_amplitude
  return my_amplitude * (-1.0 + 2.0 * my_rand());
}

void Random::insert_into(BoundingBox & bbox) const
{
  bbox.accommodate(BoundingBox::ENTIRE_DOMAIN);
}

bool Random::does_intersect(const BoundingBox & bbox) const
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

} // namespace krino
