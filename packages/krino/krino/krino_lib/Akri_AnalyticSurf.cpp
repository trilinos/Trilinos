// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AnalyticSurf.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Transformation.hpp>

#include <cmath>

namespace krino{

Cylinder::Cylinder(const std::string & n,  // surface name
                   const double e1[3],  // first endpoint of axis
                   const double e2[3],  // second endpoint of axis
                   const double r,      // radius of cylinder
                   const int sign)
  : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
    dist_sign(sign)
{
  // initialize internal description of cylinder
  radius = r;

  p1 = Vector3d(e1);
  p2 = Vector3d(e2);

  xi = p2 - p1;

  length = xi.length();

  xi.unitize();
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

  xi = p2 - p1;

  length = xi.length();

  xi.unitize();
}

BoundingBox
Cylinder::get_bounding_box()
{
  BoundingBox bbox;
  bbox.accommodate(p1);
  bbox.accommodate(p2);
  bbox.pad(radius);
  return bbox;
}

double
Cylinder::point_signed_distance(const Vector3d &x) const
{
  double D = std::numeric_limits<double>::max();

  // convert x to cylindrical coordinates

  // v is the vector from p1 to x.
  Vector3d v = x - p1;

  // Xix is the projection of v along the axis of the cylinder
  const double Xix = Dot(v,xi);

  // u is the associated vector
  Vector3d u = Xix*xi;

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

Point::Point(const std::string & n,  // surface name
             const Vector3d & coords)
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
      my_coords(coords)
{
}

double
Point::point_signed_distance(const Vector3d &x) const
{
  return (x-my_coords).length_squared();
}

BoundingBox
Point::get_bounding_box()
{
  BoundingBox bbox;
  bbox.accommodate(my_coords);
  return bbox;
}

Sphere::Sphere(const std::string & n,  // surface name
               const Vector3d & center,
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

BoundingBox
Sphere::get_bounding_box()
{
  return BoundingBox(
      Vector3d(myCenter[0]-myRadius, myCenter[1]-myRadius, myCenter[2]-myRadius),
      Vector3d(myCenter[0]+myRadius, myCenter[1]+myRadius, myCenter[2]+myRadius)
      );
}

double
Sphere::point_signed_distance(const Vector3d &x) const
{
  return (myDistSign*((x-myCenter).length() - myRadius));
}

Ellipsoid::Ellipsoid(
      const std::string & name,  // surface name
      const std::vector<double> & center,
      const std::vector<double> & semiAxes,
      const std::vector<double> & rotationVec,
      const int sign)
: SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
  mySign(sign)
{
  ThrowAssert(center.size() == 3 && semiAxes.size() == 3);
  mySemiAxesNorm = 0.0;
  for(unsigned i = 0; i < 3; ++i) {
    myCenter[i] = center[i];
    mySemiAxes[i] = semiAxes[i];
    mySemiAxesNorm += mySemiAxes[i]*mySemiAxes[i];
  }
  mySemiAxesNorm = std::sqrt(mySemiAxesNorm);

  if (!rotationVec.empty())
  {
    ThrowAssert(rotationVec.size() == 3);
    Vector3d pointRotation(-rotationVec[0], -rotationVec[1], -rotationVec[2]); // myRotation is the rotation used for the query point locations, which is the opposite of the rotation of the ellipsoid
    myRotation = std::make_unique<Quaternion>();
    myRotation->set_from_rotation_vector(pointRotation);
  }
}

BoundingBox
Ellipsoid::get_bounding_box()
{
  if (!myRotation)
    return  BoundingBox(myCenter-mySemiAxes, myCenter+mySemiAxes);

  const Vector3d min = myCenter-mySemiAxes;
  const Vector3d max = myCenter+mySemiAxes;
  BoundingBox rotatedBbox;
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(min[0],min[1],min[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(min[0],min[1],max[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(min[0],max[1],min[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(min[0],max[1],max[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(max[0],min[1],min[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(max[0],min[1],max[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(max[0],max[1],min[2])));
  rotatedBbox.accommodate(myRotation->reverse_rotate_3d_vector(Vector3d(max[0],max[1],max[2])));

  return rotatedBbox;
}

double
Ellipsoid::point_signed_distance(const Vector3d &x) const
{
  Vector3d delta = x-myCenter;
  if (myRotation)
    delta = myRotation->rotate_3d_vector(delta);
  // Not an exact distance except in the case of a sphere
  const double dx = delta[0]/mySemiAxes[0];
  const double dy = delta[1]/mySemiAxes[1];
  const double dz = delta[2]/mySemiAxes[2];
  return mySign*(mySemiAxesNorm*std::sqrt(dx*dx+dy*dy+dz*dz)-mySemiAxesNorm);
}

Plane::Plane(const std::string & n,  // surface name
               const double normal[3],
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
Plane::point_signed_distance(const Vector3d &x) const
{
  return myMultiplier*(Dot(myNormal, x) + myOffset);
}

BoundingBox
Plane::get_bounding_box()
{
  //bounding box is entire domain
  return BoundingBox(Vector3d(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()),
      Vector3d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
}

Random::Random(const unsigned long seed)
  : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign(),
    my_amplitude(1.0)
{
  my_srand(seed);
}

double
Random::point_signed_distance(const Vector3d &x) const
{
  // generate random number between -my_amplitude and my_amplitude
  return my_amplitude * (-1.0 + 2.0 * my_rand());
}

BoundingBox
Random::get_bounding_box()
{
  //bounding box is entire domain
  return BoundingBox(Vector3d(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()),
      Vector3d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
}

Analytic_Isosurface::Analytic_Isosurface()
    : SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign()
{
}

BoundingBox
Analytic_Isosurface::get_bounding_box()
{
  return BoundingBox(
      Vector3d(-1.,-1.,-1.),
      Vector3d(1.,1.,1.)
      );
}

double
Analytic_Isosurface::point_signed_distance(const Vector3d &coord) const
{
  const double x = coord[0];
  const double y = coord[1];
  const double z = coord[2];
  return 2.*y*(y*y-3.*x*x)*(1.-z*z) + std::pow(x*x+y*y,2) - (9.*z*z-1.)*(1.-z*z);
}

} // namespace krino
