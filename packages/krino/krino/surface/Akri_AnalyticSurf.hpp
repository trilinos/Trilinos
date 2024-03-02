// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_AnalyticSurf_h
#define Akri_AnalyticSurf_h

#include <Akri_BoundingBox.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_Transformation.hpp>

#include <stk_math/StkVector.hpp>
#include <Akri_String_Function_Expression.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }

namespace krino {

class Cylinder: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Cylinder(const double e1[3],  // first endpoint of axis
           const double e2[3],  // second endpoint of axis
           const double r,      // radius of cylinder
           const int sign);

  virtual ~Cylinder() {}
 
  virtual Surface_Type type() const override { return CYLINDER; }
  virtual size_t storage_size() const override { return sizeof(Cylinder); }
  
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override { bbox.accommodate(myBoundingBox); }
  bool does_intersect(const BoundingBox & bbox) const override { return bbox.intersects(myBoundingBox); }

private:
  void set_axis_and_bounding_box();
  int dist_sign;
  double radius;
  double length;

  // two end points
  stk::math::Vector3d p1, p2;

  // unit vector in direction of axis.
  stk::math::Vector3d xi;
  BoundingBox myBoundingBox;
};

class Sphere: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Sphere(const stk::math::Vector3d & center,
	 const double radius,
         const int sign = 1);

  virtual ~Sphere() {}

  virtual Surface_Type type() const override { return SPHERE; }
  virtual size_t storage_size() const override { return sizeof(Sphere); }
  
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;

private:
  int myDistSign;
  stk::math::Vector3d myCenter;
  double myRadius;
};

class Ellipsoid : public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Ellipsoid(
      const std::vector<double> & center,
      const std::vector<double> & semiAxes,
      const std::vector<double> & rotationVec,
      const int sign);

  virtual ~Ellipsoid() {}

  virtual Surface_Type type() const override { return ELLIPSOID; }
  virtual size_t storage_size() const override { return sizeof(Ellipsoid); }

  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override { bbox.accommodate(myBoundingBox); }
  bool does_intersect(const BoundingBox & bbox) const override { return bbox.intersects(myBoundingBox); }

private:
  void set_bounding_box();
  int mySign;

  //center of ellipsoid
  stk::math::Vector3d myCenter;

  stk::math::Vector3d mySemiAxes;
  double mySemiAxesNorm;
  std::unique_ptr<Quaternion> myRotation;
  BoundingBox myBoundingBox;
};

class Plane: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Plane(const stk::math::Vector3d & normal,
    const double offset) : Plane(normal, offset, 1.) {}
  Plane(const stk::math::Vector3d & normal,
    const double offset,
    const double multiplier);

  virtual ~Plane() {}

  virtual Surface_Type type() const override { return PLANE; }
  virtual size_t storage_size() const override { return sizeof(Plane); }

  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;

private:
  double myMultiplier;

  stk::math::Vector3d myNormal;       //unit normal to plane
  double myOffset;       //eq of plane is 0 = offset + normal \dot x
};

class Point: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Point(const stk::math::Vector3d & coords);

  virtual ~Point() {}

  virtual Surface_Type type() const override { return POINT; }
  virtual size_t storage_size() const override { return sizeof(Point); }
  
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;

private:
  stk::math::Vector3d my_coords;
};

class Random : public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Random(const unsigned long seed);

  virtual ~Random() {}

  virtual Surface_Type type() const override { return RANDOM; }
  virtual size_t storage_size() const override { return sizeof(Random); }
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override
  {
    if (truncation_length > 0.0) my_amplitude = truncation_length;
  }

  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;

private:
  mutable unsigned long iseed;
  double my_amplitude;

  double my_rand() const { iseed = (1664525L*iseed + 1013904223L) & 0xffffffff;
                  return( ((double) iseed) / ((double) 0xffffffff) ); }
  void my_srand(unsigned int seed) const {iseed = seed;}
};

class LevelSet_String_Function: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  LevelSet_String_Function(const std::string & expression);

  virtual ~LevelSet_String_Function() {}

  virtual Surface_Type type() const override { return STRING_FUNCTION; }
  virtual size_t storage_size() const override { return sizeof(LevelSet_String_Function); }

  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  void insert_into(BoundingBox & bbox) const override { bbox.accommodate(myBoundingBox); }
  bool does_intersect(const BoundingBox & bbox) const override { return bbox.intersects(myBoundingBox); }

  void set_bounding_box(const BoundingBox & bbox) { myBoundingBox = bbox; }

private:
  String_Function_Expression myExpression;
  BoundingBox myBoundingBox;
};

} // namespace krino

#endif // Akri_AnalyticSurf_h
