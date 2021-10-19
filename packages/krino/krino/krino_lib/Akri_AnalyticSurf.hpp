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

#include <stk_mesh/base/Field.hpp>
#include <stk_util/diag/Timer.hpp>

#include <Akri_Vec.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }

namespace krino {

class Cylinder: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Cylinder(const std::string & n,  // surface name
           const double e1[3],  // first endpoint of axis
           const double e2[3],  // second endpoint of axis
           const double r,      // radius of cylinder
           const int sign);

  virtual ~Cylinder() {}
 
  virtual Surface_Type type() const override { return CYLINDER; }
  virtual size_t storage_size() const override { return sizeof(Cylinder); }
  
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;

private:
  int dist_sign;
  double radius;
  double length;

  // two end points
  Vector3d p1, p2;

  // unit vector in direction of axis.
  Vector3d xi;
};

class Sphere: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Sphere(const std::string & n,  // surface name
	 const Vector3d & center,
	 const double radius,
         const int sign = 1);

  virtual ~Sphere() {}

  virtual Surface_Type type() const override { return SPHERE; }
  virtual size_t storage_size() const override { return sizeof(Sphere); }
  
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;

private:
  int myDistSign;
  Vector3d myCenter;
  double myRadius;
};

class Ellipsoid : public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Ellipsoid(
      const std::string & name,  // surface name
      const std::vector<double> & center,
      const std::vector<double> & semiAxes,
      const std::vector<double> & rotationVec,
      const int sign);

  virtual ~Ellipsoid() {}

  virtual Surface_Type type() const override { return ELLIPSOID; }
  virtual size_t storage_size() const override { return sizeof(Ellipsoid); }

  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;

private:
  int mySign;

  //center of ellipsoid
  Vector3d myCenter;

  Vector3d mySemiAxes;
  double mySemiAxesNorm;
  std::unique_ptr<Quaternion> myRotation;
};

class Plane: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Plane(const std::string & n,  // surface name
	 const double normal[3],
	 const double offset,
	 const double multiplier);

  virtual ~Plane() {}

  virtual Surface_Type type() const override { return PLANE; }
  virtual size_t storage_size() const override { return sizeof(Plane); }

  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;

private:
  double myMultiplier;

  Vector3d myNormal;       //unit normal to plane
  double myOffset;       //eq of plane is 0 = offset + normal \dot x
};

class Point: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Point(const std::string & n,  // surface name
        const Vector3d & coords);

  virtual ~Point() {}

  virtual Surface_Type type() const override { return POINT; }
  virtual size_t storage_size() const override { return sizeof(Point); }
  
  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;

private:
  Vector3d my_coords;
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

  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;

private:
  mutable unsigned long iseed;
  double my_amplitude;

  double my_rand() const { iseed = (1664525L*iseed + 1013904223L) & 0xffffffff;
                  return( ((double) iseed) / ((double) 0xffffffff) ); }
  void my_srand(unsigned int seed) const {iseed = seed;}
};

class Analytic_Isosurface: public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign {
public:
  Analytic_Isosurface();

  virtual ~Analytic_Isosurface() {}

  virtual Surface_Type type() const override { return SPHERE; }
  virtual size_t storage_size() const override { return sizeof(Analytic_Isosurface); }

  virtual double point_signed_distance(const Vector3d &x) const override;
  virtual BoundingBox get_bounding_box() override;
};

} // namespace krino

#endif // Akri_AnalyticSurf_h
