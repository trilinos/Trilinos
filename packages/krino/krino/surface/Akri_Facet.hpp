// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Facet_h
#define Akri_Facet_h

#include <cmath>
#include <vector>
#include <memory>

#include <Akri_BoundingBox.hpp>
#include <Akri_Segment.hpp>
#include <Akri_Triangle.hpp>
#include <stk_math/StkVector.hpp>

namespace stk { class CommBuffer; }

namespace krino {

class Facet;
class Transformation;

class Facet {
public:
  Facet() {}
  virtual ~Facet() {}

  // for debugging memory usage
  virtual size_t storage_size() const = 0;

  virtual std::ostream & put( std::ostream& os ) const = 0;
  friend std::ostream & operator << ( std::ostream &os , const Facet &f ) { return f.put(os); }

  static stk::math::Vector3d get_centroid(const Facet * facet) { return facet->centroid(); }
  static void insert_into_bounding_box(const Facet * facet, BoundingBox & bbox) { return facet->insert_into(bbox); }

  virtual const stk::math::Vector3d & facet_vertex(const int i) const = 0;
  virtual stk::math::Vector3d & facet_vertex(const int i) = 0;
  virtual double facet_area() const = 0;
  virtual stk::math::Vector3d facet_normal() const = 0;
  virtual bool degenerate() const = 0;

  virtual double point_distance_squared( const stk::math::Vector3d & x ) const = 0;
  virtual void closest_point( const stk::math::Vector3d & queryPt, stk::math::Vector3d & closestPt ) const = 0;
  virtual void closest_point( const stk::math::Vector3d & queryPt, stk::math::Vector3d & closestPt, stk::math::Vector2d & paramAtClosestPt ) const = 0;
  virtual double facet_plane_signed_distance( const stk::math::Vector3d & x ) const = 0;
  virtual int point_distance_sign( const stk::math::Vector3d & x ) const { return (facet_plane_signed_distance(x) < 0.0) ? -1 : 1; }
  virtual stk::math::Vector3d weights(const stk::math::Vector2d & parametric_coords) const = 0;
  virtual stk::math::Vector3d centroid() const = 0;

  virtual void apply_transformation(const Transformation & transformation) = 0;
  virtual void insert_into(BoundingBox & bbox) const = 0;
  virtual bool does_intersect(const BoundingBox & bbox) const = 0;
  virtual double mean_squared_edge_length() const = 0;
};

class Facet3d : public Facet {
public:
  Facet3d( const stk::math::Vector3d & x0,
           const stk::math::Vector3d & x1,
           const stk::math::Vector3d & x2 );
  Facet3d(stk::CommBuffer & b);
  Facet3d() = delete;
  virtual ~Facet3d() {}

  static constexpr int DIM=3;
  typedef CalcTriangle3<double> Calc;
  static bool is_degenerate(const std::array<stk::math::Vector3d,3> & coords) { return Calc::normal_dir(coords).zero_length(); }
  static stk::math::Vector3d get_centroid(const Facet3d * facet) { return facet->centroid(); }
  static void insert_into_bounding_box(const Facet3d * facet, BoundingBox & bbox) { return facet->insert_into(bbox); }

  virtual size_t storage_size() const override { return sizeof(Facet3d); }
  virtual std::ostream & put( std::ostream& os ) const override;

  virtual void pack_into_buffer(stk::CommBuffer & b) const;
  static void emplace_back_from_buffer( std::vector<Facet3d> & facets, stk::CommBuffer & b );

  virtual const stk::math::Vector3d & facet_vertex(const int i) const override { return myCoords[i]; }
  virtual stk::math::Vector3d & facet_vertex(const int i) override { return myCoords[i]; }
  virtual double facet_area() const override { return Calc::area(myCoords); }
  virtual stk::math::Vector3d facet_normal() const override { return Calc::normal(myCoords); }
  virtual bool degenerate() const override { return is_degenerate(myCoords); }
  virtual double point_distance_squared( const stk::math::Vector3d & x ) const override
    { return Calc::distance_squared( myCoords, x ); }
  virtual void closest_point( const stk::math::Vector3d & queryPt, stk::math::Vector3d & closestPt ) const override
    { Calc::closest_point( myCoords[0], myCoords[1], myCoords[2], queryPt, closestPt ); }
  virtual void closest_point( const stk::math::Vector3d & queryPt, stk::math::Vector3d & closestPt, stk::math::Vector2d & paramAtClosestPt ) const override
    { Calc::closest_point_and_parametric_coords( myCoords, queryPt, closestPt, paramAtClosestPt ); }
  virtual double facet_plane_signed_distance( const stk::math::Vector3d & x ) const override
    { return Dot(Calc::normal_dir(myCoords),x-facet_vertex(0)); }
  virtual stk::math::Vector3d weights(const stk::math::Vector2d & parametric_coords) const override
    { return stk::math::Vector3d(1.0-parametric_coords[0]-parametric_coords[1],parametric_coords[0],parametric_coords[1]); }
  virtual stk::math::Vector3d centroid() const override
    { return 1./3.*(myCoords[0]+myCoords[1]+myCoords[2]); }
  virtual void apply_transformation(const Transformation & transformation) override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;
  virtual double mean_squared_edge_length() const override;
private:
  std::array<stk::math::Vector3d,3> myCoords;
};

class Facet2d : public Facet {
public:
  Facet2d( const stk::math::Vector3d & x0,
               const stk::math::Vector3d & x1 );
  Facet2d(stk::CommBuffer & b);
  Facet2d() = delete;
  virtual ~Facet2d() {}

  static constexpr int DIM=2;
  typedef CalcSegment3<double> Calc;
  static stk::math::Vector3d get_centroid(const Facet2d * facet) { return facet->centroid(); }
  static void insert_into_bounding_box(const Facet2d * facet, BoundingBox & bbox) { return facet->insert_into(bbox); }

  virtual size_t storage_size() const override { return sizeof(Facet2d); }
  virtual std::ostream & put( std::ostream& os ) const override;

  virtual void pack_into_buffer(stk::CommBuffer & b) const;
  static void emplace_back_from_buffer( std::vector<Facet2d> & facets, stk::CommBuffer & b );

  virtual const stk::math::Vector3d & facet_vertex(const int i) const override { return myCoords[i]; }
  virtual stk::math::Vector3d & facet_vertex(const int i) override { return myCoords[i]; }
  virtual double facet_area() const override { return Calc::length(myCoords); }
  virtual stk::math::Vector3d facet_normal() const override
    { return (crossZ(facet_vertex(1)-facet_vertex(0))).unit_vector(); }
  virtual bool degenerate() const override { return (0.0 == Calc::length_squared(myCoords)); }
  virtual double point_distance_squared( const stk::math::Vector3d & x ) const override
    { return Calc::distance_squared(myCoords, x); }
  virtual void closest_point( const stk::math::Vector3d & queryPt, stk::math::Vector3d & closestPt ) const override
    { Calc::closest_point( myCoords, queryPt, closestPt ); }
  virtual void closest_point( const stk::math::Vector3d & queryPt, stk::math::Vector3d & closestPt, stk::math::Vector2d & paramAtClosestPt ) const override
    { Calc::closest_point_and_parametric_coord(myCoords, queryPt, closestPt, paramAtClosestPt[0]); }
  virtual double facet_plane_signed_distance( const stk::math::Vector3d & x ) const override
    { return Dot(crossZ(facet_vertex(1)-facet_vertex(0)), x-facet_vertex(0)); }
  virtual stk::math::Vector3d weights(const stk::math::Vector2d & parametric_coords) const override
    { return stk::math::Vector3d(1.0-parametric_coords[0],parametric_coords[0],0.0); }
  virtual stk::math::Vector3d centroid() const override
    { return 0.5*(myCoords[0]+myCoords[1]); }
  virtual void apply_transformation(const Transformation & transformation) override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;
  virtual double mean_squared_edge_length() const override { return Calc::length_squared(myCoords); }
private:
  std::array<stk::math::Vector3d,2> myCoords;
};

class FacetWithVelocity2d : public Facet2d {
public:
  FacetWithVelocity2d( const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const stk::math::Vector3d vel0,
    const stk::math::Vector3d vel1)
  : Facet2d(x0, x1), myVelocity{vel0, vel1} {}
  FacetWithVelocity2d( const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1)
  : FacetWithVelocity2d(x0, x1, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO) {}
  FacetWithVelocity2d(stk::CommBuffer & b);
  virtual ~FacetWithVelocity2d() {}
  void set_velocity(const std::array<stk::math::Vector3d,2> & velocity) { myVelocity = velocity; }
  const std::array<stk::math::Vector3d,2> & get_velocity() const { return myVelocity; }
  stk::math::Vector3d velocity_at_closest_point( const stk::math::Vector3d & queryPt ) const;

  virtual void pack_into_buffer(stk::CommBuffer & b) const;
  static void emplace_back_from_buffer( std::vector<FacetWithVelocity2d> & facets, stk::CommBuffer & b );
private:
  std::array<stk::math::Vector3d,2> myVelocity;
};

class FacetWithVelocity3d : public Facet3d {
public:
  FacetWithVelocity3d( const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const stk::math::Vector3d & x2,
    const stk::math::Vector3d vel0,
    const stk::math::Vector3d vel1,
    const stk::math::Vector3d vel2)
  : Facet3d(x0, x1, x2), myVelocity{vel0, vel1, vel2} {}
  FacetWithVelocity3d( const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const stk::math::Vector3d & x2)
  : FacetWithVelocity3d(x0, x1, x2, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO) {}
  FacetWithVelocity3d(stk::CommBuffer & b);
  virtual ~FacetWithVelocity3d() {}
  void set_velocity(const std::array<stk::math::Vector3d,3> & velocity) { myVelocity = velocity; }
  const std::array<stk::math::Vector3d,3> & get_velocity() const { return myVelocity; }
  stk::math::Vector3d velocity_at_closest_point( const stk::math::Vector3d & queryPt ) const;

  virtual void pack_into_buffer(stk::CommBuffer & b) const;
  static void emplace_back_from_buffer( std::vector<FacetWithVelocity3d> & facets, stk::CommBuffer & b );
private:
  std::array<stk::math::Vector3d,3> myVelocity;
};

template <class FACET>
class FacetDistanceQuery {
public:
  FacetDistanceQuery() : myFacet(nullptr), mySqrDistance(0.0) {}
  FacetDistanceQuery(const FACET & facet, const stk::math::Vector3d & queryPt) : myFacet(&facet)
  {
    myFacet->closest_point(queryPt, myClosestPt, myParamCoords);
    mySqrDistance = (queryPt-myClosestPt).length_squared();
  }

  bool empty() const { return myFacet == nullptr; }
  const FACET & facet() const { return *myFacet; }
  double distance_squared() const { return mySqrDistance; }
  const stk::math::Vector3d & closest_point() const { return myClosestPt; }
  double signed_distance(const stk::math::Vector3d & queryPt) const { return std::sqrt(mySqrDistance)*myFacet->point_distance_sign(queryPt); }
  stk::math::Vector3d closest_point_weights() const { return myFacet->weights(myParamCoords); }

private:
  const FACET* myFacet;
  stk::math::Vector3d myClosestPt;
  double mySqrDistance;
  stk::math::Vector2d myParamCoords;
};

} // namespace krino

#endif // Akri_Facet_h
