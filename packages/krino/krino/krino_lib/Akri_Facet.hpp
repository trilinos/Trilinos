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

#include <Akri_TypeDefs.hpp>
#include <Akri_Vec.hpp>
#include <Akri_Segment.hpp>
#include <Akri_Triangle.hpp>

namespace stk { class CommBuffer; }

namespace krino {

class Facet;
class Transformation;

typedef std::vector< Facet * > FacetVec;
typedef std::vector< std::unique_ptr<Facet> > FacetOwningVec;

class Facet {
public:
  Facet() {}
  virtual ~Facet() {}

  // for debugging memory usage
  virtual size_t storage_size() const = 0;

  virtual std::ostream & put( std::ostream& os ) const = 0;
  friend std::ostream & operator << ( std::ostream &os , const Facet &f ) { return f.put(os); }

  // methods for off-processor communication
  static std::unique_ptr<Facet> unpack_from_buffer( stk::CommBuffer & b ); // static method that builds facet from data in buffer for off-processor communication
  virtual void pack_into_buffer(stk::CommBuffer & b) const = 0; // pack into buffer for off-processor communication

  virtual Vector3d facet_vertex(const int i) const = 0;
  virtual double facet_area() const = 0;
  virtual Vector3d facet_normal() const = 0;
  virtual bool degenerate() const = 0;

  virtual double point_distance_squared( const Vector3d & x ) const = 0;
  virtual double point_distance_squared( const Vector3d & x, Vector2d & parametric_coords ) const = 0;
  virtual int point_distance_sign( const Vector3d & x ) const = 0;
  virtual Vector3d real_coordinates(const Vector2d & parametric_coords) const = 0;
  virtual Vector3d weights(const Vector2d & parametric_coords) const = 0;

  virtual void apply_transformation(const Transformation & transformation) = 0;

  const BoundingBox & bounding_box() const { ThrowAssert (my_bounding_box.valid()); return my_bounding_box; }
  static const BoundingBox & get_bounding_box(const Facet * facet) { return facet->bounding_box(); }

protected:
  BoundingBox my_bounding_box;
};

class Facet3d : public Facet {
public:
  Facet3d( const Vector3d & x0,
           const Vector3d & x1,
           const Vector3d & x2 );
  Facet3d() = delete;
  virtual ~Facet3d() {}

  virtual size_t storage_size() const { return sizeof(Facet3d); }
  virtual std::ostream & put( std::ostream& os ) const;

  static std::unique_ptr<Facet3d> unpack_from_buffer( stk::CommBuffer & b ); // static method that builds facet from data in buffer for off-processor communication
  virtual void pack_into_buffer(stk::CommBuffer & b) const;

  virtual Vector3d facet_vertex(const int i) const { return my_facet_tri.GetNode(i); }
  virtual double facet_area() const { return my_facet_tri.area(); }
  virtual Vector3d facet_normal() const { return my_facet_tri.normal(); }
  virtual bool degenerate() const { return my_facet_tri.normal_dir().zero_length(); }
  virtual double point_distance_squared( const Vector3d & x ) const
    { return my_facet_tri.DistanceSquared( x ); }
  virtual double point_distance_squared( const Vector3d & x, Vector2d & parametric_coords ) const
    { return my_facet_tri.DistanceSquared( x, parametric_coords ); }
  virtual int point_distance_sign( const Vector3d & x ) const
    { return (Dot(my_facet_tri.normal_dir(),x-facet_vertex(0)) < 0.0) ? -1 : 1; }
  virtual Vector3d real_coordinates(const Vector2d & parametric_coords) const
    { return my_facet_tri.ParametricToRealCoords(parametric_coords); }
  virtual Vector3d weights(const Vector2d & parametric_coords) const
    { return Vector3d(1.0-parametric_coords[0]-parametric_coords[1],parametric_coords[0],parametric_coords[1]); }
  virtual void apply_transformation(const Transformation & transformation);
private:
  Triangle3d my_facet_tri;
};

class Facet2d : public Facet {
public:
  Facet2d( const Vector3d & x0,
               const Vector3d & x1 );
  Facet2d() = delete;
  virtual ~Facet2d() {}

  virtual size_t storage_size() const { return sizeof(Facet2d); }
  virtual std::ostream & put( std::ostream& os ) const;

  static std::unique_ptr<Facet2d> unpack_from_buffer( stk::CommBuffer & b ); // static method that builds facet from data in buffer for off-processor communication
  virtual void pack_into_buffer(stk::CommBuffer & b) const;

  virtual Vector3d facet_vertex(const int i) const { return my_facet_segment.GetNode(i); }
  virtual double facet_area() const { return my_facet_segment.Length(); }
  virtual Vector3d facet_normal() const
    { return (crossZ(facet_vertex(1)-facet_vertex(0))).unit_vector(); }
  virtual bool degenerate() const { return (0.0 == my_facet_segment.Length()); }
  virtual double point_distance_squared( const Vector3d & x ) const
    { return my_facet_segment.DistanceSquared(x); }
  virtual double point_distance_squared( const Vector3d & x, Vector2d & parametric_coords ) const
    { return my_facet_segment.DistanceSquared(x, parametric_coords[0]); }
  virtual int point_distance_sign( const Vector3d & x ) const
    { return (Dot(crossZ(facet_vertex(1)-facet_vertex(0)), x-facet_vertex(0)) < 0.0) ? -1 : 1; }
  virtual Vector3d real_coordinates(const Vector2d & parametric_coords) const
    { return (1.0-parametric_coords[0])*my_facet_segment.GetNode(0) + parametric_coords[0]*my_facet_segment.GetNode(1); }
  virtual Vector3d weights(const Vector2d & parametric_coords) const
    { return Vector3d(1.0-parametric_coords[0],parametric_coords[0],0.0); }
  virtual void apply_transformation(const Transformation & transformation);
private:
  Segment3d my_facet_segment;
};

class FacetDistanceQuery {
public:
  FacetDistanceQuery() : my_facet(nullptr), my_sqr_distance(0.0) {}
  FacetDistanceQuery(const Facet & in_facet, const Vector3d & x) : my_facet(&in_facet), my_query_pt(x)
  {
    my_sqr_distance = my_facet->point_distance_squared(x, my_parametric_coords);
  }

  bool empty() const { return my_facet == nullptr; }
  const Facet & facet() const { return *my_facet; }
  double distance_squared() const { return my_sqr_distance; }
  Vector3d closest_point() const { return my_facet->real_coordinates(my_parametric_coords); }
  double signed_distance() const { return std::sqrt(my_sqr_distance)*my_facet->point_distance_sign(my_query_pt); }
  Vector3d closest_point_weights() const { return my_facet->weights(my_parametric_coords); }

private:
  const Facet* my_facet;
  Vector3d my_query_pt;
  double my_sqr_distance;
  Vector2d my_parametric_coords;
};

} // namespace krino

#endif // Akri_Facet_h
