// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Geometry_h)
#define MiniTensor_Geometry_h

#include <vector>
#include "MiniTensor_Tensor.h"

namespace minitensor {

///
/// Useful to distinguish among different finite elements.
///
namespace ELEMENT{

enum Type {
  UNKNOWN,
  SEGMENTAL,
  TRIANGULAR,
  QUADRILATERAL,
  TETRAHEDRAL,
  HEXAHEDRAL};

} //namespace ELEMENT

///
/// Length of a segment
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
length(Vector<T, N> const & p0, Vector<T, N> const & p1);

///
/// Area of a triangle
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
area(Vector<T, N> const & p0, Vector<T, N> const & p1,
     Vector<T, N> const & p2);

///
/// Area of a quadrilateral.
/// Taken from:
/// Calculation of the volume of a general hexahedron for flow predictions
/// Davies, D. E.; Salmond, D. J.
/// AIAA Journal (ISSN 0001-1452), vol. 23, June 1985, p. 954-956.
/// Their vertex naming convention: ABCD
/// Our convention:                 0123
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
area(Vector<T, N> const & p0, Vector<T, N> const & p1,
     Vector<T, N> const & p2, Vector<T, N> const & p3);

///
/// Volume of tetrahedron
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
volume(Vector<T, N> const & p0, Vector<T, N> const & p1,
       Vector<T, N> const & p2, Vector<T, N> const & p3);

///
/// Volume of hexahedron
/// Taken from:
/// Calculation of the volume of a general hexahedron for flow predictions
/// Davies, D. E.; Salmond, D. J.
/// AIAA Journal (ISSN 0001-1452), vol. 23, June 1985, p. 954-956.
/// Their vertex naming convention: ABCDEFGH
/// Our convention:                 45670123
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
volume(Vector<T, N> const & p0, Vector<T, N> const & p1,
       Vector<T, N> const & p2, Vector<T, N> const & p3,
       Vector<T, N> const & p4, Vector<T, N> const & p5,
       Vector<T, N> const & p6, Vector<T, N> const & p7);

///
/// Centroids of segment, triangle, tetrahedron, quadrilateral
/// and hexahedron.
/// For these we can just take the average of the vertices.
/// WARNING: This is not the center of mass.
///
template<typename T, Index N>
Vector<T, N>
centroid(std::vector<Vector<T, N>> const & points);

///
/// The surface normal of a face
/// Input: 3 independent nodes on the face
/// Output: normal vector
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
normal(Vector<T, N> const & p0,
       Vector<T, N> const & p1,
       Vector<T, N> const & p2);

///
/// Given 3 points p0, p1, p2 that define a plane
/// determine if point p is in the same side of the normal
/// to the plane as defined by the right hand rule.
/// If a tolrance is given, use that as criterion for minimal distance.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_normal_side(
    Vector<T, N> const & p,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    T const tolerance = 0);

///
/// Given two iterators to a container of points,
/// find the associated bounding box.
/// \param start end: define sequence of points
/// \return vectors that define the bounding box
///
template<typename T, typename I, Index N>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, N>, Vector<T, N>>
bounding_box(I start, I end);

template<typename T, typename I>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, DYNAMIC>, Vector<T, DYNAMIC>>
bounding_box(I start, I end);

///
/// Determine if a given point is inside a bounding box.
/// \param p the point
/// \param min max points defining the box
/// \return whether the point is inside
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_box(
    Vector<T, N> const & p,
    Vector<T, N> const & min,
    Vector<T, N> const & max);

///
/// Generate random point inside bounding box
/// \param min max the bounding box
/// \return p point inside box
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
random_in_box(
    Vector<T, N> const & min,
    Vector<T, N> const & max);

///
/// Given 4 points p0, p1, p2, p3 that define a tetrahedron
/// determine if point p is inside it.
/// If a tolrance is given, use that as criterion for minimal distance.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_tetrahedron(
    Vector<T, N> const & p,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3,
    T const tolerance = 0);

///
/// Given 8 points that define a hexahedron
/// determine if point p is inside it.
/// Assumption: faces are planar
/// If a tolrance is given, use that as criterion for minimal distance.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_hexahedron(
    Vector<T, N> const & p,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3,
    Vector<T, N> const & p4,
    Vector<T, N> const & p5,
    Vector<T, N> const & p6,
    Vector<T, N> const & p7,
    T const tolerance = 0);

///
/// Closest point
/// \param p the point
/// \param n vector of points to test
/// \return index to closest point
///
template<typename T, Index N>
typename std::vector<Vector<T, N>>::size_type
closest_point(Vector<T, N> const & p, std::vector<Vector<T, N>> const & n);

/// Median of a sequence defined by random
/// access iterators. Undefined for empty set.
/// \param begin end Iterators that define the sequence
/// \return median of sequence
///
template<typename T, typename Iterator>
KOKKOS_INLINE_FUNCTION
T
median(Iterator begin, Iterator end);

///
/// Given quadrilateral nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_quadrilateral(
    Vector<T, dimension_const<N, 2>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3);

///
/// Given triangle nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_triangle(
    Vector<T, dimension_const<N, 3>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2);

///
/// Given hexahedron nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_hexahedron(
    Vector<T, dimension_const<N, 3>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3,
    Vector<T, N> const & p4,
    Vector<T, N> const & p5,
    Vector<T, N> const & p6,
    Vector<T, N> const & p7);

///
/// Given tetrahedron nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_tetrahedron(
    Vector<T, dimension_const<N, 4>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3);

///
/// Given element type and nodes and a position
/// in parametric coordinates, interpolate.
/// \param element_type element type
/// \param xi position in parametric coordinates
/// \param v ... corner nodes
/// \return interpolated position
///
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_element(
    ELEMENT::Type element_type,
    Vector<T, M> & xi,
    std::vector<Vector<T, N>> const & v);

///
/// Given a vector of points, determine
/// distances between all of them.
/// \param points vector of points
/// \return distance matrix
///
template<typename T, Index N>
std::vector< std::vector<T>>
distance_matrix(std::vector<Vector<T, N>> const & points);

///
/// Given a distance matrix, determine the minimum
/// distance between two distinct points.
/// \param distances distance matrix
/// \return minimum distance
///
template<typename T>
std::vector<T>
minimum_distances(std::vector< std::vector<T>> const & distances);

///
/// Given space dimension and number of (vertex) nodes,
/// determine the type of a finite element.
///
KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type(Index const dimension, Index const number_nodes);

///
/// Spherical parametrization functor
///
template<typename T, Index N>
class SphericalParametrization
{
public:

  KOKKOS_INLINE_FUNCTION
  SphericalParametrization(Tensor4<T, N> const & A);

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 2>::value> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N> const &
  tangent_;

  T
  minimum_;

  Vector<T, 2>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 2>
  arg_maximum_;
};

///
/// Stereographic parametrization functor
///
template<typename T, Index N>
class StereographicParametrization
{
public:

  KOKKOS_INLINE_FUNCTION
  StereographicParametrization(Tensor4<T, N> const & A);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 2>::value> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value> const & parameters);

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N> const &
  tangent_;

  T
  minimum_;

  Vector<T, 2>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 2>
  arg_maximum_;
};

///
/// Projective parametrization functor
///
template<typename T, Index N>
class ProjectiveParametrization
{
public:

  ///
  /// Constructor that takes material tangent
  ///
  KOKKOS_INLINE_FUNCTION
  ProjectiveParametrization(Tensor4<T, N> const & A);

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 3>::value> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 3>::value> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N> const &
  tangent_;

  T
  minimum_;

  Vector<T, 3>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 3>
  arg_maximum_;
};

///
/// Tangent parametrization functor
///
template<typename T, Index N>
class TangentParametrization
{
public:

  ///
  /// Constructor that takes material tangent
  ///
  KOKKOS_INLINE_FUNCTION
  TangentParametrization(Tensor4<T, N> const & A);

  ///
  ///
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 2>::value> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N> const &
  tangent_;

  T
  minimum_;

  Vector<T, 2>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 2>
  arg_maximum_;
};

///
/// Cartesian parametrization functor
///
template<typename T, Index N>
class CartesianParametrization
{
public:

  KOKKOS_INLINE_FUNCTION
  CartesianParametrization(Tensor4<T, N> const & A);

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 3>::value> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 3>::value> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N> const &
  tangent_;

  T
  minimum_;

  Vector<T, 3>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 3>
  arg_maximum_;
};

///
/// Parametric grid class
///
template<typename T, Index N>
class ParametricGrid
{

public:

  ///
  /// Default constructor
  ///
  KOKKOS_INLINE_FUNCTION
  ParametricGrid() {}

  ///
  /// Constructor that defines grid limits
  /// \param lower lower limit
  /// \param upper upper limit
  /// \param points_per_dimension number of points in each dimension
  ///
  KOKKOS_INLINE_FUNCTION
  ParametricGrid(
      Vector<T, N> const & lower,
      Vector<T, N> const & upper,
      Vector<Index, N> const & points_per_dimension);

  ///
  ///
  template<typename Visitor>
  KOKKOS_INLINE_FUNCTION
  void
  traverse(Visitor & visitor) const;

private:

  Vector<T, N>
  lower_;

  Vector<T, N>
  upper_;

  Vector<Index, N>
  points_per_dimension_;

};

} // namespace minitensor

#include "MiniTensor_Geometry.i.h"
#include "MiniTensor_Geometry.t.h"

#endif // MiniTensor_Geometry_h
