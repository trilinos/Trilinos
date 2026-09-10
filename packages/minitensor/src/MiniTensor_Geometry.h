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

/// \addtogroup minitensor_geometry
/// @{

///
/// Useful to distinguish among different finite elements.
///
namespace ELEMENT{

///
/// Finite element types.
///
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
/// \param start iterator that begins the sequence of points
/// \param end iterator that ends the sequence of points
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
/// \param min lower corner of the box
/// \param max upper corner of the box
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
/// \param min lower corner of the bounding box
/// \param max upper corner of the bounding box
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
/// \param begin iterator that begins the sequence
/// \param end iterator that ends the sequence
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
/// \param p0 corner node
/// \param p1 corner node
/// \param p2 corner node
/// \param p3 corner node
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
/// \param p0 corner node
/// \param p1 corner node
/// \param p2 corner node
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
/// \param p0 corner node
/// \param p1 corner node
/// \param p2 corner node
/// \param p3 corner node
/// \param p4 corner node
/// \param p5 corner node
/// \param p6 corner node
/// \param p7 corner node
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
/// \param p0 corner node
/// \param p1 corner node
/// \param p2 corner node
/// \param p3 corner node
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

  ///
  /// Constructor that takes material tangent
  /// \param A material tangent (fourth-order tensor)
  ///
  KOKKOS_INLINE_FUNCTION
  SphericalParametrization(Tensor4<T, N> const & A);

  ///
  /// Evaluate the localization determinant for the normal defined
  /// by the spherical angles and update the extrema found so far.
  /// \param parameters spherical angles (phi, theta)
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value> const & parameters);

  ///
  /// Unit normal for the given spherical angles.
  /// \param parameters spherical angles (phi, theta)
  /// \return unit normal on the sphere
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 2>::value> const & parameters) const;

  ///
  /// Smallest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  ///
  /// Largest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  ///
  /// Parameters at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_minimum() const {return arg_minimum_;}

  ///
  /// Parameters at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_maximum() const {return arg_maximum_;}

  ///
  /// Normal at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  ///
  /// Normal at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  ///
  /// Material tangent
  ///
  Tensor4<T, N> const &
  tangent_;

  ///
  /// Minimum of the localization determinant
  ///
  T
  minimum_;

  ///
  /// Parameters of the minimum
  ///
  Vector<T, 2>
  arg_minimum_;

  ///
  /// Maximum of the localization determinant
  ///
  T
  maximum_;

  ///
  /// Parameters of the maximum
  ///
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

  ///
  /// Constructor that takes material tangent
  /// \param A material tangent (fourth-order tensor)
  ///
  KOKKOS_INLINE_FUNCTION
  StereographicParametrization(Tensor4<T, N> const & A);

  ///
  /// Unit normal from the inverse stereographic projection.
  /// \param parameters stereographic projection coordinates (x, y)
  /// \return unit normal on the sphere
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 2>::value> const & parameters) const;

  ///
  /// Evaluate the localization determinant for the normal defined
  /// by the stereographic coordinates and update the extrema found so far.
  /// \param parameters stereographic projection coordinates (x, y)
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value> const & parameters);

  ///
  /// Smallest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  ///
  /// Largest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  ///
  /// Parameters at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_minimum() const {return arg_minimum_;}

  ///
  /// Parameters at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_maximum() const {return arg_maximum_;}

  ///
  /// Normal at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  ///
  /// Normal at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  ///
  /// Material tangent
  ///
  Tensor4<T, N> const &
  tangent_;

  ///
  /// Minimum of the localization determinant
  ///
  T
  minimum_;

  ///
  /// Parameters of the minimum
  ///
  Vector<T, 2>
  arg_minimum_;

  ///
  /// Maximum of the localization determinant
  ///
  T
  maximum_;

  ///
  /// Parameters of the maximum
  ///
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

  ///
  /// Evaluate the localization determinant for the normal defined
  /// by the projective coordinates and update the extrema found so far.
  /// \param parameters homogeneous (projective) coordinates (x, y, z)
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 3>::value> const & parameters);

  ///
  /// Normalized normal for the given homogeneous coordinates.
  /// \param parameters homogeneous (projective) coordinates (x, y, z)
  /// \return unit normal, or (1, 1, 1) for the zero vector
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 3>::value> const & parameters) const;

  ///
  /// Smallest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  ///
  /// Largest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  ///
  /// Parameters at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_minimum() const {return arg_minimum_;}

  ///
  /// Parameters at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_maximum() const {return arg_maximum_;}

  ///
  /// Normal at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  ///
  /// Normal at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  ///
  /// Material tangent
  ///
  Tensor4<T, N> const &
  tangent_;

  ///
  /// Minimum of the localization determinant
  ///
  T
  minimum_;

  ///
  /// Parameters of the minimum
  ///
  Vector<T, 3>
  arg_minimum_;

  ///
  /// Maximum of the localization determinant
  ///
  T
  maximum_;

  ///
  /// Parameters of the maximum
  ///
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
  /// Evaluate the localization determinant for the normal defined
  /// by the tangent coordinates and update the extrema found so far.
  /// \param parameters tangent coordinates (x, y)
  ///
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value> const & parameters);

  ///
  /// Unit normal from the exponential map of the tangent coordinates.
  /// \param parameters tangent coordinates (x, y)
  /// \return unit normal on the sphere
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 2>::value> const & parameters) const;

  ///
  /// Smallest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  ///
  /// Largest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  ///
  /// Parameters at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_minimum() const {return arg_minimum_;}

  ///
  /// Parameters at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 2>
  get_arg_maximum() const {return arg_maximum_;}

  ///
  /// Normal at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  ///
  /// Normal at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  ///
  /// Material tangent
  ///
  Tensor4<T, N> const &
  tangent_;

  ///
  /// Minimum of the localization determinant
  ///
  T
  minimum_;

  ///
  /// Parameters of the minimum
  ///
  Vector<T, 2>
  arg_minimum_;

  ///
  /// Maximum of the localization determinant
  ///
  T
  maximum_;

  ///
  /// Parameters of the maximum
  ///
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

  ///
  /// Constructor that takes material tangent
  /// \param A material tangent (fourth-order tensor)
  ///
  KOKKOS_INLINE_FUNCTION
  CartesianParametrization(Tensor4<T, N> const & A);

  ///
  /// Evaluate the localization determinant for the normal defined
  /// by the Cartesian coordinates and update the extrema found so far.
  /// \param parameters Cartesian components (x, y, z)
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 3>::value> const & parameters);

  ///
  /// Normal formed directly from the Cartesian components.
  /// Note: not normalized.
  /// \param parameters Cartesian components (x, y, z)
  /// \return normal vector (x, y, z)
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal(Vector<T, dimension_const<N, 3>::value> const & parameters) const;

  ///
  /// Smallest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  ///
  /// Largest localization determinant found so far.
  ///
  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  ///
  /// Parameters at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_minimum() const {return arg_minimum_;}

  ///
  /// Parameters at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 3>
  get_arg_maximum() const {return arg_maximum_;}

  ///
  /// Normal at which the minimum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  ///
  /// Normal at which the maximum determinant occurs.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, N>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  ///
  /// Material tangent
  ///
  Tensor4<T, N> const &
  tangent_;

  ///
  /// Minimum of the localization determinant
  ///
  T
  minimum_;

  ///
  /// Parameters of the minimum
  ///
  Vector<T, 3>
  arg_minimum_;

  ///
  /// Maximum of the localization determinant
  ///
  T
  maximum_;

  ///
  /// Parameters of the maximum
  ///
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
  /// Traverse the grid, applying the visitor to each grid point.
  /// \param visitor functor to apply at each grid point
  ///
  template<typename Visitor>
  KOKKOS_INLINE_FUNCTION
  void
  traverse(Visitor & visitor) const;

private:

  ///
  /// Lower limits of the grid
  ///
  Vector<T, N>
  lower_;

  ///
  /// Upper limits of the grid
  ///
  Vector<T, N>
  upper_;

  ///
  /// Number of grid points in each dimension
  ///
  Vector<Index, N>
  points_per_dimension_;

};

} // namespace minitensor

namespace minitensor {

//
// Helper functions for determining the type of element
//
namespace {

KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type_1D(Index const nodes)
{
  switch (nodes) {
    case 2:   return ELEMENT::SEGMENTAL;
    default:  return ELEMENT::UNKNOWN;
  }
}

KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type_2D(Index const nodes)
{
  switch (nodes) {
    case 3:   return ELEMENT::TRIANGULAR;
    case 4:   return ELEMENT::QUADRILATERAL;
    default:  return ELEMENT::UNKNOWN;
  }
}

KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type_3D(Index const nodes)
{
  switch (nodes) {
    case 4:   return ELEMENT::TETRAHEDRAL;
    case 8:   return ELEMENT::HEXAHEDRAL;
    default:  return ELEMENT::UNKNOWN;
  }
}

} // anonymous namespace


//
//
//
KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type(Index const dimension, Index const number_nodes)
{

  ELEMENT::Type
  type = ELEMENT::UNKNOWN;

  switch (dimension) {

    case 1:
      type = find_type_1D(number_nodes);
      break;

    case 2:
      type = find_type_2D(number_nodes);
      break;

    case 3:
      type = find_type_3D(number_nodes);
      break;

    default:
      break;

  }

  if (type == ELEMENT::UNKNOWN) {
    MT_ERROR_EXIT("Unknown element type.");
  }

  return type;
}

//
// Constructor for SphericalParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
SphericalParametrization<T, N>::SphericalParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for SphericalParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
SphericalParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 2>::value> const & parameters
) const
{
  T const &
  phi = parameters(0);

  T const &
  theta = parameters(1);

  Vector<T, N>
  normal(sin(phi) * sin(theta), cos(phi), sin(phi) * cos(theta));

  return normal;
}

//
// Evaluation for SphericalParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
SphericalParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 2>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 2);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for StereographicParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
StereographicParametrization<T, N>::StereographicParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for StereographicParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
StereographicParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 2>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const
  r2 = x * x + y * y;

  Vector<T, N>
  normal(2.0 * x, 2.0 * y, r2 - 1.0);

  normal /= (r2 + 1.0);

  return normal;
}

//
// Evaluation for StereographicParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
StereographicParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 2>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 2);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for ProjectiveParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
ProjectiveParametrization<T, N>::ProjectiveParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for ProjectiveParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
ProjectiveParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 3>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const &
  z = parameters(2);

  Vector<T, N>
  normal(x, y, z);

  T const
  n = norm(normal);

  if (n > 0.0) {
    normal /= n;
  } else {
    normal = Vector<T, N>(1.0, 1.0, 1.0);
  }

  return normal;
}

//
// Evaluation for ProjectiveParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
ProjectiveParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 3>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 3);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for TangentParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
TangentParametrization<T, N>::TangentParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for TangentParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
TangentParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 2>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const
  r = std::sqrt(x * x + y * y);

  Vector<T, N>
  normal(3, Filler::ZEROS);

  if (r > 0.0) {
    normal(0) = x * sin(r) / r;
    normal(1) = y * sin(r) / r;
    normal(2) = cos(r);
  } else {
    normal(2) = 1.0;
  }

  return normal;
}

//
// Evaluation for TangentParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
TangentParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 2>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 2);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for CartesianParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
CartesianParametrization<T, N>::CartesianParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
  return;
#endif
}

//
// Normal vector for CartesianParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
CartesianParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 3>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const
  z = parameters(2);

  Vector<T, N> const
  normal(x, y, z);

  return normal;
}

//
// Evaluation for CartesianParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
CartesianParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 3>::value> const & parameters
)
{
  Vector<T, N>
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for ParametricGrid
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
ParametricGrid<T, N>::ParametricGrid(
    Vector<T, N> const & lower,
    Vector<T, N> const & upper,
    Vector<Index, N> const & points_per_dimension)
{
  assert(lower.get_dimension() == upper.get_dimension());
  assert(lower.get_dimension() == points_per_dimension.get_dimension());

  lower_ = lower;
  upper_ = upper;
  points_per_dimension_ = points_per_dimension;

  return;
}

//
// Traverse the grid and apply the visitor to each point.
//
template<typename T, Index N>
template<typename Visitor>
KOKKOS_INLINE_FUNCTION
void
ParametricGrid<T, N>::traverse(Visitor & visitor) const
{
  // Loop over the grid
  Index const
  number_parameters = lower_.get_dimension();

  LongIndex
  total_number_points = 1;

  for (Index dimension = 0; dimension < number_parameters; ++dimension) {
    total_number_points *= points_per_dimension_(dimension);
  }

  Vector<LongIndex, N> steps(number_parameters, Filler::ONES);

  for (Index dimension = 1; dimension < number_parameters; ++dimension) {
    steps(dimension) =
        steps(dimension - 1) * points_per_dimension_(dimension - 1);
  }

  Vector<Index, N>
  indices(number_parameters, Filler::ZEROS);

  Vector<T, N>
  position_in_grid(number_parameters, Filler::ZEROS);

  Vector<T, N> const
  span = upper_ - lower_;

  for (LongIndex point = 1;  point <= total_number_points; ++point) {

    //std::cout << "Indices : ";

    for (Index dimension = 0; dimension < number_parameters; ++dimension) {

/*           if ( points_per_dimension_(dimension) == 1 ) {     

             position_in_grid(dimension) = lower_(dimension);
        
          } else {
      
            position_in_grid(dimension) = indices(dimension) * span(dimension) /
                    (points_per_dimension_(dimension) - 1) + lower_(dimension);
          }
*/
        position_in_grid(dimension) = indices(dimension) * span(dimension) /
              (points_per_dimension_(dimension) - 1) + lower_(dimension);
        visitor(position_in_grid);

      //std::cout << indices(dimension) << " ";

      // Check if index needs to be increased or rolled back
      if (point % steps(dimension) == 0) {
        ++indices(dimension);
      }
      if (indices(dimension) == points_per_dimension_(dimension)) {
        indices(dimension) = 0;
      }

    }
     visitor(position_in_grid);
    //std::cout << std::endl;
    //std::cout << "Position : " << position_in_grid << std::endl;

  }

  return;
}

} // namespace minitensor
#include <iterator>

namespace minitensor {

//
// Length of a segment
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
length(Vector<T, N> const & p0, Vector<T, N> const & p1)
{
  return norm(p1 - p0);
}

//
// Area of a triangle
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
area(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2)
{
  Vector<T, N> const u = p1 - p0;
  Vector<T, N> const v = p2 - p0;

  T const area = 0.5 * norm(cross(u, v));

  return area;
}

//
// Area of a quadrilateral.
// Taken from:
// Calculation of the volume of a general hexahedron for flow predictions
// Davies, D. E.; Salmond, D. J.
// AIAA Journal (ISSN 0001-1452), vol. 23, June 1985, p. 954-956.
// Their vertex naming convention: ABCD
// Our convention:                 0123
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
area(
    Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3)
{
  Vector<T, N> const v31 = p1 - p3;
  Vector<T, N> const v02 = p2 - p0;

  T const area = 0.5 * norm(cross(v31, v02));

  return area;
}

//
// Volume of tetrahedron
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
volume(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3)
{
  Vector<T, N> const u = p1 - p0;
  Vector<T, N> const v = p2 - p0;
  Vector<T, N> const w = p3 - p0;

  T const volume =  std::abs(dot(u, cross(v, w)))/ 6.0;

  return volume;
}

//
// Volume of hexahedron
// Taken from:
// Calculation of the volume of a general hexahedron for flow predictions
// Davies, D. E.; Salmond, D. J.
// AIAA Journal (ISSN 0001-1452), vol. 23, June 1985, p. 954-956.
// Their vertex naming convention: ABCDEFGH
// Our convention:                 45670123
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
volume(
    Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3,
    Vector<T, N> const & p4, Vector<T, N> const & p5,
    Vector<T, N> const & p6, Vector<T, N> const & p7)
{
  Vector<T, N> const v24 = p4 - p2;
  Vector<T, N> const v25 = p5 - p2;
  Vector<T, N> const v20 = p0 - p2;
  Vector<T, N> const v27 = p7 - p2;

  Vector<T, N> const v75 = p5 - p7;
  Vector<T, N> const v46 = p6 - p4;
  Vector<T, N> const v50 = p0 - p5;
  Vector<T, N> const v41 = p1 - p4;
  Vector<T, N> const v07 = p7 - p0;
  Vector<T, N> const v43 = p3 - p4;

  Vector<T, N> const v26 = p6 - p2;
  Vector<T, N> const v16 = p6 - p1;
  Vector<T, N> const v21 = p1 - p2;
  Vector<T, N> const v31 = p1 - p3;
  Vector<T, N> const v23 = p3 - p2;
  Vector<T, N> const v63 = p3 - p6;

  Vector<T, N> const v7546 = cross(v75, v46);
  Vector<T, N> const v5041 = cross(v50, v41);
  Vector<T, N> const v0743 = cross(v07, v43);
  Vector<T, N> const v2616 = cross(v26, v16);
  Vector<T, N> const v2131 = cross(v21, v31);
  Vector<T, N> const v2363 = cross(v23, v63);

  T const V1 = dot(v24, v7546 + v5041 + v0743);
  T const V2 = dot(v25, v7546 + v2616);
  T const V3 = dot(v20, v5041 + v2131);
  T const V4 = dot(v27, v0743 + v2363);

  T const volume = (V1 + V2 + V3 + V4) / 12.0;

  return volume;
}

//
// Centroids of segment, triangle, tetrahedron, quadrilateral
// and hexahedron
// For these we can just take the average of the vertices.
// WARNING: This is not the center of mass.
//
template<typename T, Index N>
Vector<T, N>
centroid(std::vector<Vector<T, N>> const & points)
{
  Vector<T, N> C(points[0].get_dimension());
  C.clear();
  typedef typename std::vector<Vector<T, N>>::size_type sizeT;
  sizeT const n = points.size();

  for (sizeT i = 0; i < n; ++i) {
    C += points[i];
  }
  return C / static_cast<T>(n);
}

//
// The surface normal of a face
// Input: 3 independent nodes on the face
// Output: unit normal vector
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
normal(Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2)
{
  // Construct 2 independent vectors
  Vector<T, N> const v0 = p1 - p0;
  Vector<T, N> const v1 = p2 - p0;

  Vector<T, N> const n = unit(cross(v0, v1));

  return n;
}

//
// Given 3 points p0, p1, p2 that define a plane
// determine if point p is in the same side of the normal
// to the plane as defined by the right hand rule.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_normal_side(
    Vector<T, N> const & p,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    T const tolerance)
{
  Vector<T, N> const v0 = p1 - p0;
  Vector<T, N> const v1 = p2 - p0;
  T const h = std::min(norm(v0), norm(v1));
  Vector<T, N> const n = unit(cross(v0, v1));
  Vector<T, N> const v = p - p0;

  T const s = dot(v, n);

  if (s < -tolerance * h) return false;

  return true;
}

//
// Given two iterators to a container of points,
// find the associated bounding box.
// \param start, end: define sequence of points
// \return vectors that define the bounding box
//
template<typename T, typename I, Index N>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, N>, Vector<T, N>>
bounding_box(I start, I end)
{
  I
  it = start;

  Vector<T, N>
  min = (*it);

  Vector<T, N>
  max = min;

  Index const
  dimension = min.get_dimension();

  ++it;

  for (; it != end; ++it) {

    Vector<T, N> const &
    point = (*it);

    for (Index i = 0; i < dimension; ++i) {
      T const s = point(i);
      if (s < min(i)) min(i) = s;
      if (s > max(i)) max(i) = s;
    }

  }

  return std::make_pair(min, max);
}

///
/// Bounding box of a sequence of points.
///
template<typename T, typename I>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, DYNAMIC>, Vector<T, DYNAMIC>>
bounding_box(I start, I end)
{
  return bounding_box<T, I, DYNAMIC>(start, end);
}

//
// Determine if a given point is inside a bounding box.
// \param p the point
// \param min, max points defining the box
// \return whether the point is inside
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_box(
    Vector<T, N> const & p,
    Vector<T, N> const & min,
    Vector<T, N> const & max)
{
  Index const
  dimension = p.get_dimension();

  assert(min.get_dimension() == dimension);
  assert(max.get_dimension() == dimension);

  for (Index i = 0; i < dimension; ++i) {
    T const & s = p(i);
    if (s < min(i)) return false;
    if (s > max(i)) return false;
  }

  return true;
}

//
// Generate random point inside bounding box
// \param min, max the bounding box
// \return p point inside box
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
random_in_box(Vector<T, N> const & min, Vector<T, N> const & max)
{
  Index const
  dimension = min.get_dimension();

  assert(max.get_dimension() == dimension);

  Vector<T, N> p(dimension);

  for (Index i = 0; i < dimension; ++i) {
    p(i) = (max(i) - min(i)) * T(std::rand())/T(RAND_MAX) + min(i);
  }

  return p;
}

//
// Given 4 points p0, p1, p2, p3 that define a tetrahedron
// determine if point p is inside it.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
in_tetrahedron(
    Vector<T, N> const & p,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3,
    T const tolerance)
{
  if (in_normal_side(p, p0, p1, p2, tolerance) == false) return false;
  if (in_normal_side(p, p0, p3, p1, tolerance) == false) return false;
  if (in_normal_side(p, p1, p3, p2, tolerance) == false) return false;
  if (in_normal_side(p, p2, p3, p0, tolerance) == false) return false;

  return true;
}

//
// Given 8 points that define a hexahedron
// determine if point p is inside it.
// Assumption: faces are planar
//
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
    T const tolerance)
{
  if (in_normal_side(p, p0, p1, p2, tolerance) == false) return false;
  if (in_normal_side(p, p0, p4, p5, tolerance) == false) return false;
  if (in_normal_side(p, p1, p5, p6, tolerance) == false) return false;
  if (in_normal_side(p, p2, p6, p7, tolerance) == false) return false;
  if (in_normal_side(p, p3, p7, p4, tolerance) == false) return false;
  if (in_normal_side(p, p4, p7, p6, tolerance) == false) return false;

  return true;
}

//
// Closest point
// \param p the point
// \param n vector of points to test
// \return index to closest point
//
template<typename T, Index N>
typename std::vector<Vector<T, N>>::size_type
closest_point(Vector<T, N> const & p, std::vector<Vector<T, N>> const & n)
{
  assert(n.size() > 0);

  typename std::vector<Vector<T, N>>::size_type
  index = 0;

  Vector<T, N> const
  v0 = p - n[0];

  T
  min = norm_square(v0);

  for (typename std::vector<Vector<T, N>>::size_type i = 1;
      i < n.size();
      ++i) {

    Vector<T, N> const
    vi = p - n[i];

    T const
    s = norm_square(vi);

    if (s < min) {
      min = s;
      index = i;
    }

  }

  return index;
}

// Median of a sequence defined by random
// access iterators. Undefined for empty set.
// \param begin, end Iterators that define the sequence
// \return median of sequence
//
template<typename T, typename Iterator>
KOKKOS_INLINE_FUNCTION
T
median(Iterator begin, Iterator end)
{
  // Firewall
  if (begin == end) {
    MT_ERROR_EXIT("Median undefined for empty set.");
  }

  Index const
  size = static_cast<Index>(std::distance(begin, end));

  T
  median;

  Index const
  mid_index = size / 2;

  Iterator
  mid_iterator = begin + mid_index;
  std::partial_sort(begin, mid_iterator, end);

  if (size % 2 == 0) {

    // Even number of elements
    T const
    b = *mid_iterator;

    Iterator
    previous = mid_iterator - 1;

    T const
    a = *previous;

    median = (a + b) / 2.0;

  } else {

    // Odd number of elements
    median = *mid_iterator;

  }

  return median;
}

//
// Given quadrilateral nodes and a position
// in parametric coordinates, interpolate.
// \param xi position in parametric coordinates
// \param p0 ... corner nodes
// \return interpolated position
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_quadrilateral(
    Vector<T, dimension_const<N, 2>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3)
{

  T const
  N0 = 0.25 * (1 - xi(0)) * (1 - xi(1));

  T const
  N1 = 0.25 * (1 + xi(0)) * (1 - xi(1));

  T const
  N2 = 0.25 * (1 + xi(0)) * (1 + xi(1));

  T const
  N3 = 0.25 * (1 - xi(0)) * (1 + xi(1));

  Vector<T, N> const
  p = N0 * p0 + N1 * p1 + N2 * p2 + N3 * p3;

  return p;
}

//
// Given triangle nodes and a position
// in parametric coordinates, interpolate.
// \param xi position in parametric coordinates
// \param p0 ... corner nodes
// \return interpolated position
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_triangle(
    Vector<T, dimension_const<N, 3>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2)
{
  xi(2) = 1.0 - xi(0) - xi(1);

  Vector<T, N> const
  p = xi(0) * p0 + xi(1) * p1 + xi(2) * p2;

  return p;
}

//
// Given hexahedron nodes and a position
// in parametric coordinates, interpolate.
// \param xi position in parametric coordinates
// \param p0 ... corner nodes
// \return interpolated position
//
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
    Vector<T, N> const & p7)
{

  T const
  N0 = 0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 - xi(2));

  T const
  N1 = 0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 - xi(2));

  T const
  N2 = 0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 - xi(2));

  T const
  N3 = 0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 - xi(2));

  T const
  N4 = 0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 + xi(2));

  T const
  N5 = 0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 + xi(2));

  T const
  N6 = 0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 + xi(2));

  T const
  N7 = 0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 + xi(2));

  Vector<T, N> const
  p =
      N0 * p0 + N1 * p1 + N2 * p2 + N3 * p3 +
      N4 * p4 + N5 * p5 + N6 * p6 + N7 * p7;

  return p;
}

//
// Given tetrahedron nodes and a position
// in parametric coordinates, interpolate.
// \param xi position in parametric coordinates
// \param p0 ... corner nodes
// \return interpolated position
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_tetrahedron(
    Vector<T, dimension_const<N, 4>::value> & xi,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3)
{
  xi(3) = 1.0 - xi(0) - xi(1) - xi(2);

  Vector<T, N> const
  p = xi(0) * p0 + xi(1) * p1 + xi(2) * p2 + xi(3) * p3;

  return p;
}

//
// Given element type and nodes and a position
// in parametric coordinates, interpolate.
// \param type element type
// \param xi position in parametric coordinates
// \param v ... corner nodes
// \return interpolated position
//
template<typename T, Index M, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
interpolate_element(
    ELEMENT::Type element_type,
    Vector<T, M> &xi,
    std::vector<Vector<T, N>> const &v)
{
  Vector<T, N> p;

  switch (element_type) {

    case ELEMENT::TRIANGULAR:
      p = interpolate_triangle(xi, v[0], v[1], v[2]);
      break;

    case ELEMENT::QUADRILATERAL:
      p = interpolate_quadrilateral(xi, v[0], v[1], v[2], v[3]);
      break;

    case ELEMENT::TETRAHEDRAL:
      p = interpolate_tetrahedron(xi, v[0], v[1], v[2], v[3]);
      break;

    case ELEMENT::HEXAHEDRAL:
      p = interpolate_hexahedron(
          xi, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
      break;

    default:
      MT_ERROR_EXIT("Unknown element type in interpolation.");
      break;

  }

  return p;
}

//
// Given a vector of points, determine
// distances between all of them.
// \param vector of points
// \return distance matrix
//
template<typename T, Index N>
std::vector< std::vector<T>>
distance_matrix(std::vector<Vector<T, N>> const & points)
{
  Index const
  number_points = points.size();

  std::vector< std::vector<T>>
  distances(number_points);

  for (Index i = 0; i < number_points; ++i) {

    distances[i].resize(number_points);

    distances[i][i] = 0.0;

    for (Index j = i + 1; j < number_points; ++j) {

      T const
      distance = norm(points[i] - points[j]);

      distances[i][j] = distance;
      distances[j][i] = distance;

    }

  }

  return distances;
}

//
// Given a distance matrix, determine the minimum
// distance between two distinct points.
// \param distance matrix
// \return minimum distance
//
template<typename T>
std::vector<T>
minimum_distances(std::vector< std::vector<T>> const & distances)
{
  Index const
  number_points = distances.size();

  std::vector<T>
  minima(number_points);

  // First row
  T
  minimum = distances[0][1];

  for (Index j = 2; j < number_points; ++j) {
    minimum = std::min(minimum, distances[0][j]);
  }

  minima[0] = minimum;

  // Remaining rows
  for (Index i = 1; i < number_points; ++i) {

    minimum = distances[i][0];

    for (Index j = 1; j < number_points; ++j) {

      if (i == j) continue;

      minimum = std::min(minimum, distances[i][j]);

    }

    minima[i] = minimum;

  }

  return minima;
}

/// @}
} // namespace minitensor

#endif // MiniTensor_Geometry_h
