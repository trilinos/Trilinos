// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid2_MiniTensor_Geometry_h)
#define Intrepid2_MiniTensor_Geometry_h

#include <vector>
#include "Intrepid2_MiniTensor_Tensor.h"

namespace Intrepid2 {

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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
length(Vector<T, N, ES> const & p0, Vector<T, N, ES> const & p1);

///
/// Area of a triangle
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
area(Vector<T, N, ES> const & p0, Vector<T, N, ES> const & p1,
     Vector<T, N, ES> const & p2);

///
/// Area of a quadrilateral, assummed planar. If not planar, returns
/// the sum of the areas of the two triangles p0,p1,p2 and p0,p2,p3
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
area(Vector<T, N, ES> const & p0, Vector<T, N, ES> const & p1,
     Vector<T, N, ES> const & p2, Vector<T, N, ES> const & p3);

///
/// Volume of tetrahedron
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
volume(Vector<T, N, ES> const & p0, Vector<T, N, ES> const & p1,
       Vector<T, N, ES> const & p2, Vector<T, N, ES> const & p3);

///
/// Volume of pyramid of quadrilateral base
/// Base is assumed planar
/// Base is p0,p1,p2,p3
/// Apex is p4
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
volume(Vector<T, N, ES> const & p0, Vector<T, N, ES> const & p1,
       Vector<T, N, ES> const & p2, Vector<T, N, ES> const & p3,
       Vector<T, N, ES> const & p4);

///
/// Volume of hexahedron
/// Assumption: all faces are planar
/// Decompose into 3 pyramids
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
volume(Vector<T, N, ES> const & p0, Vector<T, N, ES> const & p1,
       Vector<T, N, ES> const & p2, Vector<T, N, ES> const & p3,
       Vector<T, N, ES> const & p4, Vector<T, N, ES> const & p5,
       Vector<T, N, ES> const & p6, Vector<T, N, ES> const & p7);

///
/// Centroids of segment, triangle, tetrahedron, quadrilateral
/// and hexahedron.
/// For these we can just take the average of the vertices.
/// WARNING: This is not the center of mass.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
centroid(std::vector<Vector<T, N, ES> > const & points);

///
/// The surface normal of a face
/// Input: 3 independent nodes on the face
/// Output: normal vector
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
normal(Vector<T, N, ES> const & p0,
       Vector<T, N, ES> const & p1,
       Vector<T, N, ES> const & p2);

///
/// Given 3 points p0, p1, p2 that define a plane
/// determine if point p is in the same side of the normal
/// to the plane as defined by the right hand rule.
/// If a tolrance is given, use that as criterion for minimal distance.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
in_normal_side(
    Vector<T, N, ES> const & p,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2,
    T const tolerance = 0);

///
/// Given two iterators to a container of points,
/// find the associated bounding box.
/// \param start end: define sequence of points
/// \return vectors that define the bounding box
///
template<typename T, typename I, Index N,  typename ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
std::pair< Vector<T, N, ES>, Vector<T, N, ES>>
bounding_box(I start, I end);

template<typename T, typename I,  typename ES=NOKOKKOS>
KOKKOS_INLINE_FUNCTION
std::pair< Vector<T, DYNAMIC, ES>, Vector<T, DYNAMIC, ES>>
bounding_box(I start, I end);

///
/// Determine if a given point is inside a bounding box.
/// \param p the point
/// \param min max points defining the box
/// \return whether the point is inside
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
in_box(
    Vector<T, N, ES> const & p,
    Vector<T, N, ES> const & min,
    Vector<T, N, ES> const & max);

///
/// Generate random point inside bounding box
/// \param min max the bounding box
/// \return p point inside box
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
random_in_box(
    Vector<T, N, ES> const & min,
    Vector<T, N, ES> const & max);

///
/// Given 4 points p0, p1, p2, p3 that define a tetrahedron
/// determine if point p is inside it.
/// If a tolrance is given, use that as criterion for minimal distance.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
in_tetrahedron(
    Vector<T, N, ES> const & p,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2,
    Vector<T, N, ES> const & p3,
    T const tolerance = 0);

///
/// Given 8 points that define a hexahedron
/// determine if point p is inside it.
/// Assumption: faces are planar
/// If a tolrance is given, use that as criterion for minimal distance.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
bool
in_hexahedron(
    Vector<T, N, ES> const & p,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2,
    Vector<T, N, ES> const & p3,
    Vector<T, N, ES> const & p4,
    Vector<T, N, ES> const & p5,
    Vector<T, N, ES> const & p6,
    Vector<T, N, ES> const & p7,
    T const tolerance = 0);

///
/// Closest point
/// \param p the point
/// \param n vector of points to test
/// \return index to closest point
///
template<typename T, Index N,  typename ES>
typename std::vector< Vector<T, N, ES> >::size_type
closest_point(Vector<T, N, ES> const & p, std::vector< Vector<T, N, ES> > const & n);

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
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
interpolate_quadrilateral(
    Vector<T, dimension_const<N, 2>::value, ES> & xi,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2,
    Vector<T, N, ES> const & p3);

///
/// Given triangle nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
interpolate_triangle(
    Vector<T, dimension_const<N, 3>::value, ES> & xi,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2);

///
/// Given hexahedron nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
interpolate_hexahedron(
    Vector<T, dimension_const<N, 3>::value, ES> & xi,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2,
    Vector<T, N, ES> const & p3,
    Vector<T, N, ES> const & p4,
    Vector<T, N, ES> const & p5,
    Vector<T, N, ES> const & p6,
    Vector<T, N, ES> const & p7);

///
/// Given tetrahedron nodes and a position
/// in parametric coordinates, interpolate.
/// \param xi position in parametric coordinates
/// \param p0 ... corner nodes
/// \return interpolated position
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
interpolate_tetrahedron(
    Vector<T, dimension_const<N, 4>::value, ES> & xi,
    Vector<T, N, ES> const & p0,
    Vector<T, N, ES> const & p1,
    Vector<T, N, ES> const & p2,
    Vector<T, N, ES> const & p3);

///
/// Given element type and nodes and a position
/// in parametric coordinates, interpolate.
/// \param element_type element type
/// \param xi position in parametric coordinates
/// \param v ... corner nodes
/// \return interpolated position
///
template<typename T, Index M, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
interpolate_element(
    ELEMENT::Type element_type,
    Vector<T, M, ES> & xi,
    std::vector< Vector<T, N, ES> > const & v);

///
/// Given a vector of points, determine
/// distances between all of them.
/// \param points vector of points
/// \return distance matrix
///
template<typename T, Index N,  typename ES>
std::vector< std::vector<T> >
distance_matrix(std::vector< Vector<T, N, ES> > const & points);

///
/// Given a distance matrix, determine the minimum
/// distance between two distinct points.
/// \param distances distance matrix
/// \return minimum distance
///
template<typename T,  typename ES>
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
template<typename T, Index N,  typename ES=NOKOKKOS>
class SphericalParametrization
{
public:

  KOKKOS_INLINE_FUNCTION
  SphericalParametrization(Tensor4<T, N, ES> const & A);

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value, ES> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal(Vector<T, dimension_const<N, 2>::value, ES> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2, ES>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2, ES>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N, ES> const &
  tangent_;

  T
  minimum_;

  Vector<T, 2, ES>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 2, ES>
  arg_maximum_;
};

///
/// Stereographic parametrization functor
///
template<typename T, Index N,  typename ES=NOKOKKOS>
class StereographicParametrization
{
public:

  KOKKOS_INLINE_FUNCTION
  StereographicParametrization(Tensor4<T, N, ES> const & A);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal(Vector<T, dimension_const<N, 2>::value, ES> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value, ES> const & parameters);

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2, ES>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2, ES>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N, ES> const &
  tangent_;

  T
  minimum_;

  Vector<T, 2, ES>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 2, ES>
  arg_maximum_;
};

///
/// Projective parametrization functor
///
template<typename T, Index N,  typename ES=NOKOKKOS>
class ProjectiveParametrization
{
public:

  ///
  /// Constructor that takes material tangent
  ///
  KOKKOS_INLINE_FUNCTION
  ProjectiveParametrization(Tensor4<T, N, ES> const & A);

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 3>::value, ES> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal(Vector<T, dimension_const<N, 3>::value, ES> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3, ES>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3, ES>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N, ES> const &
  tangent_;

  T
  minimum_;

  Vector<T, 3, ES>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 3, ES>
  arg_maximum_;
};

///
/// Tangent parametrization functor
///
template<typename T, Index N,  typename ES=NOKOKKOS>
class TangentParametrization
{
public:

  ///
  /// Constructor that takes material tangent
  ///
  KOKKOS_INLINE_FUNCTION
  TangentParametrization(Tensor4<T, N, ES> const & A);

  ///
  ///
  ///
  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 2>::value, ES> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal(Vector<T, dimension_const<N, 2>::value, ES> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2, ES>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 2, ES>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N, ES> const &
  tangent_;

  T
  minimum_;

  Vector<T, 2, ES>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 2, ES>
  arg_maximum_;
};

///
/// Cartesian parametrization functor
///
template<typename T, Index N,  typename ES=NOKOKKOS>
class CartesianParametrization
{
public:

  KOKKOS_INLINE_FUNCTION
  CartesianParametrization(Tensor4<T, N, ES> const & A);

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Vector<T, dimension_const<N, 3>::value, ES> const & parameters);

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal(Vector<T, dimension_const<N, 3>::value, ES> const & parameters) const;

  KOKKOS_INLINE_FUNCTION
  T
  get_minimum() const {return minimum_;}

  KOKKOS_INLINE_FUNCTION
  T
  get_maximum() const {return maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3, ES>
  get_arg_minimum() const {return arg_minimum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, 3, ES>
  get_arg_maximum() const {return arg_maximum_;}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_minimum() const {return get_normal(arg_minimum_);}

  KOKKOS_INLINE_FUNCTION
  Vector<T, N, ES>
  get_normal_maximum() const {return get_normal(arg_maximum_);}

private:

  Tensor4<T, N, ES> const &
  tangent_;

  T
  minimum_;

  Vector<T, 3, ES>
  arg_minimum_;

  T
  maximum_;

  Vector<T, 3, ES>
  arg_maximum_;
};

///
/// Parametric grid class
///
template<typename T, Index N,  typename ES=NOKOKKOS>
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
      Vector<T, N, ES> const & lower,
      Vector<T, N, ES> const & upper,
      Vector<Index, N, ES> const & points_per_dimension);

  ///
  ///
  template<typename Visitor>
  KOKKOS_INLINE_FUNCTION
  void
  traverse(Visitor & visitor) const;

private:

  Vector<T, N, ES>
  lower_;

  Vector<T, N, ES>
  upper_;

  Vector<Index, N, ES>
  points_per_dimension_;

};

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Geometry.i.h"
#include "Intrepid2_MiniTensor_Geometry.t.h"

#endif // Intrepid2_MiniTensor_Geometry_h
