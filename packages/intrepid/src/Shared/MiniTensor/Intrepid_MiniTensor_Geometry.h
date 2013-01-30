// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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

#if !defined(Intrepid_MiniTensor_Geometry_h)
#define Intrepid_MiniTensor_Geometry_h

#include <vector>
#include "Intrepid_MiniTensor_Tensor.h"

namespace Intrepid {

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
  }

  ///
  /// Length of a segment
  ///
  template<typename T>
  T
  length(Vector<T> const & p0, Vector<T> const & p1);

  ///
  /// Area of a triangle
  ///
  template<typename T>
  T
  area(Vector<T> const & p0, Vector<T> const & p1,
      Vector<T> const & p2);

  ///
  /// Area of a quadrilateral, assummed planar. If not planar, returns
  /// the sum of the areas of the two triangles p0,p1,p2 and p0,p2,p3
  ///
  template<typename T>
  T
  area(Vector<T> const & p0, Vector<T> const & p1,
      Vector<T> const & p2, Vector<T> const & p3);

  ///
  /// Volume of tetrahedron
  ///
  template<typename T>
  T
  volume(Vector<T> const & p0, Vector<T> const & p1,
      Vector<T> const & p2, Vector<T> const & p3);

  ///
  /// Volume of pyramid of quadrilateral base
  /// Base is assumed planar
  /// Base is p0,p1,p2,p3
  /// Apex is p4
  ///
  template<typename T>
  T
  volume(Vector<T> const & p0, Vector<T> const & p1,
      Vector<T> const & p2, Vector<T> const & p3,
      Vector<T> const & p4);

  ///
  /// Volume of hexahedron
  /// Assumption: all faces are planar
  /// Decompose into 3 pyramids
  ///
  template<typename T>
  T
  volume(Vector<T> const & p0, Vector<T> const & p1,
      Vector<T> const & p2, Vector<T> const & p3,
      Vector<T> const & p4, Vector<T> const & p5,
      Vector<T> const & p6, Vector<T> const & p7);

  ///
  /// Centroids of segment, triangle, tetrahedron, quadrilateral
  /// and hexahedron.
  /// For these we can just take the average of the vertices.
  /// WARNING: This is not the center of mass.
  ///
  template<typename T>
  Vector<T>
  centroid(std::vector<Vector<T> > const & points);

  ///
  /// The surface normal of a face
  /// Input: 3 independent nodes on the face
  /// Output: normal vector
  ///
  template<typename T>
  Vector<T>
  normal(Vector<T> const & p0,
          Vector<T> const & p1,
          Vector<T> const & p2);

  ///
  /// Given 3 points p0, p1, p2 that define a plane
  /// determine if point p is in the same side of the normal
  /// to the plane as defined by the right hand rule.
  ///
  template<typename T>
  bool
  in_normal_side(
      Vector<T> const & p,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2);

  ///
  /// Given two iterators to a container of points,
  /// find the associated bounding box.
  /// \param start, end: define sequence of points
  /// \return vectors that define the bounding box
  ///
  template<typename T, typename I>
  std::pair< Vector<T>, Vector<T> >
  bounding_box(I start, I end);

  ///
  /// Determine if a given point is inside a bounding box.
  /// \param p the point
  /// \param min, max points defining the box
  /// \return whether the point is inside
  ///
  template<typename T>
  bool
  in_box(
      Vector<T> const & p,
      Vector<T> const & min,
      Vector<T> const & max);

  ///
  /// Generate random point inside bounding box
  /// \param min, max the bounding box
  /// \return p point inside box
  ///
  template<typename T>
  Vector<T>
  random_in_box(
      Vector<T> const & min,
      Vector<T> const & max);

  ///
  /// Given 4 points p0, p1, p2, p3 that define a tetrahedron
  /// determine if point p is inside it.
  ///
  template<typename T>
  bool
  in_tetrahedron(
      Vector<T> const & p,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2,
      Vector<T> const & p3);

  ///
  /// Given 8 points that define a hexahedron
  /// determine if point p is inside it.
  /// Assumption: faces are planar
  ///
  template<typename T>
  bool
  in_hexahedron(
      Vector<T> const & p,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2,
      Vector<T> const & p3,
      Vector<T> const & p4,
      Vector<T> const & p5,
      Vector<T> const & p6,
      Vector<T> const & p7);

  ///
  /// Closest point
  /// \param p the point
  /// \param n vector of points to test
  /// \return index to closest point
  ///
  template<typename T>
  typename std::vector< Vector<T> >::size_type
  closest_point(Vector<T> const & p, std::vector< Vector<T> > const & n);

  /// Median of a sequence defined by random
  /// access iterators. Undefined for empty set.
  /// \param begin, end Iterators that define the sequence
  /// \return median of sequence
  ///
  template<typename T, typename Iterator>
  T
  median(Iterator begin, Iterator end);

  ///
  /// Given quadrilateral nodes and a position
  /// in parametric coordinates, interpolate.
  /// \param xi position in parametric coordinates
  /// \param p0 ... corner nodes
  /// \return interpolated position
  ///
  template<typename T>
  Vector<T>
  interpolate_quadrilateral(
      Vector<T> & xi,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2,
      Vector<T> const & p3);

  ///
  /// Given triangle nodes and a position
  /// in parametric coordinates, interpolate.
  /// \param xi position in parametric coordinates
  /// \param p0 ... corner nodes
  /// \return interpolated position
  ///
  template<typename T>
  Vector<T>
  interpolate_triangle(
      Vector<T> & xi,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2);

  ///
  /// Given hexahedron nodes and a position
  /// in parametric coordinates, interpolate.
  /// \param xi position in parametric coordinates
  /// \param p0 ... corner nodes
  /// \return interpolated position
  ///
  template<typename T>
  Vector<T>
  interpolate_hexahedron(
      Vector<T> & xi,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2,
      Vector<T> const & p3,
      Vector<T> const & p4,
      Vector<T> const & p5,
      Vector<T> const & p6,
      Vector<T> const & p7);

  ///
  /// Given tetrahedron nodes and a position
  /// in parametric coordinates, interpolate.
  /// \param xi position in parametric coordinates
  /// \param p0 ... corner nodes
  /// \return interpolated position
  ///
  template<typename T>
  Vector<T>
  interpolate_tetrahedron(
      Vector<T> & xi,
      Vector<T> const & p0,
      Vector<T> const & p1,
      Vector<T> const & p2,
      Vector<T> const & p3);

  ///
  /// Given element type and nodes and a position
  /// in parametric coordinates, interpolate.
  /// \param type element type
  /// \param xi position in parametric coordinates
  /// \param v ... corner nodes
  /// \return interpolated position
  ///
  template<typename T>
  Vector<T>
  interpolate_element(
      ELEMENT::Type element_type,
      Vector<T> & xi,
      std::vector< Vector<T> > const & v);

  ///
  /// Given a vector of points, determine
  /// distances between all of them.
  /// \param vector of points
  /// \return distance matrix
  ///
  template<typename T>
  std::vector< std::vector<T> >
  distance_matrix(std::vector< Vector<T> > const & points);

  ///
  /// Given a distance matrix, determine the minimum
  /// distance between two distinct points.
  /// \param distance matrix
  /// \return minimum distance
  ///
  template<typename T>
  std::vector<T>
  minimum_distances(std::vector< std::vector<T> > const & distances);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Geometry.i.h"
#include "Intrepid_MiniTensor_Geometry.t.h"

#endif // Intrepid_MiniTensor_Geometry_h
