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

#if !defined(Intrepid_MiniTensor_Geometry_t_h)
#define Intrepid_MiniTensor_Geometry_t_h

#include <iterator>

namespace Intrepid {

//
// Length of a segment
//
template<typename T, Index N>
T
length(Vector<T, N> const & p0, Vector<T, N> const & p1)
{
  Vector<T, N> const v = p1 - p0;
  return norm(v);
}

//
// Area of a triangle
//
template<typename T, Index N>
T
area(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2)
{
  Vector<T, N> const u = p1 - p0;
  Vector<T, N> const v = p2 - p0;

  T const base = norm(u);

  Vector<T, N> const i = u / base;
  Vector<T, N> const n = v - dot(v, i) * i;

  T const height = norm(n);
  T const area = 0.5 * base * height;

  return area;
}

//
// Area of a quadrilateral, assumed planar. If not planar, returns
// the sum of the areas of the two triangles p0,p1,p2 and p0,p2,p3
//
template<typename T, Index N>
T
area(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3)
{
  return area(p0, p1, p2) + area(p0, p2, p3);
}

//
// Volume of tetrahedron
//
template<typename T, Index N>
T
volume(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3)
{
  // Area of base triangle
  T const
  base = area(p0, p1, p2);

  // Height
  Vector<T, N> const u = p1 - p0;
  Vector<T, N> const v = p2 - p0;
  Vector<T, N> const w = p3 - p0;

  Vector<T, N> const i = u / norm(u);
  Vector<T, N> const j = v / norm(v);

  Vector<T, N> const n = w - dot(w, i) * i - dot(w, j) * j;

  T const height = norm(n);

  // Volume
  T const volume = base * height / 3.0;

  return volume;
}

//
// Volume of pyramid of quadrilateral base
// Base is assumed planar
// Base is p0,p1,p2,p3
// Apex is p4
//
template<typename T, Index N>
T
volume(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3,
    Vector<T, N> const & p4)
{
  // Area of base quadrilateral
  T const
  base = area(p0, p1, p2, p3);

  // Height
  Vector<T, N> const u = p1 - p0;
  Vector<T, N> const v = p2 - p0;
  Vector<T, N> const w = p4 - p0;

  Vector<T, N> const i = u / norm(u);
  Vector<T, N> const j = v / norm(v);

  Vector<T, N> const n = w - dot(w, i) * i - dot(w, j) * j;

  T const height = norm(n);

  // Volume
  T const volume = base * height / 3.0;

  return volume;
}

//
// Volume of hexahedron
// Assumption: all faces are planar
// Decompose into 3 pyramids
//
template<typename T, Index N>
T
volume(Vector<T, N> const & p0, Vector<T, N> const & p1,
    Vector<T, N> const & p2, Vector<T, N> const & p3,
    Vector<T, N> const & p4, Vector<T, N> const & p5,
    Vector<T, N> const & p6, Vector<T, N> const & p7)
{
  // 1st pyramid
  T const V1 = volume(p4, p7, p6, p5, p0);

  // 2nd pyramid
  T const V2 = volume(p3, p2, p6, p7, p0);

  // 3rd pyramid
  T const V3 = volume(p1, p5, p6, p2, p0);

  return V1 + V2 + V3;
}

//
// Centroids of segment, triangle, tetrahedron, quadrilateral
// and hexahedron
// For these we can just take the average of the vertices.
// WARNING: This is not the center of mass.
//
template<typename T, Index N>
Vector<T, N>
centroid(std::vector<Vector<T, N> > const & points)
{
  Vector<T, N> C(points[0].get_dimension());
  C.clear();
  typedef typename std::vector<Vector<T, N> >::size_type sizeT;
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
std::pair< Vector<T, N>, Vector<T, N> >
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

template<typename T, typename I>
std::pair< Vector<T, DYNAMIC>, Vector<T, DYNAMIC> >
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
bool
in_tetrahedron(
    Vector<T, N> const & p,
    Vector<T, N> const & p0,
    Vector<T, N> const & p1,
    Vector<T, N> const & p2,
    Vector<T, N> const & p3,
    T const tolerance)
{
  if (in_normal_side(p, p0, p1, p2, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p0, p3, p1, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p1, p3, p2, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p2, p3, p0, tolerance) == false) {

    return false;

  }

  return true;
}

//
// Given 8 points that define a hexahedron
// determine if point p is inside it.
// Assumption: faces are planar
//
template<typename T, Index N>
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
  if (in_normal_side(p, p0, p1, p2, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p0, p4, p5, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p1, p5, p6, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p2, p6, p7, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p3, p7, p4, tolerance) == false) {

    return false;

  } else if (in_normal_side(p, p4, p7, p6, tolerance) == false) {

    return false;

  }

  return true;
}

//
// Closest point
// \param p the point
// \param n vector of points to test
// \return index to closest point
//
template<typename T, Index N>
typename std::vector< Vector<T, N> >::size_type
closest_point(Vector<T, N> const & p, std::vector< Vector<T, N> > const & n)
{
  assert(n.size() > 0);

  typename std::vector< Vector<T, N> >::size_type
  index = 0;

  Vector<T, N> const
  v0 = p - n[0];

  T
  min = norm_square(v0);

  for (typename std::vector< Vector<T, N> >::size_type i = 1;
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
T
median(Iterator begin, Iterator end)
{
  // Firewall
  if (begin == end) {
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
    std::cerr << std::endl;
    std::cerr << "Median undefined for empty set." << std::endl;
    exit(1);
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

///
/// Given element type and nodes and a position
/// in parametric coordinates, interpolate.
/// \param type element type
/// \param xi position in parametric coordinates
/// \param v ... corner nodes
/// \return interpolated position
///
template<typename T, Index M, Index N>
Vector<T, N>
interpolate_element(
    ELEMENT::Type element_type,
    Vector<T, M> & xi,
    std::vector< Vector<T, N> > const & v)
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
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Unknown element type in interpolation.";
      std::cerr << std::endl;
      exit(1);
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
std::vector< std::vector<T> >
distance_matrix(std::vector< Vector<T, N> > const & points)
{
  Index const
  number_points = points.size();

  std::vector< std::vector<T> >
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
minimum_distances(std::vector< std::vector<T> > const & distances)
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

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Geometry_t_h
