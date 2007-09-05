// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
    \brief Header file for auxiliary classes providing linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/

#ifndef INTREPID_REALSPACE_HPP
#define INTREPID_REALSPACE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"

namespace Intrepid {

/** \class Intrepid::Point
    \brief Implements points (vectors) in 1D, 2D, and 3D Euclidean space, with 
  basic vector operations. The main purpose of this class is to provide
  basic algebraic functionality needed by the reconstruction methods.

  A Point object models a point/vector in 1-3D Euclidean space. It stores the 
  Cartesian coordinates of the point/vector and information about its type. The
  PointType is defined in Intrepid::Intrepid_Types. The default type of a point
  is AMBIENT, and indicates that the point instance models a point in the
  "physical" domain. The other available type is REFERENCE. It indicates that
  the point belongs to a cell in a "reference" coordinate frame. Such cells are
  used with FEM reconstructions. These reconstruction methods are based on
  pullbacks and need to be able to distinguish between points in ambient vs. 
  points in reference coordinates. 
  
  The class overloads several operators to provide basic operations on vectors:

  \arg ^ computes a cross (vector) product of two Point objects
  \arg + computes a sum of two Point objects
  \arg * computes a dot (inner) product of two Point objects
  \arg * computes a product of a scalar with a Point object
  \arg etc.

  I/O operators are also overloaded, and may be useful for debugging.  
*/
template<class Scalar>
class Point {
private:
  
  std::vector<Scalar> data;
  PointType point_type;

public:

  /** \brief Default destructor.
  */
  ~Point() {};

  /** \brief Default constructor.
  */
  Point() {
    data.resize(0); 
    point_type = AMBIENT;
  };

  /** \brief Copy constructor.
  */
  Point(const Point& right);

  /** \brief Create a 1-3D zero point (vector).
  */
  Point(const int dim, 
        const PointType point_type_ = AMBIENT);

  /** \brief Create a 1D point based on one Scalar coordinate.
  */
  Point(const Scalar x, 
        const PointType point_type_ = AMBIENT);

  /** \brief Create a 2D point based on two Scalar coordinates.
  */
  Point(const Scalar x, 
        const Scalar y, 
        const PointType point_type_ = AMBIENT);

  /** \brief Create a 2D point based on two int coordinates. This constructor
    helps to avoid undesired behavior when the user specifies two int coordinates.
    Without this constructor, the first int coordinate is interpreted as a 
    pointer and the wrong ctor (based on a Scalar pointer and dimension) is 
    invoked.
    */
  Point(const int x, 
        const int y, 
        const PointType point_type_ = AMBIENT);

  /** \brief Create a 3D point based on three Scalar coordinates.
  */
  Point(const Scalar x, 
        const Scalar y, 
        const Scalar z, 
        const PointType point_type_ = AMBIENT);

  /** \brief Create a 1-3D point based on a Scalar pointer and dimension.
  */
  Point(const Scalar* dataptr, 
        const int dim, 
        const PointType point_type_ = AMBIENT);

  /** \brief Return dimension of the point.
  */
  int getDim() const;
  
  /** \brief Returns the point type (reference or ambient)
    */
  PointType getType() const;
  
  /** \brief Returns a string with the name of the point type.
    */
  const char* getName() const;

  /** \brief Returns const reference to data member.
  */
  const std::vector<Scalar> & getData() const;
  
  /** \brief Load 1-3D point with new coordinate data, do not change its type.
    */
  void setData(const Scalar* dataptr, 
               const int dim);
  
  /** \brief Changes the type (AMBIENT or REFERENCE) of Point instance.
    */
  void setType(const PointType target_type);

  /** \brief Checks if a point with type = REFERENCE belongs to a reference cell
    */
  ErrorCode inRefCell(const CellType cell_type) const;
  
  /** \brief   Overloaded [] operator.
    */
  const Scalar & operator [] (const int coordinate_id_) const;
  
  
  /** \brief Returns Euclidian distance from <var>*this</var> to <var>v2</var>.
  */
  Scalar distance(const Point& v2) const;

  /** \brief Assignment operator <var>*this = right</var>.
  */
  Point& operator  = (const Point& right);

  /** \brief Cross (vector) product <var>*this = *this</var> X <var>right</var>.
  */
  Point& operator ^= (const Point& right);

  /** \brief Addition <var>*this = *this + right</var>.
  */
  Point& operator += (const Point& right);

  /** \brief Subtraction <var>*this = *this - right</var>.
  */
  Point& operator -= (const Point& right);

  /** \brief Scalar multiplication <var>*this = left</var> * <var>(*this)</var>.
  */
  Point& operator *= (const Scalar left);
}; // class Point

/** \relates Point
    Cross product.
*/
template<class Scalar>
const Point<Scalar> operator ^ (const Point<Scalar>& left, const Point<Scalar>& right);

/** \relates Point
    Addition.
*/
template<class Scalar>
const Point<Scalar> operator + (const Point<Scalar>& left, const Point<Scalar>& right);

/** \relates Point
    Subtraction.
*/
template<class Scalar>
const Point<Scalar> operator - (const Point<Scalar>& left, const Point<Scalar>& right);

/** \relates Point
    Scalar multiplication.
*/
template<class Scalar>
const Point<Scalar> operator * (const Scalar left, const Point<Scalar>& right);

/** \relates Point
    Dot product.
*/
template<class Scalar>
const Scalar operator * (const Point<Scalar>& left, const Point<Scalar>& right);

/** \relates Point
    Outputs a formated stream representing a Point. For debugging purposes.
*/
template<class Scalar>
std::ostream& operator << (std::ostream& os, const Point<Scalar>& point);
// end function declarations related to class Point



/** \class Intrepid::LinearMap
    \brief Linear operators in 1D, 2D, and 3D Euclidean space.

    A LinearMap models a linear operator in 1-3D Euclidean space (a SQUARE matrix).
    The main purpose of this class is to provide basic algebraic functionality
    needed by the interpolation/reconstruction method.
    The class provides several methods and overloaded operators for this purpose:

    \arg getElement extracts an entry
    \arg getColumn extracts a column
    \arg getRow extracts a row
    \arg * computes product of LinearMap with another LinearMap or a Point object
    \arg etc.
 */
template<class Scalar>
class LinearMap{
private:

    std::vector< std::vector<Scalar> > elements_;
public:

  /** \brief Default destructor.
  */
  ~LinearMap() {};

  /** \brief Default constructor.
  */
  LinearMap() {elements_.resize(0);};

  /** \brief Copy constructor.
  */
  LinearMap(const LinearMap& right);

  /** \brief Create <var>dim</var>-by-<var>dim</var> zero matrix.
  */
  LinearMap(const int dim);

  /** \brief Create <var>dim</var>-by-<var>dim</var> matrix, and fill it
    with Scalars referenced by <var>dataptr</var>. It is assumed that the Scalars
    are ordered by row, i.e., if dim = 3, dataptr contains
    
            (a, b, c, d ,e, f, g, h, i), then
    
                        | a b c |
            LinearMap = | d e f |
                        | g h i |
    
  */
  LinearMap(const Scalar* dataptr, const int dim);


  /** \brief Assignment operator <var>*this = right</var>.
  */
  LinearMap& operator = (const LinearMap& right);

  
  /** \brief Fills an existing <var>dim</var>-by-<var>dim</var> matrix
    with Scalars referenced by <var>dataptr</var>. It is assumed that the scalares are
    ordered by row.
    */
  void setData(const Scalar* dataptr, const int dim);
  
  
  /** \brief Return matrix dimension.
  */
  int getDim() const;

  /** \brief Return Scalar element in <var>rowID</var>-th row,  <var>colID</var>-th column.
  */
  Scalar getElement(int rowID, int colID) const;

  /** \brief Return <var>colID</var>-th column as Point.
  */
  Point<Scalar> getColumn(int colID) const;

  /** \brief Return <var>rowID</var>-th row as Point.
  */
  Point<Scalar> getRow(int rowID) const;

  /** \brief Compute determinant.
  */
  Scalar Det() const;

  /** \brief Returns transposed matrix.
  */
  LinearMap getTranspose() const;

  /** \brief In-place transpose.
      This is a mutator, i.e. it will change the data of <var>*this</var>.
  */
  void Transpose(); // mutator

  /** \brief Returns inverse matrix.
  */
  LinearMap getInverse() const;

  /** \brief In-place inverse.
      This is a mutator, i.e. it will change the data of <var>*this</var>.
  */
  void Invert(); // mutator
}; // class LinearMap

/** \relates LinearMap
    Overloaded matrix-vector multiplication.
*/
template<class Scalar>
const Point<Scalar> operator * (const LinearMap<Scalar>& mat, const Point<Scalar>& vec);

/** \relates LinearMap
    Overloaded matrix-matrix multiplication.
*/
template<class Scalar>
const LinearMap<Scalar> operator * (const LinearMap<Scalar>& lmat, const LinearMap<Scalar>& rmat);

/** \relates LinearMap
    Outputs a formated stream representing a LinearMap (matrix). For debugging purposes.
*/
template<class Scalar>
std::ostream& operator << (std::ostream& os, const LinearMap<Scalar>& matrix);
// end function declarations related to class LinearMap

}  // namespace Intrepid

// include templated definitions
#include <Intrepid_RealSpaceDef.hpp>

#endif
