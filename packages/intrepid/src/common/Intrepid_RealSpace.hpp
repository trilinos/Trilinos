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

/** \file   Intrepid_RealSpace.hpp
    \brief  Header file for classes providing basic linear algebra functionality in 1D, 2D and 3D.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_REALSPACE_HPP
#define INTREPID_REALSPACE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::Point
  \brief Implementation of a point in 1-, 2-, and 3-dimensional Euclidean space with Cartesian coordinates. 

  Provides basic vector space operations. A Point object stores the Cartesian coordinates
  and the frame kind (physical or reference) of the point. The frame kind is of the enumerated type EFrame
  defined in Intrepid::Intrepid_Types. The default frame is FRAME_PHYSICAL which means that the Point
  instance represents a point in the "physical" space. FRAME_REFERENCE indicates that the point is
  in a "reference" cell, i.e., one of the standard cells used in FEM reconstructions. 
  
  Overloaded vector space operations:
  
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
    
    /** \brief Array with coordinates of the point
    */
    Teuchos::Array<Scalar> data_;
    
    /**\brief The frame of the Point object. Indicates if the Point is in the reference or physical space
    */
    EFrame frameKind_;
    
  public:
      
    /** \brief Default destructor.
    */
    ~Point() {};
    
    /** \brief Default constructor.
    */
    Point() {
      data_.resize(0); 
      frameKind_ = FRAME_PHYSICAL;
    };
    
    /** \brief Copy constructor.
    */
    Point(const Point& right);
    
    /** \brief Create a 1-3D Point (vector) and initialize by zero.
    */
    Point(const int dim, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    /** \brief Create a 1D Point based on one Scalar coordinate.
    */
    Point(const Scalar x, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    /** \brief Create a 2D Point based on two Scalar coordinates.
    */
    Point(const Scalar x, 
          const Scalar y, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    /** \brief Create a 2D Point based on two int coordinates.

        This constructor helps to avoid undesired behavior when the user specifies
        two int coordinates. Without this constructor, the first int coordinate is interpreted
        as a pointer and the wrong ctor, based on a Scalar pointer and dimension, is invoked.
    */
    Point(const int x, 
          const int y, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    /** \brief Create a 3D Point based on three Scalar coordinates.
    */
    Point(const Scalar x, 
          const Scalar y, 
          const Scalar z, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    /** \brief Create a 1-3D Point based on a Scalar pointer and dimension.
      */
    Point(const Scalar* dataPtr, 
          const int dim, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    
    /** \brief Create a 1-3D Point based on an Array.
      */
    Point(const Teuchos::Array<Scalar>& dataArray, 
          const EFrame frameKind = FRAME_PHYSICAL);
    
    /** \brief Return dimension of the Point.
    */
    int getDim() const;
    
    /** \brief Returns the frame kind (reference or physical) of a Point
    */
    EFrame getFrameKind() const;
    
    /** \brief Returns a string with the name of the frame kind.
    */
    std::string getFrameName() const;
    
    /** \brief Returns const reference to data member containing the Point coordinates.
    */
    const Teuchos::Array<Scalar> & getCoordinates() const;
    
    /** \brief Load 1-3D Point with new coordinate data, do not change its frame kind.
    */
    void setCoordinates(const Scalar* dataPtr,
                        const int dim);
    
    /** \brief Changes the frame kind (PHYSICAL or REFERENCE) of the Point instance.
    */
    void setFrameKind(const EFrame newFrameKind);
    
    /** \brief Returns Euclidian distance from <var>*this</var> to <var>endPoint</var>.
    */
    Scalar distance(const Point& endPoint) const;
    
    /** \brief Returns the norm (1,2 or infinity) of the Point.
    */
    Scalar norm(ENorm normType) const;
    
    /** \brief   Overloaded [] operator.
    */
    const Scalar & operator [] (const int coordinateId) const;
    
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
}; // end class Point
  
/////////////////// Function declarations related to class Point:

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

////////////////// end function declarations related to class Point
  
  
/** \class Intrepid::Matrix
    \brief Implementation of a matrix in 1-, 2-, and 3-dimensional Euclidean space.

    Provides methods and overloaded operators for basic matrix algebra:
    
    \arg getElement extracts an entry
    \arg getColumn extracts a column as a Point object
    \arg getRow extracts a row as a Point object
    \arg * product of a Scalar and a Matrix
    \arg * product of a Point and a Matrix (right Matrix multiply) 
    \arg * product of a Matrix and a Point (left Matrix multiply) 
    \arg * product of two Matrix objects
    \arg etc.
*/
template<class Scalar>
class Matrix{
  private:
      Teuchos::Array<Scalar> elements_;
      int dim_;

  public:
      
      /** \brief Default destructor.
      */
      ~Matrix() {};
      
      /** \brief Default constructor.
      */
      Matrix() {elements_.resize(0);};
      
      /** \brief Copy constructor.
      */
      Matrix(const Matrix& right);
      
      /** \brief Create <var>dim</var>-by-<var>dim</var> zero matrix.
      */
      Matrix(const int dim);
      
      /** \brief Create <var>dim</var>-by-<var>dim</var> matrix, and fill it with Scalars  
                 referenced by <var>dataPtr</var>.

          It is assumed that the Scalars are ordered by row, i.e., for dim = 3, if <var>dataPtr</var> contains
          \code
          (a, b, c, d ,e, f, g, h, i),
          \endcode
          then
          \code
                   | a b c |
          Matrix = | d e f |
                   | g h i |
          \endcode
      */
      Matrix(const Scalar* dataPtr, const int dim);
      
      
      /** \brief Fills an existing <var>dim</var>-by-<var>dim</var> matrix with Scalars
                 referenced by <var>dataPtr</var>.

          It is assumed that the Scalars are ordered by row, i.e., for dim = 3, if <var>dataPtr</var> contains
          \code
          (a, b, c, d ,e, f, g, h, i),
          \endcode
          then
          \code
                   | a b c |
          Matrix = | d e f |
                   | g h i |
          \endcode
      */
      void setElements(const Scalar* dataPtr, const int dim);
      
      /** \brief Resize an existing matrix to a <var>dim</var>-by-<var>dim</var> zero matrix.
        
        \param dim   [in]  - desired matrix dimension
        */
      void resize(const int dim);
      
      
      /** \brief Return matrix dimension.
      */
      int getDim() const;
      
      /** \brief Return Scalar element in <var>rowID</var>-th row, <var>colID</var>-th column.
      */
      Scalar getElement(int rowID, int colID) const;
      
      /** \brief Return <var>colID</var>-th column as Point.
      */
      Point<Scalar> getColumn(int colID) const;
      
      /** \brief Return <var>rowID</var>-th row as Point.
      */
      Point<Scalar> getRow(int rowID) const;
            
      /** \brief Returns transposed matrix.
      */
      Matrix getTranspose() const;
      
      /** \brief Returns inverse matrix.
      */
      Matrix getInverse() const;
      
      /** \brief In-place transpose. This is a mutator, i.e. it will change the data of <var>*this</var>.
      */
      void transpose(); // mutator

      /** \brief Compute determinant.
      */
      Scalar det() const;
      
      /** \brief Returns 1-, Infinity- or Frobenius norm of the Matrix object.
      */
      Scalar norm(ENorm normType) const;
      
      
      /** \brief In-place inverse. This is a mutator, i.e. it will change the data of <var>*this</var>.
      */
      void invert(); // mutator
      
      
      /** \brief In-place multiply by scalar. This is a mutator, i.e. it will change the 
        data of <var>*this</var>.
        */
      void scalarMultiply(const Scalar& scalar); // mutator
      
      
      /** \brief Matrix-vector left multiply. An array, thought of as a column vector, is
        multiplied on the left by the Matrix.  This method allows to perform matrix-vector
        multiplication without creating a temporary Point object. It operates unchecked
        on Scalar pointers. 
        
        \param matVec       [out]         - matrix-vector product
        \param vec          [in]          - the vector argument
      */
      void multiplyLeft(Scalar*        matVec, 
                        const Scalar*  vec) const;


      /** \brief Matrix-vector left multiply. An array, thought of as a column vector, is
        multiplied on the left by the Matrix.  This method allows to perform matrix-vector
        multiplication without creating a temporary Point object. 
        
        \param matVec       [out]         - matrix-vector product
        \param vec          [in]          - the vector argument
      */
      void multiplyLeft(Teuchos::Array<Scalar>&       matVec, 
                        const Teuchos::Array<Scalar>& vec) const;
      
      
      /** \brief Returns a scalar representing product of Matrix row with the specified <var>rowId</var>
        and a vector whose components are stored contiguously in a FieldContainer, starting at
        position <var>vecBegin</var>. This method allows to perform matrix-vector multiplications
        without creating temporary Point objects and should be used when performance is of importance.
        
        \warning This method expects that the vector has the same dimension as the Matrix object
        whose row is being multiplied. 
        
        \param rowVec        [out]          -  product of Matrix row and the vector
        \param rowId          [in]          - order of the row that is being multiplied by the vector
        \param vec            [in]          - FieldContainer storing the vector
        \param vecBegin       [in]          - enumeration of the first element of the vector
      */      
      void rowMultiply(Scalar &                       rowVec,
                       const int                      rowId,
                       const FieldContainer<Scalar>&  vec,
                       const int                      vecBegin)  const;
      
      
      /** \brief   Overloaded () operator. Allows to access Matrix elements using (rowId,colId).
        */
      const Scalar& operator () (const int rowId, const int colId) const;
      
      
      /** \brief Assignment operator <var>*this = right</var>.
      */
      Matrix& operator = (const Matrix& right);
      
      
      /** \brief Addition <var>*this = *this + right</var>.
      */
      Matrix& operator += (const Matrix& right);
      
      
      /** \brief Subtraction <var>*this = *this - right</var>.
      */
      Matrix& operator -= (const Matrix& right);
      
      
      /** \brief Scalar multiplication <var>*this = left</var> * <var>(*this)</var>.
      */
      Matrix& operator *= (const Scalar left);
      
    }; // class Matrix


/////////////////// Function declarations related to class Matrix:

/** \relates Matrix
    Overloaded matrix addition.
*/
template<class Scalar>
const Matrix<Scalar> operator + (const Matrix<Scalar>& left, 
                                 const Matrix<Scalar>& right);

/** \relates Matrix
    Overloaded matrix subtraction.
*/
template<class Scalar>
const Matrix<Scalar> operator - (const Matrix<Scalar>& left, 
                                 const Matrix<Scalar>& right);
  
/** \relates Matrix
    Overloaded scalar-matrix multiplication.
*/
template<class Scalar>
const Matrix<Scalar> operator * (const Scalar sca, const Matrix<Scalar>& mat);

/** \relates Matrix
Overloaded matrix-vector multiplication. Matrix is multiplied on the right by a Point
object thought of as a column vector.
*/
template<class Scalar>
const Point<Scalar> operator * (const Matrix<Scalar>& mat, const Point<Scalar>& vec);



/** \relates Matrix
    Overloaded vector-matrix multiplication. Matrix is multiplied on the left by a Point
    object thought of as a row vector.
*/
template<class Scalar>
const Point<Scalar> operator * (const Point<Scalar>& vec, const Matrix<Scalar>& mat);
  
/** \relates Matrix
    Overloaded matrix-matrix multiplication.
*/
template<class Scalar>
const Matrix<Scalar> operator * (const Matrix<Scalar>& lmat, const Matrix<Scalar>& rmat);
  
/** \relates Matrix
    Outputs a formated stream representing a Matrix (matrix). For debugging purposes.
*/
template<class Scalar>
std::ostream& operator << (std::ostream& os, const Matrix<Scalar>& matrix);

///////////// end function declarations related to class Matrix
  
} // end namespace Intrepid

// include templated definitions
#include <Intrepid_RealSpaceDef.hpp>

#endif
