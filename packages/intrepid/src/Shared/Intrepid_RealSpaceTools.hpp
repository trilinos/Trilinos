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

/** \file   Intrepid_RealSpaceTools.hpp
    \brief  Header file for classes providing basic linear algebra functionality in 1D, 2D and 3D.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_REALSPACETOOLS_HPP
#define INTREPID_REALSPACETOOLS_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
/** \class Intrepid::RealSpaceTools
    \brief Implementation of basic linear algebra functionality in Euclidean space.
*/
template<class Scalar>
class RealSpaceTools {
    
  public:

    /** \brief Computes absolute value of contiguous input data  <b><var>inArray</var></b>
               of size <b><var>size</var></b>.

        \param absArray   [out]  - output data
        \param inArray     [in]  - input data
        \param size        [in]  - size
    */
    static void absval(Scalar* absArray, const Scalar* inArray, const int size);


    /** \brief Computes absolute value of contiguous data  <b><var>inoutAbsArray</var></b>
               of size <b><var>size</var></b> in place.

        \param inoutAbsArray  [in/out]  - input/output data
        \param size               [in]  - size
    */
    static void absval(Scalar* inoutArray, const int size);


    /** \brief Computes absolute value of an array.

        \param outArray   [out]  - output array
        \param inArray     [in]  - input array

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>absArray</var></b>) == rank(<b><var>inArray</var></b>)
               \li dimensions(<b><var>absArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<class ArrayScalar>
    static void absval(ArrayScalar & absArray, const ArrayScalar & inArray);


    /** \brief Computes, in place, absolute value of an array.

        \param inoutAbsArray  [in/out]  - input/output array
    */
    template<class ArrayScalar>
    static void absval(ArrayScalar & inoutAbsArray);


    /** \brief Computes norm (1, 2, infinity) of the vector <b><var>inVec</var></b>
               of size <b><var>dim</var></b>.

        \param inVec       [in]  - vector
        \param dim         [in]  - vector dimension
        \param normType    [in]  - norm type
    */
    static Scalar vectorNorm(const Scalar* inVec, const int dim, const ENorm normType);


    /** \brief Computes norm (1, 2, infinity) of a single vector stored in
               an array of rank 1.

        \param inVec     [in]  - array representing a single vector
        \param normType  [in]  - norm type

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inVec</var></b>) == 1
    */
    template<class VecArray>
    static Scalar vectorNorm(const VecArray & inVec, const ENorm normType);


    /** \brief Computes norms (1, 2, infinity) of vectors stored in a
               array of total rank 2 (array of vectors), indexed by (i0, D),
               or 3 (array of arrays of vectors), indexed by (i0, i1, D).

        \param normArray  [out]  - norm array indexed by (i0) or (i0, i1)
        \param inVecs      [in]  - array of vectors indexed by (i0, D) or (i0, i1, D)
        \param normType    [in]  - norm type

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>normArray</var></b>) == rank(<b><var>inVecs</var></b>) - 1
               \li rank(<b><var>inVecs</var></b>) == 2 or 3
               \li dimensions i0, i1 of <b><var>normArray</var></b> and <b><var>inVecs</var></b> must agree
    */
    template<class NormArray, class VecArray>
    static void vectorNorm(NormArray & normArray, const VecArray & inVecs, const ENorm normType);


    /** \brief Computes transpose of the square matrix <b><var>inMat</var></b>
               of size <b><var>dim</var></b> by <b><var>dim</var></b>.

        \param transposeMat  [out] - matrix transpose
        \param inMat          [in] - matrix
        \param dim            [in] - matrix dimension
    */
    static void transpose(Scalar* transposeMat, const Scalar* inMat, const int dim);


    /** \brief Computes transposes of square matrices stored in
               an array of total rank 2 (single matrix), indexed by (D, D),
               3 (array of matrices), indexed by (i0, D, D),
               or 4 (array of arrays of matrices), indexed by (i0, i1, D, D).

        \param transposeMats  [out]  - array of transposes indexed by (D, D), (i0, D, D) or (i0, i1, D, D)
        \param inMats          [in]  - array of matrices indexed by (D, D), (i0, D, D) or (i0, i1, D, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>transposeMats</var></b>) == rank(<b><var>inMats</var></b>)
               \li rank(<b><var>inMats</var></b>) == 3 or 4
               \li dimensions(<b><var>transposeMats</var></b>) == dimensions(<b><var>inMats</var></b>)
               \li matrices must be square
    */
    template<class ArrayScalar>
    static void transpose(ArrayScalar & transposeMats, const ArrayScalar & inMats);


    /** \brief Computes inverse of the square matrix <b><var>inMat</var></b>
               of size <b><var>dim</var></b> by <b><var>dim</var></b>.

        \param inverseMat  [out] - matrix inverse
        \param inMat        [in] - matrix
        \param dim          [in] - matrix dimension
    */
    static void inverse(Scalar* inverseMat, const Scalar* inMat, const int dim);


    /** \brief Computes inverses of nonsingular matrices stored in
               an array of total rank 2 (single matrix), indexed by (D, D),
               3 (array of matrices), indexed by (i0, D, D),
               or 4 (array of arrays of matrices), indexed by (i0, i1, D, D).

        \param inverseMats  [out]  - array of inverses indexed by (D, D), (i0, D, D) or (i0, i1, D, D)
        \param inMats        [in]  - array of matrices indexed by (D, D), (i0, D, D) or (i0, i1, D, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inverseMats</var></b>) == rank(<b><var>inMats</var></b>)
               \li rank(<b><var>inMats</var></b>) == 3 or 4
               \li dimensions(<b><var>inverseMats</var></b>) == dimensions(<b><var>inMats</var></b>)
               \li matrices must be square
               \li matrix dimensions are limited to 1, 2, and 3
    */
    template<class MatArray>
    static void inverse(MatArray & inverseMats, const MatArray & inMats);


    /** \brief Computes determinant of the square matrix <b><var>inMat</var></b>
               of size <b><var>dim</var></b> by <b><var>dim</var></b>.

        \param inMat  [in]  - matrix
        \param dim    [in]  - matrix dimension
    */
    static Scalar det(const Scalar* inMat, const int dim);


    /** \brief Computes determinant of a single square matrix stored in
               an array of rank 2.

        \param inMat  [in]  - array representing a single matrix, indexed by (D, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inMats</var></b>) == 2
               \li matrix dimension is limited to 1, 2, and 3
    */
    template<class ArrayScalar>
    static Scalar det(const ArrayScalar & inMat);


    /** \brief Computes determinants of matrices stored in
               an array of total rank 3 (array of matrices),
               indexed by (i0, D, D), or 4 (array of arrays of matrices),
               indexed by (i0, i1, D, D).

        \param detArray  [out]  - array of determinants indexed by (i0) or (i0, i1)
        \param inMats     [in]  - array of matrices indexed by (i0, D, D) or (i0, i1, D, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>detArray</var></b>) == rank(<b><var>inMats</var></b>) - 2 
               \li rank(<b><var>inMats</var></b>) == 3 or 4
               \li dimensions i0, i1 of <b><var>detArray</var></b> and <b><var>inMats</var></b> must agree
               \li matrix dimensions are limited to 1, 2, and 3
    */
    template<class DetArray, class MatArray>
    static void det(DetArray & detArray, const MatArray & inMats);


    /** \brief Adds contiguous data <b><var>inArray1</var></b> and <b><var>inArray2</var></b>
               of size <b><var>size</var></b>:\n
               <b><var>sumArray</var></b> = <b><var>inArray1</var></b> + <b><var>inArray2</var></b>.

        \param sumArray  [out]  - sum
        \param inArray1   [in]  - first summand
        \param inArray2   [in]  - second summand
        \param size       [in]  - size of input/output data
    */
    static void add(Scalar* sumArray, const Scalar* inArray1, const Scalar* inArray2, const int size);


    /** \brief Adds, in place, contiguous data <b><var>inArray</var></b> into
               <b><var>inoutSumArray</var></b> of size <b><var>size</var></b>:\n
               <b><var>inoutSumArray</var></b> = <b><var>inoutSumArray</var></b> + <b><var>inArray</var></b>.

        \param inoutSumArray  [in/out]  - sum / first summand
        \param inArray            [in]  - second summand
        \param size               [in]  - size of input/output data
    */
    static void add(Scalar* inoutSumArray, const Scalar* inArray, const int size);


    /** \brief Adds arrays <b><var>inArray1</var></b> and <b><var>inArray2</var></b>:\n
               <b><var>sumArray</var></b> = <b><var>inArray1</var></b> + <b><var>inArray2</var></b>.

        \param sumArray  [out]  - sum
        \param inArray1   [in]  - first summand
        \param inArray2   [in]  - second summand

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>sumArray</var></b>) == rank(<b><var>inArray1</var></b>) == rank(<b><var>inArray2</var></b>)
               \li dimensions(<b><var>sumArray</var></b>) == dimensions(<b><var>inArray1</var></b>) == dimensions(<b><var>inArray2</var></b>)
    */
    template<class ArrayScalar>
    static void add(ArrayScalar & sumArray, const ArrayScalar & inArray1, const ArrayScalar & inArray2);


    /** \brief Adds, in place, <b><var>inArray</var></b> into <b><var>inoutSumArray</var></b>:\n
               <b><var>inoutSumArray</var></b> = <b><var>inoutSumArray</var></b> + <b><var>inArray</var></b>.

        \param inoutSumArray  [in/out]  - sum/first summand
        \param inArray            [in]  - second summand

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inoutSumArray</var></b>) == rank(<b><var>inArray</var></b>)
               \li dimensions(<b><var>inoutSumArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<class ArrayScalar>
    static void add(ArrayScalar & inoutSumArray, const ArrayScalar & inArray);


    /** \brief Subtracts contiguous data <b><var>inArray2</var></b> from <b><var>inArray1</var></b>
               of size <b><var>size</var></b>:\n
               <b><var>diffArray</var></b> = <b><var>inArray1</var></b> - <b><var>inArray2</var></b>.

        \param diffArray  [out]  - difference
        \param inArray1    [in]  - minuend
        \param inArray2    [in]  - subtrahend
        \param size        [in]  - size of input/output data
    */
    static void subtract(Scalar* diffArray, const Scalar* inArray1, const Scalar* inArray2, const int size);


    /** \brief Subtracts, in place, contiguous data <b><var>inArray</var></b> from
               <b><var>inoutDiffArray</var></b> of size <b><var>size</var></b>:\n
               <b><var>inoutDiffArray</var></b> = <b><var>inoutDiffArray</var></b> - <b><var>inArray</var></b>.

        \param inoutDiffArray  [in/out]  - difference/minuend
        \param inArray             [in]  - subtrahend
        \param size                [in]  - size of input/output data
    */
    static void subtract(Scalar* inoutDiffArray, const Scalar* inArray, const int size);


    /** \brief Subtracts <b><var>inArray2</var></b> from <b><var>inArray1</var></b>:\n
               <b><var>diffArray</var></b> = <b><var>inArray1</var></b> - <b><var>inArray2</var></b>.

        \param diffArray  [out]  - difference
        \param inArray1    [in]  - minuend
        \param inArray2    [in]  - subtrahend

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>sumArray</var></b>) == rank(<b><var>inArray1</var></b>) == rank(<b><var>inArray2</var></b>)
               \li dimensions(<b><var>sumArray</var></b>) == dimensions(<b><var>inArray1</var></b>) == dimensions(<b><var>inArray2</var></b>)
    */
    template<class ArrayScalar>
    static void subtract(ArrayScalar & diffArray, const ArrayScalar & inArray1, const ArrayScalar & inArray2);


    /** \brief Subtracts, in place, <b><var>inArray</var></b> from <b><var>inoutDiffArray</var></b>:\n
               <b><var>inoutDiffArray</var></b> = <b><var>inoutDiffArray</var></b> - <b><var>inArray</var></b>.

        \param inoutDiffArray  [in/out]  - difference/minuend
        \param inArray             [in]  - subtrahend

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inoutDiffArray</var></b>) == rank(<b><var>inArray</var></b>)
               \li dimensions(<b><var>inoutDiffArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<class ArrayScalar>
    static void subtract(ArrayScalar & inoutDiffArray, const ArrayScalar & inArray);


    /** \brief Multiplies contiguous data <b><var>inArray</var></b> of size
               <b><var>size</var></b> by a scalar (componentwise):\n
               <b><var>scaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inArray</var></b>.

        \param scaledArray  [out]  - scaled array
        \param inArray       [in]  - input array
        \param size          [in]  - size of the input array
        \param scalar        [in]  - multiplier
    */
    static void scale(Scalar* scaledArray, const Scalar* inArray, const int size, const Scalar scalar);


    /** \brief Multiplies, in place, contiguous data <b><var>inoutScaledArray</var></b> of size
               <b><var>size</var></b> by a scalar (componentwise):\n
               <b><var>inoutScaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inoutScaledArray</var></b>.

        \param inoutScaledArray  [in/out]  - input/scaled array
        \param size                  [in]  - size of array
        \param scalar                [in]  - multiplier
    */
    static void scale(Scalar* inoutScaledArray, const int size, const Scalar scalar);


    /** \brief Multiplies array <b><var>inArray</var></b> by the scalar <b><var>scalar</var></b> (componentwise):\n
               <b><var>scaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inArray</var></b>.

        \param scaledArray  [out]  - scaled array
        \param inArray       [in]  - input array
        \param scalar        [in]  - multiplier

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>scaledArray</var></b>) == rank(<b><var>inArray</var></b>)
               \li dimensions(<b><var>scaledArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<class ArrayScalar>
    static void scale(ArrayScalar & scaledArray, const ArrayScalar & inArray, const Scalar scalar);


    /** \brief Multiplies, in place, array <b><var>inoutScaledArray</var></b> by the scalar <b><var>scalar</var></b> (componentwise):\n
               <b><var>inoutScaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inoutScaledArray</var></b>.

        \param inoutScaledArray  [in/out]  - input/output array
        \param scalar                [in]  - multiplier
    */
    template<class ArrayScalar>
    static void scale(ArrayScalar & inoutScaledArray, const Scalar scalar);


    /** \brief Computes dot product of contiguous data <b><var>inArray1</var></b> and <b><var>inArray2</var></b>
               of size <b><var>size</var></b>.

        \param inArray1    [in]  - first array
        \param inArray2    [in]  - second array
        \param size        [in]  - size of input arrays
    */
    static Scalar dot(const Scalar* inArray1, const Scalar* inArray2, const int size);


    /** \brief Computes dot product of two vectors stored in
               arrays of rank 1.

        \param inVec1    [in]  - first vector
        \param inVec2    [in]  - second vector

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inVec1</var></b>) == rank(<b><var>inVec2</var></b>) == 1
               \li <b><var>inVec1</var></b> and <b><var>inVec2</var></b> have same dimension
    */
    template<class VecArray>
    static Scalar dot(const VecArray & inVec1, const VecArray & inVec2);


    /** \brief Computes dot product of vectors stored in an
               array of total rank 2 (array of vectors), indexed by (i0, D),
               or 3 (array of arrays of vectors), indexed by (i0, i1, D).

        \param dotArray  [out]  - dot product array indexed by (i0) or (i0, i1)
        \param inVecs1    [in]  - first array of vectors indexed by (i0, D) or (i0, i1, D)
        \param inVecs2    [in]  - second array of vectors indexed by (i0, D) or (i0, i1, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>dotArray</var></b>) == rank(<b><var>inVecs1</var></b>) - 1 == rank(<b><var>inVecs2</var></b>) - 1
               \li rank(<b><var>inVecs1</var></b>) == 2 or 3
               \li dimensions i0, i1 of <b><var>dotArray</var></b> and <b><var>inVecs1</var></b> / <b><var>inVecs2</var></b> must agree
    */
    template<class DotArray, class VecArray>
    static void dot(DotArray & dotArray, const VecArray & inVecs1, const VecArray & inVecs2);


    /** \brief Matrix-vector left multiply using contiguous data:\n
               <b><var>matVec</var></b> = <b><var>inMat</var></b> * <b><var>inVec</var></b>.

               A single "column" vector <b><var>inVec</var></b> of size <b><var>dim</var></b> is
               multiplied on the left by a square matrix <b><var>inMat</var></b> of size
               <b><var>dim</var></b> by <b><var>dim</var></b>.
        
        \param matVec  [out]  - matrix-vector product
        \param inMat    [in]  - the matrix argument
        \param inVec    [in]  - the vector argument
        \param dim      [in]  - matrix/vector dimension
    */
    static void matvec(Scalar* matVec, const Scalar* inMat, const Scalar* inVec, const int dim);


    /** \brief Matrix-vector left multiply using multidimensional arrays:\n
               <b><var>matVec</var></b> = <b><var>inMat</var></b> * <b><var>inVec</var></b>.

               An array (rank 1, 2 or 3) of "column" vectors, indexed by
               (D), (i0, D) or (i0, i1, D), is multiplied on the left by an
               array (rank 2, 3 or 4) of square matrices, indexed by (D, D),
               (i0, D, D) or (i0, i1, D, D).
        
        \param matVec  [out]  - matrix-vector product indexed by (D), (i0, D) or (i0, i1, D)
        \param inMat    [in]  - the matrix argument indexed by (D, D), (i0, D, D) or (i0, i1, D, D)
        \param inVec    [in]  - the vector argument indexed by (D), (i0, D) or (i0, i1, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>matVec</var></b>) == rank(<b><var>inVec</var></b>) == rank(<b><var>inMat</var></b>) - 1
               \li dimensions(<b><var>matVec</var></b>) == dimensions(<b><var>inVec</var></b>)
               \li matrix and vector dimensions D, i0 and i1 must agree
               \li matrices are square
    */
    template<class MatArray, class VecArray>
    static void matvec(VecArray & matVecs, const MatArray & inMats, const VecArray & inVecs);

    
    /* \brief Vector product using multidimensional arrays:\n
              <b><var>vecProd</var></b> = <b><var>inVecLeft</var></b> x <b><var>inVecRight</var></b>
      
               Vector multiplication of two "column" vectors stored in arrays (rank 1, 2, or 3)
               indexed by (D), (i0, D) or (i0, i1, D). 
      
       \param vecProd   [in]  - vector product indexed by (D), (i0, D) or (i0, i1, D)
       \param inLeft    [in]  - left vector argument indexed by (D), (i0, D) or (i0, i1, D)
       \param inRight   [in]  - right vector argument indexed by (D), (i0, D) or (i0, i1, D)
      
       \todo Need to decide on how to handle vecprod in 2D: is the result a vector, i.e., 
      there's dimension D or a scalar?
      */
    template<class VecArrayOut, class VecArrayIn>
    static void vecprod(VecArrayOut &       vecProd, 
                        const VecArrayIn &  inLeft, 
                        const VecArrayIn &  inRight);
    
    
    
}; // class RealSpaceTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_RealSpaceToolsDef.hpp>

#endif
