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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
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
#include "Teuchos_Assert.hpp"


#include <Intrepid_Rank.hpp>
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
    template<class ArrayAbs, class ArrayIn>
    static void absval(ArrayAbs & absArray, const ArrayIn & inArray);


    /** \brief Computes, in place, absolute value of an array.

        \param inoutAbsArray  [in/out]  - input/output array
    */
    template<class ArrayInOut>
    static void absval(ArrayInOut & inoutAbsArray);


    /** \brief Computes norm (1, 2, infinity) of the vector <b><var>inVec</var></b>
               of size <b><var>dim</var></b>.

        \param inVec       [in]  - vector
        \param dim         [in]  - vector dimension
        \param normType    [in]  - norm type
    */
    static Scalar vectorNorm(const Scalar* inVec, const size_t dim, const ENorm normType);


    /** \brief Computes norm (1, 2, infinity) of a single vector stored in
               an array of rank 1.

        \param inVec     [in]  - array representing a single vector
        \param normType  [in]  - norm type

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inVec</var></b>) == 1
    */
    template<class ArrayIn>
    static Scalar vectorNorm(const ArrayIn & inVec, const ENorm normType);


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
    template<class ArrayNorm, class ArrayIn>
    static void vectorNorm(ArrayNorm & normArray, const ArrayIn & inVecs, const ENorm normType);

  /*  template<class ArrayNorm, class ArrayIn>
    static void vectorNormTemp(ArrayNorm & normArray, const ArrayIn & inVecs, const ENorm normType);
*/
    /** \brief Computes transpose of the square matrix <b><var>inMat</var></b>
               of size <b><var>dim</var></b> by <b><var>dim</var></b>.

        \param transposeMat  [out] - matrix transpose
        \param inMat          [in] - matrix
        \param dim            [in] - matrix dimension
    */
    static void transpose(Scalar* transposeMat, const Scalar* inMat, const size_t dim);
    
  /*  template<class ArrayTranspose, class ArrayIn>
	static void transpose(ArrayTranspose transposeMat, const ArrayIn inMat, const size_t dim);*/
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
    template<class ArrayTranspose, class ArrayIn>
    static void transpose(ArrayTranspose & transposeMats, const ArrayIn & inMats);

   /* template<class ArrayTranspose, class ArrayIn>
    static void transposeTemp(ArrayTranspose & transposeMats, const ArrayIn & inMats);*/
    /** \brief Computes inverse of the square matrix <b><var>inMat</var></b>
               of size <b><var>dim</var></b> by <b><var>dim</var></b>.

        \param inverseMat  [out] - matrix inverse
        \param inMat        [in] - matrix
        \param dim          [in] - matrix dimension
    */
    static void inverse(Scalar* inverseMat, const Scalar* inMat, const size_t dim);


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
    template<class ArrayInverse, class ArrayIn>
    static void inverse(ArrayInverse & inverseMats, const ArrayIn & inMats);
    

    static Scalar det(const Scalar* inMat, const size_t dim);


    /** \brief Computes determinant of a single square matrix stored in
               an array of rank 2.

        \param inMat  [in]  - array representing a single matrix, indexed by (D, D)

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inMats</var></b>) == 2
               \li matrix dimension is limited to 1, 2, and 3
    */
    template<class ArrayIn>
    static Scalar det(const ArrayIn & inMat);


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
    template<class ArrayDet, class ArrayIn>
    static void det(ArrayDet & detArray, const ArrayIn & inMats);
    
    template<class ArrayDet, class ArrayIn, int matRank>
	struct detTempSpec;
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
    template<class ArraySum, class ArrayIn1, class ArrayIn2>
    static void add(ArraySum & sumArray, const ArrayIn1 & inArray1, const ArrayIn2 & inArray2);


    /** \brief Adds, in place, <b><var>inArray</var></b> into <b><var>inoutSumArray</var></b>:\n
               <b><var>inoutSumArray</var></b> = <b><var>inoutSumArray</var></b> + <b><var>inArray</var></b>.

        \param inoutSumArray  [in/out]  - sum/first summand
        \param inArray            [in]  - second summand

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inoutSumArray</var></b>) == rank(<b><var>inArray</var></b>)
               \li dimensions(<b><var>inoutSumArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<class ArraySum, class ArrayIn>
    static void add(ArraySum & inoutSumArray, const ArrayIn & inArray);


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
    template<class ArrayDiff, class ArrayIn1, class ArrayIn2>
    static void subtract(ArrayDiff & diffArray, const ArrayIn1 & inArray1, const ArrayIn2 & inArray2);


    /** \brief Subtracts, in place, <b><var>inArray</var></b> from <b><var>inoutDiffArray</var></b>:\n
               <b><var>inoutDiffArray</var></b> = <b><var>inoutDiffArray</var></b> - <b><var>inArray</var></b>.

        \param inoutDiffArray  [in/out]  - difference/minuend
        \param inArray             [in]  - subtrahend

        \note  Requirements (checked at runtime, in debug mode): \n
               \li rank(<b><var>inoutDiffArray</var></b>) == rank(<b><var>inArray</var></b>)
               \li dimensions(<b><var>inoutDiffArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<class ArrayDiff, class ArrayIn>
    static void subtract(ArrayDiff & inoutDiffArray, const ArrayIn & inArray);

    template<class ArrayDiff, class ArrayIn>
    static void subtractTemp(ArrayDiff & inoutDiffArray, const ArrayIn & inArray);
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
    template<class ArrayScaled, class ArrayIn>
    static void scale(ArrayScaled & scaledArray, const ArrayIn & inArray, const Scalar scalar);


    /** \brief Multiplies, in place, array <b><var>inoutScaledArray</var></b> by the scalar <b><var>scalar</var></b> (componentwise):\n
               <b><var>inoutScaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inoutScaledArray</var></b>.

        \param inoutScaledArray  [in/out]  - input/output array
        \param scalar                [in]  - multiplier
    */
    template<class ArrayScaled>
    static void scale(ArrayScaled & inoutScaledArray, const Scalar scalar);


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
    template<class ArrayVec1, class ArrayVec2>
    static Scalar dot(const ArrayVec1 & inVec1, const ArrayVec2 & inVec2);


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
    template<class ArrayDot, class ArrayVec1, class ArrayVec2>
    static void dot(ArrayDot & dotArray, const ArrayVec1 & inVecs1, const ArrayVec2 & inVecs2);


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
    static void matvec(Scalar* matVec, const Scalar* inMat, const Scalar* inVec, const size_t dim);


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
    template<class ArrayMatVec, class ArrayMat, class ArrayVec>
    static void matvec(ArrayMatVec & matVecs, const ArrayMat & inMats, const ArrayVec & inVecs);

    
    /** \brief Vector product using multidimensional arrays:\n
              <b><var>vecProd</var></b> = <b><var>inVecLeft</var></b> x <b><var>inVecRight</var></b>
      
               Vector multiplication of two "column" vectors stored in arrays (rank 1, 2, or 3)
               indexed by (D), (i0, D) or (i0, i1, D). 
      
        \param vecProd   [in]  - vector product indexed by (D), (i0, D) or (i0, i1, D)
        \param inLeft    [in]  - left vector argument indexed by (D), (i0, D) or (i0, i1, D)
        \param inRight   [in]  - right vector argument indexed by (D), (i0, D) or (i0, i1, D)
      
        \todo Need to decide on how to handle vecprod in 2D: is the result a vector, i.e., 
      there's dimension D or a scalar?
      */
    template<class ArrayVecProd, class ArrayIn1, class ArrayIn2>
    static void vecprod(ArrayVecProd & vecProd, const ArrayIn1 & inLeft, const ArrayIn2 & inRight);
    
   /* template<class ArrayVecProd, class ArrayIn1, class ArrayIn2>
    static void vecprodTemp(ArrayVecProd & vecProd, const ArrayIn1 & inLeft, const ArrayIn2 & inRight); */   
    
}; // class RealSpaceTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_RealSpaceToolsDef.hpp>

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

