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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_RealSpaceTools.hpp
    \brief  Header file for classes providing basic linear algebra functionality in 1D, 2D and 3D.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_REALSPACETOOLS_HPP__
#define __INTREPID2_REALSPACETOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {
  
  /** \class Intrepid2::RealSpaceTools
      \brief Implementation of basic linear algebra functionality in Euclidean space.

      Note:
      - Compiled on devices (KOKKOS_INLINE_FUNCTION)
      - Callable on devices and Kokkos functor (no temporary allocation)
      - parallel_for inside of intrepid functions is dictated by the provided execution space
      - With Kokkos::Serial, functions (that already contain parallel_for) can be nested in
        Kokkos functors
      - When a function is decorated with KOKKOS_INLINE_FUNCTION, remove
        Teuchos testings and std::vectors
      - For norm and det computations, it choose the value_type from the first input container.
      - As the name says, the class assume that input containers has "real type" and no complex type.
      - Functions that are designed for small problems are not parallelized e.g., dot, norm and det.
        These functions has non-void return types and expected to be used for small problems.
        However, high level (working on cell-level) functions are parallelized.
  */
  template<typename ExecSpaceType>
  class RealSpaceTools {
  public:

    /** \brief Computes absolute value of an array.

        \param outArray   [out]  - output array
        \param inArray     [in]  - input array

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>absArray</var></b>) == rank(<b><var>inArray</var></b>)
        \li dimensions(<b><var>absArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<typename absArrayValueType, class ...absArrayProperties, 
             typename inArrayValueType,  class ...inArrayProperties>
    KOKKOS_INLINE_FUNCTION
    static void 
    absval( /**/  Kokkos::DynRankView<absArrayValueType,absArrayProperties...> absArray, 
            const Kokkos::DynRankView<inArrayValueType, inArrayProperties...>   inArray );


    /** \brief Computes, in place, absolute value of an array.

        \param inoutAbsArray  [in/out]  - input/output array
    */
    template<typename inoutArrayValueType, class ...inoutArrayProperties>
    KOKKOS_INLINE_FUNCTION
    static void 
    absval( Kokkos::DynRankView<inoutArrayValueType,inoutArrayProperties...> inoutArray );

    /** \brief Computes norm (1, 2, infinity) of a single vector stored in
        an array of rank 1.

        \param inVec     [in]  - array representing a single vector
        \param normType  [in]  - norm type

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>inVec</var></b>) == 1
    */
    template<typename inVecValueType, class ...inVecProperties>
    KOKKOS_INLINE_FUNCTION
    static inVecValueType
    vectorNorm( const Kokkos::DynRankView<inVecValueType,inVecProperties...> inVec, 
                const ENorm normType );

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
    template<typename normArrayValueType, class ...normArrayProperties,
             typename inVecValueType,     class ...inVecProperties> 
    KOKKOS_INLINE_FUNCTION
    static void 
    vectorNorm( /**/  Kokkos::DynRankView<normArrayValueType,normArrayProperties...> normArray, 
                const Kokkos::DynRankView<inVecValueType,    inVecProperties...>     inVecs, 
                const ENorm normType );
    
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
    template<typename transposeMatValueType, class ...transposeMatProperties,
             typename inMatValueType,        class ...inMatProperties> 
    KOKKOS_INLINE_FUNCTION
    static void 
    transpose( /**/  Kokkos::DynRankView<transposeMatValueType,transposeMatProperties...> transposeMats, 
               const Kokkos::DynRankView<inMatValueType,       inMatProperties...>        inMats );
    
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
    template<typename inverseMatValueType, class ...inverseMatProperties, 
             typename inMatValueType,      class ...inMatProperties>
    KOKKOS_INLINE_FUNCTION
    static void 
    inverse( /**/  Kokkos::DynRankView<inverseMatValueType,inverseMatProperties...> inverseMats, 
             const Kokkos::DynRankView<inMatValueType,     inMatProperties...>      inMats );
    
    /** \brief Computes determinant of a single square matrix stored in
        an array of rank 2.

        \param inMat  [in]  - array representing a single matrix, indexed by (D, D)

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>inMats</var></b>) == 2
        \li matrix dimension is limited to 1, 2, and 3
    */
    template<typename inMatValueType, class ...inMatProperties>
    KOKKOS_INLINE_FUNCTION
    static inMatValueType
    det( const Kokkos::DynRankView<inMatValueType,inMatProperties...> inMat );

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
    template<typename detArrayValueType, class ...detArrayProperties,
             typename inMatValueType,    class ...inMatProperties>
    KOKKOS_INLINE_FUNCTION
    static void 
    det( /**/  Kokkos::DynRankView<detArrayValueType,detArrayProperties...> detArray, 
         const Kokkos::DynRankView<inMatValueType,   inMatProperties...>    inMats );
    
    /** \brief Adds arrays <b><var>inArray1</var></b> and <b><var>inArray2</var></b>:\n
        <b><var>sumArray</var></b> = <b><var>inArray1</var></b> + <b><var>inArray2</var></b>.

        \param sumArray  [out]  - sum
        \param inArray1   [in]  - first summand
        \param inArray2   [in]  - second summand

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>sumArray</var></b>) == rank(<b><var>inArray1</var></b>) == rank(<b><var>inArray2</var></b>)
        \li dimensions(<b><var>sumArray</var></b>) == dimensions(<b><var>inArray1</var></b>) == dimensions(<b><var>inArray2</var></b>)
    */
    template<typename sumArrayValueType, class ...sumArrayProperties,
             typename inArray1ValueType, class ...inArray1Properties,
             typename inArray2ValueType, class ...inArray2Properties>
    KOKKOS_INLINE_FUNCTION
    static void
    add( /**/  Kokkos::DynRankView<sumArrayValueType,sumArrayProperties...> sumArray, 
         const Kokkos::DynRankView<inArray1ValueType,inArray1Properties...> inArray1, 
         const Kokkos::DynRankView<inArray2ValueType,inArray2Properties...> inArray2 );

    /** \brief Adds, in place, <b><var>inArray</var></b> into <b><var>inoutSumArray</var></b>:\n
        <b><var>inoutSumArray</var></b> = <b><var>inoutSumArray</var></b> + <b><var>inArray</var></b>.

        \param inoutSumArray  [in/out]  - sum/first summand
        \param inArray            [in]  - second summand

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>inoutSumArray</var></b>) == rank(<b><var>inArray</var></b>)
        \li dimensions(<b><var>inoutSumArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<typename inoutSumArrayValueType, class ...inoutSumArrayProperties,
             typename inArrayValueType,       class ...inArrayProperties>
    KOKKOS_INLINE_FUNCTION
    static void
    add( /**/  Kokkos::DynRankView<inoutSumArrayValueType,inoutSumArrayProperties...> inoutSumArray, 
         const Kokkos::DynRankView<inArrayValueType,      inArrayProperties...>       inArray );

    /** \brief Subtracts <b><var>inArray2</var></b> from <b><var>inArray1</var></b>:\n
        <b><var>diffArray</var></b> = <b><var>inArray1</var></b> - <b><var>inArray2</var></b>.

        \param diffArray  [out]  - difference
        \param inArray1    [in]  - minuend
        \param inArray2    [in]  - subtrahend

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>sumArray</var></b>) == rank(<b><var>inArray1</var></b>) == rank(<b><var>inArray2</var></b>)
        \li dimensions(<b><var>sumArray</var></b>) == dimensions(<b><var>inArray1</var></b>) == dimensions(<b><var>inArray2</var></b>)
    */
    template<typename diffArrayValueType, class ...diffArrayProperties,
             typename inArray1ValueType,  class ...inArray1Properties,
             typename inArray2ValueType,  class ...inArray2Properties>
    KOKKOS_INLINE_FUNCTION
    static void
    subtract( /**/  Kokkos::DynRankView<diffArrayValueType,diffArrayProperties...> diffArray, 
              const Kokkos::DynRankView<inArray1ValueType, inArray1Properties...>  inArray1, 
              const Kokkos::DynRankView<inArray2ValueType, inArray2Properties...>  inArray2 );

    /** \brief Subtracts, in place, <b><var>inArray</var></b> from <b><var>inoutDiffArray</var></b>:\n
        <b><var>inoutDiffArray</var></b> = <b><var>inoutDiffArray</var></b> - <b><var>inArray</var></b>.

        \param inoutDiffArray  [in/out]  - difference/minuend
        \param inArray             [in]  - subtrahend

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>inoutDiffArray</var></b>) == rank(<b><var>inArray</var></b>)
        \li dimensions(<b><var>inoutDiffArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<typename inoutDiffArrayValueType, class ...inoutDiffArrayProperties,
             typename inArrayValueType,        class ...inArrayProperties>
    KOKKOS_INLINE_FUNCTION
    static void
    subtract( /**/  Kokkos::DynRankView<inoutDiffArrayValueType,inoutDiffArrayProperties...> diffArray,
              const Kokkos::DynRankView<inArrayValueType,       inArrayProperties...>        inArray );

    /** \brief Multiplies array <b><var>inArray</var></b> by the scalar <b><var>scalar</var></b> (componentwise):\n
        <b><var>scaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inArray</var></b>.

        \param scaledArray  [out]  - scaled array
        \param inArray       [in]  - input array
        \param alpha         [in]  - multiplier

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>scaledArray</var></b>) == rank(<b><var>inArray</var></b>)
        \li dimensions(<b><var>scaledArray</var></b>) == dimensions(<b><var>inArray</var></b>)
    */
    template<typename ValueType,
             typename scaledArrayValueType, class ...scaledArrayProperties,
             typename inArrayValueType,     class ...inArrayProperties>
    KOKKOS_INLINE_FUNCTION
    static void
    scale( /**/  Kokkos::DynRankView<scaledArrayValueType,scaledArrayProperties...> scaledArray,
           const Kokkos::DynRankView<inArrayValueType,    inArrayProperties...>     inArray,
           const ValueType alpha );

    /** \brief Multiplies, in place, array <b><var>inoutScaledArray</var></b> by the scalar <b><var>scalar</var></b> (componentwise):\n
        <b><var>inoutScaledArray</var></b> = <b><var>scalar</var></b> * <b><var>inoutScaledArray</var></b>.

        \param inoutScaledArray  [in/out]  - input/output array
        \param alpha                [in]  - multiplier
    */
    template<typename ValueType,
             typename inoutScaledArrayValueType, class ...inoutScaledArrayProperties>
    KOKKOS_INLINE_FUNCTION
    static void
    scale( /**/  Kokkos::DynRankView<inoutScaledArrayValueType,inoutScaledArrayProperties...> inoutScaledArray,
           const ValueType alpha );
    
    /** \brief Computes dot product of two vectors stored in
        arrays of rank 1.

        \param inVec1    [in]  - first vector
        \param inVec2    [in]  - second vector

        \note  Requirements (checked at runtime, in debug mode): \n
        \li rank(<b><var>inVec1</var></b>) == rank(<b><var>inVec2</var></b>) == 1
        \li <b><var>inVec1</var></b> and <b><var>inVec2</var></b> have same dimension
    */
    template<typename inVec1ValueType, class ...inVec1Properties,
             typename inVec2ValueType, class ...inVec2Properties>
    KOKKOS_INLINE_FUNCTION
    static inVec1ValueType
    dot( const Kokkos::DynRankView<inVec1ValueType,inVec1Properties...> inVec1, 
         const Kokkos::DynRankView<inVec2ValueType,inVec2Properties...> inVec2 );

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
    template<typename dotArrayValueType, class ...dotArrayProperties,
             typename inVec1ValueType,   class ...inVec1Properties,
             typename inVec2ValueType,   class ...inVec2Properties>
    KOKKOS_INLINE_FUNCTION
    static void
    dot( /**/  Kokkos::DynRankView<dotArrayValueType,dotArrayProperties...> dotArray, 
         const Kokkos::DynRankView<inVec1ValueType,  inVec1Properties...>   inVecs1, 
         const Kokkos::DynRankView<inVec2ValueType,  inVec2Properties...>   inVecs2 );

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
    template<typename matVecValueType, class ...matVecProperties,
             typename inMatValueType,  class ...inMatProperties,
             typename inVecValueType,  class ...inVecProperties>
    KOKKOS_INLINE_FUNCTION
    static void
    matvec( /**/  Kokkos::DynRankView<matVecValueType,matVecProperties...> matVecs,
            const Kokkos::DynRankView<inMatValueType, inMatProperties...>  inMats,
            const Kokkos::DynRankView<inVecValueType, inVecProperties...>  inVecs );
    
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
    template<typename vecProdValueType, class ...vecProdProperties,
             typename inLeftValueType,  class ...inLeftProperties,
             typename inRightValueType, class ...inRightProperties>
    KOKKOS_INLINE_FUNCTION
    static void
    vecprod( /**/  Kokkos::DynRankView<vecProdValueType,vecProdProperties...> vecProd,
             const Kokkos::DynRankView<inLeftValueType, inLeftProperties...>  inLeft,
             const Kokkos::DynRankView<inRightValueType,inRightProperties...> inRight );
    
  }; // class RealSpaceTools

} // end namespace Intrepid2

// include templated definitions
#include "Intrepid2_RealSpaceToolsDef.hpp"

#endif
