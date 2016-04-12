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

/** \file   Intrepid_Utils.hpp
    \brief  Intrepid utilities.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_UTILS_HPP__
#define __INTREPID2_UTILS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"

namespace Intrepid2 {
  
#define INTREPID2_TEST_FOR_ABORT(test, msg)                             \
  if (test) {                                                           \
    fprintf(stderr, "[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    fprintf(stderr, "            Test that evaluated to true: %s\n", #test); \
    fprintf(stderr, "            %s \n", msg);                          \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }

  // check the first error only
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#define INTREPID2_TEST_FOR_DEBUG_ABORT(test, info, msg)                 \
  if (!info && test) {                                                  \
    fprintf(stderr, "[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    fprintf(stderr, "            Test that evaluated to true: %s\n", #test); \
    fprintf(stderr, "            %s \n", msg);                          \
    info = true;                                                        \
  }
#else  
#define INTREPID2_TEST_FOR_DEBUG_ABORT(test, info, msg)                 \
  if (!info && test) {                                                  \
    fprintf(stderr, "[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    fprintf(stderr, "            Test that evaluated to true: %s\n", #test); \
    fprintf(stderr, "            %s \n", msg);                          \
    info = true ;                                                       \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }
#endif

  class Util {
  public:

    template<typename IdxType, typename DimType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void 
    unrollIndex(IdxType &i, IdxType &j, 
                const DimType dim0,
                const IterType iter) {
      j = iter/dim0;
      i = iter%dim0;
    }
    
    template<typename IdxType, typename DimType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void 
    unrollIndex(IdxType &i, IdxType &j, IdxType &k, 
                const DimType dim0,
                const DimType dim1,
                const IterType iter) {
      DimType tmp;
      unrollIndex(tmp, k, dim0*dim1, iter);
      unrollIndex(  i, j, dim0,      tmp);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T min(const T a, const T b) {
      return (a < b ? a : b);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T max(const T a, const T b) {
      return (a > b ? a : b);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T abs(const T a) {
      return (a > 0 ? a : -a);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T real(const T a) {
      return a;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T imag(const T a) {
      return 0;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T conj(const T a) {
      return a;
    }
    
    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static void swap(T &a, T &b) {
      T c(a); a = b; b = c;
    }
  };

  
  // /***************************************************************************************************
  //  **                                                                                               **
  //  **  Declarations of non-templated utility functions for order and cardinality of operators       **
  //  **                                                                                               **
  //  ***************************************************************************************************/

  // /** \brief  Returns the rank of fields in a function space of the specified type.

  //     Field rank is defined as the number of indices needed to specify function value and
  //     equals 0, 1,or 2 for scalars, vectors and tensors, respectively. The scalar field
  //     spaces in Intrepid are FUNCTION_SPACE_HGRAD and FUNCTION_SPACE_HVOL. The vector field
  //     spaces are FUNCTION_SPACE_HCURL, FUNCTION_SPACE_HDIV and FUNCTION_SPACE_VECTOR_HGRAD.
  //     FUNCTION_SPACE_TENSOR_HGRAD contains rank-2 tensors.

  //     \param  spaceType         [in]     -  function space type
  //     \return rank of the fields in the specified function space
  // */
  // KOKKOS_INLINE_FUNCTION
  // ordinal_type getFieldRank(const EFunctionSpace spaceType);

  // /** \brief  Returns rank of an operator.

  //     When an operator acts on a field of a certain rank, the result can be a field with the
  //     same or a different rank. Operator rank is defined the difference between the ranks of
  //     the output field and the input field:
  //     \verbatim
  //     Rank(OPERATOR) = Rank(OPERATOR(FIELD)) - Rank(FIELD)
  //     \endverbatim
  //     Therefore, operator rank allows us to figure out the rank of the result:
  //     \verbatim
  //     Rank(OPERATOR(FIELD)) = Rank(FIELD) + Rank(OPERATOR)
  //     \endverbatim
  //     and provides means to size properly arrays for output results. The following table
  //     summarizes operator ranks (~ denotes undefined, below slash means 3D).
  //     By default, in 1D any operator other than VALUE has rank 1, i.e., GRAD, CURL and DIV
  //     reduce to d/dx and Dk are the higher-order derivatives d^k/dx^k. Only scalar functions
  //     are allowed in 1D.

  //     \verbatim
  //     |========|======|============================|=========|==========|==========|==========|
  //     | field  | rank |  FUNCTION_SPACE_[type]     |  VALUE  | GRAD, Dk |   CURL   |    DIV   |
  //     |--------|------|----------------------------|---------|----------|----------|----------|
  //     | scalar |   0  |  HGRAD, HVOL               |    0    |     1    | 3-dim/~  |     ~    |
  //     | vector |   1  |  HCURL, HDIV, VECTOR_HGRAD |    0    |     1    | dim - 3  |    -1    |
  //     | tensor |   2  |  TENSOR_HGRAD              |    0    |     1    | dim - 3  |    -1    |
  //     |--------|------|----------------------------|---------|----------|----------|----------|
  //     |   1D   |   0  |  HGRAD, HVOL only          |    0    |     1    |     1    |     1    |
  //     |=======================================================================================|
  //     \endverbatim

  //     \param  spaceType        [in]    - function space type
  //     \param  operatorType     [in]    - the operator acting on the specified function space
  //     \param  spaceDim         [in]    - spatial dimension
  //     \return rank of the operator as defined in the table
  // */
  // KOKKOS_INLINE_FUNCTION
  // ordinal_type getOperatorRank(const EFunctionSpace spaceType,
  //                              const EOperator      operatorType,
  //                              const ordinal_type   spaceDim);

  // /** \brief  Returns order of an operator.

  //     \param  operatorType       [in]    - type of the operator whose order we want to know
  //     \return result ranges from 0 (for OPERATOR_VALUE) to 10 (OPERATOR_D10)
  // */
  // KOKKOS_INLINE_FUNCTION
  // ordinal_type getOperatorOrder(const EOperator operatorType);

  // /** \brief  Returns the ordinal of a partial derivative of order k based on the multiplicities of
  //     the partials dx, dy, and dz.

  //     By default, any implementation of Intrepid2::Basis method returns partials of order k
  //     (specified by OPERATOR_Dk) as a multiset ordered by the lexicographical order of the
  //     partial derivatives multiplicities. For example, the 10 derivatives of order 3 in 3D
  //     are enumerated as:
  //     \verbatim
  //     D3={(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),(0,3,0),(0,2,1),(0,1,2),(0,0,3)}
  //     \endverbatim
  //     The enumeration formula for this lexicographical order is
  //     <table>
  //     <tr> <td>\f$i(xMult)            = 0\f$                            </td> <td>in 1D (only 1 derivative)</td> </tr>
  //     <tr> <td>\f$i(xMult,yMult)      =yMult\f$                         </td> <td>in 2D</td> </tr>
  //     <tr> <td>\f$i(xMult,yMult,zMult)=zMult+\sum_{r = 0}^{k-xMult} r\f$</td> <td>in 3D</td> </tr>
  //     </table>
  //     where the order k of Dk is implicitly defined by xMult + yMult + zMult. Space dimension is
  //     implicitly defined by the default values of the multiplicities of y and z derivatives.

  //     \param  xMult            [in]    - multiplicity of dx
  //     \param  yMult            [in]    - multiplicity of dy (default = -1)
  //     \param  zMult            [in]    - multiplicity of dz (default = -1)
  //     \return the ordinal of partial derivative of order k as function of dx, dy, dz multiplicities
  // */
  // KOKKOS_INLINE_FUNCTION
  // ordinal_type getDkEnumeration(const ordinal_type xMult,
  //                               const ordinal_type yMult = -1,
  //                               const ordinal_type zMult = -1);

  // /** \brief  Returns multiplicities of dx, dy, and dz based on the enumeration of the partial
  //     derivative, its order and the space dimension. Inverse of the getDkEnumeration() method.

  //     \param  partialMult      [out]    - array with the multiplicities f dx, dy and dz
  //     \param  derivativeEnum   [in]     - enumeration of the partial derivative
  //     \param  operatorType     [in]     - k-th partial derivative Dk
  //     \param  spaceDim         [in]     - space dimension
  // */
  // template<typename OrdinalArrayType>
  // KOKKOS_INLINE_FUNCTION
  // void getDkMultiplicities(OrdinalArrayType   partialMult,
  //                          const ordinal_type derivativeEnum,
  //                          const EOperator    operatorType,
  //                          const ordinal_type spaceDim);

  // /** \brief  Returns cardinality of Dk, i.e., the number of all derivatives of order k.

  //     The set of all partial derivatives of order k is isomorphic to the set of all multisets
  //     of cardinality k with elements taken from the sets {x}, {x,y}, and {x,y,z} in 1D, 2D,
  //     and 3D respectively. For example, the partial derivative
  //     \f$\displaystyle D\{1,2,1\}f = \frac{d^4 f}{dx dy^2 dz}\f$  maps to the multiset
  //     \f$\{x, y, z\}\f$ with multiplicities \f$\{1,2,1\}\f$. The number of all such multisets
  //     is given by the binomial coefficient
  //     \f[       \begin{pmatrix} spaceDim + k - 1 \\ spaceDim - 1 \end{pmatrix}              \f]
  //     Therefore:
  //     \li     in 1D: cardinality = 1\n
  //     \li     in 2D: cardinality = k + 1\n
  //     \li     in 3D: cardinality = (k + 1)*(k + 2)/2

  //     \param  operatorType     [in]     - k-th derivative operator Dk
  //     \param  spaceDim         [in]     - space dimension
  //     \return the number of all partial derivatives of order k
  // */
  // KOKKOS_INLINE_FUNCTION
  // ordinal_type getDkCardinality(const EOperator    operatorType,
  //                               const ordinal_type spaceDim);


  // /***************************************************************************************************
  //  **                                                                                               **
  //  **                      Declarations of helper functions for the basis class                     **
  //  **                                                                                               **
  //  ***************************************************************************************************/

  // /** \brief  Fills <var>ordinalToTag_</var> and <var>tagToOrdinal_</var> by basis-specific tag data

  //     \param  tagToOrdinal     [out]  - Lookup table for the DoF's ordinal by its tag
  //     \param  ordinalToTag     [out]  - Lookup table for the DoF's tag by its ordinal
  //     \param  tags             [in]   - a set of basis-dependent tags in flat (rank-1) array format.
  //     \param  basisCard        [in]   - cardinality of the basis
  //     \param  tagSize          [in]   - number of fields in a DoF tag
  //     \param  posScDim         [in]   - position in the tag, counting from 0, of the subcell dim
  //     \param  posScOrd         [in]   - position in the tag, counting from 0, of the subcell ordinal
  //     \param  posDfOrd         [in]   - position in the tag, counting from 0, of DoF ordinal relative to the subcell

  //     Note:
  //     - Compiled on device (KOKKOS_INLINE_FUNCTION) 
  //     - NOT Callable in Kokkos functors (tagToOrdinal and ordinalToTag are created)
  //     - This function should be placed in Basis
  // */
  // template<class ...tagToOrdinalProperties,
  //          class ...ordinalToTagProperties,
  //          class ...tagProperties>
  // KOKKOS_INLINE_FUNCTION
  // void setOrdinalTagData(Kokkos::View<ordinal_type***,  tagToOrdinalProperties...> &tagToOrdinal,
  //                        Kokkos::View<ordinal_type*[4], ordinalToTagProperties...> &ordinalToTag,
  //                        const Kokkos::View<ordinal_type*, tagProperties...> tags,
  //                        const ordinal_type  basisCard,
  //                        const ordinal_type  tagSize,
  //                        const ordinal_type  posScDim,
  //                        const ordinal_type  posScOrd,
  //                        const ordinal_type  posDfOrd);


} // end namespace Intrepid2

#endif
