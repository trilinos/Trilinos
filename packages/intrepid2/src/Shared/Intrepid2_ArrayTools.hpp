// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ArrayTools.hpp
    \brief  Header file for Intrepid2::ArrayTools class providing utilities for array operations
    \author Created by P. Bochev and D. Ridzal,
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_HPP__
#define __INTREPID2_ARRAYTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_Kernels.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::ArrayTools
      \brief Utility class that provides methods for higher-order algebraic
      manipulation of user-defined arrays, such as tensor contractions.
      For low-order operations, see Intrepid2::RealSpaceTools.

      Note:
      - Compiled on devices (KOKKOS_INLINE_FUNCTION)
      - Callable on devices and Kokkos functor (no temporary allocation)
      - parallel_for inside of intrepid functions is dictated by the provided execution space
      - With Kokkos::Serial, functions (that already contain parallel_for) can be nested in
        Kokkos functors
      - When a function is decorated with KOKKOS_INLINE_FUNCTION, remove
        Teuchos testings and std::vectors
  */

  template<typename DeviceType>
  class ArrayTools {
    using ExecSpaceType = typename DeviceType::execution_space;
  public:

    /** \brief Contracts the "point" dimension P of two rank-3 containers with
               dimensions (C,L,P) and (C,R,P), and returns the result in a 
               rank-3 container with dimensions (C,L,R).

        \code
          C - num. integration domains       dim0 in both input containers
          L - num. "left" fields             dim1 in "left" container
          R - num. "right" fields            dim1 in "right" container
          P - num. integration points        dim2 in both input containers
        \endcode

        \param  outputFields   [out] - Output (product) fields array.
        \param  leftFields      [in] - Left input array.
        \param  rightFields     [in] - Right input array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename leftFieldValueType,   class ...leftFieldProperties,
             typename rightFieldValueType,  class ...rightFieldProperties>
    static void
    contractFieldFieldScalar(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                              const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                              const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                              const bool sumInto = false );

    /** \brief Contracts the "point" and "space" dimensions P and D1 of two rank-4
        containers with dimensions (C,L,P,D1) and (C,R,P,D1), and returns the
        result in a rank-3 container with dimensions (C,L,R).

        For a fixed index "C", (C,L,R) represents a rectangular L X R matrix
        where L and R may be different.
        \code
        C - num. integration domains       dim0 in both input containers
        L - num. "left" fields             dim1 in "left" container
        R - num. "right" fields            dim1 in "right" container
        P - num. integration points        dim2 in both input containers
        D1- vector dimension               dim3 in both input containers
        \endcode

        \param  outputFields   [out] - Output array.
        \param  leftFields      [in] - Left input array.
        \param  rightFields     [in] - Right input array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename leftFieldValueType,   class ...leftFieldProperties,
             typename rightFieldValueType,  class ...rightFieldProperties>
    static void
    contractFieldFieldVector(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                              const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                              const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                              const bool sumInto = false );

    /** \brief Contracts the "point" and "space" dimensions P, D1, and D2 of two rank-5
        containers with dimensions (C,L,P,D1,D2) and (C,R,P,D1,D2), and returns
        the result in a rank-3 container with dimensions (C,L,R).

        For a fixed index "C", (C,L,R) represents a rectangular L X R matrix
        where L and R may be different.
        \code
        C - num. integration domains       dim0 in both input containers
        L - num. "left" fields             dim1 in "left" container
        R - num. "right" fields            dim1 in "right" container
        P - num. integration points        dim2 in both input containers
        D1- vector dimension               dim3 in both input containers
        D2- 2nd tensor dimension           dim4 in both input containers
        \endcode

        \param  outputFields   [out] - Output array.
        \param  leftFields      [in] - Left input array.
        \param  rightFields     [in] - Right input array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename leftFieldValueType,   class ...leftFieldProperties,
             typename rightFieldValueType,  class ...rightFieldProperties>
    static void
    contractFieldFieldTensor(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                              const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                              const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                              const bool sumInto = false );

    /** \brief Contracts the "point" dimensions P of a rank-3 containers and
        a rank-2 container with dimensions (C,F,P) and (C,P), respectively,
        and returns the result in a rank-2 container with dimensions (C,F).

        For a fixed index "C", (C,F) represents a (column) vector of length F.
        \code
        C - num. integration domains       dim0 in both input containers
        F - num. fields                    dim1 in fields input container
        P - num. integration points        dim2 in fields input container and dim1 in scalar data container
        \endcode

        \param  outputFields   [out] - Output fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    contractDataFieldScalar(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>  outputFields,
                             const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>    inputData,
                             const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>   inputFields,
                             const bool sumInto = false );

    /** \brief Contracts the "point" and "space" dimensions P and D of a rank-4 container and
        a rank-3 container with dimensions (C,F,P,D) and (C,P,D), respectively,
        and returns the result in a rank-2 container with dimensions (C,F).

        For a fixed index "C", (C,F) represents a (column) vector of length F.
        \code
        C - num. integration domains                dim0 in both input containers
        F - num. fields                             dim1 in fields input container
        P - num. integration points                 dim2 in fields input container and dim1 in vector data container
        D - spatial (vector) dimension index        dim3 in fields input container and dim2 in vector data container
        \endcode

        \param  outputFields   [out] - Output fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    contractDataFieldVector(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                             const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                             const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                             const bool sumInto = false );

    /** \brief Contracts the "point" and "space" dimensions P, D1 and D2 of a rank-5 container and
        a rank-4 container with dimensions (C,F,P,D1,D2) and (C,P,D1,D2), respectively,
        and returns the result in a rank-2 container with dimensions (C,F).

        For a fixed index "C", (C,F) represents a (column) vector of length F.
        \code
        C  - num. integration domains                       dim0 in both input containers
        F  - num. fields                                    dim1 in fields input container
        P  - num. integration points                        dim2 in fields input container and dim1 in tensor data container
        D1 - first spatial (tensor) dimension index         dim3 in fields input container and dim2 in tensor data container
        D2 - second spatial (tensor) dimension index        dim4 in fields input container and dim3 in tensor data container
        \endcode

        \param  outputFields   [out] - Output fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    contractDataFieldTensor(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                             const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                             const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                             const bool sumInto = false);

    /** \brief Contracts the "point" dimensions P of rank-2 containers
        with dimensions (C,P), and returns the result in a rank-1 container
        with dimensions (C).

        \code
        C - num. integration domains       dim0 in both input containers
        P - num. integration points        dim1 in both input containers
        \endcode

        \param  outputData     [out] - Output data array.
        \param  inputDataLeft   [in] - Left data input array.
        \param  inputDataRight  [in] - Right data input array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    contractDataDataScalar(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                            const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                            const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                            const bool sumInto = false );

    /** \brief Contracts the "point" and "space" dimensions P and D of rank-3 containers
        with dimensions (C,P,D) and returns the result in a rank-1 container with dimensions (C).

        \code
        C - num. integration domains                dim0 in both input containers
        P - num. integration points                 dim1 in both input containers
        D - spatial (vector) dimension index        dim2 in both input containers
        \endcode

        \param  outputData     [out] - Output data array.
        \param  inputDataLeft   [in] - Left data input array.
        \param  inputDataRight  [in] - Right data input array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    contractDataDataVector(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                            const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                            const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                            const bool sumInto = false );

    /** \brief Contracts the "point" and "space" dimensions P, D1 and D2 of rank-4 containers
        with dimensions (C,P,D1,D2) and returns the result in a rank-1 container with dimensions (C).

        \code
        C - num. integration domains                     dim0 in both input containers
        P - num. integration points                      dim1 in both input containers
        D1 - first spatial (tensor) dimension index      dim2 in both input containers
        D2 - second spatial (tensor) dimension index     dim3 in both input containers
        \endcode

        \param  outputData     [out] - Output data array.
        \param  inputDataLeft   [in] - Left data input array.
        \param  inputDataRight  [in] - Right data input array.
        \param  sumInto         [in] - If TRUE, sum into given output array,
        otherwise overwrite it. Default: FALSE.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    contractDataDataTensor(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                            const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                            const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                            const bool sumInto = false );

    /** \brief There are two use cases:
        (1) multiplies a rank-3, 4, or 5 container \a <b>inputFields</b> with dimensions (C,F,P),
        (C,F,P,D1) or (C,F,P,D1,D2), representing the values of a set of scalar, vector
        or tensor fields, by the values in a rank-2 container \a <b>inputData</b> indexed by (C,P),
        representing the values of scalar data, OR
        (2) multiplies a rank-2, 3, or 4 container \a <b>inputFields</b> with dimensions (F,P),
        (F,P,D1) or (F,P,D1,D2), representing the values of a scalar, vector or a
        tensor field, by the values in a rank-2 container \a <b>inputData</b> indexed by (C,P),
        representing the values of scalar data;
        the output value container \a <b>outputFields</b> is indexed by (C,F,P), (C,F,P,D1)
        or (C,F,P,D1,D2), regardless of which of the two use cases is considered.

        \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \note   The argument <var><b>inputFields</b></var> can be changed!
        This enables in-place multiplication.

        \param  outputFields   [out] - Output (product) fields array.
        \param  inputData       [in] - Data (multiplying) array.
        \param  inputFields     [in] - Input (being multiplied) fields array.
        \param  reciprocal      [in] - If TRUE, <b>divides</b> input fields by the data
        (instead of multiplying). Default: FALSE.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    scalarMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                             const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                             const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                             const bool reciprocal = false );

    /** \brief There are two use cases:
        (1) multiplies a rank-2, 3, or 4 container \a <b>inputDataRight</b> with dimensions (C,P),
        (C,P,D1) or (C,P,D1,D2), representing the values of a set of scalar, vector
        or tensor data, by the values in a rank-2 container \a <b>inputDataLeft</b> indexed by (C,P),
        representing the values of scalar data, OR
        (2) multiplies a rank-1, 2, or 3 container \a <b>inputDataRight</b> with dimensions (P),
        (P,D1) or (P,D1,D2), representing the values of scalar, vector or
        tensor data, by the values in a rank-2 container \a <b>inputDataLeft</b> indexed by (C,P),
        representing the values of scalar data;
        the output value container \a <b>outputData</b> is indexed by (C,P), (C,P,D1) or (C,P,D1,D2),
        regardless of which of the two use cases is considered.

        \code
        C  - num. integration domains
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \note   The arguments <var><b>inputDataLeft</b></var>, <var><b>inputDataRight</b></var> can be changed!
        This enables in-place multiplication.

        \param  outputData      [out] - Output data array.
        \param  inputDataLeft    [in] - Left (multiplying) data array.
        \param  inputDataRight   [in] - Right (being multiplied) data array.
        \param  reciprocal       [in] - If TRUE, <b>divides</b> input fields by the data
        (instead of multiplying). Default: FALSE.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    scalarMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                            const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                            const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                            const bool reciprocal = false );

    /** \brief There are two use cases:
        (1) dot product of a rank-3, 4 or 5 container \a <b>inputFields</b> with dimensions (C,F,P)
        (C,F,P,D1) or (C,F,P,D1,D2), representing the values of a set of scalar, vector
        or tensor fields, by the values in a rank-2, 3 or 4 container \a <b>inputData</b> indexed by
        (C,P), (C,P,D1), or (C,P,D1,D2) representing the values of scalar, vector or
        tensor data, OR
        (2) dot product of a rank-2, 3 or 4 container \a <b>inputFields</b> with dimensions (F,P),
        (F,P,D1) or (F,P,D1,D2), representing the values of a scalar, vector or tensor
        field, by the values in a rank-2 container \a <b>inputData</b> indexed by (C,P), (C,P,D1) or
        (C,P,D1,D2), representing the values of scalar, vector or tensor data;
        the output value container \a <b>outputFields</b> is indexed by (C,F,P),
        regardless of which of the two use cases is considered.

        For input fields containers without a dimension index, this operation reduces to
        scalar multiplication.
        \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode>

        \param  outputFields    [out] - Output (dot product) field array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Field array.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    dotMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                          const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields );

    /** \brief There are two use cases:
        (1) dot product of a rank-2, 3 or 4 container \a <b>inputDataRight</b> with dimensions (C,P)
        (C,P,D1) or (C,P,D1,D2), representing the values of a scalar, vector or a
        tensor set of data, by the values in a rank-2, 3 or 4 container \a <b>inputDataLeft</b> indexed by
        (C,P), (C,P,D1), or (C,P,D1,D2) representing the values of scalar, vector or
        tensor data, OR
        (2) dot product of a rank-2, 3 or 4 container \a <b>inputDataRight</b> with dimensions (P),
        (P,D1) or (P,D1,D2), representing the values of scalar, vector or tensor
        data, by the values in a rank-2 container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D1) or
        (C,P,D1,D2), representing the values of scalar, vector, or tensor data;
        the output value container \a <b>outputData</b> is indexed by (C,P),
        regardless of which of the two use cases is considered.

        For input fields containers without a dimension index, this operation reduces to
        scalar multiplication.
        \code
        C  - num. integration domains
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \param  outputData      [out] - Output (dot product) data array.
        \param  inputDataLeft    [in] - Left input data array.
        \param  inputDataRight   [in] - Right input data array.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    dotMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                         const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight );


    /** \brief There are two use cases:
        (1) cross product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
        representing the values of a set of vector fields, on the left by the values in a rank-3
        container \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
        (2) cross product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
        representing the values of a vector field, on the left by the values in a rank-3 container
        \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data;
        the output value container \a <b>outputFields</b> is indexed by (C,F,P,D) in 3D (vector output)
        and by (C,F,P) in 2D (scalar output), regardless of which of the two use cases is considered.

        \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D  - spatial dimension of vector data and vector fields
        \endcode

        \param  outputFields   [out] - Output (cross product) fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    crossProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields );

    /** \brief There are two use cases:
        (1) cross product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
        representing the values of a set of vector data, on the left by the values in a rank-3
        container \a <b>inputDataLeft</b> indexed by (C,P,D) representing the values of vector data, OR
        (2) cross product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
        representing the values of vector data, on the left by the values in a rank-3 container
        \a <b>inputDataLeft</b> indexed by (C,P,D), representing the values of vector data;
        the output value container \a <b>outputData</b> is indexed by (C,P,D) in 3D (vector output) and by
        (C,P) in 2D (scalar output), regardless of which of the two use cases is considered.

        \code
        C  - num. integration domains
        P  - num. integration points
        D  - spatial dimension of vector data and vector fields
        \endcode

        \param  outputData      [out] - Output (cross product) data array.
        \param  inputDataLeft    [in] - Left input data array.
        \param  inputDataRight   [in] - Right input data array.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    crossProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight );

    /** \brief There are two use cases:
        (1) outer product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
        representing the values of a set of vector fields, on the left by the values in a rank-3
        container \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
        (2) outer product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
        representing the values of a vector field, on the left by the values in a rank-3 container
        \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data;
        the output value container \a <b>outputFields</b> is indexed by (C,F,P,D,D),
        regardless of which of the two use cases is considered.

        \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \param  outputFields   [out] - Output (outer product) fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    outerProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields );

    /** \brief There are two use cases:
        (1) outer product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
        representing the values of a set of vector data, on the left by the values in a rank-3
        container \a <b>inputDataLeft</b> indexed by (C,P,D) representing the values of vector data, OR
        (2) outer product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
        representing the values of vector data, on the left by the values in a rank-3 container
        \a <b>inputDataLeft</b> indexed by (C,P,D), representing the values of vector data;
        the output value container \a <b>outputData</b> is indexed by (C,P,D,D),
        regardless of which of the two use cases is considered.

        \code
        C  - num. integration domains
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \param  outputData      [out] - Output (outer product) data array.
        \param  inputDataLeft    [in] - Left input data array.
        \param  inputDataRight   [in] - Right input data array.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValuetype,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    outerProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValuetype, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight );

    /** \brief There are two use cases:
        (1) matrix-vector product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
        representing the values of a set of vector fields, on the left by the values in a rank-2, 3, or 4
        container \a <b>inputData</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively,
        representing the values of tensor data, OR
        (2) matrix-vector product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
        representing the values of a vector field, on the left by the values in a rank-2, 3, or 4
        container \a <b>inputData</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively,
        representing the values of tensor data; the output value container \a <b>outputFields</b> is
        indexed by (C,F,P,D), regardless of which of the two use cases is considered.

        \remarks
        The rank of <b>inputData</b> implicitly defines the type of tensor data:
        \li rank = 2 corresponds to a constant diagonal tensor \f$ diag(a,\ldots,a) \f$
        \li rank = 3 corresponds to a nonconstant diagonal tensor \f$ diag(a_1,\ldots,a_d) \f$
        \li rank = 4 corresponds to a full tensor \f$ \{a_{ij}\}\f$

        \note  It is assumed that all tensors are square!

        \note  The method is defined for spatial dimensions D = 1, 2, 3

        \code
        C    - num. integration domains
        F    - num. fields
        P    - num. integration points
        D    - spatial dimension
        D1*  - first tensor dimensions, equals the spatial dimension D
        D2** - second tensor dimension, equals the spatial dimension D
        \endcode

        \param  outputFields   [out] - Output (matrix-vector product) fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
        \param  transpose       [in] - If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    matvecProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                            const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                            const char transpose = 'N');

    /** \brief There are two use cases:
               (1) matrix-vector product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
               representing the values of a set of vector data, on the left by the values in a rank-2, 3, or 4 
               container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively, 
               representing the values of tensor data, OR
               (2) matrix-vector product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
               representing the values of vector data, on the left by the values in a rank-2, 3, or 4 
               container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively, 
               representing the values of tensor data; the output value container \a <b>outputData</b> 
               is indexed by (C,P,D), regardless of which of the two use cases is considered.
      
        \remarks 
              The rank of <b>inputDataLeft</b> implicitly defines the type of tensor data:
              \li rank = 2 corresponds to a constant diagonal tensor \f$ diag(a,\ldots,a) \f$
              \li rank = 3 corresponds to a nonconstant diagonal tensor \f$ diag(a_1,\ldots,a_d) \f$
              \li rank = 4 corresponds to a full tensor \f$ \{a_{ij}\}\f$  
      
        \note   It is assumed that all tensors are square!
      
        \code
          C    - num. integration domains
          P    - num. integration points
          D    - spatial dimension
          D1*  - first tensor dimensions, equals the spatial dimension D
          D2** - second tensor dimension, equals the spatial dimension D
        \endcode

        \param  outputData      [out] - Output (matrix-vector product) data array.
        \param  inputDataLeft    [in] - Left input data array.
        \param  inputDataRight   [in] - Right input data array.
        \param  transpose        [in] - If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    matvecProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>    outputData,
                           const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...> inputDataLeft,
                           const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...>   inputDataRight,
                           const char transpose = 'N');

    /** \brief There are two use cases:
               (1) matrix-matrix product of a rank-5 container \a <b>inputFields</b> with dimensions (C,F,P,D1,D2),
               representing the values of a set of tensor fields, on the left by the values in a rank-2, 3, or 4
               container \a <b>inputData</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively,
               representing the values of tensor data, OR
               (2) matrix-matrix product of a rank-4 container \a <b>inputFields</b> with dimensions (F,P,D1,D2),
               representing the values of a tensor field, on the left by the values in a rank-2, 3, or 4
               container \a <b>inputData</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively,
               representing the values of tensor data; the output value container \a <b>outputFields</b> is
               indexed by (C,F,P,D1,D2), regardless of which of the two use cases is considered.

        \remarks
               The rank of <b>inputData</b> implicitly defines the type of tensor data:
               \li rank = 2 corresponds to a constant diagonal tensor \f$ diag(a,\ldots,a) \f$
               \li rank = 3 corresponds to a nonconstant diagonal tensor \f$ diag(a_1,\ldots,a_d) \f$
               \li rank = 4 corresponds to a full tensor \f$ \{a_{ij}\}\f$

        \note  It is assumed that all tensors are square!

        \note  The method is defined for spatial dimensions D = 1, 2, 3

        \code
          C    - num. integration domains
          F    - num. fields
          P    - num. integration points
          D1*  - first spatial (tensor) dimension index
          D2** - second spatial (tensor) dimension index
        \endcode

        \param  outputFields   [out] - Output (matrix-matrix product) fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
        \param  transpose       [in] - If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputDataValueType,   class ...inputDataProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    matmatProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                            const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                            const char transpose = 'N' );
    
    /** \brief There are two use cases:
        (1) matrix-matrix product of a rank-4 container \a <b>inputDataRight</b> with dimensions (C,P,D1,D2),
        representing the values of a set of tensor data, on the left by the values in a rank-2, 3, or 4
        container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively,
        representing the values of tensor data, OR
        (2) matrix-matrix product of a rank-3 container \a <b>inputDataRight</b> with dimensions (P,D1,D2),
        representing the values of tensor data, on the left by the values in a rank-2, 3, or 4
        container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D1) or (C,P,D1,D2), respectively,
        representing the values of tensor data; the output value container \a <b>outputData</b>
        is indexed by (C,P,D1,D2), regardless of which of the two use cases is considered.

        \remarks
        The rank of <b>inputData</b> implicitly defines the type of tensor data:
        \li rank = 2 corresponds to a constant diagonal tensor \f$ diag(a,\ldots,a) \f$
        \li rank = 3 corresponds to a nonconstant diagonal tensor \f$ diag(a_1,\ldots,a_d) \f$
        \li rank = 4 corresponds to a full tensor \f$ \{a_{ij}\}\f$

        \note  It is assumed that all tensors are square!

        \note  The method is defined for spatial dimensions D = 1, 2, 3

        \code
        C    - num. integration domains
        P    - num. integration points
        D1*  - first spatial (tensor) dimension index
        D2** - second spatial (tensor) dimension index
        \endcode

        \param  outputData      [out] - Output (matrix-vector product) data array.
        \param  inputDataLeft    [in] - Left input data array.
        \param  inputDataRight   [in] - Right input data array.
        \param  transpose        [in] - If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'.
    */
    template<typename outputDataValueType,     class ...outputDataProperties,
             typename inputDataLeftValueType,  class ...inputDataLeftProperties,
             typename inputDataRightValueType, class ...inputDataRightProperties>
    static void
    matmatProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                           const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                           const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                           const char transpose = 'N' );

    /** \brief Replicates a rank-2, 3, or 4 container with dimensions (F,P),
        (F,P,D1) or (F,P,D1,D2), representing the values of a scalar, vector or a
        tensor field, into an output value container of size (C,F,P),
        (C,F,P,D1) or (C,F,P,D1,D2).

        \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \param  outputFields   [out] - Output fields array.
        \param  inputFields     [in] - Input fields array.
    */
    template<typename outputFieldValueType, class ...outputFieldProperties,
             typename inputFieldValueType,  class ...inputFieldProperties>
    static void
    cloneFields(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                 const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields );

    /** \brief Replicates a rank-2, 3, or 4 container with dimensions (F,P),
        (F,P,D1) or (F,P,D1,D2), representing the values of a scalar, vector or a
        tensor field, into an output value container of size (C,F,P),
        (C,F,P,D1) or (C,F,P,D1,D2).

        \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D1 - first spatial (tensor) dimension index
        D2 - second spatial (tensor) dimension index
        \endcode

        \param  outputData   [out] - Output data array.
        \param  inputData     [in] - Input data array.
    */
    template<typename outputDataValueType, class ...outputDataProperties,
             typename inputDataValueType,  class ...inputDataProperties>
    static void
    cloneData(       Kokkos::DynRankView<outputDataValueType,outputDataProperties...> outputData,
               const Kokkos::DynRankView<inputDataValueType, inputDataProperties...> inputData );

    // =====================================================================================
    // Internal universal implementations
    //
    //
  private:

    class Internal {
    public:

      template<typename outputFieldValueType, class ...outputFieldProperties,
               typename leftFieldValueType,   class ...leftFieldProperties,
               typename rightFieldValueType,  class ...rightFieldProperties>
      static void
      contractFieldField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                          const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                          const bool sumInto );
      
      template<typename outputFieldValueType, class ...outputFieldProperties,
               typename inputDataValueType,   class ...inputDataProperties,
               typename inputFieldValuetype,  class ...inputFieldProperties>
      static void
      contractDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>      outputFields,
                         const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>        inputData,
                         const Kokkos::DynRankView<inputFieldValuetype, inputFieldProperties...> inputFields,
                         const bool sumInto );

      template<typename outputDataValueType,     class ...outputDataProperties,
               typename inputDataLeftValueType,  class ...inputDataLeftProperties,
               typename inputDataRightValueType, class ...inputDataRightProperties>
      static void
      contractDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>          outputData,
                        const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                        const bool sumInto );
      
      template<typename outputValueType,     class ...outputProperties,
               typename leftInputValueType,  class ...leftInputProperties,
               typename rightInputValueType, class ...rightInputProperties>
      static void
      dotMultiply(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                   const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                   const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                   const bool hasField );
      
      template<typename outputValueType,     class ...outputProperties,
               typename leftInputValueType,  class ...leftInputProperties,
               typename rightInputValueType, class ...rightInputProperties>
      static void
      crossProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                    const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                    const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                    const bool hasField );
      
      template<typename outputValueType,     class ...outputProperties,
               typename leftInputValueType,  class ...leftInputProperties,
               typename rightInputValueType, class ...rightInputProperties>
      static void
      outerProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                    const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                    const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                    const bool hasField );
      
      template<typename outputValueType,     class ...outputProperties,
               typename leftInputValueType,  class ...leftInputProperties,
               typename rightInputValueType, class ...rightInputProperties>
      static void
      matvecProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                     const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                     const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                     const bool hasField,
                     const bool isTranspose );
      
      template<typename outputValueType,     class ...outputProperties,
               typename leftInputValueType,  class ...leftInputProperties,
               typename rightInputValueType, class ...rightInputProperties>
      static void
      matmatProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                     const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                     const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                     const bool hasField,
                     const bool isTranspose );
    };

  }; // end class ArrayTools

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_ArrayToolsDefContractions.hpp>
#include <Intrepid2_ArrayToolsDefScalar.hpp>
#include <Intrepid2_ArrayToolsDefDot.hpp>
#include <Intrepid2_ArrayToolsDefTensor.hpp>
#include <Intrepid2_ArrayToolsDefCloneScale.hpp>

#endif
