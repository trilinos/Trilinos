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

/** \file   Intrepid_ArrayTools.hpp
    \brief  Header file for utility class to provide array tools,
            such as tensor contractions, etc.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_ARRAYTOOLS_HPP
#define INTREPID_ARRAYTOOLS_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::ArrayTools
      \brief Utility class that provides methods for higher-order algebraic
             manipulation of user-defined arrays, such as tensor contractions.
             For low-order operations, see Intrepid::RealSpaceTools.
  */
  class ArrayTools {
  public:

    /** \brief Contracts the "point" dimension P of two rank-3 containers with
               dimensions (C,L,P) and (C,R,P), and returns the result in a
               rank-3 container with dimensions (C,L,R).

               For a fixed index "C", (C,L,R) represents a rectangular L X R matrix
               where L and R may be different.
        \code
          C - num. integration domains       dim0 in both input containers
          L - num. "left" fields             dim1 in "left" container
          R - num. "right" fields            dim1 in "right" container
          P - num. integration points        dim2 in both input containers
        \endcode

        \param  outputFields   [out] - Output array.
        \param  leftFields      [in] - Left input array.
        \param  rightFields     [in] - Right input array.
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
    static void contractFieldFieldScalar(ArrayOutFields &            outputFields,
                                         const ArrayInFieldsLeft &   leftFields,
                                         const ArrayInFieldsRight &  rightFields,
                                         const ECompEngine           compEngine,
                                         const bool                  sumInto = false);


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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
    static void contractFieldFieldVector(ArrayOutFields &            outputFields,
                                         const ArrayInFieldsLeft &   leftFields,
                                         const ArrayInFieldsRight &  rightFields,
                                         const ECompEngine           compEngine,
                                         const bool                  sumInto = false);

    
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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
    static void contractFieldFieldTensor(ArrayOutFields &            outputFields,
                                         const ArrayInFieldsLeft &   leftFields,
                                         const ArrayInFieldsRight &  rightFields,
                                         const ECompEngine           compEngine,
                                         const bool                  sumInto = false);
    
    
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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void contractDataFieldScalar(ArrayOutFields &       outputFields,
                                        const ArrayInData &    inputData,
                                        const ArrayInFields &  inputFields,
                                        const ECompEngine      compEngine,
                                        const bool             sumInto = false);


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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void contractDataFieldVector(ArrayOutFields &       outputFields,
                                        const ArrayInData &    inputData,
                                        const ArrayInFields &  inputFields,
                                        const ECompEngine      compEngine,
                                        const bool             sumInto = false);

    
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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void contractDataFieldTensor(ArrayOutFields &       outputFields,
                                        const ArrayInData &    inputData,
                                        const ArrayInFields &  inputFields,
                                        const ECompEngine      compEngine,
                                        const bool             sumInto = false);


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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void contractDataDataScalar(ArrayOutData &            outputData,
                                       const ArrayInDataLeft &   inputDataLeft,
                                       const ArrayInDataRight &  inputDataRight,
                                       const ECompEngine         compEngine,
                                       const bool                sumInto = false);


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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void contractDataDataVector(ArrayOutData &            outputData,
                                       const ArrayInDataLeft &   inputDataLeft,
                                       const ArrayInDataRight &  inputDataRight,
                                       const ECompEngine         compEngine,
                                       const bool                sumInto = false);


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
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE. 
    */
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void contractDataDataTensor(ArrayOutData &            outputData,
                                       const ArrayInDataLeft &   inputDataLeft,
                                       const ArrayInDataRight &  inputDataRight,
                                       const ECompEngine         compEngine,
                                       const bool                sumInto = false);


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
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void scalarMultiplyDataField(ArrayOutFields &     outputFields,
                                        const ArrayInData &  inputData,
                                        ArrayInFields &      inputFields,
                                        const bool           reciprocal = false);


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
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void scalarMultiplyDataData(ArrayOutData &           outputData,
                                       ArrayInDataLeft &        inputDataLeft,
                                       ArrayInDataRight &       inputDataRight,
                                       const bool               reciprocal = false);


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
        \endcode

        \param  outputFields   [out] - Output (dot product) fields array.
        \param  inputData       [in] - Data array.
        \param  inputFields     [in] - Input fields array.
    */
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void dotMultiplyDataField(ArrayOutFields &       outputFields,
                                     const ArrayInData &    inputData,
                                     const ArrayInFields &  inputFields);


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
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void dotMultiplyDataData(ArrayOutData &            outputData,
                                    const ArrayInDataLeft &   inputDataLeft,
                                    const ArrayInDataRight &  inputDataRight);


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
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void crossProductDataField(ArrayOutFields &       outputFields,
                                      const ArrayInData &    inputData,
                                      const ArrayInFields &  inputFields);


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
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void crossProductDataData(ArrayOutData &            outputData,
                                     const ArrayInDataLeft &   inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight);


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
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void outerProductDataField(ArrayOutFields &       outputFields,
                                      const ArrayInData &    inputData,
                                      const ArrayInFields &  inputFields);


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
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void outerProductDataData(ArrayOutData &            outputData,
                                     const ArrayInDataLeft &   inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight);


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
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void matvecProductDataField(ArrayOutFields &       outputFields,
                                       const ArrayInData &    inputData,
                                       const ArrayInFields &  inputFields,
                                       const char             transpose = 'N');

    

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
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void matvecProductDataData(ArrayOutData &            outputData,
                                      const ArrayInDataLeft &   inputDataLeft,
                                      const ArrayInDataRight &  inputDataRight,
                                      const char                transpose = 'N');
    
    
    
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
    template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
    static void matmatProductDataField(ArrayOutFields &       outputFields,
                                       const ArrayInData &    inputData,
                                       const ArrayInFields &  inputFields,
                                       const char             transpose = 'N');

    

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
    template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
    static void matmatProductDataData(ArrayOutData &            outputData,
                                      const ArrayInDataLeft &   inputDataLeft,
                                      const ArrayInDataRight &  inputDataRight,
                                      const char                transpose = 'N');


    
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
    template<class Scalar, class ArrayOutFields, class ArrayInFields>
    static void cloneFields(ArrayOutFields &       outputFields,
                            const ArrayInFields &  inputFields);


    /** \brief Multiplies a rank-2, 3, or 4 container with dimensions (F,P),
               (F,P,D1) or (F,P,D1,D2), representing the values of a scalar, vector or a
               tensor field, F-componentwise with a scalar container indexed by (C,F),
               and stores the result in an output value container of size (C,F,P),
               (C,F,P,D1) or (C,F,P,D1,D2).

        \code
          C  - num. integration domains               
          F  - num. fields                            
          P  - num. integration points                
          D1 - first spatial (tensor) dimension index 
          D2 - second spatial (tensor) dimension index
        \endcode

        \param  outputFields   [out] - Output fields array.
        \param  inputFactors    [in] - Input field factors array.
        \param  inputFields     [in] - Input fields array.
    */
    template<class Scalar, class ArrayOutFields, class ArrayInFactors, class ArrayInFields>
    static void cloneScaleFields(ArrayOutFields &        outputFields,
                                 const ArrayInFactors &  inputFactors,
                                 const ArrayInFields &   inputFields);


    /** \brief Multiplies, in place, a rank-2, 3, or 4 container with dimensions (C,F,P),
               (C,F,P,D1) or (C,F,P,D1,D2), representing the values of a scalar, vector or a
               tensor field, F-componentwise with a scalar container indexed by (C,F).

        \code
          C  - num. integration domains               
          F  - num. fields                            
          P  - num. integration points                
          D1 - first spatial (tensor) dimension index 
          D2 - second spatial (tensor) dimension index
        \endcode

        \param  inoutFields    [in/out] - Input / output fields array.
        \param  inputFactors       [in] - Scaling field factors array.
    */
    template<class Scalar, class ArrayInOutFields, class ArrayInFactors>
    static void scaleFields(ArrayInOutFields &      inoutFields,
                            const ArrayInFactors &  inputFactors);


  }; // end class ArrayTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_ArrayToolsDefContractions.hpp>
#include <Intrepid_ArrayToolsDefScalar.hpp>
#include <Intrepid_ArrayToolsDefDot.hpp>
#include <Intrepid_ArrayToolsDefTensor.hpp>
#include <Intrepid_ArrayToolsDefCloneScale.hpp>

#endif
