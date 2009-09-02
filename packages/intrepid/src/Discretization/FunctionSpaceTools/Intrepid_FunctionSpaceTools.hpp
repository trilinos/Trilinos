// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copytest (2007) Sandia Corporation
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

/** \file   Intrepid_FunctionSpaceTools.hpp
    \brief  Header file for the Intrepid::FunctionSpaceTools class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_FUNCTIONSPACETOOLS_HPP
#define INTREPID_FUNCTIONSPACETOOLS_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"


namespace Intrepid {

/** \class Intrepid::FunctionSpaceTools
    \brief Defines expert-level interfaces for the evaluation of functions
           and operators in physical space (supported for FE, FV, and FD methods)
           and FE reference space; in addition, provides several function
           transformation utilities.
*/
class FunctionSpaceTools {
  public:
  /** \brief Transformation of a (scalar) value field in the H-grad space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn>
  static void HGRADtransformVALUE(ArrayTypeOut       & outVals,
                                  const ArrayTypeIn  & inVals);

  /** \brief Transformation of a gradient field in the H-grad space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |   D  |         space dim    |  0 <= D < spatial dimension                      |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
  static void HGRADtransformGRAD(ArrayTypeOut       & outVals,
                                 const ArrayTypeJac & jacobianInverse,
                                 const ArrayTypeIn  & inVals,
                                 const char           transpose = 'T');

  /** \brief Transformation of a (vector) value field in the H-curl space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of native basis                   |
    |   P  |         point        |  0 <= P < num. integration points                |
    |   D  |         space dim    |  0 <= D < spatial dimension                      |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
  static void HCURLtransformVALUE(ArrayTypeOut        & outVals,
                                  const ArrayTypeJac  & jacobianInverse,
                                  const ArrayTypeIn   & inVals,
                                  const char            transpose = 'T');

  /** \brief Transformation of a curl field in the H-curl space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |   D  |         space dim    |  0 <= D < spatial dimension                      |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
  static void HCURLtransformCURL(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');

  /** \brief Transformation of a (vector) value field in the H-div space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |   D  |         space dim    |  0 <= D < spatial dimension                      |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
  static void HDIVtransformVALUE(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');

  /** \brief Transformation of a divergence field in the H-div space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
  static void HDIVtransformDIV(ArrayTypeOut        & outVals,
                               const ArrayTypeDet  & jacobianDet,
                               const ArrayTypeIn   & inVals);

  /** \brief Transformation of a (scalar) value field in the H-vol space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P).

             Math here ...
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
  static void HVOLtransformVALUE(ArrayTypeOut        & outVals,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals);

  /** \brief   Contracts \a <b>leftValues</b> and \a <b>rightValues</b> arrays on the
               point and possibly space dimensions and stores the result in \a <b>outputValues</b>;
               this is a generic, high-level integration routine that calls either
               FunctionSpaceTools::operatorIntegral, or FunctionSpaceTools::functionalIntegral,
               or FunctionSpaceTools::dataIntegral methods, depending on the rank of the
               \a <b>outputValues</b> array.

        \param  outputValues   [out] - Output array.
        \param  leftValues      [in] - Left input array.
        \param  rightValues     [in] - Right input array.
        \param  compEngine      [in] - Computational engine.
        \param  sumInto         [in] - If TRUE, sum into given output array,
                                       otherwise overwrite it. Default: FALSE.
  */
  template<class Scalar, class ArrayOut, class ArrayInLeft, class ArrayInRight>
  static void integrate(ArrayOut            & outputValues,
                        const ArrayInLeft   & leftValues,
                        const ArrayInRight  & rightValues,
                        const ECompEngine     compEngine,
                        const bool            sumInto = false);

  /** \brief   Contracts the point (and space) dimensions P (and D1 and D2) of
               two rank-3, 4, or 5 containers with dimensions (C,L,P) and (C,R,P),
               or (C,L,P,D1) and (C,R,P,D1), or (C,L,P,D1,D2) and (C,R,P,D1,D2),
               and returns the result in a rank-3 container with dimensions (C,L,R).

               For a fixed index "C", (C,L,R) represents a rectangular L X R matrix
               where L and R may be different.
        \code
          C - num. integration domains       dim0 in both input containers
          L - num. "left" fields             dim1 in "left" container
          R - num. "right" fields            dim1 in "right" container
          P - num. integration points        dim2 in both input containers
          D1- vector (1st tensor) dimension  dim3 in both input containers
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
  static void operatorIntegral(ArrayOutFields &            outputFields,
                               const ArrayInFieldsLeft &   leftFields,
                               const ArrayInFieldsRight &  rightFields,
                               const ECompEngine           compEngine,
                               const bool                  sumInto = false);

 /** \brief    Contracts the point (and space) dimensions P (and D1 and D2) of a
               rank-3, 4, or 5 container and a rank-2, 3, or 4 container, respectively,
               with dimensions (C,F,P) and (C,P), or (C,F,P,D1) and (C,P,D1),
               or (C,F,P,D1,D2) and (C,P,D1,D2), respectively, and returns the result
               in a rank-2 container with dimensions (C,F).

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
  static void functionalIntegral(ArrayOutFields &       outputFields,
                                 const ArrayInData &    inputData,
                                 const ArrayInFields &  inputFields,
                                 const ECompEngine      compEngine,
                                 const bool             sumInto = false);

  /** \brief   Contracts the point (and space) dimensions P (and D1 and D2) of two
               rank-2, 3, or 4 containers with dimensions (C,P), or (C,P,D1), or
               (C,P,D1,D2), respectively, and returns the result in a rank-1 container
               with dimensions (C).

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
  static void dataIntegral(ArrayOutData &            outputData,
                           const ArrayInDataLeft &   inputDataLeft,
                           const ArrayInDataRight &  inputDataRight,
                           const ECompEngine         compEngine,
                           const bool                sumInto = false);


  template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
  static void computeMeasure(ArrayOut             & outVals,
                             const ArrayDet       & inDet,
                             const ArrayWeights   & inWeights);

  template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
  static void computeCellMeasure(ArrayOut             & outVals,
                                 const ArrayDet       & inDet,
                                 const ArrayWeights   & inWeights);

  template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  static void computeFaceMeasure(ArrayOut                   & outVals,
                                 const ArrayJac             & inJac,
                                 const ArrayWeights         & inWeights,
                                 const int                    whichFace,
                                 const shards::CellTopology & parentCell);

  template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  static void computeEdgeMeasure(ArrayOut                   & outVals,
                                 const ArrayJac             & inJac,
                                 const ArrayWeights         & inWeights,
                                 const int                    whichEdge,
                                 const shards::CellTopology & parentCell);

  template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
  static void multiplyMeasure(ArrayTypeOut             & outVals,
                              const ArrayTypeMeasure   & inMeasure,
                              const ArrayTypeIn        & inVals);

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
                                      ArrayInData &        inputData,
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

  /** \brief There are four use cases:
             (1) cross product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
             representing the values of a set of vector fields, on the left by the values in a rank-3
             container \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
             (2) cross product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
             representing the values of a vector field, on the left by the values in a rank-3 container
             \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
             (3) outer product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
             representing the values of a set of vector fields, on the left by the values in a rank-3
             container \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
             (4) outer product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
             representing the values of a vector field, on the left by the values in a rank-3 container
             \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data;
             for cross products, the output value container \a <b>outputFields</b> is indexed by
             (C,F,P,D) in 3D (vector output) and by (C,F,P) in 2D (scalar output);
             for outer products, the output value container \a <b>outputFields</b> is indexed by (C,F,P,D,D).

      \code
        C  - num. integration domains
        F  - num. fields
        P  - num. integration points
        D  - spatial dimension, must be 2 or 3
      \endcode

      \param  outputFields   [out] - Output (cross or outer product) fields array.
      \param  inputData       [in] - Data array.
      \param  inputFields     [in] - Input fields array.
  */
  template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
  static void vectorMultiplyDataField(ArrayOutFields &       outputFields,
                                      const ArrayInData &    inputData,
                                      const ArrayInFields &  inputFields);

  /** \brief There are four use cases:
             (1) cross product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
             representing the values of a set of vector data, on the left by the values in a rank-3
             container \a <b>inputDataLeft</b> indexed by (C,P,D) representing the values of vector data, OR
             (2) cross product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
             representing the values of vector data, on the left by the values in a rank-3 container
             \a <b>inputDataLeft</b> indexed by (C,P,D), representing the values of vector data, OR
             (3) outer product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
             representing the values of a set of vector data, on the left by the values in a rank-3
             container \a <b>inputDataLeft</b> indexed by (C,P,D) representing the values of vector data, OR
             (4) outer product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
             representing the values of vector data, on the left by the values in a rank-3 container
             \a <b>inputDataLeft</b> indexed by (C,P,D), representing the values of vector data;
             for cross products, the output value container \a <b>outputData</b> is indexed by
             (C,P,D) in 3D (vector output) and by (C,P) in 2D (scalar output);
             for outer products, the output value container \a <b>outputData</b> is indexed by (C,P,D,D).

      \code
        C  - num. integration domains
        P  - num. integration points
        D  - spatial dimension, must be 2 or 3
      \endcode

      \param  outputData      [out] - Output (cross or outer product) data array.
      \param  inputDataLeft    [in] - Left input data array.
      \param  inputDataRight   [in] - Right input data array.
  */
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  static void vectorMultiplyDataData(ArrayOutData &            outputData,
                                     const ArrayInDataLeft &   inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight);

  /** \brief There are four use cases:
             (1) matrix-vector product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
             representing the values of a set of vector fields, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             (2) matrix-vector product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
             representing the values of a vector field, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             (3) matrix-matrix product of a rank-5 container \a <b>inputFields</b> with dimensions (C,F,P,D,D),
             representing the values of a set of tensor fields, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             (4) matrix-matrix product of a rank-4 container \a <b>inputFields</b> with dimensions (F,P,D,D),
             representing the values of a tensor field, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data;
             for matrix-vector products, the output value container \a <b>outputFields</b> is
             indexed by (C,F,P,D);
             for matrix-matrix products the output value container \a <b>outputFields</b> is
             indexed by (C,F,P,D,D).

      \remarks
             The rank of \a <b>inputData</b> implicitly defines the type of tensor data:
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
      \endcode

      \param  outputFields   [out] - Output (matrix-vector or matrix-matrix product) fields array.
      \param  inputData       [in] - Data array.
      \param  inputFields     [in] - Input fields array.
      \param  transpose       [in] - If 'T', use transposed left data tensor; if 'N', no transpose. Default: 'N'.
  */
  template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
  static void tensorMultiplyDataField(ArrayOutFields &       outputFields,
                                      const ArrayInData &    inputData,
                                      const ArrayInFields &  inputFields,
                                      const char             transpose = 'N');

  /** \brief There are four use cases:
             (1) matrix-vector product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
             representing the values of a set of vector data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             (2) matrix-vector product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
             representing the values of vector data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             (3) matrix-matrix product of a rank-4 container \a <b>inputDataRight</b> with dimensions (C,P,D,D),
             representing the values of a set of tensor data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             (4) matrix-matrix product of a rank-3 container \a <b>inputDataRight</b> with dimensions (P,D,D),
             representing the values of tensor data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data;
             for matrix-vector products, the output value container \a <b>outputData</b>
             is indexed by (C,P,D);
             for matrix-matrix products, the output value container \a <b>outputData</b>
             is indexed by (C,P,D1,D2).

      \remarks
            The rank of <b>inputDataLeft</b> implicitly defines the type of tensor data:
            \li rank = 2 corresponds to a constant diagonal tensor \f$ diag(a,\ldots,a) \f$
            \li rank = 3 corresponds to a nonconstant diagonal tensor \f$ diag(a_1,\ldots,a_d) \f$
            \li rank = 4 corresponds to a full tensor \f$ \{a_{ij}\}\f$

      \note  It is assumed that all tensors are square!

      \note  The method is defined for spatial dimensions D = 1, 2, 3

      \code
        C    - num. integration domains
        P    - num. integration points
        D    - spatial dimension
      \endcode

      \param  outputData      [out] - Output (matrix-vector product) data array.
      \param  inputDataLeft    [in] - Left input data array.
      \param  inputDataRight   [in] - Right input data array.
      \param  transpose        [in] - If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'.
  */
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  static void tensorMultiplyDataData(ArrayOutData &            outputData,
                                     const ArrayInDataLeft &   inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight,
                                     const char                transpose = 'N');


  /** \brief Applies left (row) signs, stored in the user-provided container
             <var><b>fieldSigns</b></var> and indexed by (C,L), to the operator
             <var><b>inoutFields</b></var> indexed by (C,L,R).

             Math here ...
    \code
    |------|----------------------|----------------------------------------------------|
    |      |         Index        |                   Dimension                        |
    |------|----------------------|----------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains                 |
    |   L  | num. "left" fields   |  0 <= L < dim. of the left basis (in inoutFields)  |
    |   R  | num. "right" fields  |  0 <= R < dim. of the right basis (in inoutFields) |
    |------|----------------------|----------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyLeftFieldSigns(ArrayTypeInOut        & inoutOperator,
                                  const ArrayTypeSign   & fieldSigns);

  /** \brief Applies right (column) signs, stored in the user-provided container
             <var><b>fieldSigns</b></var> and indexed by (C,R), to the operator
             <var><b>inoutFields</b></var> indexed by (C,L,R).

             Math here ...
    \code
    |------|----------------------|----------------------------------------------------|
    |      |         Index        |                   Dimension                        |
    |------|----------------------|----------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains                 |
    |   L  | num. "left" fields   |  0 <= L < dim. of the left basis (in inoutFields)  |
    |   R  | num. "right" fields  |  0 <= R < dim. of the right basis (in inoutFields) |
    |------|----------------------|----------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyRightFieldSigns(ArrayTypeInOut        & inoutOperator,
                                   const ArrayTypeSign   & fieldSigns);

  /** \brief Applies field signs, stored in the user-provided container
             <var><b>fieldSigns</b></var> and indexed by (C,F), to the functional
             <var><b>inoutFields</b></var> indexed by (C,F).

             Math here ...
    \code
    |------|----------------------|----------------------------------------------------|
    |      |         Index        |                   Dimension                        |
    |------|----------------------|----------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains                 |
    |   F  |      num. fields     |  0 <= F < dim. of the basis (in inoutFields)       |
    |------|----------------------|----------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyFieldSigns(ArrayTypeInOut        & inoutFunctional,
                              const ArrayTypeSign   & fieldSigns);

  
};  // end FunctionSpaceTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FunctionSpaceToolsDef.hpp>

#endif
