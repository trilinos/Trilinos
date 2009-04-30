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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeSign, class ArrayTypeIn>
  static void HCURLtransformVALUE(ArrayTypeOut        & outVals,
                                  const ArrayTypeJac  & jacobianInverse,
                                  const ArrayTypeSign & fieldSigns,
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeSign, class ArrayTypeIn>
  static void HCURLtransformCURL(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeSign & fieldSigns,
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeSign, class ArrayTypeIn>
  static void HDIVtransformVALUE(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeSign & fieldSigns,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');

  /** \brief Transformation of a divergence field in the H-div space, defined at points on a
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeSign, class ArrayTypeIn>
  static void HDIVtransformDIV(ArrayTypeOut        & outVals,
                               const ArrayTypeDet  & jacobianDet,
                               const ArrayTypeSign & fieldSigns,
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

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData, class ArrayTypeIn>
  static void multiplyData(ArrayTypeOut             & outVals,
                           const ArrayTypeData      & inData,
                           const ArrayTypeIn        & inVals,
                           const char               transpose = 'N');

  template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
  static void multiplyMeasure(ArrayTypeOut             & outVals,
                              const ArrayTypeMeasure   & inMeasure,
                              const ArrayTypeIn        & inVals);


  // updated stuff



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


  template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
  static void computeMeasure(ArrayOut             & outVals,
                             const ArrayDet       & inDet,
                             const ArrayWeights   & inWeights);






  /////////////////////////////////// new stuff

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



  
};  // end FunctionSpaceTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FunctionSpaceToolsDef.hpp>

#endif
