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

  template<class Scalar, class ArrayTypeOut, class ArrayTypeInLeft, class ArrayTypeInRight>
  static void integrate(ArrayTypeOut            & outputValues,
                        const ArrayTypeInLeft   & leftValues,
                        const ArrayTypeInRight  & rightValues,
                        const ECompEngine         compEngine);

  template<class Scalar, class ArrayTypeOut, class ArrayTypeWeights, class ArrayTypeDet>
  static void computeMeasure(ArrayTypeOut             & outVals,
                             const ArrayTypeWeights   & inWeights,
                             const ArrayTypeDet       & inDet);
  
};  // end FunctionSpaceTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FunctionSpaceToolsDef.hpp>

#endif
