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

/** \file   Intrepid_FunctionSpaceToolsInPlace.hpp
    \brief  Header file for the Intrepid::FunctionSpaceToolsInPlace class.
    \author Created by R. Kirby
*/

#ifndef INTREPID_FUNCTIONSPACETOOLSINPLACE_HPP
#define INTREPID_FUNCTIONSPACETOOLSINPLACE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"


namespace Intrepid {

/** \class Intrepid::FunctionSpaceToolsInPlace
    \brief Defines expert-level interfaces for the evaluation of functions
           and operators in physical space (supported for FE, FV, and FD methods)
           and FE reference space; in addition, provides several function
           transformation utilities.

           The functionality here largely mirrors that in Intrepid::FunctionSpaceTools,
	   except that the input data is overwrriten with the result of the particular
	   transformation.  This can reduce intermediate storage when the input values
	   are not required by later calculations.

	   A new feature compared to Intrepid::FunctionSpaceTools is that of "dual"
	   transforms that are useful in alternate workflow patterns.  At the innermost
	   loop nest, we may have a computation of the form (T_1 D_1 u).(T_2 D_2 v), and
	   we would like to rewrite this as (T_2^t T_1 D_1 u).(D_2 v).  The transposes
	   of each of these transformations are supplied in this routine.
*/
class FunctionSpaceToolsInPlace {
  public:
  /** \brief Transformation of a (scalar) value field in the H-grad space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P).
   
             Computes pullback of \e HGRAD functions \f$\Phi^*(\widehat{u}_f) = \widehat{u}_f\circ F^{-1}_{c} \f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the values of the function set \f$\{\widehat{u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p) = \widehat{u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
             outVals(c,f,p) 
                = \widehat{u}_f\circ F^{-1}_{c}(x_{c,p}) 
                = \widehat{u}_f(\widehat{x}_p) =  inVals(f,p) \qquad 0\le c < C \,,
      \f]
            i.e., it simply replicates the values in the user-provided container to every cell. 
            See Section \ref sec_pullbacks for more details about pullbacks. 
    
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < number of fields per cell to transform |
    |   P  |         point        |  0 <= P < num. integration points                |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayType>
  static void HGRADtransformVALUE(ArrayType       & inOutVals );

  /** \brief Since there is no matrix involved, this is the same transformation 
      as HGRADtransformVALUE */
  template<class Scalar, class ArrayType>
  static void HGRADtransformVALUEDual(ArrayType       & inOutVals );

  /** \brief Transformation of a gradient field in the H-grad space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Computes pullback of gradients of \e HGRAD functions 
             \f$\Phi^*(\nabla\widehat{u}_f) = \left((DF_c)^{-{\sf T}}\cdot\nabla\widehat{u}_f\right)\circ F^{-1}_{c}\f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the gradients of the function set \f$\{\widehat{u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p,*) = \nabla\widehat{u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
             outVals(c,f,p,*) 
                  = \left((DF_c)^{-{\sf T}}\cdot\nabla\widehat{u}_f\right)\circ F^{-1}_{c}(x_{c,p}) 
                  = (DF_c)^{-{\sf T}}(\widehat{x}_p)\cdot\nabla\widehat{u}_f(\widehat{x}_p) \qquad 0\le c < C \,.
      \f]
             See Section \ref sec_pullbacks for more details about pullbacks.
  
    \code
    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                   Dimension                      |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   F  |         field        |  0 <= F < number of fields per cell              |
    |   P  |         point        |  0 <= P < num. integration points                |
    |   D  |         space dim    |  0 <= D < spatial dimension                      |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayType, class ArrayTypeJac>
  static void HGRADtransformGRAD(ArrayType          & inOutVals,
                                 const ArrayTypeJac & jacobianInverse,
                                 const char           transpose = 'T');

  /** \brief Applies the transpose of the HGRADtransformGRAD to the data */

  template<class Scalar, class ArrayType, class ArrayTypeJac>
  static void HGRADtransformGRADDual(ArrayType       & inOutVals,
				     const ArrayTypeJac & jacobianInverse,
				     const char           transpose = 'T');

  /** \brief Transformation of a (vector) value field in the H-curl space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Computes pullback of \e HCURL functions 
             \f$\Phi^*(\widehat{\bf u}_f) = \left((DF_c)^{-{\sf T}}\cdot\widehat{\bf u}_f\right)\circ F^{-1}_{c}\f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the values of the vector function set \f$\{\widehat{\bf u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p,*) = \widehat{\bf u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
              outVals(c,f,p,*) 
                = \left((DF_c)^{-{\sf T}}\cdot\widehat{\bf u}_f\right)\circ F^{-1}_{c}(x_{c,p}) 
                = (DF_c)^{-{\sf T}}(\widehat{x}_p)\cdot\widehat{\bf u}_f(\widehat{x}_p) \qquad 0\le c < C \,.
      \f]
            See Section \ref sec_pullbacks for more details about pullbacks.
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
  template<class Scalar, class ArrayType, class ArrayTypeJac>
  static void HCURLtransformVALUE(ArrayType        & inOutVals,
                                  const ArrayTypeJac  & jacobianInverse,
                                  const char            transpose = 'T');

  /** \brief Applies the dual of the HCURLtransformVALUE  transformation */
  template<class Scalar, class ArrayType, class ArrayTypeJac>
  static void HCURLtransformVALUEDual(ArrayType           & outVals,
				      const ArrayTypeJac  & jacobianInverse,
				      const char            transpose = 'T');

  /** \brief Transformation of a curl field in the H-curl space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Computes pullback of curls of \e HCURL functions 
             \f$\Phi^*(\widehat{\bf u}_f) = \left(J^{-1}_{c} DF_{c}\cdot\nabla\times\widehat{\bf u}_f\right)\circ F^{-1}_{c}\f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the curls of the vector function set \f$\{\widehat{\bf u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p,*) = \nabla\times\widehat{\bf u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
             outVals(c,f,p,*) 
                = \left(J^{-1}_{c} DF_{c}\cdot\nabla\times\widehat{\bf u}_f\right)\circ F^{-1}_{c}(x_{c,p}) 
                = J^{-1}_{c}(\widehat{x}_p) DF_{c}(\widehat{x}_p)\cdot\nabla\times\widehat{\bf u}_f(\widehat{x}_p)
              \qquad 0\le c < C \,.
      \f]
             See Section \ref sec_pullbacks for more details about pullbacks.
    
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
  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  static void HCURLtransformCURL(ArrayType           & inOutVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const char            transpose = 'N');

  /** \brief Applies the dual of the HCURLtransformCURL transformation */
  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  static void HCURLtransformCURLDual(ArrayType           & outVals,
				     const ArrayTypeJac  & jacobian,
				     const ArrayTypeDet  & jacobianDet,
				     const char            transpose = 'N');

  /** \brief Transformation of a (vector) value field in the H-div space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P,D), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P,D).

             Computes pullback of \e HDIV functions 
             \f$\Phi^*(\widehat{\bf u}_f) = \left(J^{-1}_{c} DF_{c}\cdot\widehat{\bf u}_f\right)\circ F^{-1}_{c} \f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the values of the vector function set \f$\{\widehat{\bf u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p,*) = \widehat{\bf u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
             outVals(c,f,p,*) 
              = \left(J^{-1}_{c} DF_{c}\cdot \widehat{\bf u}_f\right)\circ F^{-1}_{c}(x_{c,p}) 
              = J^{-1}_{c}(\widehat{x}_p) DF_{c}(\widehat{x}_p)\cdot\widehat{\bf u}_f(\widehat{x}_p)
              \qquad 0\le c < C \,.
      \f]
             See Section \ref sec_pullbacks for more details about pullbacks.
    
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
  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  static void HDIVtransformVALUE(ArrayType           & inOutVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const char            transpose = 'N');

  /** \brief Applies the dual of HDIVtransformVALUE */
  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  static void HDIVtransformVALUEDual(ArrayType           & outVals,
				     const ArrayTypeJac  & jacobian,
				     const ArrayTypeDet  & jacobianDet,
				     const char            transpose = 'N');

  /** \brief Transformation of a divergence field in the H-div space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P).

             Computes pullback of the divergence of \e HDIV functions 
             \f$\Phi^*(\widehat{\bf u}_f) = \left(J^{-1}_{c}\nabla\cdot\widehat{\bf u}_{f}\right) \circ F^{-1}_{c} \f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the divergencies of the vector function set \f$\{\widehat{\bf u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p) = \nabla\cdot\widehat{\bf u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
             outVals(c,f,p,*) 
                = \left(J^{-1}_{c}\nabla\cdot\widehat{\bf u}_{f}\right) \circ F^{-1}_{c} (x_{c,p}) 
                = J^{-1}_{c}(\widehat{x}_p) \nabla\cdot\widehat{\bf u}_{f} (\widehat{x}_p)
                \qquad 0\le c < C \,.
      \f]
             See Section \ref sec_pullbacks for more details about pullbacks.
    
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
  template<class Scalar, class ArrayType, class ArrayTypeDet>
  static void HDIVtransformDIV(ArrayType           & inOutVals,
                               const ArrayTypeDet  & jacobianDet);

  /** \brief Applies the dual of HDIVtransformDIV, which is the same */
  template<class Scalar, class ArrayType, class ArrayTypeDet>
  static void HDIVtransformDIVDual(ArrayType           & inOutVals,
				   const ArrayTypeDet  & jacobianDet);

  /** \brief Transformation of a (scalar) value field in the H-vol space, defined at points on a
             reference cell, stored in the user-provided container <var><b>inVals</b></var>
             and indexed by (F,P), into the output container <var><b>outVals</b></var>,
             defined on cells in physical space and indexed by (C,F,P).

             Computes pullback of \e HVOL functions 
             \f$\Phi^*(\widehat{u}_f) = \left(J^{-1}_{c}\widehat{u}_{f}\right) \circ F^{-1}_{c} \f$ 
             for points in one or more physical cells that are images of a given set of points in the reference cell:
      \f[
             \{ x_{c,p} \}_{p=0}^P = \{ F_{c} (\widehat{x}_p) \}_{p=0}^{P}\qquad 0\le c < C \,.
      \f]     
             In this case \f$ F^{-1}_{c}(x_{c,p}) = \widehat{x}_p \f$ and the user-provided container
             should contain the values of the functions in the set \f$\{\widehat{\bf u}_f\}_{f=0}^{F}\f$ at the 
             reference points:
      \f[
             inVals(f,p) = \widehat{u}_f(\widehat{x}_p) \,.
      \f]
             The method returns   
      \f[
             outVals(c,f,p,*) 
                = \left(J^{-1}_{c}\widehat{u}_{f}\right) \circ F^{-1}_{c} (x_{c,p}) 
                = J^{-1}_{c}(\widehat{x}_p) \widehat{u}_{f} (\widehat{x}_p)
                \qquad 0\le c < C \,.
      \f]
             See Section \ref sec_pullbacks for more details about pullbacks.
    
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
  template<class Scalar, class ArrayType, class ArrayTypeDet>
  static void HVOLtransformVALUE(ArrayType           & inOutVals,
                                 const ArrayTypeDet  & jacobianDet);

  /** \brief Applies the dual of HVOLtransformVALUE */
  template<class Scalar, class ArrayType, class ArrayTypeDet>
  static void HVOLtransformVALUEDual(ArrayType           & inOutVals,
				     const ArrayTypeDet  & jacobianDet);
  
  template<class Scalar, class ArrayType, class ArrayTypeMeasure>
  static void multiplyMeasure(ArrayType                & inOutVals,
                              const ArrayTypeMeasure   & inMeasure);


};  // end FunctionSpaceToolsInPlace

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FunctionSpaceToolsInPlaceDef.hpp>

#endif






















#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

