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

#include <Intrepid_Rank.hpp>


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
    |   F  |         field        |  0 <= F < dim. of the basis                      |
    |   P  |         point        |  0 <= P < num. integration points                |
    |------|----------------------|--------------------------------------------------|
    \endcode
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn>
  static void HGRADtransformVALUE(ArrayTypeOut       & outVals,
                                  const ArrayTypeIn  & inVals);
/*
  template<class Scalar, class ArrayTypeOut, class ArrayTypeIn>
  static void HGRADtransformVALUETemp(ArrayTypeOut       & outVals,
                                  const ArrayTypeIn  & inVals);*/
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
                                 /*
    template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
  static void HGRADtransformGRADTemp(ArrayTypeOut       & outVals,
                                 const ArrayTypeJac & jacobianInverse,
                                 const ArrayTypeIn  & inVals,
                                 const char           transpose = 'T');
*/
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
  static void HCURLtransformVALUE(ArrayTypeOut        & outVals,
                                  const ArrayTypeJac  & jacobianInverse,
                                  const ArrayTypeIn   & inVals,
                                  const char            transpose = 'T');
                                  /*
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
  static void HCURLtransformVALUETemp(ArrayTypeOut        & outVals,
                                  const ArrayTypeJac  & jacobianInverse,
                                  const ArrayTypeIn   & inVals,
                                  const char            transpose = 'T');*/
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
  static void HCURLtransformCURL(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');
                                 /*
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
  static void HCURLtransformCURLTemp(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');*/
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
  static void HDIVtransformVALUE(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');
                                 /*
  template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
  static void HDIVtransformVALUETemp(ArrayTypeOut        & outVals,
                                 const ArrayTypeJac  & jacobian,
                                 const ArrayTypeDet  & jacobianDet,
                                 const ArrayTypeIn   & inVals,
                                 const char            transpose = 'N');*/
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
  template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
  static void HDIVtransformDIV(ArrayTypeOut        & outVals,
                               const ArrayTypeDet  & jacobianDet,
                               const ArrayTypeIn   & inVals);
                               /*
  template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
  static void HDIVtransformDIVTemp(ArrayTypeOut        & outVals,
                               const ArrayTypeDet  & jacobianDet,
                               const ArrayTypeIn   & inVals);                               
*/
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
  template<class Scalar>
  static void integrate(Intrepid::FieldContainer<Scalar>            & outputValues,
                                   const Intrepid::FieldContainer<Scalar>   & leftValues,
                                   const Intrepid::FieldContainer<Scalar>  & rightValues,
                                   const ECompEngine           compEngine,
                                   const bool            sumInto = false);

  template<class Scalar, class ArrayOut, class ArrayInLeft, class ArrayInRight>
  static void integrate(ArrayOut            & outputValues,
                        const ArrayInLeft   & leftValues,
                        const ArrayInRight  & rightValues,
                        const ECompEngine     compEngine,
                        const bool            sumInto = false);
    
  /*  template<class Scalar, class ArrayOut, class ArrayInLeft, class ArrayInRight>
    static void integrateTemp(ArrayOut            & outputValues,
                          const ArrayInLeft   & leftValues,
                          const ArrayInRight  & rightValues,
                          const ECompEngine     compEngine,
                          const bool            sumInto = false);
                          */
    template<class Scalar, class ArrayOut, class ArrayInLeft, class ArrayInRight,int leftrank,int outrank>
    struct integrateTempSpec;

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
/* template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
  static void operatorIntegralTemp(ArrayOutFields &            outputFields,
                               const ArrayInFieldsLeft &   leftFields,
                               const ArrayInFieldsRight &  rightFields,
                               const ECompEngine           compEngine,
                               const bool                  sumInto = false);*/
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
 /*   template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
  static void functionalIntegralTemp(ArrayOutFields &       outputFields,
                                 const ArrayInData &    inputData,
                                 const ArrayInFields &  inputFields,
                                 const ECompEngine      compEngine,
                                 const bool             sumInto = false);
*/
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
   /* template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  static void dataIntegralTemp(ArrayOutData &            outputData,
                           const ArrayInDataLeft &   inputDataLeft,
                           const ArrayInDataRight &  inputDataRight,
                           const ECompEngine         compEngine,
                           const bool                sumInto = false);
*/
  /** \brief   Returns the weighted integration measures \a <b>outVals</b> with dimensions
               (C,P) used for the computation of cell integrals, by multiplying absolute values 
               of the user-provided cell Jacobian determinants \a <b>inDet</b> with dimensions (C,P) 
               with the user-provided integration weights \a <b>inWeights</b> with dimensions (P).

               Returns a rank-2 array (C, P) array such that
        \f[
               \mbox{outVals}(c,p)   = |\mbox{det}(DF_{c}(\widehat{x}_p))|\omega_{p} \,,
        \f]
               where \f$\{(\widehat{x}_p,\omega_p)\}\f$ is a cubature rule defined on a reference cell
               (a set of integration points and their associated weights; see
               Intrepid::Cubature::getCubature for getting cubature rules on reference cells). 
        \warning 
               The user is responsible for providing input arrays with consistent data: the determinants
               in \a <b>inDet</b> should be evaluated at integration points on the <b>reference cell</b> 
               corresponding to the weights in \a <b>inWeights</b>.
 
        \remark
               See Intrepid::CellTools::setJacobian for computation of \e DF and 
               Intrepid::CellTools::setJacobianDet for computation of its determinant.

        \code
          C - num. integration domains                     dim0 in all containers
          P - num. integration points                      dim1 in all containers
        \endcode

        \param  outVals     [out] - Output array with weighted cell measures.
        \param  inDet        [in] - Input array containing determinants of cell Jacobians.
        \param  inWeights    [in] - Input integration weights.
  */
  template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
  static void computeCellMeasure(ArrayOut             & outVals,
                                 const ArrayDet       & inDet,
                                 const ArrayWeights   & inWeights);

  /*template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
  static void computeCellMeasureTemp(ArrayOut             & outVals,
                                 const ArrayDet       & inDet,
                                 const ArrayWeights   & inWeights);*/
  /** \brief   Returns the weighted integration measures \a <b>outVals</b> with dimensions
               (C,P) used for the computation of face integrals, based on the provided
               cell Jacobian array \a <b>inJac</b> with dimensions (C,P,D,D) and the
               provided integration weights \a <b>inWeights</b> with dimensions (P). 

               Returns a rank-2 array (C, P) array such that
      \f[
               \mbox{outVals}(c,p)   = 
                \left\|\frac{\partial\Phi_c(\widehat{x}_p)}{\partial u}\times 
                       \frac{\partial\Phi_c(\widehat{x}_p)}{\partial v}\right\|\omega_{p} \,,
      \f]
               where: 
      \li      \f$\{(\widehat{x}_p,\omega_p)\}\f$ is a cubature rule defined on \b reference 
               \b face \f$\widehat{\mathcal{F}}\f$, with ordinal \e whichFace relative to the specified parent reference cell;
      \li      \f$ \Phi_c : R \mapsto \mathcal{F} \f$ is parameterization of the physical face
               corresponding to \f$\widehat{\mathcal{F}}\f$; see Section \ref sec_cell_topology_subcell_map.
    
      \warning 
               The user is responsible for providing input arrays with consistent data: the Jacobians
               in \a <b>inJac</b> should be evaluated at integration points on the <b>reference face</b>
               corresponding to the weights in \a <b>inWeights</b>.
    
      \remark 
              Cubature rules on reference faces are defined by a two-step process:
      \li     A cubature rule is defined on the parametrization domain \e R of the face 
              (\e R is the standard 2-simplex {(0,0),(1,0),(0,1)} or the standard 2-cube [-1,1] X [-1,1]).
      \li     The points are mapped to a reference face using Intrepid::CellTools::mapToReferenceSubcell

      \remark
               See Intrepid::CellTools::setJacobian for computation of \e DF and 
               Intrepid::CellTools::setJacobianDet for computation of its determinant.
    
        \code
          C - num. integration domains                     dim0 in all input containers
          P - num. integration points                      dim1 in all input containers
          D - spatial dimension                            dim2 and dim3 in Jacobian container
        \endcode

        \param  outVals     [out] - Output array with weighted face measures.
        \param  inJac        [in] - Input array containing cell Jacobians.
        \param  inWeights    [in] - Input integration weights.
        \param  whichFace    [in] - Index of the face subcell relative to the parent cell; defines the domain of integration.
        \param  parentCell   [in] - Parent cell topology.
  */
  template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  static void computeFaceMeasure(ArrayOut                   & outVals,
                                 const ArrayJac             & inJac,
                                 const ArrayWeights         & inWeights,
                                 const int                    whichFace,
                                 const shards::CellTopology & parentCell);

/*  template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  static void computeFaceMeasureTemp(ArrayOut                   & outVals,
                                 const ArrayJac             & inJac,
                                 const ArrayWeights         & inWeights,
                                 const int                    whichFace,
                                 const shards::CellTopology & parentCell);*/
  /** \brief   Returns the weighted integration measures \a <b>outVals</b> with dimensions
               (C,P) used for the computation of edge integrals, based on the provided
               cell Jacobian array \a <b>inJac</b> with dimensions (C,P,D,D) and the
               provided integration weights \a <b>inWeights</b> with dimensions (P). 

               Returns a rank-2 array (C, P) array such that
      \f[
               \mbox{outVals}(c,p)   = 
                    \left\|\frac{d \Phi_c(\widehat{x}_p)}{d s}\right\|\omega_{p} \,,
      \f]
               where: 
      \li      \f$\{(\widehat{x}_p,\omega_p)\}\f$ is a cubature rule defined on \b reference 
               \b edge \f$\widehat{\mathcal{E}}\f$, with ordinal \e whichEdge relative to the specified parent reference cell;
      \li      \f$ \Phi_c : R \mapsto \mathcal{E} \f$ is parameterization of the physical edge
               corresponding to \f$\widehat{\mathcal{E}}\f$; see Section \ref sec_cell_topology_subcell_map.
    
      \warning 
               The user is responsible for providing input arrays with consistent data: the Jacobians
               in \a <b>inJac</b> should be evaluated at integration points on the <b>reference edge</b>
               corresponding to the weights in \a <b>inWeights</b>.
    
      \remark 
               Cubature rules on reference edges are defined by a two-step process:
      \li      A cubature rule is defined on the parametrization domain \e R = [-1,1] of the edge. 
      \li      The points are mapped to a reference edge using Intrepid::CellTools::mapToReferenceSubcell
    
      \remark
               See Intrepid::CellTools::setJacobian for computation of \e DF and 
               Intrepid::CellTools::setJacobianDet for computation of its determinant.

        \code
          C - num. integration domains                     dim0 in all input containers
          P - num. integration points                      dim1 in all input containers
          D - spatial dimension                            dim2 and dim3 in Jacobian container
        \endcode

        \param  outVals     [out] - Output array with weighted edge measures.
        \param  inJac        [in] - Input array containing cell Jacobians.
        \param  inWeights    [in] - Input integration weights.
        \param  whichEdge    [in] - Index of the edge subcell relative to the parent cell; defines the domain of integration.
        \param  parentCell   [in] - Parent cell topology.
  */
  template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  static void computeEdgeMeasure(ArrayOut                   & outVals,
                                 const ArrayJac             & inJac,
                                 const ArrayWeights         & inWeights,
                                 const int                    whichEdge,
                                 const shards::CellTopology & parentCell);

/*  template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  static void computeEdgeMeasureTemp(ArrayOut                   & outVals,
                                 const ArrayJac             & inJac,
                                 const ArrayWeights         & inWeights,
                                 const int                    whichEdge,
                                 const shards::CellTopology & parentCell);*/
  /** \brief   Multiplies fields \a <b>inVals</b> by weighted measures \a <b>inMeasure</b> and
               returns the field array \a <b>outVals</b>; this is a simple redirection to the call
               FunctionSpaceTools::scalarMultiplyDataField.

        \param  outVals     [out] - Output array with scaled field values.
        \param  inMeasure    [in] - Input array containing weighted measures.
        \param  inVals       [in] - Input fields.
  */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
  static void multiplyMeasure(ArrayTypeOut             & outVals,
                              const ArrayTypeMeasure   & inMeasure,
                              const ArrayTypeIn        & inVals);
 /*   template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
  static void multiplyMeasureTemp(ArrayTypeOut             & outVals,
                              const ArrayTypeMeasure   & inMeasure,
                              const ArrayTypeIn        & inVals);*/
  /** \brief Scalar multiplication of data and fields; please read the description below.
             
             There are two use cases:
             \li
             multiplies a rank-3, 4, or 5 container \a <b>inputFields</b> with dimensions (C,F,P),
             (C,F,P,D1) or (C,F,P,D1,D2), representing the values of a set of scalar, vector
             or tensor fields, by the values in a rank-2 container \a <b>inputData</b> indexed by (C,P),
             representing the values of scalar data, OR
             \li
             multiplies a rank-2, 3, or 4 container \a <b>inputFields</b> with dimensions (F,P),
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

  /** \brief Scalar multiplication of data and data; please read the description below.

             There are two use cases:
             \li
             multiplies a rank-2, 3, or 4 container \a <b>inputDataRight</b> with dimensions (C,P),
             (C,P,D1) or (C,P,D1,D2), representing the values of a set of scalar, vector
             or tensor data, by the values in a rank-2 container \a <b>inputDataLeft</b> indexed by (C,P),
             representing the values of scalar data, OR
             \li
             multiplies a rank-1, 2, or 3 container \a <b>inputDataRight</b> with dimensions (P),
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

  /** \brief Dot product of data and fields; please read the description below.
             
             There are two use cases:
             \li
             dot product of a rank-3, 4 or 5 container \a <b>inputFields</b> with dimensions (C,F,P)
             (C,F,P,D1) or (C,F,P,D1,D2), representing the values of a set of scalar, vector
             or tensor fields, by the values in a rank-2, 3 or 4 container \a <b>inputData</b> indexed by
             (C,P), (C,P,D1), or (C,P,D1,D2) representing the values of scalar, vector or
             tensor data, OR
             \li
             dot product of a rank-2, 3 or 4 container \a <b>inputFields</b> with dimensions (F,P),
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

  /** \brief Dot product of data and data; please read the description below.

             There are two use cases:
             \li
             dot product of a rank-2, 3 or 4 container \a <b>inputDataRight</b> with dimensions (C,P)
             (C,P,D1) or (C,P,D1,D2), representing the values of a scalar, vector or a
             tensor set of data, by the values in a rank-2, 3 or 4 container \a <b>inputDataLeft</b> indexed by
             (C,P), (C,P,D1), or (C,P,D1,D2) representing the values of scalar, vector or
             tensor data, OR
             \li
             dot product of a rank-2, 3 or 4 container \a <b>inputDataRight</b> with dimensions (P),
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

  /** \brief Cross or outer product of data and fields; please read the description below.

             There are four use cases:
             \li
             cross product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
             representing the values of a set of vector fields, on the left by the values in a rank-3
             container \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
             \li
             cross product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
             representing the values of a vector field, on the left by the values in a rank-3 container
             \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
             \li
             outer product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
             representing the values of a set of vector fields, on the left by the values in a rank-3
             container \a <b>inputData</b> indexed by (C,P,D), representing the values of vector data, OR
             \li
             outer product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
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

  /** \brief Cross or outer product of data and data; please read the description below.

             There are four use cases:
             \li
             cross product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
             representing the values of a set of vector data, on the left by the values in a rank-3
             container \a <b>inputDataLeft</b> indexed by (C,P,D) representing the values of vector data, OR
             \li
             cross product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
             representing the values of vector data, on the left by the values in a rank-3 container
             \a <b>inputDataLeft</b> indexed by (C,P,D), representing the values of vector data, OR
             \li
             outer product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
             representing the values of a set of vector data, on the left by the values in a rank-3
             container \a <b>inputDataLeft</b> indexed by (C,P,D) representing the values of vector data, OR
             \li
             outer product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
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

  /** \brief Matrix-vector or matrix-matrix product of data and fields; please read the description below.

             There are four use cases:
             \li
             matrix-vector product of a rank-4 container \a <b>inputFields</b> with dimensions (C,F,P,D),
             representing the values of a set of vector fields, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             \li
             matrix-vector product of a rank-3 container \a <b>inputFields</b> with dimensions (F,P,D),
             representing the values of a vector field, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             \li
             matrix-matrix product of a rank-5 container \a <b>inputFields</b> with dimensions (C,F,P,D,D),
             representing the values of a set of tensor fields, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputData</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             \li
             matrix-matrix product of a rank-4 container \a <b>inputFields</b> with dimensions (F,P,D,D),
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

  /** \brief Matrix-vector or matrix-matrix product of data and data; please read the description below.

             There are four use cases:
             \li
             matrix-vector product of a rank-3 container \a <b>inputDataRight</b> with dimensions (C,P,D),
             representing the values of a set of vector data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             \li
             matrix-vector product of a rank-2 container \a <b>inputDataRight</b> with dimensions (P,D),
             representing the values of vector data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             \li
             matrix-matrix product of a rank-4 container \a <b>inputDataRight</b> with dimensions (C,P,D,D),
             representing the values of a set of tensor data, on the left by the values in a rank-2, 3, or 4
             container \a <b>inputDataLeft</b> indexed by (C,P), (C,P,D) or (C,P,D,D), respectively,
             representing the values of tensor data, OR
             \li
             matrix-matrix product of a rank-3 container \a <b>inputDataRight</b> with dimensions (P,D,D),
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
                                     
/* template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  static void tensorMultiplyDataDataTemp(ArrayOutData &            outputData,
                                     const ArrayInDataLeft &   inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight,
                                     const char                transpose = 'N');
          */                           
   template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight,int outvalRank>
	struct tensorMultiplyDataDataTempSpec;
  /** \brief Applies left (row) signs, stored in the user-provided container
             <var><b>fieldSigns</b></var> and indexed by (C,L), to the operator
             <var><b>inoutOperator</b></var> indexed by (C,L,R).

             Mathematically, this method computes the matrix-matrix product
      \f[
             \mathbf{K}^{c} = \mbox{diag}(\sigma^c_0,\ldots,\sigma^c_{L-1}) \mathbf{K}^c 
      \f]
             where \f$\mathbf{K}^{c} \in \mathbf{R}^{L\times R}\f$ is array of matrices  indexed by 
             cell number \e c and stored in the rank-3 array \e inoutOperator, and 
             \f$\{\sigma^c_l\}_{l=0}^{L-1}\f$  is array of left field signs indexed by cell number \e c
             and stored in the rank-2 container \e fieldSigns;  
             see Section \ref sec_pullbacks for discussion of field signs. This operation is 
             required for operators generated by \e HCURL and \e HDIV-conforming vector-valued 
             finite element basis functions; see Sections \ref sec_pullbacks and Section 
             \ref sec_ops for applications of this method.
    
      \code
        C    - num. integration domains
        L    - num. left fields
        R    - num. right fields
      \endcode

      \param  inoutOperator [in/out] - Input / output operator array.
      \param  fieldSigns        [in] - Left field signs.
  */
  template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyLeftFieldSigns(ArrayTypeInOut        & inoutOperator,
                                  const ArrayTypeSign   & fieldSigns);

  /** \brief Applies right (column) signs, stored in the user-provided container
             <var><b>fieldSigns</b></var> and indexed by (C,R), to the operator
             <var><b>inoutOperator</b></var> indexed by (C,L,R).

             Mathematically, this method computes the matrix-matrix product
      \f[
             \mathbf{K}^{c} = \mathbf{K}^c \mbox{diag}(\sigma^c_0,\ldots,\sigma^c_{R-1})
      \f]
             where \f$\mathbf{K}^{c} \in \mathbf{R}^{L\times R}\f$ is array of matrices indexed by 
             cell number \e c and stored in the rank-3 container \e inoutOperator, and 
             \f$\{\sigma^c_r\}_{r=0}^{R-1}\f$ is array of right field signs indexed by cell number \e c
             and stored in the rank-2 container \e fieldSigns;  
             see Section \ref sec_pullbacks for discussion of field signs. This operation is 
             required for operators generated by \e HCURL and \e HDIV-conforming vector-valued 
             finite element basis functions; see Sections \ref sec_pullbacks and Section 
             \ref sec_ops for applications of this method.
    
      \code
        C    - num. integration domains
        L    - num. left fields
        R    - num. right fields
      \endcode

      \param  inoutOperator [in/out] - Input / output operator array.
      \param  fieldSigns        [in] - Right field signs.
  */
  template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyRightFieldSigns(ArrayTypeInOut        & inoutOperator,
                                   const ArrayTypeSign   & fieldSigns);

  /** \brief Applies field signs, stored in the user-provided container
             <var><b>fieldSigns</b></var> and indexed by (C,F), to the function
             <var><b>inoutFunction</b></var> indexed by (C,F), (C,F,P),
             (C,F,P,D1) or (C,F,P,D1,D2).

             Returns
      \f[    
             \mbox{inoutFunction}(c,f,*) = \mbox{fieldSigns}(c,f)*\mbox{inoutFunction}(c,f,*)
      \f]
             See Section \ref sec_pullbacks for discussion of field signs. 

      \code
        C    - num. integration domains
        F    - num. fields
        P    - num. integration points
        D1   - spatial dimension
        D2   - spatial dimension
      \endcode

      \param  inoutFunction [in/out] - Input / output function array.
      \param  fieldSigns        [in] - Right field signs.
  */
  template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyFieldSigns(ArrayTypeInOut        & inoutFunction,
                              const ArrayTypeSign   & fieldSigns);
 /* template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
  static void applyFieldSignsTemp(ArrayTypeInOut        & inoutFunction,
                              const ArrayTypeSign   & fieldSigns);
*/
  /** \brief Computes point values \a <b>outPointVals</b> of a discrete function
             specified by the basis \a <b>inFields</b> and coefficients
             \a <b>inCoeffs</b>.

             The array \a <b>inFields</b> with dimensions (C,F,P), (C,F,P,D1),
             or (C,F,P,D1,D2) represents the signed, transformed field (basis) values at
             points in REFERENCE frame; the \a <b>outPointVals</b> array with
             dimensions (C,P), (C,P,D1), or (C,P,D1,D2), respectively, represents
             values of a discrete function at points in PHYSICAL frame.
             The array \a <b>inCoeffs</b> dimensioned (C,F) supplies the coefficients
             for the field (basis) array.
   
             Returns rank-2,3 or 4 array such that
      \f[
             outPointValues(c,p,*) = \sum_{f=0}^{F-1} \sigma_{c,f} u_{c,f}(x_p)
      \f]
             where \f$\{u_{c,f}\}_{f=0}^{F-1} \f$ is scalar, vector or tensor valued finite element
             basis defined on physical cell \f$\mathcal{C}\f$ and \f$\{\sigma_{c,f}\}_{f=0}^{F-1} \f$
             are the field signs of the basis functions; see Section \ref sec_pullbacks. 
             This method implements the last step in a four step process; please see Section
             \ref sec_evaluate for details about the first three steps that prepare the 
             necessary data for this method. 

      \code
        C    - num. integration domains
        F    - num. fields
        P    - num. integration points
        D1   - spatial dimension
        D2   - spatial dimension
      \endcode

      \param  outPointVals [out] - Output point values of a discrete function.
      \param  inCoeffs      [in] - Coefficients associated with the fields (basis) array.
      \param  inFields      [in] - Field (basis) values.
  */
  template<class Scalar, class ArrayOutPointVals, class ArrayInCoeffs, class ArrayInFields>
  static void evaluate(ArrayOutPointVals     & outPointVals,
                       const ArrayInCoeffs   & inCoeffs,
                       const ArrayInFields   & inFields);
  
};  // end FunctionSpaceTools

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FunctionSpaceToolsDef.hpp>

#endif

/***************************************************************************************************
 **                                                                                               **
 **                           D O C U M E N T A T I O N   P A G E S                               **
 **                                                                                               **
 **************************************************************************************************/

/**
 \page    function_space_tools_page                 Function space tools
 
 <b>Table of contents </b>
 \li \ref sec_fst_overview 
 \li \ref sec_pullbacks
 \li \ref sec_measure
 \li \ref sec_evaluate
 
 \section sec_fst_overview                          Overview

 Intrepid::FunctionSpaceTools is a stateless class of \e expert \e methods for operations on finite
 element subspaces of \f$H(grad,\Omega)\f$, \f$H(curl,\Omega)\f$, \f$H(div,\Omega)\f$ and \f$L^2(\Omega)\f$.
 In Intrepid these spaces are referred to as \e HGRAD, \e HCURL, \e HDIV and \e HVOL. There are four 
 basic groups of methods:
 
 - Transformation methods provide implementation of pullbacks for \e HGRAD, \e HCURL, \e HDIV and \e HVOL 
   finite element functions. Thease are essentialy the "change of variables rules" needed to transform 
   values of basis functions and their derivatives defined on a reference element \f$\widehat{{\mathcal C}}\f$ 
   to a physical element \f${\mathcal C}\f$. See Section \ref sec_pullbacks for details
 - Measure computation methods implement the volume, surface and line measures required for computation
   of integrals in the physical frame by changing variables to reference frame. See Section \ref sec_measure
   for details.
 - Integration methods implement the algebraic operations to compute ubiquitous integrals of finite element 
   functions: integrals arising in bilinear forms and linear functionals.
 - Methods for algebraic and vector-algebraic operations on multi-dimensional arrays with finite element
   function values. These methods are used to prepare multidimensional arrays with data and finite
   element function values for the integration routines. They also include evaluation methods to compute
   finite element function values at some given points in physical frame; see Section \ref sec_evaluate.
 
 
 \section sec_pullbacks                             Pullbacks
 
 Notation in this section follows the standard definition of a finite element space by Ciarlet; see
 <var> The Finite Element Method for Elliptic Problems, Classics in Applied Mathematics, SIAM, 2002. </var>
 Given a reference cell \f$\{\widehat{{\mathcal C}},\widehat{P},\widehat{\Lambda}\}\f$ with a basis  
 \f$\{\widehat{u}_i\}_{i=0}^n\f$, the basis \f$\{{u}_i\}_{i=0}^n\f$ of  \f$\{{\mathcal C},P,\Lambda\}\f$ is defined 
 as follows:
 \f[
      u_i = \sigma_i \Phi^*(\widehat{u}_i), \qquad i=1,\ldots,n \,.
 \f]  
 In this formula \f$\{\sigma_i\}_{i=0}^n\f$, where \f$\sigma_i = \pm 1\f$, are the \e field \e signs, 
 and \f$\Phi^*\f$ is the \e pullback ("change of variables") transformation. For scalar spaces
 such as \e HGRAD and \e HVOL the field signs are always equal to 1 and can be disregarded. For vector
 field spaces such as \e HCURL or \e HDIV, the field sign of a basis function can be +1 or -1, 
 depending on the orientation of the physical edge or face, associated with the basis function.
  
 The actual form of the pullback depends on which one of the four function spaces \e HGRAD, \e HCURL, 
 \e HDIV and \e HVOL is being approximated and is computed as follows. Let \f$F_{\mathcal C}\f$ 
 denote the reference-to-physical map (see Section \ref sec_cell_topology_ref_map);
 \f$DF_{\mathcal C}\f$ is its Jacobian (see Section \ref sec_cell_topology_ref_map_DF) and 
 \f$J_{\mathcal C} = \det(DF_{\mathcal C})\f$. Then,
 \f[
    \begin{array}{ll}
      \Phi^*_G : HGRAD(\widehat{{\mathcal C}}) \mapsto HGRAD({\mathcal C})&
      \qquad \Phi^*_G(\widehat{u}) = \widehat{u}\circ F^{-1}_{\mathcal C} \\[2ex]
      \Phi^*_C : HCURL(\widehat{{\mathcal C}}) \mapsto HCURL({\mathcal C})&
      \qquad \Phi^*_C(\widehat{\bf u}) = \left((DF_{\mathcal C})^{-{\sf T}}\cdot\widehat{\bf u}\right)\circ F^{-1}_{\mathcal C} \\[2ex]
      \Phi^*_D : HDIV(\widehat{{\mathcal C}}) \mapsto HDIV({\mathcal C})&
      \qquad \Phi^*_D(\widehat{\bf u}) = \left(J^{-1}_{\mathcal C} DF_{\mathcal C}\cdot\widehat{\bf u}\right)\circ F^{-1}_{\mathcal C} 
      \\[2ex]
      \Phi^*_S : HVOL(\widehat{{\mathcal C}}) \mapsto HVOL({\mathcal C})&
      \qquad \Phi^*_S(\widehat{u}) = \left(J^{-1}_{\mathcal C} \widehat{u}\right) \circ F^{-1}_{\mathcal C} \,.
    \end{array}
 \f]
 Intrepid supports pullbacks only for cell topologies that have reference cells; see 
 \ref cell_topology_ref_cells.
 
 
 \section sec_measure                             Measure
 
 In Intrepid integrals of finite element functions over cells, 2-subcells (faces) and 1-subcells (edges) 
 are computed by change of variables to reference frame and require three  different kinds of measures. 
 
 -# The integral of a scalar function over a cell \f${\mathcal C}\f$
      \f[
          \int_{{\mathcal C}} f(x) dx = \int_{\widehat{{\mathcal C}}} f(F(\widehat{x})) |J | d\widehat{x}
      \f]
      requires the volume measure defined by the determinant of the Jacobian. This measure is computed 
      by Intrepid::FunctionSpaceTools::computeCellMeasure
 -# The integral of a scalar function over 2-subcell \f$\mathcal{F}\f$
      \f[
          \int_{\mathcal{F}} f(x) dx = \int_{R} f(\Phi(u,v)) 
          \left\|\frac{\partial\Phi}{\partial u}\times \frac{\partial\Phi}{\partial v}\right\| du\,dv
      \f]
      requires the surface measure defined by the norm of the vector product of the surface tangents. This   
      measure is computed by Intrepid::FunctionSpaceTools::computeFaceMeasure. In this formula \e R is the parametrization 
      domain for the 2-subcell; see Section \ref sec_cell_topology_subcell_map for details.
 -# The integral of a scalar function over a 1-subcell \f$\mathcal{E}\f$
      \f[
          \int_{\mathcal{E}} f(x) dx = \int_{R} f(\Phi(s)) \|\Phi'\| ds
      \f]
      requires the arc measure defined by the norm of the arc tangent vector. This measure is computed 
      by Intrepid::FunctionSpaceTools::computeEdgeMeasure. In this formula \e R is the parametrization 
      domain for the 1-subcell; see Section \ref sec_cell_topology_subcell_map for details.
 
 
 \section sec_evaluate                          Evaluation of finite element fields
 
 To make this example more specific, assume curl-conforming finite element spaces.
 Suppose that we have a physical cell \f$\{{\mathcal C},P,\Lambda\}\f$ with a basis
 \f$\{{\bf u}_i\}_{i=0}^n\f$. A finite element function on this cell is defined by a set of \e n
 coefficients \f$\{c_i\}_{i=0}^n\f$:
 \f[
      {\bf u}^h(x) = \sum_{i=0}^n c_i {\bf u}_i(x) \,.
 \f]
 
 From Section \ref sec_pullbacks it follows that
 \f[
      {\bf u}^h(x) = \sum_{i=0}^n c_i \sigma_i
               \left((DF_{\mathcal C})^{-{\sf T}}\cdot\widehat{\bf u}_i\right)\circ 
                     F^{-1}_{\mathcal C}(x) 
             = \sum_{i=0}^n c_i \sigma_i 
                    (DF_{\mathcal C}(\widehat{x}))^{-{\sf T}}\cdot\widehat{\bf u}_i(\widehat{x})\,,
 \f]
 where \f$ \widehat{x} = F^{-1}_{\mathcal C}(x) \in \widehat{\mathcal C} \f$ is the pre-image 
 of \e x in the reference cell. 
 
 Consequently, evaluation of finite element functions at a given set of points 
 \f$\{x_p\}_{p=0}^P \subset {\mathcal C}\f$ comprises of the following four steps:
 
 -#   Application of the inverse map \f$F^{-1}_{\mathcal C}\f$ to obtain the pre-images
      \f$\{\widehat{x}_p\}_{p=0}^P\f$ of the evaluation points in the reference cell 
      \f$\widehat{\mathcal{C}}\f$; see Intrepid::CellTools::mapToReferenceFrame
 -#   Evaluation of the appropriate reference basis set \f$\{\widehat{\bf u}_i\}_{i=1}^n\f$
      at the pre-image set \f$\{\widehat{x}_p\}_{p=0}^P\f$; see Intrepid::Basis::getValues
 -#   Application of the appropriate transformation and field signs. In our example the finite
      element space is curl-conforming and the appropriate transformation is implemented in
      Intrepid::FunctionSpaceTools::HCURLtransformVALUE. Application of the signs to the
      transformed functions is done by Intrepid::FunctionSpaceTools::applyFieldSigns.
 -#   The final step is to compute the sum of the transformed and signed basis function values
      multiplied by the coefficients of the finite element function using 
      Intrepid::FunctionSpaceTools::evaluate.
 

 Evaluation of adimssible derivatives of finite element functions is completely analogous 
 and follows the same four steps. Evaluation of scalar finite element functions is simpler
 because application of the signes can be skipped for these functions. 
 
 
 
 \section sec_ops                          Evaluation of finite element operators and functionals
 
 Assume the same setting as in Section \ref sec_evaluate. A finite element operator defined 
 by the finite element basis on the physical cell \f$\mathcal{C}\f$ is a matrix
 \f[
      \mathbf{K}^{\mathcal{C}}_{i,j} = \int_{\mathcal C} {\mathcal L}_L {\bf u}_i(x)\, {\mathcal L}_R {\bf u}_j(x) \, dx \,.
 \f]
 where \f${\mathcal L}_L\f$ and \f${\mathcal L}_R \f$ are \e left and \e right operators acting on the basis
 functions. Typically, when the left and the right basis functions are from the same finite
 element basis (as in this example), the left and right operators are the same. If they are set
 to \e VALUE we get a mass matrix; if they are set to an admissible differential operator we get
 a stiffnesss matrix. Assume again that the basis is curl-conforming and the operators are 
 set to \e VALUE. Using the basis definition from Section \ref sec_pullbacks we have that
 \f[
     \mathbf{K}^{\mathcal{C}}_{i,j} = \int_{\widehat{\mathcal C}} \sigma_i \sigma_j
     (DF_{\mathcal C}(\widehat{x}))^{-{\sf T}}\cdot\widehat{\bf u}_i(\widehat{x})\cdot
     (DF_{\mathcal C}(\widehat{x}))^{-{\sf T}}\cdot\widehat{\bf u}_i(\widehat{x})\,d\widehat{x}
 \f]
 It follows that 
 \f[
   \mathbf{K}^{\mathcal{C}}_{i,j} = 
   \mbox{diag}(\sigma_0,\ldots,\sigma_n)\widehat{\mathbf{K}}^{\mathcal{C}}\mbox{diag}(\sigma_0,\ldots,\sigma_n)
 \f]
where 
 \f[ 
   \widehat{\mathbf{K}}^{\mathcal{C}}_{i,j} = \int_{\widehat{\mathcal C}}
   (DF_{\mathcal C}(\widehat{x}))^{-{\sf T}}\cdot\widehat{\bf u}_i(\widehat{x})\cdot
   (DF_{\mathcal C}(\widehat{x}))^{-{\sf T}}\cdot\widehat{\bf u}_i(\widehat{x})\,d\widehat{x}
 \f]
 is the raw cell operator matrix. The methods Intrepid::FunctionSpaceTools::applyLeftFieldSigns and
 Intrepid::FunctionSpaceTools::applyRightFieldSigns apply the left and right diagonal sign matrices to
 the raw cell operator.

 
 A finite element operator defined by the finite element basis on the physical cell is a vector
 \f[
   \mathbf{f}^{\mathcal{C}}_{i} = \int_{\mathcal C} f(x) {\mathcal L}_R u_i(x) \, dx \,.
 \f]
 Assuming again operator \e VALUE and using the same arguments as above, we see that
 \f[
   \mathbf{f}^{\mathcal{C}} = 
      \mbox{diag}(\sigma_0,\ldots,\sigma_n)\widehat{\mathbf{f}}^{\mathcal{C}}\,,
 \f]
 where 
 \f[
   \widehat{\mathbf{f}}^{\mathcal{C}} = \int_{\widehat{\mathcal C}}
         \mathbf{f}\circ F_{\mathcal C}(\widehat{x}) 
         (DF_{\mathcal C}(\widehat{x}))^{-{\sf T}}\cdot\widehat{\bf u}_i(\widehat{x})\,d\widehat{x}
 \f]
 is the raw cell functional.
*/





























#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

