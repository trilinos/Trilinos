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

/** \file   Intrepid_AdaptiveSparseGridInterface.hpp
    \brief  Header file for the Intrepid::AdaptiveSparseGridInterface class.
    \author Created by D. Kouri and D. Ridzal.
*/

#ifndef INTREPID_ADAPTIVESPARSEGRIDINTERFACE_HPP
#define INTREPID_ADAPTIVESPARSEGRIDINTERFACE_HPP

#include "Intrepid_CubatureTensorSorted.hpp"

namespace Intrepid {

/** \class Intrepid::AdaptiveSparseGrid
    \brief Builds general adaptive sparse grid rules (Gerstner and Griebel) 
           using the 1D cubature rules in the Intrepid::CubatureLineSorted 
	   class.
*/

template<class Scalar, class UserVector>
class AdaptiveSparseGridInterface {
protected:  
  /** \brief The dimension of the integration domain.
  */
  int    dimension_;
 
  /** \brief The maximum allowable level of quadrature.
  */
  int    maxLevel_;

  /** \brief The initial contribution to the integral.
  */
  Scalar initialDiff_;

  /** \brief Whether or not to normalize the weights.
  */
  bool   isNormalized_;

  /** \brief The user defined 1D cubature rules.
  */
  std::vector<EIntrepidBurkardt> rule1D_;
  
  /** \brief The user defined 1D growth rules.
  */
  std::vector<EIntrepidGrowth>   growth1D_;

public:  
  /** \brief Destructor.
  */
  virtual ~AdaptiveSparseGridInterface() {};

  /** \brief Constructor starts with index [1,...,1].     

      \param dimension     [in]  - Dimension of integration domain.
      \param rule1D        [in]  - 1D cubature rule names.
      \param growth1D      [in]  - 1D cubature growth rules.
      \param maxlevel      [in]  - Maximum level of sparse grid.
      \param isNormalized  [in]  - Flag whether to normalize cubature weights.
  */
  AdaptiveSparseGridInterface(
		     int dimension,
		     std::vector<EIntrepidBurkardt> rule1D,
		     std::vector<EIntrepidGrowth> growth1D,
		     int maxLevel,
                     bool isNormalized);

  /** \brief Evaluate the integrand function.
      
      \param output   [out] - Output of integrand evaluation.
      \param input    [in]  - Evaluation points.
  */
  virtual void eval_integrand(
		     UserVector & output, 
		     std::vector<Scalar> & input) = 0;

  /** \brief Evaluate the cubature rule.
            
      \param output   [out] - Output of cubature evaluation.
      \param cubRule  [in]  - Cubature rule.
  */
  virtual void eval_cubature(
		     UserVector & output, 
		     CubatureTensorSorted<Scalar> & cubRule);

  /** \brief User defined error indicator function. 
           
      \param input [in] - Reduction of high dimensional integral value to a Scalar.
  */
  virtual Scalar error_indicator(
		     UserVector & input) = 0;

  /** \brief User defined test for maximum level of cubature.

      \param index [in] - Multi-index of cubature levels.
  */
  virtual bool max_level(
		     std::vector<int> index);

  /** \brief Compute initial quantities for sparse grid adaptation
      \param output [out] - Vector of outputs.
  */
  void init(UserVector & output);

  /** \brief Return user defined 1D quadrature rules.
      \param rule1D [out] - 1D quadrature rules names.
  */
  void getRule(std::vector<EIntrepidBurkardt> & rule1D);

  /** \brief Return user defined 1D growth rules.
      \param growth1D [out] - 1D quadrature growth rules.
  */
  void getGrowth(std::vector<EIntrepidGrowth> & growth1D);

  /** \brief Return dimension of integration domain.
  */
  int getDimension();
 
  /** \brief Return initial error indicator.
  */
  Scalar getInitialDiff();

  /** \brief Return whether or not cubature weights are normalized.
  */
  bool isNormalized();
};

} // End Intrepid namespace

// include templated definitions
#include <Intrepid_AdaptiveSparseGridInterfaceDef.hpp>

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

