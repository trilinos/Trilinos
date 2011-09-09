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

/** \file   Intrepid_CubatureLineSorted.hpp
    \brief  Header file for the Intrepid::CubatureLineSorted class.
    \author Created by D. Kouri and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_LINESORTED_HPP
#define INTREPID_CUBATURE_LINESORTED_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_BurkardtRules.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace Intrepid {

enum EIntrepidGrowth {
  GROWTH_DEFAULT = 0,
  GROWTH_SLOWLIN,
  GROWTH_SLOWLINODD,
  GROWTH_MODLIN,
  GROWTH_SLOWEXP,
  GROWTH_MODEXP,
  GROWTH_FULLEXP
};

/** \class Intrepid::CubatureLineSorted
    \brief Utilizes cubature (integration) rules contained in the library sandia_rules
           (John Burkardt, Scientific Computing, Florida State University) within Intrepid.
	   
	   This class contains ten rules:
	   \li Gauss-Chebyshev Type 1 - integration domain is [-1,+1] with weight w(x) = 1/sqrt(1-x^2)
	   \li Gauss-Chebyshev Type 2 - integration domain is [-1,+1] with weight w(x) = sqrt(1-x^2)
	   \li Clenshaw-Curtis        - integration domain is [-1,+1] with weight w(x) = 1
	   \li Fejer Type 2           - integration domain is [-1,+1] with weight w(x) = 1
	   \li Gauss-Legendre         - integration domain is [-1,+1] with weight w(x) = 1
	   \li Gauss-Patterson        - integration domain is [-1,+1] with weight w(x) = 1
	   \li Trapezoidal            - integration domain is [-1,+1] with weight w(x) = 1
	   \li Gauss-Hermite          - integration domain is (-oo,+oo) with weight w(x) = exp(-x^2)
	   \li Hermite-Genz-Keister   - integration domain is (-oo,+oo) with weight w(x) = exp(-x^2)
	   \li Gauss-Laguerre         - integration domain is [0,+oo) with weight w(x) = exp(-x)
*/

template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint> 
class CubatureLineSorted : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight> {

private: 
  
  /** \brief Contains points of this cubature rule.
  */
  std::map<Scalar,int> points_; // keys = nodes, map values = location of weights

  /** \brief Contains points of this cubature rule.
  */
  std::vector<Scalar> weights_;

  /** \brief Contains the number of nodes for this cubature rule.
  */
  int numPoints_; 
  
  /** \brief The degree of polynomials that are integrated
             exactly by this cubature rule.
  */
  int degree_;

  /** \brief Type of integration points.
  */
  EIntrepidBurkardt rule_type_;
 
  /** \brief Cubature name.
  */
  static const char *cubature_name_;

public:

  ~CubatureLineSorted() {}
  
  /** \brief Constructor.

      \param degree       [in]   - The degree of polynomials that are integrated
                                   exactly by this cubature rule. Default: 0.
  */
  CubatureLineSorted(int degree = 0, EIntrepidBurkardt rule = BURK_CLENSHAWCURTIS, bool isNormalized = false);

  /** \brief Constructor.

      \param numPoints    [in]   - The number of cubature points. Default: 0.
  */
  CubatureLineSorted(EIntrepidBurkardt rule = BURK_CLENSHAWCURTIS, int numPoints = 0, bool isNormalized = false);

  CubatureLineSorted(std::vector<Scalar> & points, std::vector<Scalar> & weights);

  // Access Operators  
  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints   [out]   - Array containing the cubature points.
      \param cubWeights  [out]   - Array of corresponding cubature weights.
  */
  void getCubature(ArrayPoint  & cubPoints,
		   ArrayWeight & cubWeights) const;

  /** \brief Returns the number of cubature points.
  */
  int getNumPoints() const;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has size 1.
  */
  void getAccuracy(std::vector<int> & accuracy) const;

  /** \brief Returns dimension of domain of integration.
  */  
  int getDimension() const;

  /** \brief Returns cubature name.
  */
  const char* getName() const;

  /** \brief Get a specific node described by the iterator location.
  */
  Scalar getNode(typename std::map<Scalar,int>::iterator it);

  /** \brief Get a specific weight described by the integer location.
  */
  Scalar getWeight(int weight);

  /** \brief Get a specific weight described by the corresponding node.
  */
  Scalar getWeight(Scalar point);

  /** \brief Initiate iterator at the beginning of data.
  */
  typename std::map<Scalar,int>::iterator begin(void);

  /** \brief Initiate iterator at the end of data.
  */
  typename std::map<Scalar,int>::iterator end(void);

  /** \brief Replace CubatureLineSorted values with 
             "this = alpha1*this+alpha2*cubRule2".
  */
  void update(Scalar alpha2, CubatureLineSorted<Scalar> & cubRule2, Scalar alpha1); 
};

int growthRule1D(int index, EIntrepidGrowth growth, EIntrepidBurkardt rule);

} // Intrepid Namespace

// include templated definitions
#include <Intrepid_CubatureLineSortedDef.hpp>

#endif
