// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   ROL_BurkardtRules.hpp
    \brief  Header file for integration rules provided by John Burkardt.
	    <A HREF="http://people.sc.fsu.edu/~jburkardt/cpp_src/sandia_rules/sandia_rules.html">
	    <\A>
    \author Created by D. Kouri and D. Ridzal.
 */

#ifndef ROL_BURKARDTRULES_HPP
#define ROL_BURKARDTRULES_HPP

#include "ROL_Types.hpp"
#include "Teuchos_Assert.hpp"
//#include <fstream>
//#include <string>

namespace ROL {
 
  /** \enum  ROL::EROLBurkardt
      \brief Enumeration of integration rules provided by John Burkardt.
  */
  enum EROLBurkardt {
    BURK_CHEBYSHEV1 = 0,
    BURK_CHEBYSHEV2,
    BURK_CLENSHAWCURTIS,
    BURK_FEJER2,
    BURK_LEGENDRE,
    BURK_PATTERSON,
    BURK_TRAPEZOIDAL,
    BURK_HERMITE,
    BURK_GENZKEISTER,
    BURK_LAGUERRE,
    BURK_LAST
  };

  inline std::string EROLBurkardtToString(EROLBurkardt rule) {
    std::string retString;
    switch(rule) {
      case BURK_CHEBYSHEV1:       retString = "Gauss-Chebyshev Type 1";    break;
      case BURK_CHEBYSHEV2:       retString = "Gauss-Chebyshev Type 2";    break;
      case BURK_CLENSHAWCURTIS:   retString = "Clenshaw-Curtis";           break;
      case BURK_FEJER2:           retString = "Fejer Type 2";              break;
      case BURK_LEGENDRE:         retString = "Gauss-Legendre";            break;
      case BURK_PATTERSON:        retString = "Gauss-Patterson";           break;
      case BURK_TRAPEZOIDAL:      retString = "Trapezoidal Rule";          break;
      case BURK_HERMITE:          retString = "Gauss-Hermite";             break;
      case BURK_GENZKEISTER:      retString = "Hermite-Genz-Keister";      break;
      case BURK_LAGUERRE:         retString = "Gauss-Laguerre";            break;
      default:                    retString = "INVALID EROLBurkardt";
    }
    return retString;
  }

  inline EROLBurkardt & operator++(EROLBurkardt &type) {
    return type = static_cast<EROLBurkardt>(type+1);
  }

  inline EROLBurkardt operator++(EROLBurkardt &type, int) {
    EROLBurkardt oldval = type;
    ++type;
    return oldval;
  }
  
  /** \class ROL::ROLBurkardtRules
      \brief Providing integration rules, 
             created by John Burkardt, Scientific Computing, 
             Florida State University,
             modified and redistributed by D. Kouri.

             See <A HREF="http://people.sc.fsu.edu/~jburkardt/cpp_src/sandia_rules/sandia_rules.html">
	     <\A>
  */
  class ROLBurkardtRules {
      
  public:

    /* HELPER FUNCTIONS */
    template<class Real> 
    static void imtqlx ( int n, Real d[], Real e[], Real z[] );
    template<class Real> static Real r8_epsilon( Real one );
    template<class Real> static Real r8_sign( Real x );

    /* COMPUTE CHEBYSHEV TYPE 1 NODES AND WEIGHTS                      */
    /* Integrates functions on [-1,1] weighted by w(x) = 1/sqrt(1-x^2) */ 
    /** \brief  Gauss-Chebyshev of Type 1; returns points and weights.
    */
    template<class Real> 
    static void chebyshev1_compute ( int order, Real x[], Real w[] );
    /** \brief  Gauss-Chebyshev of Type 1; returns points.
    */
    template<class Real> 
    static void chebyshev1_compute_points ( int order, Real x[] );
    /** \brief  Gauss-Chebyshev of Type 1; returns weights.
    */
    template<class Real> 
    static void chebyshev1_compute_weights ( int order, Real w[] );

    /* COMPUTE CHEBYSHEV TYPE 2 NODES AND WEIGHTS                      */
    /* Integrates functions on [-1,1] weighted by w(x) = sqrt(1-x^2)   */
    /** \brief  Gauss-Chebyshev of Type 2; returns points and weights.
    */
    template<class Real> 
    static void chebyshev2_compute ( int order, Real x[], Real w[] );
    /** \brief  Gauss-Chebyshev of Type 2; returns points.
    */
    template<class Real> 
    static void chebyshev2_compute_points ( int order, Real x[] );
    /** \brief  Gauss-Chebyshev of Type 2; returns weights.
    */
    template<class Real> 
    static void chebyshev2_compute_weights ( int order, Real w[] );
    
    /* COMPUTE CLENSHAW CURTIS NODES AND WEIGHTS                       */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */    
    /** \brief  Clenshaw-Curtis; returns points and weights.
    */
    template<class Real> 
    static void clenshaw_curtis_compute ( int order, Real x[], Real w[] );
    /** \brief  Clenshaw-Curtis; returns points.
    */    
    template<class Real> 
    static void clenshaw_curtis_compute_points ( int order, Real x[] );
    /** \brief  Clenshaw-Curtis; returns weights.
    */    
    template<class Real> 
    static void clenshaw_curtis_compute_weights ( int order, Real w[] );
    
    /* COMPUTE FEJER TYPE 2 NODES AND WEIGHTS                          */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */    
    /** \brief  Fejer type 2; returns points and weights.
    */
    template<class Real> 
    static void fejer2_compute ( int order, Real x[], Real w[] );       
    /** \brief  Fejer type 2; returns points.
    */ 
    template<class Real> 
    static void fejer2_compute_points ( int order, Real x[] );       
    /** \brief  Fejer type 2; returns weights.
    */ 
    template<class Real> 
    static void fejer2_compute_weights ( int order, Real w[] );
    
    /* COMPUTE GAUSS HERMITE NODES AND WEIGHTS                         */
    /* Integrates functions on (-oo,oo) weighted by w(x) = exp(-x^2)   */    
    /** \brief  Gauss-Hermite; returns points and weights.
    */
    template<class Real> 
    static void hermite_compute ( int order, Real x[], Real w[] );    
    /** \brief  Gauss-Hermite; returns points.
    */
    template<class Real> 
    static void hermite_compute_points ( int order, Real x[] );   
    /** \brief  Gauss-Hermite; returns weights.
    */
    template<class Real> 
    static void hermite_compute_weights ( int order, Real w[] );
    
    /** \brief  Gauss-Hermite; returns points and weights.
    */
    template<class Real> 
    static void hermite_lookup ( int n, Real x[], Real w[] );    
    /** \brief  Gauss-Hermite; returns points.
    */
    template<class Real> 
    static void hermite_lookup_points ( int n, Real x[] );    
    /** \brief  Gauss-Hermite; returns weights.
    */
    template<class Real> 
    static void hermite_lookup_weights ( int n, Real w[] );

    /* COMPUTE GENZ KEISTER NODES AND WEIGHTS                          */
    /* Integrates functions on (-oo,oo) weighted by w(x) = exp(-x^2)   */    
    /** \brief  Hermite-Genz-Keister; returns points and weights.
    */
    template<class Real> 
    static void hermite_genz_keister_lookup ( int n, Real x[], Real w[] );
    /** \brief  Hermite-Genz-Keister; returns points.
    */
    template<class Real> 
    static void hermite_genz_keister_lookup_points ( int n, Real x[] ); 
    /** \brief  Hermite-Genz-Keister; returns weights.
    */
    template<class Real> 
    static void hermite_genz_keister_lookup_weights ( int n, Real w[] );

    /* COMPUTE GAUSS LAGUERRE NODES AND WEIGHTS                        */
    /* Integrates functons on [0,oo) weighted by w(x) = exp(-x)        */   
    /** \brief  Gauss-Laguerre; returns points and weights.
    */
    template<class Real> 
    static void laguerre_compute ( int n, Real x[], Real w[] );   
    /** \brief  Gauss-Laguerre; returns points.
    */
    template<class Real> 
    static void laguerre_compute_points ( int order, Real x[] ); 
    /** \brief  Gauss-Laguerre; returns weights.
    */
    template<class Real> 
    static void laguerre_compute_weights ( int order, Real w[] );
    
    /** \brief  Gauss-Laguerre; returns points and weights.
    */
    template<class Real> 
    static void laguerre_lookup ( int n, Real x[], Real w[] );    
    /** \brief  Gauss-Laguerre; returns points.
    */
    template<class Real> 
    static void laguerre_lookup_points ( int n, Real x[] );   
    /** \brief  Gauss-Laguerre; returns weights.
    */
    template<class Real> 
    static void laguerre_lookup_weights ( int n, Real w[] );

    /* COMPUTE GAUSS LEGENDRE NODES AND WEIGHTS                        */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */    
    /** \brief  Gauss-Legendre; returns points and weights.
    */
    template<class Real> 
    static void legendre_compute ( int n, Real x[], Real w[] );   
    /** \brief  Gauss-Legendre; returns points.
    */
    template<class Real> 
    static void legendre_compute_points ( int order, Real x[] );  
    /** \brief  Gauss-Legendre; returns weights.
    */
    template<class Real> 
    static void legendre_compute_weights ( int order, Real w[] );
    
    /** \brief  Gauss-Legendre; returns points and weights.
    */
    template<class Real> 
    static void legendre_lookup ( int n, Real x[], Real w[] );    
    /** \brief  Gauss-Legendre; returns points.
    */
    template<class Real> 
    static void legendre_lookup_points ( int n, Real x[] ); 
    /** \brief  Gauss-Legendre; returns weights.
    */
    template<class Real> 
    static void legendre_lookup_weights ( int n, Real w[] );

    /* COMPUTE GAUSS PATTERSON NODES AND WEIGHTS                       */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */  
    /** \brief  Gauss-Patterson; returns points and weights.
    */
    template<class Real> 
    static void patterson_lookup ( int n, Real x[], Real w[] );  
    /** \brief  Gauss-Patterson; returns points.
    */
    template<class Real> 
    static void patterson_lookup_points ( int n, Real x[] );   
    /** \brief  Gauss-Patterson; returns weights.
    */
    template<class Real> 
    static void patterson_lookup_weights ( int n, Real w[] );

    /* COMPUTE TRAPEZOIDAL RULE NODES AND WEIGHTS                      */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */  
    /** \brief  Trapezoidal rule; returns points and weights.
    */
    template<class Real> 
    static void trapezoidal_compute ( int n, Real x[], Real w[] ); 
    /** \brief  Trapezoidal rule; returns points.
    */
    template<class Real> 
    static void trapezoidal_compute_points ( int order, Real x[] );  
    /** \brief  Trapezoidal rule; returns weights.
    */
    template<class Real> 
    static void trapezoidal_compute_weights ( int order, Real w[] );
       
  }; // class ROLBurkardtRules

} // end ROL namespace

// include templated definitions
#include "ROL_BurkardtRulesDef.hpp"
 
#endif
