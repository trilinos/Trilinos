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

/** \file   Intrepid_BurkardtRules.hpp
    \brief  Header file for integration rules provided by John Burkardt.
	    <A HREF="http://people.sc.fsu.edu/~jburkardt/cpp_src/sandia_rules/sandia_rules.html">
	    <\A>
    \author Created by D. Kouri and D. Ridzal.
 */

#ifndef INTREPID_BURKARDTRULES_HPP
#define INTREPID_BURKARDTRULES_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_TestForException.hpp"
//#include <fstream>
//#include <string>

namespace Intrepid {
 
  /** \enum  Intrepid::EIntrepidBurkardt
      \brief Enumeration of integration rules provided by John Burkardt.
  */
  enum  EIntrepidBurkardt {
    BURK_CHEBYSHEV1 = 0,
    BURK_CHEBYSHEV2,
    BURK_CLENSHAWCURTIS,
    BURK_FEJER2,
    BURK_LEGENDRE,
    BURK_PATTERSON,
    BURK_TRAPEZOIDAL,
    BURK_HERMITE,
    BURK_GENZKEISTER,
    BURK_LAGUERRE
  };

  inline std::string EIntrepidBurkardtToString(EIntrepidBurkardt rule) {
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
      default:                    retString = "INVALID EIntrepidBurkardt";
    }
    return retString;
  }

  inline EIntrepidBurkardt & operator++(EIntrepidBurkardt &type) {
    return type = static_cast<EIntrepidBurkardt>(type+1);
  }

  inline EIntrepidBurkardt operator++(EIntrepidBurkardt &type, int) {
    EIntrepidBurkardt oldval = type;
    ++type;
    return oldval;
  }
  
  /** \class Intrepid::IntrepidBurkardtRules
      \brief Providing integration rules, 
             created by John Burkardt, Scientific Computing, 
             Florida State University,
             modified and redistributed by D. Kouri.

             See <A HREF="http://people.sc.fsu.edu/~jburkardt/cpp_src/sandia_rules/sandia_rules.html">
	     <\A>
  */
  class IntrepidBurkardtRules {
      
  public:

    /* HELPER FUNCTIONS */
    template<class Scalar> 
    static void imtqlx ( int n, Scalar d[], Scalar e[], Scalar z[] );
    template<class Scalar> static Scalar r8_epsilon( Scalar one );
    template<class Scalar> static Scalar r8_sign( Scalar x );

    /* COMPUTE CHEBYSHEV TYPE 1 NODES AND WEIGHTS                      */
    /* Integrates functions on [-1,1] weighted by w(x) = 1/sqrt(1-x^2) */ 
    /** \brief  Gauss-Chebyshev of Type 1; returns points and weights.
    */
    template<class Scalar> 
    static void chebyshev1_compute ( int order, Scalar x[], Scalar w[] );
    /** \brief  Gauss-Chebyshev of Type 1; returns points.
    */
    template<class Scalar> 
    static void chebyshev1_compute_points ( int order, Scalar x[] );
    /** \brief  Gauss-Chebyshev of Type 1; returns weights.
    */
    template<class Scalar> 
    static void chebyshev1_compute_weights ( int order, Scalar w[] );

    /* COMPUTE CHEBYSHEV TYPE 2 NODES AND WEIGHTS                      */
    /* Integrates functions on [-1,1] weighted by w(x) = sqrt(1-x^2)   */
    /** \brief  Gauss-Chebyshev of Type 2; returns points and weights.
    */
    template<class Scalar> 
    static void chebyshev2_compute ( int order, Scalar x[], Scalar w[] );
    /** \brief  Gauss-Chebyshev of Type 2; returns points.
    */
    template<class Scalar> 
    static void chebyshev2_compute_points ( int order, Scalar x[] );
    /** \brief  Gauss-Chebyshev of Type 2; returns weights.
    */
    template<class Scalar> 
    static void chebyshev2_compute_weights ( int order, Scalar w[] );
    
    /* COMPUTE CLENSHAW CURTIS NODES AND WEIGHTS                       */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */    
    /** \brief  Clenshaw-Curtis; returns points and weights.
    */
    template<class Scalar> 
    static void clenshaw_curtis_compute ( int order, Scalar x[], Scalar w[] );
    /** \brief  Clenshaw-Curtis; returns points.
    */    
    template<class Scalar> 
    static void clenshaw_curtis_compute_points ( int order, Scalar x[] );
    /** \brief  Clenshaw-Curtis; returns weights.
    */    
    template<class Scalar> 
    static void clenshaw_curtis_compute_weights ( int order, Scalar w[] );
    
    /* COMPUTE FEJER TYPE 2 NODES AND WEIGHTS                          */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */    
    /** \brief  Fejer type 2; returns points and weights.
    */
    template<class Scalar> 
    static void fejer2_compute ( int order, Scalar x[], Scalar w[] );       
    /** \brief  Fejer type 2; returns points.
    */ 
    template<class Scalar> 
    static void fejer2_compute_points ( int order, Scalar x[] );       
    /** \brief  Fejer type 2; returns weights.
    */ 
    template<class Scalar> 
    static void fejer2_compute_weights ( int order, Scalar w[] );
    
    /* COMPUTE GAUSS HERMITE NODES AND WEIGHTS                         */
    /* Integrates functions on (-oo,oo) weighted by w(x) = exp(-x^2)   */    
    /** \brief  Gauss-Hermite; returns points and weights.
    */
    template<class Scalar> 
    static void hermite_compute ( int order, Scalar x[], Scalar w[] );    
    /** \brief  Gauss-Hermite; returns points.
    */
    template<class Scalar> 
    static void hermite_compute_points ( int order, Scalar x[] );   
    /** \brief  Gauss-Hermite; returns weights.
    */
    template<class Scalar> 
    static void hermite_compute_weights ( int order, Scalar w[] );
    
    /** \brief  Gauss-Hermite; returns points and weights.
    */
    template<class Scalar> 
    static void hermite_lookup ( int n, Scalar x[], Scalar w[] );    
    /** \brief  Gauss-Hermite; returns points.
    */
    template<class Scalar> 
    static void hermite_lookup_points ( int n, Scalar x[] );    
    /** \brief  Gauss-Hermite; returns weights.
    */
    template<class Scalar> 
    static void hermite_lookup_weights ( int n, Scalar w[] );

    /* COMPUTE GENZ KEISTER NODES AND WEIGHTS                          */
    /* Integrates functions on (-oo,oo) weighted by w(x) = exp(-x^2)   */    
    /** \brief  Hermite-Genz-Keister; returns points and weights.
    */
    template<class Scalar> 
    static void hermite_genz_keister_lookup ( int n, Scalar x[], Scalar w[] );
    /** \brief  Hermite-Genz-Keister; returns points.
    */
    template<class Scalar> 
    static void hermite_genz_keister_lookup_points ( int n, Scalar x[] ); 
    /** \brief  Hermite-Genz-Keister; returns weights.
    */
    template<class Scalar> 
    static void hermite_genz_keister_lookup_weights ( int n, Scalar w[] );

    /* COMPUTE GAUSS LAGUERRE NODES AND WEIGHTS                        */
    /* Integrates functons on [0,oo) weighted by w(x) = exp(-x)        */   
    /** \brief  Gauss-Laguerre; returns points and weights.
    */
    template<class Scalar> 
    static void laguerre_compute ( int n, Scalar x[], Scalar w[] );   
    /** \brief  Gauss-Laguerre; returns points.
    */
    template<class Scalar> 
    static void laguerre_compute_points ( int order, Scalar x[] ); 
    /** \brief  Gauss-Laguerre; returns weights.
    */
    template<class Scalar> 
    static void laguerre_compute_weights ( int order, Scalar w[] );
    
    /** \brief  Gauss-Laguerre; returns points and weights.
    */
    template<class Scalar> 
    static void laguerre_lookup ( int n, Scalar x[], Scalar w[] );    
    /** \brief  Gauss-Laguerre; returns points.
    */
    template<class Scalar> 
    static void laguerre_lookup_points ( int n, Scalar x[] );   
    /** \brief  Gauss-Laguerre; returns weights.
    */
    template<class Scalar> 
    static void laguerre_lookup_weights ( int n, Scalar w[] );

    /* COMPUTE GAUSS LEGENDRE NODES AND WEIGHTS                        */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */    
    /** \brief  Gauss-Legendre; returns points and weights.
    */
    template<class Scalar> 
    static void legendre_compute ( int n, Scalar x[], Scalar w[] );   
    /** \brief  Gauss-Legendre; returns points.
    */
    template<class Scalar> 
    static void legendre_compute_points ( int order, Scalar x[] );  
    /** \brief  Gauss-Legendre; returns weights.
    */
    template<class Scalar> 
    static void legendre_compute_weights ( int order, Scalar w[] );
    
    /** \brief  Gauss-Legendre; returns points and weights.
    */
    template<class Scalar> 
    static void legendre_lookup ( int n, Scalar x[], Scalar w[] );    
    /** \brief  Gauss-Legendre; returns points.
    */
    template<class Scalar> 
    static void legendre_lookup_points ( int n, Scalar x[] ); 
    /** \brief  Gauss-Legendre; returns weights.
    */
    template<class Scalar> 
    static void legendre_lookup_weights ( int n, Scalar w[] );

    /* COMPUTE GAUSS PATTERSON NODES AND WEIGHTS                       */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */  
    /** \brief  Gauss-Patterson; returns points and weights.
    */
    template<class Scalar> 
    static void patterson_lookup ( int n, Scalar x[], Scalar w[] );  
    /** \brief  Gauss-Patterson; returns points.
    */
    template<class Scalar> 
    static void patterson_lookup_points ( int n, Scalar x[] );   
    /** \brief  Gauss-Patterson; returns weights.
    */
    template<class Scalar> 
    static void patterson_lookup_weights ( int n, Scalar w[] );

    /* COMPUTE TRAPEZOIDAL RULE NODES AND WEIGHTS                      */
    /* Integrates functions on [-1,1] weighted by w(x) = 1             */  
    /** \brief  Trapezoidal rule; returns points and weights.
    */
    template<class Scalar> 
    static void trapezoidal_compute ( int n, Scalar x[], Scalar w[] ); 
    /** \brief  Trapezoidal rule; returns points.
    */
    template<class Scalar> 
    static void trapezoidal_compute_points ( int order, Scalar x[] );  
    /** \brief  Trapezoidal rule; returns weights.
    */
    template<class Scalar> 
    static void trapezoidal_compute_weights ( int order, Scalar w[] );
       
  }; // class IntrepidBurkardtRules

} // end Intrepid namespace

// include templated definitions
#include "Intrepid_BurkardtRulesDef.hpp"
 
#endif
