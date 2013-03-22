// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_POLYNOMIAL_DECL_HPP
#define TEUCHOS_POLYNOMIAL_DECL_HPP

#include "Teuchos_Describable.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_PolynomialTraits.hpp"

namespace Teuchos {

  //! Lightweight container class to represent a simple polynomial.
  /*!
   * This class represents a simple polynomial of the form:
   * \f[
   *     x(t) = \sum_{i=0}^d x_i t^i
   * \f]
   * where \f$d\f$ is the degree of the polynomial and \f$t\f$ is a scalar.
   * The coefficients \f$x_i\f$, \f$i=0,\dots,d\f$ can be scalars, vectors,
   * operators, etc., any type that can implement Teuchos::PolynomialTraits.
   * A template specialization of Teuchos::PolynomialTraits must be provided
   * for the coefficient type, however Teuchos::PolynomailTraits does provide 
   * a default template definition for scalars.
   *
   * This class provides methods for creating a polynomial of some degree, 
   * setting and retreiving coefficients, and evaluating the polynomial and
   * its derivative at some value \f$t\f$.
   */
  template <typename CoeffT>
  class Polynomial : virtual public Teuchos::Describable {
  public:

    //! Typename of coefficients
    typedef CoeffT coeff_type;

    //! Typename of scalars
    typedef typename Teuchos::PolynomialTraits<coeff_type>::scalar_type scalar_type;

    //! Create a polynomial of degree \c deg
    /*!
     * If \c reserve > \c deg, a polynomial of degree \c deg will be created,
     * but space for \c reserve + 1 coefficients will be created to allow 
     * future increases in the degree of the polynomial to be more efficient.
     * A clone of \c cloneCoeff will be placed at each coefficient.
     */
    Polynomial(unsigned int deg, const CoeffT& cloneCoeff, 
	       unsigned int reserve = 0);

    //! Create a polynomial of degree \c deg without cloning.  DANGEROUS!
    /*!
     * In this version of the constructor, no clone object is provided,
     * and therefore no storage will be created for each coefficient.
     * In this case, setCoefficientPtr() below should be used to set
     * each coefficient pointer to a valid cofficient.  This constructor exists
     * to be able to create an efficient "view" of another polynomial.
     */
    Polynomial(unsigned int deg, unsigned int reserve = 0);

    //! Destructor
    ~Polynomial();

    //! Return degree of polynomial
    unsigned int degree() const { return d; }

    //! Set degree of polynomial to \c deg
    void setDegree(unsigned int deg);

    //! Return ref-count pointer to coefficient \c i
    Teuchos::RCP<CoeffT>
    getCoefficient(unsigned int i);

    //! Return ref-count pointer to constant coefficient \c i
    Teuchos::RCP<const CoeffT>
    getCoefficient(unsigned int i) const;

    //! Set coefficient \c i to \c c
    void setCoefficient(unsigned int i, const CoeffT& v);

    //! Set pointer for coefficient \c i to \c c_ptr.  DANGEROUS!
    /*!
     * Directly set the coefficient pointer to c_ptr.  This method should
     * be used with care since future calls to setCoefficient(i,c) will
     * also modify the coefficient pointed to \c c_ptr.  However, in certain
     * situations it is necessary to do this for efficiency.
     */
    void setCoefficientPtr(unsigned int i, 
			   const Teuchos::RCP<CoeffT>& c_ptr);

    //! Evaluate polynomial and possibly its derivative at time \c t
    /*!
     * The value of the polynomial at \c t is computed and stored in \c *x. 
     * If \c xdot is not \c NULL, the derivative with respect to t is 
     * evaluated and stored in \c *xdot.
     *
     * Horner's method is used to efficiently evaluate the polynomial
     * and its derivative.
     */
    void evaluate(typename Teuchos::Polynomial<CoeffT>::scalar_type  t, 
		  CoeffT* x, CoeffT* xdot = NULL) const;

  private:

    //! Prohibit copying
    Polynomial(const Polynomial&);

    //! Prohibit copying
    Polynomial& operator=(const Polynomial&);

  protected:

    //! Degree of polynomial
    unsigned int d; 

    //! Size of polynomial (may be > d)
    unsigned int sz;

    //! Vector of polynomial coefficients
    /*!
     * \c coeff[i] corresponds to the degree \c i term, \c i=0,...,d
     */
    std::vector< Teuchos::RCP<CoeffT> > coeff;

  }; // class Polynomial

} // end namespace Teuchos

#endif  // TEUCHOS_POLYNOMIAL_DECL_HPP
