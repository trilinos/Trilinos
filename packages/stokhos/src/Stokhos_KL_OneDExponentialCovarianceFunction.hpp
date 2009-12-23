// $Id: Stokhos_Quadrature.hpp,v 1.4 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_Quadrature.hpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_KL_ONE_D_EXPONENTIAL_COVARIANCE_FUNCTION_HPP
#define STOKHOS_KL_ONE_D_EXPONENTIAL_COVARIANCE_FUNCTION_HPP

#include <string>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Stokhos_KL_OneDEigenPair.hpp"

namespace Stokhos {

  namespace KL {

    /*! 
     * \brief Class representing an exponential covariance function and its %KL 
     * eigevalues/eigenfunctions.
     */
    /*!
     * This class provides the exponential covariance function 
     * \f[
     *     \mbox{cov}(x,x') = \exp(-|x-x'|/L).
     * \f]
     * The corresponding eigenfunctions can be shown to be
     * \f$A_n\cos(\omega_n(x-\beta))\f$ and 
     * \f$B_n\sin(\omega^\ast_n(x-\beta))\f$ for
     * \f$x\in[a,b]\f$ where \f$\omega_n\f$ and \f$\omega^\ast_n\f$
     * satisfy
     * \f[
     *     1 - L\omega_n\tan(\omega_n\alpha) = 0
     * \f]
     * and
     * \f[
     *     L\omega^\ast_n + \tan(\omega^\ast_n\alpha) = 0
     * \f]
     * respectively, where \f$\alpha=(b-a)/2\f$ and \f$\beta=(b+a)/2\f$.  Then
     * \f[
     *  A_n = \frac{1}{\left(\int_a^b\cos^2(\omega_n(x-\beta)) dx\right)^{1/2}} 
     *      = \frac{1}{\sqrt{\alpha + \frac{\sin(2\omega_n\alpha)}{2\omega_n}}},
     * \f]
     * \f[
     *  B_n 
     *    = \frac{1}{\left(\int_a^b\sin^2(\omega^\ast_n(x-\beta)) dx\right)^{1/2}} 
     *    = \frac{1}{\sqrt{\alpha - \frac{\sin(2\omega_n^\ast\alpha)}{2\omega^\ast_n}}}
     * \f]
     * and the corresponding eigenvalues are given by
     * \f[
     *     \lambda_n = \frac{2L}{(L\omega_n)^2 + 1}
     * \f]
     * and
     * \f[
     *     \lambda^\ast_n = \frac{2L}{(L\omega^\ast_n)^2 + 1}.
     * \f]
     * It is straightforward to show that for each \f$n\f$, 
     * \f$\omega^\ast_n < \omega_n\f$, and thus 
     * \f$\lambda_n < \lambda^\ast_n\f$.  Hence when sorted on decreasing
     * eigenvalue, the eigenfunctions alternate between cosine and sine.
     * See "Stochastic Finite Elements" by Ghanem and Spanos for a complete 
     * description of how to derive these relationships.
     *
     * For a given value of \f$M\f$, the code works by computing the \f$M\f$
     * eigenfunctions using a bisection root solver to compute the frequencies 
     * \f$\omega_n\f$ and \f$\omega^\ast_n\f$.  
     *
     * Data for the root solver is passed through a Teuchos::ParameterList,
     * which accepts the following parameters:
     * <ul>
     *   <li> "Bound Peturbation Size" -- [value_type] (1.0e-6)
     *        Perturbation away from \f$i\pi/2\f$ for bounding the 
     *        frequencies \f$\omega\f$ in the bisection algorithm
     *   <li> "Nonlinear Solver Tolerance" -- [value_type] (1.0e-10)
     *        Tolerance for bisection nonlinear solver for computing frequencies
     *   <li> "Maximum Nonlinear Solver Iterations" -- [int] (100)
     *        Maximum number of nonlinear solver iterations for computing
     *        the frequencies.
     * </ul>
     */
    template <typename value_type>
    class OneDExponentialCovarianceFunction {
    public:

      //! Constructor
      OneDExponentialCovarianceFunction(int M, 
					const value_type& a,
					const value_type& b,
					const value_type& L,
					const std::string& dim_name,
					Teuchos::ParameterList& solverParams);

      //! Destructor
      ~OneDExponentialCovarianceFunction() {}

      //! Evaluate covariance
      value_type evaluateCovariance(const value_type& x,
				    const value_type& xp) const;

      //! Get eigenpairs
      const Teuchos::Array< OneDEigenPair<value_type> >& getEigenPairs() const;

    private:

      //! Prohibit copying
      OneDExponentialCovarianceFunction(const OneDExponentialCovarianceFunction&);
      
      //! Prohibit copying
      OneDExponentialCovarianceFunction& operator=(const OneDExponentialCovarianceFunction&);

    protected:

      typedef typename Teuchos::ScalarTraits<value_type>::magnitudeType magnitude_type;

      //! Correlation length
      value_type L;

      //! Eigenpairs
      Teuchos::Array< OneDEigenPair<value_type> > eig_pair;

      //! A basic root finder based on Newton's method
      template <class Func>
      value_type newton(const Func& func, const value_type& a, 
			const value_type& b, magnitude_type tol, 
			int max_num_its);
      
      //! A basic root finder based on bisection
      template <class Func>
      value_type bisection(const Func& func, const value_type& a, 
			   const value_type& b, magnitude_type tol, 
			   int max_num_its);

      /*! 
       * \brief Nonlinear function whose roots define eigenvalues for sin() 
       * eigenfunction
       */
      struct EigFuncSin {
	const value_type& alpha;
	EigFuncSin(const value_type& alpha_) : alpha(alpha_) {}
	value_type eval(const value_type& u) const { 
	  return alpha*std::tan(u) + u; }
	value_type deriv(const value_type& u) const { 
	  return alpha/(std::cos(u)*std::cos(u)) + 1.0; }
      };

      /*! 
       * \brief Nonlinear function whose roots define eigenvalues for cos() 
       * eigenfunction
       */
      struct EigFuncCos {
	const value_type& alpha;
	EigFuncCos(const value_type& alpha_) : alpha(alpha_) {}
	value_type eval(const value_type& u) const { 
	  return alpha - u*std::tan(u); }
	value_type deriv(const value_type& u) const { 
	  return -std::tan(u) - u/(std::cos(u)*std::cos(u)); }
      };

    }; // class OneDExponentialCovarianceFunction

  } // namespace KL

} // namespace Stokhos

// Include template definitions
#include "Stokhos_KL_OneDExponentialCovarianceFunctionImp.hpp"

#endif // STOKHOS_KL_ONE_D_EXPONENTIAL_COVARIANCE_FUNCTION_HPP
