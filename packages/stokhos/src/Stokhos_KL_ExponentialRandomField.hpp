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

#ifndef STOKHOS_KL_EXPONENTIAL_RANDOM_FIELD_HPP
#define STOKHOS_KL_EXPONENTIAL_RANDOM_FIELD_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_PromotionTraits.hpp"

#include "Stokhos_KL_ProductEigenPair.hpp"

namespace Stokhos {

  namespace KL {

    //! Class representing a %KL expansion of an exponential random field
    /*!
     * This class provides a means for evaluating a random field 
     * \f$a(x,\theta)\f$, \f$x\in D\f$, \f$\theta\in\Omega\f$ through the 
     * %KL expansion
     * \f[
     *     a(x,\theta) \approx a_0 + 
     *       \sigma\sum_{k=1}^M \sqrt{\lambda_k}b_k(x)\xi_k(\theta)
     * \f]
     * where \f$D\f$ is a \f$d\f$-dimensional hyper-rectangle, for the case 
     * when the covariance function of \f$a\f$ is exponential:
     * \f[
     *     \mbox{cov}(x,x') = \sigma\exp(-|x_1-x_1'|/L_1-...-|x_d-x_d'|/L_d).
     * \f]
     * In this case, the covariance function and domain factor into a product
     * 1-dimensional covariance functions over 1-dimensional domains, and thus
     * the eigenfunctions \f$b_k\f$ and eigenvalues \f$\lambda_k\f$ factor into
     * a product of corresponding 1-dimensional eigenfunctions and values. 
     * These are computed by the OneDExponentialCovarianceFunction class
     * For a given value of \f$M\f$, the code works by computing the \f$M\f$
     * eigenfunctions for each direction using this class.  
     * Then all possible tensor products of these
     * one-dimensional eigenfunctions are formed, with corresponding 
     * eigenvalue given by the product of the one-dimensional eigenvalues.
     * These eigenvalues are then sorted in increasing order, and the first
     * \f$M\f$ are chosen as the \f$M\f$ %KL eigenpairs.  Then given values
     * for the random variables \f$\xi_k\f$, the class provides a routine
     * for evaluating a realization of the random field.
     *
     * The idea for this approach was provided by Chris Miller.
     *
     * All data is passed into this class through a Teuchos::ParameterList,
     * which accepts the following parameters:
     * <ul>
     *   <li> "Number of KL Terms" -- [int] (Required) 
     *        Number \f$M\f$ of %KL terms
     *   <li> "Domain Upper Bounds" -- [Teuchos::Array<value_type>] (Required)
     *        Domain upper bounds \f$b_i\f$ for each dimension \f$i\f$
     *   <li> "Domain Lower Bounds" -- [Teuchos::Array<value_type>] (Required)
     *        Domain lower bounds \f$a_i\f$ for each dimension \f$i\f$
     *   <li> "Correlation Lengths" -- [Teuchos::Array<value_type>[ (Required)
     *        Correlation length \f$L_i\f$ for each dimension \f$i\f$
     *   <li> "Mean" -- [value_type] (Required)
     *        Mean \f$a_0\f$ of the random field
     *   <li> "Standard Deviation" -- [value_type] (Required)
     *        Standard devation \f$\sigma\f$ of the random field
     * </ul>
     * Additionally parameters for the OneDExponentialCovarianceFunction are
     * also accepted.
     */
    template <typename value_type>
    class ExponentialRandomField {
    public:

      //! Constructor
      ExponentialRandomField(Teuchos::ParameterList& solverParams);

      //! Destructor 
      ~ExponentialRandomField() {}

      //! Return spatial dimension of the field
      int spatialDimension() const { return dim; }

      //! Return stochastic dimension of the field
      int stochasticDimension() const { return num_KL; }

      //! Evaluate random field at a point
      template <typename rvar_type>
      typename Teuchos::PromotionTraits<rvar_type, value_type>::promote
      evaluate(const Teuchos::Array<value_type>& point,
	       const Teuchos::Array<rvar_type>& random_variables) const;

      //! Get eigenpairs
      const Teuchos::Array< ProductEigenPair<value_type> >&
      getEigenPairs() const;
      
      //! Print %KL expansion
      void print(std::ostream& os) const;
      
    private:
      
      //! Prohibit copying
      ExponentialRandomField(const ExponentialRandomField&);
      
      //! Prohibit copying
      ExponentialRandomField& operator=(const ExponentialRandomField&);

    protected:
      
      //! Number of %KL terms
      int num_KL;
      
      //! Dimension of expansion
      int dim;
      
      //! Domain upper bounds
      Teuchos::Array<value_type> domain_upper_bound;

      //! Domeain lower bounds
      Teuchos::Array<value_type> domain_lower_bound;
      
      //! Correlation lengths
      Teuchos::Array<value_type> correlation_length;
      
      //! Mean of random field
      value_type mean;
      
      //! Standard deviation of random field
      value_type std_dev;
      
      //! Product eigenfunctions
      Teuchos::Array< ProductEigenPair<value_type> > product_eig_pairs;
      
    }; // class ExponentialRandomField

  } // namespace KL

} // namespace Stokhos

// Include template definitions
#include "Stokhos_KL_ExponentialRandomFieldImp.hpp"

#endif // STOKHOS_KL_EXPONENTIAL_RANDOM_FIELD_HPP
