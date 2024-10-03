// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_KL_EXPONENTIAL_RANDOM_FIELD_HPP
#define STOKHOS_KL_EXPONENTIAL_RANDOM_FIELD_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_PromotionTraits.hpp"
#include "Kokkos_Core.hpp"

#include "Stokhos_KL_OneDExponentialEigenPair.hpp"
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
    template <typename value_type,
              typename execution_space = Kokkos::DefaultExecutionSpace>
    class ExponentialRandomField {
    public:

      typedef ExponentialOneDEigenFunction<value_type> one_d_eigen_func_type;
      typedef OneDEigenPair<one_d_eigen_func_type> one_d_eigen_pair_type;
      typedef ProductEigenPair<one_d_eigen_func_type,execution_space> product_eigen_pair_type;
      typedef Kokkos::View<one_d_eigen_func_type**,execution_space> eigen_func_array_type;
      typedef Kokkos::View<value_type*,execution_space> eigen_value_array_type;

      //! Default constructor
      ExponentialRandomField() : num_KL(0), mean(0), std_dev(0) {}

      //! Constructor
      ExponentialRandomField(Teuchos::ParameterList& solverParams);

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~ExponentialRandomField() {}

      //! Return spatial dimension of the field
      KOKKOS_INLINE_FUNCTION
      int spatialDimension() const { return dim; }

      //! Return stochastic dimension of the field
      KOKKOS_INLINE_FUNCTION
      int stochasticDimension() const { return num_KL; }

      //! Evaluate random field at a point
      template <typename point_type, typename rv_type>
      KOKKOS_INLINE_FUNCTION
      typename Teuchos::PromotionTraits<typename rv_type::value_type,
                                        value_type>::promote
      evaluate(const point_type& point,
               const rv_type& random_variables) const;

      //! Evaluate mean of random field at a point
      template <typename point_type>
      KOKKOS_INLINE_FUNCTION
      value_type evaluate_mean(const point_type& point) const { return mean; }

      //! Evaluate standard deviation of random field at a point
      template <typename point_type>
      KOKKOS_INLINE_FUNCTION
      typename Teuchos::PromotionTraits<typename point_type::value_type,
                                        value_type>::promote
      evaluate_standard_deviation(const point_type& point) const;

      //! Evaluate given eigenfunction at a point
      template <typename point_type>
      KOKKOS_INLINE_FUNCTION
      typename Teuchos::PromotionTraits<typename point_type::value_type, value_type>::promote
      evaluate_eigenfunction(const point_type& point, int i) const;

      //! Return eigenvalue
      value_type
      KOKKOS_INLINE_FUNCTION
      eigenvalue(int i) const { return product_eigen_values(i); }

      //! Print %KL expansion
      void print(std::ostream& os) const;

    protected:

      //! Number of %KL terms
      int num_KL;

      //! Dimension of expansion
      int dim;

      //! Mean of random field
      value_type mean;

      //! Standard deviation of random field
      value_type std_dev;

      //! Product eigenfunctions
      eigen_func_array_type product_eigen_funcs;

      //! Product eigenvalues
      eigen_value_array_type product_eigen_values;

    }; // class ExponentialRandomField

  } // namespace KL

} // namespace Stokhos

// Include template definitions
#include "Stokhos_KL_ExponentialRandomFieldImp.hpp"

#endif // STOKHOS_KL_EXPONENTIAL_RANDOM_FIELD_HPP
