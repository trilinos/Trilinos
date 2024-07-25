// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DERIVBASIS_HPP
#define STOKHOS_DERIVBASIS_HPP

#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {
  
  /*! 
   * \brief Abstract base class for multivariate orthogonal polynomials
   * that support computing double and triple products involving derivatives
   * of the basis polynomials.
   */
  template <typename ordinal_type, typename value_type>
  class DerivBasis : 
    public virtual OrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    DerivBasis() {};

    //! Destructor
    virtual ~DerivBasis() {};

    /*! 
     * \brief Compute triple product tensor 
     * \f$D_{ijk} = \langle\Psi_i\Psi_j D_v\Psi_k\rangle\f$ where 
     * \f$D_v\Psi_k\f$ represents the derivative of \f$\Psi_k\f$ in the 
     * direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the derived class implementation.
     */
    virtual 
    Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > 
    computeDerivTripleProductTensor(
      const Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij,
      const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const = 0;

    /*! 
     * \brief Compute double product tensor 
     * \f$B_{ij} = \langle \Psi_i D_v\Psi_j\rangle\f$ where \f$D_v\Psi_j\f$
     * represents the derivative of \f$\Psi_j\f$ in the direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the derived class implementation.
     */
    virtual 
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > 
    computeDerivDoubleProductTensor() const = 0;

  private:

    // Prohibit copying
    DerivBasis(const DerivBasis&);

    // Prohibit Assignment
    DerivBasis& operator=(const DerivBasis& b);

  }; // class DerivBasis

} // Namespace Stokhos

#endif // STOKHOS_DERIVBASIS
