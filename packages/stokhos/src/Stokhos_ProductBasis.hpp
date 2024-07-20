// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_PRODUCTBASIS_HPP
#define STOKHOS_PRODUCTBASIS_HPP

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  /*! 
   * \brief Abstract base class for multivariate orthogonal polynomials
   * generated from tensor products of univariate polynomials.
   */
  /*!
   * * The multivariate polynomials are given by 
   * \f[
   *     \Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)
   * \f]
   * where \f$d\f$ is the dimension of the basis.  This class adds methods 
   * for indexing the multivariate polynomial and getting the coordinate bases.
   */
  template <typename ordinal_type, typename value_type>
  class ProductBasis : 
    public virtual OrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    ProductBasis() {};

    //! Destructor
    virtual ~ProductBasis() {};

    //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual const MultiIndex<ordinal_type>& term(ordinal_type i) const = 0;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type index(const MultiIndex<ordinal_type>& term) const = 0;

    //! Return array of coordinate bases
    /*!
     * Array is of size dimension().
     */
    virtual 
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
							   value_type> > > 
    getCoordinateBases() const = 0;

    //! Return maximum order allowable for each coordinate basis
    virtual MultiIndex<ordinal_type> getMaxOrders() const = 0;

  private:

    // Prohibit copying
    ProductBasis(const ProductBasis&);

    // Prohibit Assignment
    ProductBasis& operator=(const ProductBasis& b);

  }; // class ProductBasis

} // Namespace Stokhos

#endif // STOKHOS_PRODUCTBASIS
