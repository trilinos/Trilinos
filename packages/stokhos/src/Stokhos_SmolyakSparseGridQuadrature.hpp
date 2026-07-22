// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SMOLYAK_SPARSE_GRID_QUADRATURE_HPP
#define STOKHOS_SMOLYAK_SPARSE_GRID_QUADRATURE_HPP

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_SmolyakPseudoSpectralOperator.hpp"
#include "Teuchos_RCP.hpp"

namespace Stokhos {

  /*!
   * \brief Defines quadrature for a tensor product basis by Smolyak sparse
   * grids.
   */
  /*!
   * This class generates the sparse grids using the
   * SmolyakPseudoSpectralOperator and doesn't rely on Dakota.
   */
  template < typename ordinal_type, typename value_type,
             typename point_compare_type =
             typename DefaultPointCompare<ordinal_type,value_type>::type >
  class SmolyakSparseGridQuadrature :
    public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param product_basis product basis
     * \param index_set     index set defining growth levels
     * \param point_compare comparison functor used in ordering points
     */
    template <typename index_set_type>
    SmolyakSparseGridQuadrature(
      const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,
      const index_set_type& index_set,
      const value_type duplicate_tol = 1.0e-12,
      const point_compare_type& point_compare = point_compare_type());

    //! Destructor
    virtual ~SmolyakSparseGridQuadrature() {}

    //! Get number of quadrature points
    virtual ordinal_type size() const { return quad_weights.size(); }

    //! Get quadrature points
    /*!
     * Array is dimensioned Q-by-d where Q is the number of quadrature
     * points and d is the dimension of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> >&
    getQuadPoints() const;

    //! Get quadrature weights
    /*!
     * Array is of size Q where Q is the number of quadrature points.
     */
    virtual const Teuchos::Array<value_type>&
    getQuadWeights() const;

    //! Get values of basis at quadrature points
    /*!
     * Array is dimensioned Q-by-P where Q is the number of quadrature
     * points and P is the size of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> > &
    getBasisAtQuadPoints() const;

    //! Print quadrature data
    virtual std::ostream& print(std::ostream& os) const;

  private:

    // Prohibit copying
    SmolyakSparseGridQuadrature(const SmolyakSparseGridQuadrature&);

    // Prohibit Assignment
    SmolyakSparseGridQuadrature& operator=(const SmolyakSparseGridQuadrature& b);

  protected:

    //! Quadrature points
    Teuchos::Array< Teuchos::Array<value_type> > quad_points;

    //! Quadrature weights
    Teuchos::Array<value_type> quad_weights;

    //! Quadrature values
    Teuchos::Array< Teuchos::Array<value_type> >  quad_values;

  }; // class SmolyakSparseGridQuadrature

} // namespace Stokhos

// Include template definitions
#include "Stokhos_SmolyakSparseGridQuadratureImp.hpp"

#endif //STOKHOS_SMOLYAK_SPARSE_GRID_QUADRATURE_HPP
