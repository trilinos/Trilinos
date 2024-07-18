// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TENSORPRODUCTQUADRATURE
#define STOKHOS_TENSORPRODUCTQUADRATURE

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Teuchos_RCP.hpp"

namespace Stokhos {

  /*! 
   * \brief Defines quadrature for a tensor product basis by tensor products of
   * 1-D quadrature rules.
   */
  template <typename ordinal_type, typename value_type>
  class TensorProductQuadrature : public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param product_basis product basis
     * The order of the quadrature is \f$2*p\f$, where \f$p\f$ is the order
     * of the basis.
     */
    TensorProductQuadrature(const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis);

    //! Variable order constructor
    /*!
     * \param product_basis product basis
     * \param quad_order order of quadrature to use
     */
    TensorProductQuadrature(const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis, const ordinal_type& quad_order);

    //! Destructor
    virtual ~TensorProductQuadrature() {}

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
    TensorProductQuadrature(const TensorProductQuadrature&);

    // Prohibit Assignment
    TensorProductQuadrature& operator=(const TensorProductQuadrature& b);

  protected:

    //! Quadrature points
    Teuchos::Array< Teuchos::Array<value_type> > quad_points;

    //! Quadrature weights
    Teuchos::Array<value_type> quad_weights;

    //! Quadrature values
    Teuchos::Array< Teuchos::Array<value_type> >  quad_values;

  }; // class TensorProductQuadrature

} // namespace Stokhos

// Include template definitions
#include "Stokhos_TensorProductQuadratureImp.hpp"

#endif // STOKHOS_TENSORPRODUCTQUADRATURE
