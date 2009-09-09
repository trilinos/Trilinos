// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_TENSORPRODUCTQUADRATURE
#define STOKHOS_TENSORPRODUCTQUADRATURE

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Teuchos_RCP.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class TensorProductQuadrature : public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    TensorProductQuadrature(const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis);

    //! Variable order Constructor
    TensorProductQuadrature(const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis, const ordinal_type& quad_order);

    //! Destructor
    virtual ~TensorProductQuadrature() {}

    //! Get quadrature points
    virtual const Teuchos::Array< Teuchos::Array<value_type> >& 
    getQuadPoints() const;

    //! Get quadrature weights
    virtual const Teuchos::Array<value_type>& 
    getQuadWeights() const;

    //! Get values of basis at quadrature points
    virtual const Teuchos::Array< Teuchos::Array<value_type> > & 
    getBasisAtQuadPoints() const;

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
