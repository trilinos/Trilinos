// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_USERDEFINEDQUADRATURE
#define STOKHOS_USERDEFINEDQUADRATURE

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Teuchos_RCP.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class UserDefinedQuadrature : public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    UserDefinedQuadrature(
     const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
     const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
     const Teuchos::RCP<const Teuchos::Array<value_type> >& weights);

    //! Constructor
    UserDefinedQuadrature(
     const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
     const Teuchos::RCP<const Teuchos::Array<value_type> >& weights,
     const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& values);

    //! Destructor
    virtual ~UserDefinedQuadrature() {}

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
    UserDefinedQuadrature(const UserDefinedQuadrature&);

    // Prohibit Assignment
    UserDefinedQuadrature& operator=(const UserDefinedQuadrature& b);

  protected:

    //! Quadrature points
    Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > > quad_points;

    //! Quadrature weights
    Teuchos::RCP<const Teuchos::Array<value_type> > quad_weights;

    //! Quadrature values
    Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > > quad_values;

  }; // class UserDefinedQuadrature

} // namespace Stokhos

// Include template definitions
#include "Stokhos_UserDefinedQuadratureImp.hpp"

#endif // STOKHOS_USERDEFINEDQUADRATURE
