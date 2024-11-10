// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

    //! Get number of quadrature points
    virtual ordinal_type size() const { return quad_weights->size(); }

    //! Get quadrature points
    virtual const Teuchos::Array< Teuchos::Array<value_type> >& 
    getQuadPoints() const;

    //! Get quadrature weights
    virtual const Teuchos::Array<value_type>& 
    getQuadWeights() const;

    //! Get values of basis at quadrature points
    virtual const Teuchos::Array< Teuchos::Array<value_type> > & 
    getBasisAtQuadPoints() const;

    //! Print quadrature data
    virtual std::ostream& print(std::ostream& os) const;

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
