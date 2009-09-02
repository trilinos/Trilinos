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

#ifndef STOKHOS_ONEDORTHOGPOLYBASISBASE_HPP
#define STOKHOS_ONEDORTHOGPOLYBASISBASE_HPP

#include <string>
#include "Stokhos_OneDOrthogPolyBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Partial concrete implementation of OneDOrthogPolyBasis suitable for
   * many classes of polynomials.
   */
  template <typename ordinal_type, typename value_type>
  class OneDOrthogPolyBasisBase : 
    public OneDOrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Destructor
    ~OneDOrthogPolyBasisBase();

    //! Return order of basis
    virtual ordinal_type order() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Compute norm squared of each basis element
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Compute norm squared of ith element
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    virtual Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> > getTripleProductTensor() const;

    //! Compute derivative double product tensor
    virtual Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> > getDerivDoubleProductTensor() const;

    //! Project product of two basis polynomials into this basis
    virtual void projectProduct(ordinal_type i, ordinal_type j,
                                Teuchos::Array<value_type>& coeffs) const;

    //! Evaluate basis polynomial at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    virtual void print(std::ostream& os) const;

    //! Return name of basis
    virtual const std::string& getName() const;

  protected:

    //! Constructor
    OneDOrthogPolyBasisBase(const std::string& name, ordinal_type p);

  private:

    // Prohibit copying
    OneDOrthogPolyBasisBase(const OneDOrthogPolyBasisBase&);

    // Prohibit Assignment
    OneDOrthogPolyBasisBase& operator=(const OneDOrthogPolyBasisBase& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Order of basis
    ordinal_type p;

    //! Basis polynomials
    Teuchos::Array< Polynomial<value_type> > basis;

    //! double-sized basis polynomials for accurately projecting products
    Teuchos::Array< Polynomial<value_type> > double_basis;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Triple product tensor
    mutable Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Cijk;

    //! Deriv double product tensor
    mutable Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij;

  }; // class OrthogPolyBasisBase

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_OneDOrthogPolyBasisBaseImp.hpp"

#endif
