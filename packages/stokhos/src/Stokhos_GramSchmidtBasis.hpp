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

#ifndef STOKHOS_GRAMSCHMIDTBASIS_HPP
#define STOKHOS_GRAMSCHMIDTBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class GramSchmidtBasis : 
    public OrthogPolyBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    GramSchmidtBasis(
     const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
     const Teuchos::Array< Teuchos::Array<value_type> >& points,
     const Teuchos::Array<value_type>& weights,
     const value_type& sparse_tol = 1.0e-5);

    //! Destructor
    virtual ~GramSchmidtBasis();

    //! Return order of basis
    ordinal_type order() const;

    //! Return dimension of basis
    ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Compute norm squared of each basis element
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Compute norm squared of ith element
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    virtual Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> > getTripleProductTensor() const;

    virtual Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    getLowOrderTripleProductTensor(ordinal_type order) const;

    //! Evaluate basis polynomial at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const Teuchos::Array<value_type>& point,
			       Teuchos::Array<value_type>& basis_vals) const;

    //! Print basis
    virtual void print(std::ostream& os) const;

    //! Return name of basis
    virtual const std::string& getName() const;

    //! Transform coefficients from original basis to this basis
    void transformCoeffs(const value_type *in, value_type *out) const;

  private:

    // Prohibit copying
    GramSchmidtBasis(const GramSchmidtBasis&);

    // Prohibit Assignment
    GramSchmidtBasis& operator=(const GramSchmidtBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Original basis (not orthogonal w.r.t. inner product)
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> > basis;

    Teuchos::Array<value_type> weights;
    Teuchos::Array< Teuchos::Array<value_type> > basis_values;

    //! Tolerance for computing sparse Cijk
    value_type sparse_tol;

    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Matrix storing gram-schmidt coefficients
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> gs_mat;

    //! Triple product 3 tensor
    mutable Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk;

    //! Temporary array for basis evaluation
    mutable Teuchos::Array<value_type> basis_vals_tmp;

  }; // class GramSchmidtBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_GramSchmidtBasisImp.hpp"

#endif
