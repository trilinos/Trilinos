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

#ifndef STOKHOS_COMPLETEPOLYNOMIALBASIS_HPP
#define STOKHOS_COMPLETEPOLYNOMIALBASIS_HPP

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class CompletePolynomialBasis : 
    public OrthogPolyBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    CompletePolynomialBasis(const std::vector< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > >& bases,
			    const value_type& sparse_tol = 1.0e-12,
			    const Teuchos::RCP< std::vector<value_type> >& deriv_coeffs = Teuchos::null);

    //! Destructor
    virtual ~CompletePolynomialBasis();

    //! Return order of basis
    ordinal_type order() const;

    //! Return dimension of basis
    ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Compute norm squared of each basis element
    virtual const std::vector<value_type>& norm_squared() const;

    //! Compute norm squared of ith element
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    virtual Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> > getTripleProductTensor() const;

    //! Compute derivative triple product tensor
    virtual Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> > getDerivTripleProductTensor() const;

    //! Compute derivative double product tensor
    virtual Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> > getDerivDoubleProductTensor() const;

    //! Project product of basis polynomials i and j onto this basis
    virtual void projectProduct(ordinal_type i, ordinal_type j, std::vector<value_type>& coeffs) const;

    //! Project derivative of basis polynomial into this basis
    virtual void projectDerivative(ordinal_type i, 
                                   std::vector<value_type>& coeffs) const;

    //! Evaluate basis polynomial at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point
    virtual const std::vector<value_type>& 
    evaluateBases(const std::vector<value_type>& point) const;

    //! Print basis
    virtual void print(std::ostream& os) const;

    //! Get term
    virtual std::vector<ordinal_type> getTerm(ordinal_type i) const;

    //! Get index
    virtual ordinal_type 
    getIndex(const std::vector<ordinal_type>& term) const;

    //! Return name of basis
    virtual const std::string& getName() const;

    //! Return coordinate bases
    const std::vector< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& getCoordinateBases() const;

  protected:

    /*! 
     * \brief Computes the number of terms in an expansion of dimension \c dim
     * and order \c order.
     */
    /*!
     * Returns (order+dim)!/(order!*dim!)
     */
    ordinal_type compute_num_terms(ordinal_type dim, ordinal_type order) const;

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension
     */
    void compute_terms();

    /*!
     * \brief Compute basis index given the orders for each basis
     * dimension.
     */
    ordinal_type compute_index(const std::vector<ordinal_type>& terms) const;

  private:

    // Prohibit copying
    CompletePolynomialBasis(const CompletePolynomialBasis&);

    // Prohibit Assignment
    CompletePolynomialBasis& operator=(const CompletePolynomialBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Array of bases
    std::vector< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > > bases;

    //! Tolerance for computing sparse Cijk
    value_type sparse_tol;

    //! Coefficients for derivative
    Teuchos::RCP< std::vector<value_type> > deriv_coeffs;

    //! Norms
    std::vector<value_type> norms;

    //! 2-D array of basis terms
    std::vector< std::vector<ordinal_type> > terms;

    //! Array of Triple products for computing product projections
    std::vector< Teuchos::RCP<const Dense3Tensor<ordinal_type,value_type> > > Cijk_1d;

    //! Array of double products for computing derivative projections
    std::vector< Teuchos::RCP<const Teuchos::SerialDenseMatrix<ordinal_type,value_type> > > Bij_1d;

    //! Triple product 3 tensor
    mutable Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk;

    //! Derivative triple product 3 tensor
    mutable Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Dijk;

    //! Derivative double product 2 tensor
    mutable Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij;

    //! Array to hold basis evaluations
    mutable std::vector<value_type> basis_pts;

    //! Temporary array used in basis evaluation
    mutable std::vector< std::vector<value_type> > basis_eval_tmp;

  }; // class CompletePolynomialBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_CompletePolynomialBasisImp.hpp"

#endif
