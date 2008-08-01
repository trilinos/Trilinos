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
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_TripleProduct.hpp"

namespace Stokhos {

  template <typename T>
  class CompletePolynomialBasis : public OrthogPolyBasis<T> {
  public:

    //! Typename of values
    typedef typename OrthogPolyBasis<T>::value_type value_type;

    //! Constructor
    CompletePolynomialBasis(
	   const std::vector< Teuchos::RCP<const OrthogPolyBasis<T> > >& bases,
	   const std::vector<T>& deriv_coeffs);

    //! Destructor
    ~CompletePolynomialBasis();

    //! Return order of basis
    unsigned int order() const;

    //! Return dimension of basis
    unsigned int dimension() const;

    //! Return total size of basis
    virtual unsigned int size() const;

    //! Compute norm squared of each basis element
    virtual const std::vector<T>& norm_squared() const;

    //! Project a polynomial into this basis
    virtual void projectPoly(const Polynomial<T>& poly, 
			     std::vector<T>& coeffs) const;

    //! Project product of two basis polynomials into this basis
    virtual void projectProduct(unsigned int i, unsigned int j,
				std::vector<T>& coeffs) const;

    //! Project derivative of basis polynomial into this basis
    /*!
     * Derivative here is defined to be:
     * D = a_1D_1 + ... + a_d D_d
     * Where a_1,...,a_d are the values in \c deriv_coeffs supplied in the
     * constructor.  All a_1,...,a_d should be non-zero.
     */
    virtual void projectDerivative(unsigned int i, 
				   std::vector<T>& coeffs) const;

    //! Write polynomial in standard basis
    virtual Polynomial<T> toStandardBasis(const T coeffs[], 
					  unsigned int n) const;

    //! Evaluate basis polynomial at zero
    virtual T evaluateZero(unsigned int i) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const std::vector<T>& point,
			       std::vector<T>& basis_pts) const;
    //! Print basis
    virtual void print(std::ostream& os) const;

    //! Get term
    virtual std::vector<unsigned int> getTerm(unsigned int i) const;

    //! Get index
    virtual unsigned int 
    getIndex(const std::vector<unsigned int>& term) const;

    //! Return name of basis
    virtual const std::string& getName() const;

  protected:

    /*! 
     * \brief Computes the number of terms in an expansion of dimension \c dim
     * and order \c order.
     */
    /*!
     * Returns (order+dim)!/(order!*dim!)
     */
    unsigned int compute_num_terms(unsigned int dim, unsigned int order) const;

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension
     */
    void compute_terms();

    /*!
     * \brief Compute basis index given the orders for each basis
     * dimension.
     */
    unsigned int compute_index(const std::vector<unsigned int>& terms) const;

  private:

    // Prohibit copying
    CompletePolynomialBasis(const CompletePolynomialBasis&);

    // Prohibit Assignment
    CompletePolynomialBasis& operator=(const CompletePolynomialBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Total order of basis
    unsigned int p;

    //! Total dimension of basis
    unsigned int d;

    //! Total size of basis
    unsigned int sz;

    //! Array of bases
    std::vector< Teuchos::RCP<const OrthogPolyBasis<T> > > bases;

    //! Coefficients for derivative
    std::vector<T> deriv_coeffs;

    //! Norms
    std::vector<T> norms;

    //! 2-D array of basis terms
    std::vector< std::vector<unsigned int> > terms;

    //! Array of Triple products for computing product projections
    std::vector< Teuchos::RCP<TripleProduct< OrthogPolyBasis<T> > > > Cijk;

  }; // class CompletePolynomialBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_CompletePolynomialBasisImp.hpp"

#endif
