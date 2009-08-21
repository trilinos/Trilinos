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

#ifndef STOKHOS_STIELTJESPCEBASIS_HPP
#define STOKHOS_STIELTJESPCEBASIS_HPP

#include <string>
#include <vector>
#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class StieltjesPCEBasis : 
    public OneDOrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    StieltjesPCEBasis(
       ordinal_type p,
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
       const Stokhos::OrthogPolyBasis<ordinal_type, value_type>& pce_basis,
       const Stokhos::Quadrature<ordinal_type, value_type>& quad,
       bool use_pce_quad_points);

    //! Destructor
    ~StieltjesPCEBasis();

    //! Return order of basis
    virtual ordinal_type order() const;

    //! Return dimension of basis
    virtual ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Compute norm squared of each basis element
    virtual const std::vector<value_type>& norm_squared() const;

    //! Compute norm squared of ith element
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    virtual Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> > getTripleProductTensor() const;

    //! Compute derivative double product tensor
    virtual Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> > getDerivDoubleProductTensor() const;

    //! Project a polynomial into this basis
    virtual void projectPoly(const Polynomial<value_type>& poly, 
                             std::vector<value_type>& coeffs) const;

    //! Project product of two basis polynomials into this basis
    virtual void projectProduct(ordinal_type i, ordinal_type j,
                                std::vector<value_type>& coeffs) const;

    //! Project derivative of basis polynomial into this basis
    virtual void projectDerivative(ordinal_type i, 
                                   std::vector<value_type>& coeffs) const;

    //! Write polynomial in standard basis
    virtual Polynomial<value_type> toStandardBasis(const value_type coeffs[], 
						   ordinal_type n) const;

    //! Evaluate basis polynomial at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const value_type& point,
                               std::vector<value_type>& basis_pts) const;

    //! Print basis
    virtual void print(std::ostream& os) const;

    //! Return name of basis
    virtual const std::string& getName() const;

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  std::vector<value_type>& points,
		  std::vector<value_type>& weights,
		  std::vector< std::vector<value_type> >& values) const;

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const;

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const;

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const;

    void transformCoeffsFromStieltjes(const value_type *in, value_type *out) const;

  protected:

    //! Compute 3-term recurrence using Stieljtes procedure
    void stieltjes(ordinal_type nstart,
		   ordinal_type nfinish,
		   const std::vector<value_type>& weights,
		   const std::vector<value_type>& points,
		   std::vector<value_type>& a,
		   std::vector<value_type>& b,
		   std::vector<value_type>& nrm,
		   std::vector< std::vector<value_type> >& phi_vals) const;

    //! Compute \int \phi^2_k(t) d\lambda(t) and \int t\phi^2_k(t) d\lambda(t)
    void integrateBasisSquared(ordinal_type k, 
			       const std::vector<value_type>& a,
			       const std::vector<value_type>& b,
			       const std::vector<value_type>& weights,
			       const std::vector<value_type>& points,
			       std::vector< std::vector<value_type> >& phi_vals,
			       value_type& val1, value_type& val2) const;

    //! Evaluate polynomials via 3-term recurrence
    void evaluateRecurrence(ordinal_type k,
			    const std::vector<value_type>& a,
			    const std::vector<value_type>& b,
			    const std::vector<value_type>& points,
			    std::vector< std::vector<value_type> >& values) const;

  private:

    // Prohibit copying
    StieltjesPCEBasis(const StieltjesPCEBasis&);

    // Prohibit Assignment
    StieltjesPCEBasis& operator=(const StieltjesPCEBasis& b);

  protected:

    //! Name
    std::string name;

    //! Order
    ordinal_type p;

    //! Polynomial norms
    mutable std::vector<value_type> norms;

    //! 3-term recurrence alpha coefficients
    mutable std::vector<value_type> alpha;

    //! 3-term recurrence beta coefficients
    mutable std::vector<value_type> beta;

    //! Values of PCE at quadrature points
    mutable std::vector<value_type> pce_vals;

    //! PCE quadrature weights
    std::vector<value_type> pce_weights;

    //! Values of generated polynomials at PCE quadrature points
    mutable std::vector< std::vector<value_type> > phi_vals;

    //! Triple product tensor
    mutable Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Cijk;

    //! Use underlying pce's quadrature data
    bool use_pce_quad_points;

    //! Matrix mapping coefficients in Stieltjes basis back to original basis
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> fromStieltjesMat;

  }; // class StieltjesPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_StieltjesPCEBasisImp.hpp"

#endif
