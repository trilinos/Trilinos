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

#ifndef STOKHOS_RECURRENCEBASIS_HPP
#define STOKHOS_RECURRENCEBASIS_HPP

#include "Stokhos_OneDOrthogPolyBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Partial concrete implementation of OneDOrthogPolyBasis suitable for
   * many classes of polynomials.
   */
  template <typename ordinal_type, typename value_type>
  class RecurrenceBasis : 
    public OneDOrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Destructor
    virtual ~RecurrenceBasis();

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

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const value_type& point,
                               Teuchos::Array<value_type>& basis_pts) const;

    //! Evaluate basis polynomial given by order at given point 
    virtual value_type evaluate(const value_type& point, ordinal_type order) const;

    //! Print basis
    virtual void print(std::ostream& os) const;

    //! Return name of basis
    virtual const std::string& getName() const;

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

    //! Get recurrence coefficients
    virtual void getRecurrenceCoefficients(Teuchos::Array<value_type>& alpha,
					   Teuchos::Array<value_type>& beta,
					   Teuchos::Array<value_type>& delta,
					   Teuchos::Array<value_type>& gamma) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBasesAndDerivatives(const value_type& point,
					     Teuchos::Array<value_type>& vals,
					     Teuchos::Array<value_type>& derivs) const;

  protected:

    //! Constructor
    RecurrenceBasis(const std::string& name, ordinal_type p, bool normalize);

    //! Setup after computing recurrence coefficients
    void setup();

    //! Compute recurrence coefficients
    virtual void 
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta) const = 0;

  private:

    // Prohibit copying
    RecurrenceBasis(const RecurrenceBasis&);

    // Prohibit Assignment
    RecurrenceBasis& operator=(const RecurrenceBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Order of basis
    ordinal_type p;

    //! Normalize basis
    bool normalize;

    //! Recurrence alpha coefficients
    Teuchos::Array<value_type> alpha;

    //! Recurrence alpha coefficients
    Teuchos::Array<value_type> beta;

    //! Recurrence alpha coefficients
    Teuchos::Array<value_type> delta;

    //! Recurrence gamma coefficients
    Teuchos::Array<value_type> gamma;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Triple product tensor
    mutable Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Cijk;

    //! Deriv double product tensor
    mutable Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij;

  }; // class RecurrenceBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_RecurrenceBasisImp.hpp"

#endif
