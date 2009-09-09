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

#ifndef STOKHOS_HERMITEBASIS_HPP
#define STOKHOS_HERMITEBASIS_HPP

#include "Stokhos_RecurrenceBasis.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class HermiteBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:
    
    //! Constructor
    HermiteBasis(ordinal_type p, bool normalize = false);
    
    //! Destructor
    ~HermiteBasis();

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const { return 5; }

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const { 
      return 0.5/std::sqrt(std::atan(1.0)); // 1/sqrt(pi)
    }

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const { return std::sqrt(2.0); }

  protected:

    //! Compute recurrence coefficients
    virtual void 
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta) const;

  private:

    // Prohibit copying
    HermiteBasis(const HermiteBasis&);

    // Prohibit Assignment
    HermiteBasis& operator=(const HermiteBasis& b);
    
  }; // class HermiteBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_HermiteBasisImp.hpp"

#endif
