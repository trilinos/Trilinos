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

#ifndef STOKHOS_ORTHOGPOLYEXPANSION_HPP
#define STOKHOS_ORTHOGPOLYEXPANSION_HPP

#include <cmath>
#include <algorithm>	// for std::min and std::max

#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  //! Abstract base class for orthogonal polynomial-based expansions
  template <typename ordinal_type, typename value_type> 
  class OrthogPolyExpansion {
  public:

    //! Typename of TripleProduct tensor
    typedef Sparse3Tensor<ordinal_type, value_type> tp_type;

    //! Constructor
    OrthogPolyExpansion() {}

    //! Desturctor
    virtual ~OrthogPolyExpansion() {}

    //! Get expansion size
    virtual ordinal_type size() const = 0;

    //! Get basis
    virtual Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >
    getBasis() const = 0;
 
    // Operations
    virtual void unaryMinus(OrthogPolyApprox<ordinal_type, value_type>& c, 
			    const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;

    virtual void plusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			   const value_type& x) = 0;
    virtual void minusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			    const value_type& x) = 0;
    virtual void timesEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			    const value_type& x) = 0;
    virtual void divideEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			     const value_type& x) = 0;

    virtual void plusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			   const OrthogPolyApprox<ordinal_type, value_type>& x) = 0;
    virtual void minusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			    const OrthogPolyApprox<ordinal_type, value_type>& x) = 0;
    virtual void timesEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			    const OrthogPolyApprox<ordinal_type, value_type>& x) = 0;
    virtual void divideEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
			     const OrthogPolyApprox<ordinal_type, value_type>& x) = 0;
    
    virtual void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a, 
		      const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const value_type& a, 
		      const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a, 
		      const value_type& b) = 0;
    virtual void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a,
		       const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const value_type& a, 
		       const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a, 
		       const value_type& b) = 0;
    virtual void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a, 
		       const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const value_type& a, 
		       const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a, 
		       const value_type& b) = 0;
    virtual void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
			const OrthogPolyApprox<ordinal_type, value_type>& a, 
			const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
			const value_type& a, 
			const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
			const OrthogPolyApprox<ordinal_type, value_type>& a, 
			const value_type& b) = 0;

    virtual void exp(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void log(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void log10(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void sqrt(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a, 
		     const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const value_type& a, 
		     const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a, 
		     const value_type& b) = 0;
    virtual void cos(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void sin(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void tan(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void cosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void sinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void tanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void acos(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void asin(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void atan(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
//     virtual void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
// 	       const OrthogPolyApprox<ordinal_type, value_type>& a,
// 	       const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
//     virtual void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
// 	       const T& a, 
// 	       const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
//     virtual void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
// 	       const OrthogPolyApprox<ordinal_type, value_type>& a, 
// 	       const T& b) = 0;
    virtual void acosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void asinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void atanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
		       const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void abs(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void fabs(OrthogPolyApprox<ordinal_type, value_type>& c, 
		      const OrthogPolyApprox<ordinal_type, value_type>& a) = 0;
    virtual void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a,
		     const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const value_type& a, 
		     const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a, 
		     const value_type& b) = 0;
    virtual void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a,
		     const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const value_type& a, 
		     const OrthogPolyApprox<ordinal_type, value_type>& b) = 0;
    virtual void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const OrthogPolyApprox<ordinal_type, value_type>& a, 
		     const value_type& b) = 0;

  private:

    // Prohibit copying
    OrthogPolyExpansion(const OrthogPolyExpansion&);

    // Prohibit Assignment
    OrthogPolyExpansion& operator=(const OrthogPolyExpansion& b);
    
  }; // class OrthogPolyExpansion

} // namespace Stokhos

#endif // STOKHOS_ORTHOGPOLYEXPANSION_HPP
