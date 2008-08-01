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

#ifndef STOKHOS_TAYORTHOGPOLYEXPANSION_HPP
#define STOKHOS_TAYORTHOGPOLYEXPANSION_HPP

#include <vector>
#include <cmath>
#include <algorithm>	// for std::min and std::max

#include "Teuchos_RCP.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_TripleProduct.hpp"
#include "Stokhos_TayTripleProduct.hpp"

namespace Stokhos {

  //! Functions for Hermite polynomial chaos expansions
  template <typename T> 
  class TayOrthogPolyExpansion {
  public:
    
    //! Typename of values
    typedef T value_type;

    //! Typename of TripleProduct tensor
    typedef TayTripleProduct< OrthogPolyBasis<T> >  tp_type;
    //typedef TripleProduct< OrthogPolyBasis<T> > tp_type;

    //! Constructor
    TayOrthogPolyExpansion(const Teuchos::RCP<const OrthogPolyBasis<T> >& basis);

    //! Get expansion size
    unsigned int size() const { return sz; }

    //! Get basis
    const OrthogPolyBasis<T>& getBasis() const {return *basis; }

    //! Get triple product
    const tp_type& getTripleProduct() const { return Cijk; }
 
    // Operations
    void unaryMinus(OrthogPolyApprox<value_type>& c, 
		    const OrthogPolyApprox<value_type>& a);

    void plusEqual(OrthogPolyApprox<value_type>& c, const value_type& x);
    void minusEqual(OrthogPolyApprox<value_type>& c, const value_type& x);
    void timesEqual(OrthogPolyApprox<value_type>& c, const value_type& x);
    void divideEqual(OrthogPolyApprox<value_type>& c, const value_type& x);

    void plusEqual(OrthogPolyApprox<value_type>& c, 
		   const OrthogPolyApprox<value_type>& x);
    void minusEqual(OrthogPolyApprox<value_type>& c, 
		    const OrthogPolyApprox<value_type>& x);
    void timesEqual(OrthogPolyApprox<value_type>& c, 
		    const OrthogPolyApprox<value_type>& x);
    void divideEqual(OrthogPolyApprox<value_type>& c, 
		     const OrthogPolyApprox<value_type>& x);

    void plus(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a, 
	      const OrthogPolyApprox<value_type>& b);
    void plus(OrthogPolyApprox<value_type>& c, 
	      const value_type& a, 
	      const OrthogPolyApprox<value_type>& b);
    void plus(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a, 
	      const value_type& b);
    void minus(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a,
	       const OrthogPolyApprox<value_type>& b);
    void minus(OrthogPolyApprox<value_type>& c, 
	       const value_type& a, 
	       const OrthogPolyApprox<value_type>& b);
    void minus(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a, 
	       const value_type& b);
    void times(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a, 
	       const OrthogPolyApprox<value_type>& b);
    void times(OrthogPolyApprox<value_type>& c, 
	       const value_type& a, 
	       const OrthogPolyApprox<value_type>& b);
    void times(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a, 
	       const value_type& b);
    void divide(OrthogPolyApprox<value_type>& c, 
		const OrthogPolyApprox<value_type>& a, 
		const OrthogPolyApprox<value_type>& b);
    void divide(OrthogPolyApprox<value_type>& c, 
		const value_type& a, 
		const OrthogPolyApprox<value_type>& b);
    void divide(OrthogPolyApprox<value_type>& c, 
		const OrthogPolyApprox<value_type>& a, 
		const value_type& b);

    void exp(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a);
    void log(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a);
    void log10(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a);
    void sqrt(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void pow(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a, 
	     const OrthogPolyApprox<value_type>& b);
    void pow(OrthogPolyApprox<value_type>& c, 
	     const value_type& a, 
	     const OrthogPolyApprox<value_type>& b);
    void pow(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a, 
	     const value_type& b);
    void cos(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a);
    void sin(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a);
    void tan(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a);
    void cosh(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void sinh(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void tanh(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void acos(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void asin(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void atan(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
//     void atan2(OrthogPolyApprox<value_type>& c, 
// 	       const OrthogPolyApprox<value_type>& a,
// 	       const OrthogPolyApprox<value_type>& b);
//     void atan2(OrthogPolyApprox<value_type>& c, 
// 	       const T& a, 
// 	       const OrthogPolyApprox<value_type>& b);
//     void atan2(OrthogPolyApprox<value_type>& c, 
// 	       const OrthogPolyApprox<value_type>& a, 
// 	       const T& b);
    void acosh(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a);
    void asinh(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a);
    void atanh(OrthogPolyApprox<value_type>& c, 
	       const OrthogPolyApprox<value_type>& a);
    void abs(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a);
    void fabs(OrthogPolyApprox<value_type>& c, 
	      const OrthogPolyApprox<value_type>& a);
    void max(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a,
	     const OrthogPolyApprox<value_type>& b);
    void max(OrthogPolyApprox<value_type>& c, 
	     const value_type& a, 
	     const OrthogPolyApprox<value_type>& b);
    void max(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a, 
	     const value_type& b);
    void min(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a,
	     const OrthogPolyApprox<value_type>& b);
    void min(OrthogPolyApprox<value_type>& c, 
	     const value_type& a, 
	     const OrthogPolyApprox<value_type>& b);
    void min(OrthogPolyApprox<value_type>& c, 
	     const OrthogPolyApprox<value_type>& a, 
	     const value_type& b);

  private:

    // Prohibit copying
    TayOrthogPolyExpansion(const TayOrthogPolyExpansion&);

    // Prohibit Assignment
    TayOrthogPolyExpansion& operator=(const TayOrthogPolyExpansion& b);

  protected:

    //! Basis
    Teuchos::RCP<const OrthogPolyBasis<T> > basis;

    //! Triple-product tensor
    tp_type Cijk;

    //! Order
    int order;
    
    //! Dimension
    int dim;

    //! Total size
    unsigned int sz;

    //! Tolerance for Taylor method
    double rtol;

  protected:
    
  }; // class TayOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_TayOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_TAYORTHOGPOLYEXPANSION_HPP
