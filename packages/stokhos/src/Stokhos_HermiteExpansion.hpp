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

#ifndef STOKHOS_HERMITEEXPANSION_HPP
#define STOKHOS_HERMITEEXPANSION_HPP

#include <vector>
#include <cmath>
#include <algorithm>	// for std::min and std::max

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Stokhos_HermitePoly.hpp"
#include "Stokhos_TripleProduct.hpp"

namespace Stokhos {

  //! Functions for Hermite polynomial chaos expansions
  template <typename T, typename BasisT> 
  class HermiteExpansion {
  public:
    
    //! Typename of values
    typedef T value_type;

    //! Typename of basis
    typedef BasisT basis_type;

    //! Constructor
    HermiteExpansion(unsigned int d);

    //! Resize workspace
    void resize(unsigned int d);

    //! Get expansion degree
    unsigned int degree() const { return sz-1; }

    //! Get basis
    const basis_type& getBasis() const {return Cijk.getBasis(); }
 
    // Operations
    void unaryMinus(HermitePoly<T>& c, const HermitePoly<T>& a);

    void plusEqual(HermitePoly<T>& c, const T& x);
    void minusEqual(HermitePoly<T>& c, const T& x);
    void timesEqual(HermitePoly<T>& c, const T& x);
    void divideEqual(HermitePoly<T>& c, const T& x);

    void plusEqual(HermitePoly<T>& c, const HermitePoly<T>& x);
    void minusEqual(HermitePoly<T>& c, const HermitePoly<T>& x);
    void timesEqual(HermitePoly<T>& c, const HermitePoly<T>& x);
    void divideEqual(HermitePoly<T>& c, const HermitePoly<T>& x);

    void plus(HermitePoly<T>& c, const HermitePoly<T>& a, 
	      const HermitePoly<T>& b);
    void plus(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void plus(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);
    void minus(HermitePoly<T>& c, const HermitePoly<T>& a,
	       const HermitePoly<T>& b);
    void minus(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void minus(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);
    void times(HermitePoly<T>& c, const HermitePoly<T>& a, 
	       const HermitePoly<T>& b);
    void times(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void times(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);
    void divide(HermitePoly<T>& c, const HermitePoly<T>& a, 
		const HermitePoly<T>& b);
    void divide(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void divide(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);

    void exp(HermitePoly<T>& c, const HermitePoly<T>& a);
    void log(HermitePoly<T>& c, const HermitePoly<T>& a);
    void log10(HermitePoly<T>& c, const HermitePoly<T>& a);
    void sqrt(HermitePoly<T>& c, const HermitePoly<T>& a);
    void pow(HermitePoly<T>& c, const HermitePoly<T>& a, 
	     const HermitePoly<T>& b);
    void pow(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void pow(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);
    void sincos(HermitePoly<T>& s, HermitePoly<T>& c, const HermitePoly<T>& a);
    void cos(HermitePoly<T>& c, const HermitePoly<T>& a);
    void sin(HermitePoly<T>& c, const HermitePoly<T>& a);
    void tan(HermitePoly<T>& c, const HermitePoly<T>& a);
    void sinhcosh(HermitePoly<T>& s, HermitePoly<T>& c, 
		  const HermitePoly<T>& a);
    void cosh(HermitePoly<T>& c, const HermitePoly<T>& a);
    void sinh(HermitePoly<T>& c, const HermitePoly<T>& a);
    void tanh(HermitePoly<T>& c, const HermitePoly<T>& a);
    template <typename OpT> 
    void quad(const OpT& quad_func, HermitePoly<T>& c, const HermitePoly<T>& a,
	      const HermitePoly<T>& b);
    void acos(HermitePoly<T>& c, const HermitePoly<T>& a);
    void asin(HermitePoly<T>& c, const HermitePoly<T>& a);
    void atan(HermitePoly<T>& c, const HermitePoly<T>& a);
    void atan2(HermitePoly<T>& c, const HermitePoly<T>& a,
	       const HermitePoly<T>& b);
    void atan2(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void atan2(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);
    void acosh(HermitePoly<T>& c, const HermitePoly<T>& a);
    void asinh(HermitePoly<T>& c, const HermitePoly<T>& a);
    void atanh(HermitePoly<T>& c, const HermitePoly<T>& a);
    void abs(HermitePoly<T>& c, const HermitePoly<T>& a);
    void fabs(HermitePoly<T>& c, const HermitePoly<T>& a);
    void max(HermitePoly<T>& c, const HermitePoly<T>& a,
	     const HermitePoly<T>& b);
    void max(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void max(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);
    void min(HermitePoly<T>& c, const HermitePoly<T>& a,
	     const HermitePoly<T>& b);
    void min(HermitePoly<T>& c, const T& a, const HermitePoly<T>& b);
    void min(HermitePoly<T>& c, const HermitePoly<T>& a, const T& b);

  protected:

    //! Ordinal type
    typedef int ordinal_type;

    //! Typename of matrix
    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> matrix_type;

    //! Typename of TripleProduct tensor
    typedef TripleProduct<BasisT> tp_type;

    //! Workspace size
    unsigned int sz;
    
    //! Matrix
    matrix_type A;

    //! RHS
    matrix_type B;

    //! Pivot array
    std::vector<ordinal_type> piv;

    //! Triple-product tensor
    tp_type Cijk;
    
    //! LAPACK wrappers
    Teuchos::LAPACK<ordinal_type,value_type> lapack;

  protected:

    //! Solve linear system
    ordinal_type solve(ordinal_type s, ordinal_type nrhs);
    
  }; // class Hermite

} // namespace Stokhos

#include "Stokhos_HermiteExpansionImp.hpp"

#endif // STOKHOS_HERMITEEXPANSION_HPP
