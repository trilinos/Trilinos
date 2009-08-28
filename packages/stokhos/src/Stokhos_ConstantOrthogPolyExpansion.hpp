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

#ifndef STOKHOS_CONSTANTORTHOGPOLYEXPANSION_HPP
#define STOKHOS_CONSTANTORTHOGPOLYEXPANSION_HPP

#include "Stokhos_OrthogPolyExpansion.hpp"

namespace Stokhos {

  //! Orthogonal polynomial expansion class for constant (size 1) expansions
  /*!
   * This is used primarily by the Sacado overloaded operators to provide
   * an expansion for constant expressions, which simplifies the logic of
   * the overloaded operators signficantly.
   */
  template <typename ordinal_type, typename value_type> 
  class ConstantOrthogPolyExpansion : 
    public OrthogPolyExpansion<ordinal_type, value_type> {
  public:

    //! Constructor
    ConstantOrthogPolyExpansion();

    //! Destructor
    virtual ~ConstantOrthogPolyExpansion() {}

    //! Get expansion size
    ordinal_type size() const { return 1; }

    //! Get basis
    Teuchos::RCP< const OrthogPolyBasis<ordinal_type, value_type> > 
    getBasis() const {return Teuchos::null; }
 
    // Operations
    void unaryMinus(OrthogPolyApprox<ordinal_type, value_type>& c, 
                    const OrthogPolyApprox<ordinal_type, value_type>& a);

    void plusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		   const value_type& x);
    void minusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		    const value_type& x);
    void timesEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		    const value_type& x);
    void divideEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const value_type& x);

    void plusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                   const OrthogPolyApprox<ordinal_type, value_type>& x);
    void minusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                    const OrthogPolyApprox<ordinal_type, value_type>& x);
    void timesEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                    const OrthogPolyApprox<ordinal_type, value_type>& x);
    void divideEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                     const OrthogPolyApprox<ordinal_type, value_type>& x);

    void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a, 
              const OrthogPolyApprox<ordinal_type, value_type>& b);
    void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const value_type& a, 
              const OrthogPolyApprox<ordinal_type, value_type>& b);
    void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a, 
              const value_type& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a,
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const value_type& b);
    void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const value_type& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type>& a, 
                const OrthogPolyApprox<ordinal_type, value_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
                const value_type& a, 
                const OrthogPolyApprox<ordinal_type, value_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type>& a, 
                const value_type& b);

    void exp(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void log(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void log10(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sqrt(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const value_type& b);
    void cos(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sin(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void tan(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void cosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void tanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void acos(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void asin(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atan(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a,
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const value_type& b);
    void acosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void asinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void abs(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void fabs(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a,
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const value_type& b);
    void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a,
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const value_type& b);

  private:

    // Prohibit copying
    ConstantOrthogPolyExpansion(const ConstantOrthogPolyExpansion&);

    // Prohibit Assignment
    ConstantOrthogPolyExpansion& operator=(const ConstantOrthogPolyExpansion& b);
    
  }; // class ConstantOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_ConstantOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_CONSTANTORTHOGPOLYEXPANSION_HPP
