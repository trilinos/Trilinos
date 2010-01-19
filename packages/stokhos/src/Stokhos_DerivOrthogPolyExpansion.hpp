// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef STOKHOS_DERIVORTHOGPOLYEXPANSION_HPP
#define STOKHOS_DERIVORTHOGPOLYEXPANSION_HPP

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_DerivBasis.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Stokhos {

  //! Othogonal polynomial expansions based on derivative calculations
  template <typename ordinal_type, typename value_type> 
  class DerivOrthogPolyExpansion : public OrthogPolyExpansion<ordinal_type, value_type> {
  public:

    //! Constructor
    DerivOrthogPolyExpansion(
      const Teuchos::RCP<const DerivBasis<ordinal_type, value_type> >& basis,
      const Teuchos::RCP<const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
      const Teuchos::RCP<const Stokhos::Dense3Tensor<ordinal_type, value_type> >& Dijk);

    //! Destructor
    virtual ~DerivOrthogPolyExpansion() {}

    //! Get expansion size
    ordinal_type size() const { return sz; }

    //! Get basis
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >
    getBasis() const {return basis; }
 
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
    void sincos(OrthogPolyApprox<ordinal_type, value_type>& s, 
		OrthogPolyApprox<ordinal_type, value_type>& c, 
		const OrthogPolyApprox<ordinal_type, value_type>& a);
    void cos(OrthogPolyApprox<ordinal_type, value_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sin(OrthogPolyApprox<ordinal_type, value_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type>& a);
    void tan(OrthogPolyApprox<ordinal_type, value_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sinhcosh(OrthogPolyApprox<ordinal_type, value_type>& s, 
		  OrthogPolyApprox<ordinal_type, value_type>& c, 
		  const OrthogPolyApprox<ordinal_type, value_type>& a);
    void cosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a);
    void tanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a);
    template <typename OpT> 
    void quad(const OpT& quad_func, 
	      OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a,
	      const OrthogPolyApprox<ordinal_type, value_type>& b);
    void acos(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a);
    void asin(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atan(OrthogPolyApprox<ordinal_type, value_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type>& a);
//     void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
// 	       const OrthogPolyApprox<ordinal_type, value_type>& a,
// 	       const OrthogPolyApprox<ordinal_type, value_type>& b);
//     void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
// 	       const T& a, 
// 	       const OrthogPolyApprox<ordinal_type, value_type>& b);
//     void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
// 	       const OrthogPolyApprox<ordinal_type, value_type>& a, 
// 	       const T& b);
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
    void derivative(OrthogPolyApprox<ordinal_type, value_type>& c, 
		    const OrthogPolyApprox<ordinal_type, value_type>& a);

  private:

    // Prohibit copying
    DerivOrthogPolyExpansion(const DerivOrthogPolyExpansion&);

    // Prohibit Assignment
    DerivOrthogPolyExpansion& operator=(const DerivOrthogPolyExpansion& b);

  protected:

    //! Basis
    Teuchos::RCP< const Stokhos::DerivBasis<ordinal_type, value_type> > basis;

    //! Derivative double-product tensor
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij;

    //! Triple-product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk;

    //! Derivative Triple-product tensor
    Teuchos::RCP<const Stokhos::Dense3Tensor<ordinal_type, value_type> > Dijk;

    //! Workspace size
    ordinal_type sz;
    
    //! Matrix
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> A;

    //! RHS
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> B;

    //! Pivot array
    Teuchos::Array<ordinal_type> piv;
    
    //! LAPACK wrappers
    Teuchos::LAPACK<ordinal_type,value_type> lapack;

  protected:

    //! Solve linear system
    ordinal_type solve(ordinal_type s, ordinal_type nrhs);

    struct acos_quad_func { 
      value_type operator() (const value_type& a) const { 
	return std::acos(a); 
      } 
    };

    struct asin_quad_func { 
      value_type operator() (const value_type& a) const { 
	return std::asin(a); 
      } 
    };

    struct atan_quad_func { 
      value_type operator() (const value_type& a) const { 
	return std::atan(a); 
      } 
    };

    struct acosh_quad_func { 
      value_type operator() (const value_type & a) const { 
	return std::log(a+std::sqrt(a*a-value_type(1.0))); 
      }
    };

    struct asinh_quad_func { 
      value_type operator() (const value_type& a) const { 
	return std::log(a+std::sqrt(a*a+value_type(1.0))); 
      }
    };

    struct atanh_quad_func { 
      value_type operator() (const value_type& a) const { 
	return 0.5*std::log((value_type(1.0)+a)/(value_type(1.0)-a)); 
      } 
    };
    
  }; // class DerivOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_DerivOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_DERIVORTHOGPOLYEXPANSION_HPP
