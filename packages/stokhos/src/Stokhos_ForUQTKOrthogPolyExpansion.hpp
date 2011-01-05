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

#ifndef STOKHOS_FORUQTKORTHOGPOLYEXPANSION_HPP
#define STOKHOS_FORUQTKORTHOGPOLYEXPANSION_HPP

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_FORUQTK

#include "Stokhos_OrthogPolyExpansionBase.hpp"

namespace Stokhos {

  /*! 
   * Orthogonal polynomial expansions based on methods provided by UQTK 
   * (fortran version).  These include Taylor series and time integration
   * calculations for transcendental operations.
   */
  template <typename ordinal_type, typename value_type> 
  class ForUQTKOrthogPolyExpansion : 
    public OrthogPolyExpansionBase<ordinal_type, value_type, 
				   Stokhos::StandardStorage<ordinal_type, value_type> > {
  public:

    typedef Stokhos::StandardStorage<ordinal_type, value_type> node_type;

    enum EXPANSION_METHOD {
      TAYLOR,
      INTEGRATION
    };

    //! Constructor
    ForUQTKOrthogPolyExpansion(
      const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
			       EXPANSION_METHOD method = TAYLOR,
			       value_type rtol = 1.0e-12);

    //! Destructor
    virtual ~ForUQTKOrthogPolyExpansion() {}
 
    // Operations
    void timesEqual(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
		    const value_type& x);
    void divideEqual(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
		     const value_type& x);

    void timesEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);
    void divideEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);

   
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const value_type& a, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
	       const value_type& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
		const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
		const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
		const value_type& a, 
		const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
		const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
		const value_type& b);

    void exp(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void log(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void log10(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void sqrt(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void pow(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const value_type& a, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
	     const value_type& b);
    void cos(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void sin(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void tan(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	     const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void cosh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void sinh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void tanh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void acos(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void asin(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void atan(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void atan2(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const value_type& a, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
	       const value_type& b);
    void acosh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void asinh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void atanh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	       const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);

  private:

    // Prohibit copying
    ForUQTKOrthogPolyExpansion(const ForUQTKOrthogPolyExpansion&);

    // Prohibit Assignment
    ForUQTKOrthogPolyExpansion& operator=(const ForUQTKOrthogPolyExpansion& b);

  protected:

    //! Order
    int order;
    
    //! Dimension
    int dim;

    //! Total size
    ordinal_type sz;

    //! Tolerance for Taylor method
    double rtol;

    //! Expansion method
    EXPANSION_METHOD method;
    
  }; // class ForUQTKOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_ForUQTKOrthogPolyExpansionImp.hpp"

#endif // HAVE_STOKHOS_FORUQTK

#endif // STOKHOS_FORUQTKORTHOGPOLYEXPANSION_HPP
