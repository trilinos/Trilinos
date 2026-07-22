// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    void cbrt(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
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
