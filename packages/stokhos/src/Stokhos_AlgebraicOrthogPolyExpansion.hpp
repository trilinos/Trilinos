// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_ALGEBRAICORTHOGPOLYEXPANSION_HPP
#define STOKHOS_ALGEBRAICORTHOGPOLYEXPANSION_HPP

#include "Stokhos_OrthogPolyExpansionBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Stokhos {

  //! Orthogonal polynomial expansions limited to algebraic operations
  /*!
   * Implements +, -, *, / (with constant denominator), and other operations
   * with constant arguments.  Useful for problems that only involve
   * algebraic operations, or constant non-algebraic operations in that it
   * doesn't require creation of additional data structures (like sparse
   * quadrature grids) that may take significant computational time.
   */
  template <typename ordinal_type, typename value_type> 
  class AlgebraicOrthogPolyExpansion : 
    public OrthogPolyExpansionBase<ordinal_type, value_type, 
				   Stokhos::StandardStorage<ordinal_type, value_type> > {
  public:

    typedef Stokhos::StandardStorage<ordinal_type, value_type> node_type;

    //! Constructor
    AlgebraicOrthogPolyExpansion(
      const Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor
    virtual ~AlgebraicOrthogPolyExpansion() {}
 
    // Operations
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
    AlgebraicOrthogPolyExpansion(const AlgebraicOrthogPolyExpansion&);

    // Prohibit Assignment
    AlgebraicOrthogPolyExpansion& operator=(const AlgebraicOrthogPolyExpansion& b);

  protected:

    //! Short-hand for Cijk
    typedef typename OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::Cijk_type Cijk_type;
    
  }; // class AlgebraicOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_AlgebraicOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_ALGEBRAICORTHOGPOLYEXPANSION_HPP
