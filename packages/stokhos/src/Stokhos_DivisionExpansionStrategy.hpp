// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
#define STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP

#include "Stokhos_StandardStorage.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"

namespace Stokhos {

  //! Strategy interface for computing PCE of a/b
  template <typename ordinal_type, typename value_type, typename node_type> 
  class DivisionExpansionStrategy {
  public:

    //! Constructor
    DivisionExpansionStrategy() {}

    //! Destructor
    virtual ~DivisionExpansionStrategy() {}
 
    // Division operation:  c = \alpha*(a/b) + beta*c
    virtual void divide(
      Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
      const value_type& alpha,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
      const value_type& beta) = 0;

  private:

    // Prohibit copying
    DivisionExpansionStrategy(const DivisionExpansionStrategy&);

    // Prohibit Assignment
    DivisionExpansionStrategy& operator=(const DivisionExpansionStrategy& b);
    
  }; // class DivisionExpansionStrategy

} // namespace Stokhos

#endif // STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
