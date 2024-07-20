// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MEAN_BASED_DIVISION_EXPANSION_STRATEGY_HPP
#define STOKHOS_MEAN_BASED_DIVISION_EXPANSION_STRATEGY_HPP

#include "Stokhos_DivisionExpansionStrategy.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Stokhos {

  //! Strategy interface for computing PCE of a/b using only b[0]
  /*!
   * Such a strategy is only useful when the division occurs in a preconditioner
   */
  template <typename ordinal_type, typename value_type, typename node_type> 
  class MeanBasedDivisionExpansionStrategy :
    public DivisionExpansionStrategy<ordinal_type,value_type,node_type> {
  public:

    //! Constructor
    MeanBasedDivisionExpansionStrategy() {}

    //! Destructor
    virtual ~MeanBasedDivisionExpansionStrategy() {}
 
    // Division operation:  c = \alpha*(a/b) + beta*c
    virtual void divide(
      Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
      const value_type& alpha,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
      const value_type& beta);

  private:

    // Prohibit copying
    MeanBasedDivisionExpansionStrategy(
      const MeanBasedDivisionExpansionStrategy&);

    // Prohibit Assignment
    MeanBasedDivisionExpansionStrategy& operator=(
      const MeanBasedDivisionExpansionStrategy& b);
    
  }; // class DivisionExpansionStrategy

} // namespace Stokhos

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::MeanBasedDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
       const value_type& alpha,
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
       const value_type& beta)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::MeanBasedDivisionStrategy::divide()");
#endif
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = alpha*ca[i]/cb[0] + beta*cc[i];
}

#endif // STOKHOS_MEAN_BASED_DIVISION_EXPANSION_STRATEGY_HPP
