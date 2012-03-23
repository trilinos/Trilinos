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
