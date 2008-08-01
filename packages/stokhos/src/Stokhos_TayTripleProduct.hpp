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
// Questions? Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_TAYTRIPLEPRODUCT_HPP
#define STOKHOS_TAYTRIPLEPRODUCT_HPP

namespace Stokhos {

  //! 3-tensor that stores C_{ijk} = < \Psi_i \Psi_j \Psi_k >
  template <typename BasisT>
  class TayTripleProduct {
  public:

    typedef BasisT basis_type;
    
    typedef typename BasisT::value_type value_type;
    
    //! Constructor
    TayTripleProduct();
    
    //! Destructor
    ~TayTripleProduct();

    //! Return number of non-zero's in Cijk for a given k
    unsigned int num_values(unsigned int k) const;
      
    //! Get value (i,j,k)
    void triple_value(unsigned int k, unsigned int l, 
		      unsigned int& i, unsigned int& j, value_type& c) const;
    
    //! Get norm-squared
    value_type norm_squared(unsigned int i) const;

  private:

    // Prohibit copying
    TayTripleProduct(const TayTripleProduct&);

    // Prohibit Assignment
    TayTripleProduct& operator=(const TayTripleProduct& b);

  }; // class Triple Product

} // namespace Stokhos

// Include template definitions
#include "Stokhos_TayTripleProductImp.hpp"

#endif // STOKHOS_TAYTRIPLEPRODUCT_HPP
