// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_TRIPLEPRODUCT_HPP
#define SACADO_PCE_TRIPLEPRODUCT_HPP

#include <vector>

namespace Sacado {

  namespace PCE {

    //! 3-tensor that stores C_{ijk} = < \Psi_i \Psi_j \Psi_k >
    template <typename BasisT>
    class TripleProduct {
    public:

      typedef BasisT basis_type;

      typedef typename BasisT::value_type value_type;

      //! Constructor
      TripleProduct(unsigned int degree);

      //! Copy constructor
      TripleProduct(const TripleProduct& tp);
      
      //! Destructor
      ~TripleProduct();

      //! Assignment
      TripleProduct& operator=(const TripleProduct& tp);
      
      //! Get value (i,j,k)
      const value_type& value(unsigned int i, unsigned int j, 
			      unsigned int k) const;
      
      //! Get norm-squared
      const value_type& norm_squared(unsigned int i) const;

      //! Return size
      unsigned int size() const { return l; }

      //! Return basis
      const BasisT& getBasis() const { return basis; }

      //! Resize to new dimension
      void resize(unsigned int degree);

    protected:

      //! Compute values
      void compute();
      
    protected:

      //! Size of each dimension
      unsigned int l;

      //! Basis
      BasisT basis;

      //! Cijk data
      std::vector<value_type> Cijk;

    }; // class Triple Product

  } // namespace PCE

} // namespace Sacado

// Include template definitions
#include "Sacado_PCE_TripleProductImp.hpp"

#endif // SACADO_PCE_TRIPLEPRODUCT_HPP
