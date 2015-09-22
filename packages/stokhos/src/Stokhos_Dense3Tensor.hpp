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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_DENSE3TENSOR_HPP
#define STOKHOS_DENSE3TENSOR_HPP

#include <ostream>

#include "Teuchos_Array.hpp"

namespace Stokhos {

  /*! 
   * \brief Data structure storing a dense 3-tensor C(i,j,k).
   */
  template <typename ordinal_type, typename value_type>
  class Dense3Tensor {
  public:
    
    //! Constructor
    Dense3Tensor(ordinal_type sz);
    
    //! Destructor
    ~Dense3Tensor();

    //! Return size
    ordinal_type size() const;
      
    //! Get value (i,j,k)
    const value_type& operator() (ordinal_type i, ordinal_type j, 
				  ordinal_type k) const;

    //! Get value (i,j,k)
    value_type& operator() (ordinal_type i, ordinal_type j, ordinal_type k);

    //! Return number of non-zero's in Cijk for a given k
    ordinal_type num_values(ordinal_type k) const;
      
    //! Get value (i,j,k) using sparse access
    void value(ordinal_type k, ordinal_type l, ordinal_type& i, ordinal_type& j,
	       value_type& c) const;
    
    //! Print tensor
    void print(std::ostream& os) const;

  private:

    // Prohibit copying
    Dense3Tensor(const Dense3Tensor&);

    // Prohibit Assignment
    Dense3Tensor& operator=(const Dense3Tensor& b);
    
  protected:

    //! Size of each dimension
    ordinal_type l;
    
    //! Dense tensor array
    Teuchos::Array<value_type> Cijk_values;

  }; // class Dense3Tensor

  template <typename ordinal_type, typename value_type>
  std::ostream& 
  operator << (std::ostream& os, 
	       const Dense3Tensor<ordinal_type, value_type>& Cijk) {
    Cijk.print(os);
    return os;
  }

} // namespace Stokhos

// Include template definitions
#include "Stokhos_Dense3TensorImp.hpp"

#endif // STOKHOS_DENSE3TENSOR_HPP
