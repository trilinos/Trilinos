// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
