// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_OPERATOR_HPP
#define STOKHOS_OPERATOR_HPP

#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {
    
  template <typename ordinal_type, typename value_type>
  class Operator {      
  public:

    //! Constructor
    Operator() {} 
  
    //! Destructor
    virtual ~Operator() {}

    //! Returns the result of a Operator inverse applied to a Teuchos::SerialDenseMatrix Input in Result.
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type m) const = 0;
    
  }; // class Operator
  
} // namespace Stokhos

#endif // STOKHOS_OPERATOR_HPP
