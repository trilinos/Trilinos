// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_BLOCKPRECONDITIONER_HPP
#define STOKHOS_BLOCKPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class BlockPreconditioner : 
    public Stokhos::Operator<ordinal_type, value_type> {
  public:

    //! Constructor 
    BlockPreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & K, 
      const ordinal_type p, const ordinal_type m);
  
    //! Destructor
    virtual ~BlockPreconditioner(); 
    
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type m) const;
   
  protected:

    ordinal_type facto(ordinal_type n) const;
    ordinal_type siz (ordinal_type n, ordinal_type m) const;

    const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & K;
    const ordinal_type p;
    const ordinal_type m;

  }; // class BlockPreconditioner

} // namespace Stokhos

#include "Stokhos_BlockPreconditionerImp.hpp"

#endif // STOKHOS_BLOCKPRECONDITIONER_HPP

