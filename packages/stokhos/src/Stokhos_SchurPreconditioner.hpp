/*
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
*/



#ifndef STOKHOS_SCHURPRECONDITIONER_HPP
#define STOKHOS_SCHURPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class SchurPreconditioner : 
    public Stokhos::Operator<ordinal_type, value_type> {
  public:

    //! Constructor 
    SchurPreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & K, 
      const ordinal_type p, const ordinal_type m, const ordinal_type diag);
    
    //! Destructor
    virtual ~SchurPreconditioner(); 
    
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type prec_iters) const;
      
  protected:
    ordinal_type fact(ordinal_type n) const;
    ordinal_type size(ordinal_type n, ordinal_type m) const;

    const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & K;
    const ordinal_type p;
    const ordinal_type m;
    const ordinal_type diag;

  }; // class SchurPreconditioner

} // namespace Stokhos

#include "Stokhos_SchurPreconditionerImp.hpp"

#endif // STOKHOS_SCHURPRECONDITIONER_HPP

