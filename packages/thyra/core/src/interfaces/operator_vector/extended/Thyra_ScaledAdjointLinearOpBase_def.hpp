// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_HPP
#define THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_HPP

#include "Thyra_ScaledAdjointLinearOpBase_decl.hpp"
#include "Thyra_LinearOpBase.hpp"


template<class Scalar>
void Thyra::unwrap(
  const LinearOpBase<Scalar> &Op,
  Scalar *scalar,
  EOpTransp *transp,
  const LinearOpBase<Scalar>* *origOp
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( scalar==NULL );
  TEST_FOR_EXCEPT( transp==NULL );
  TEST_FOR_EXCEPT( origOp==NULL );
#endif
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  const ScaledAdjointLinearOpBase<Scalar>
    *saOp = dynamic_cast<const ScaledAdjointLinearOpBase<Scalar>*>(&Op);
  if(saOp) {
    *scalar = saOp->overallScalar();
    *transp = saOp->overallTransp();
    *origOp = &*saOp->getOrigOp();
  }
  else {
    *scalar = ST::one();
    *transp = NOTRANS;
    *origOp = &Op;
  }
}


template<class Scalar>
void Thyra::unwrap(
  const RCP<const LinearOpBase<Scalar> > &Op,
  Scalar *scalar,
  EOpTransp *transp,
  RCP<const LinearOpBase<Scalar> > *origOp
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( scalar==NULL );
  TEST_FOR_EXCEPT( transp==NULL );
  TEST_FOR_EXCEPT( origOp==NULL );
#endif
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  RCP<const ScaledAdjointLinearOpBase<Scalar> >
    saOp = Teuchos::rcp_dynamic_cast<const ScaledAdjointLinearOpBase<Scalar> >(Op);
  if(saOp.get()) {
    *scalar = saOp->overallScalar();
    *transp = saOp->overallTransp();
    *origOp = saOp->getOrigOp();
  }
  else {
    *scalar = ST::one();
    *transp = NOTRANS;
    *origOp = Op;
  }
}


//
// Explicit instant macro
//

#define THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_INSTANT(SCALAR) \
  \
  template void unwrap( \
    const LinearOpBase<SCALAR > &Op, \
    SCALAR  *scalar, \
    EOpTransp *transp, \
    const LinearOpBase<SCALAR >* *origOp \
    ); \
   \
  template void unwrap( \
    const RCP<const LinearOpBase<SCALAR > > &Op, \
    SCALAR  *scalar, \
    EOpTransp *transp, \
    RCP<const LinearOpBase<SCALAR > > *origOp \
    );


#endif	// THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_HPP
