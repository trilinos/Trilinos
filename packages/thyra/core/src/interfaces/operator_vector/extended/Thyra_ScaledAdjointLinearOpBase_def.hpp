// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  TEUCHOS_TEST_FOR_EXCEPT( scalar==NULL );
  TEUCHOS_TEST_FOR_EXCEPT( transp==NULL );
  TEUCHOS_TEST_FOR_EXCEPT( origOp==NULL );
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
  TEUCHOS_TEST_FOR_EXCEPT( scalar==NULL );
  TEUCHOS_TEST_FOR_EXCEPT( transp==NULL );
  TEUCHOS_TEST_FOR_EXCEPT( origOp==NULL );
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
