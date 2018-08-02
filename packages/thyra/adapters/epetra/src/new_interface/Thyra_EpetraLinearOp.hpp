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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_EPETRA_LINEAR_OP_HPP
#define THYRA_EPETRA_LINEAR_OP_HPP

// Not directly needed in this file, but this way we made the
// macro HAVE_THYRA_EPETRA_REFACTOR available to files that include
// this header. This way, they do not need to include the config.h
// header manually. That's nice, because in the future we may deprecate
// and then remove the old interface, making the config.h file pointless.
// If that happens, we may remove it, and at that point all files including
// it would have to be updated. This was, only the adapters headers need to
// be updated.
#include "ThyraEpetraAdapters_config.h"

#include "Thyra_EpetraTypes.hpp"

#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

class Epetra_Operator;

namespace Thyra {

class EpetraVectorSpace;

/** \brief Concrete LinearOpBase subclass for Epetra_Operator.
 *
 * \todo Finish Documentation
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraLinearOp
  : virtual public LinearOpDefaultBase<double>,
    virtual public ScaledLinearOpBase<double>,
    virtual public RowStatLinearOpBase<double>,
    virtual public EpetraLinearOpBase
{
public:

  /** \name Constructors/initializers. */
  //@{

  /** \brief Construct to uninitialized. */
  EpetraLinearOp();

  /** \brief Initialize. */
  void initialize(const RCP<Epetra_Operator> &epetraOperator,
    const EOpTransp opTrans = NOTRANS,
    const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
    const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED);

  void initialize(
    const RCP<const VectorSpaceBase<double>> &rangeSpace,
    const RCP<const VectorSpaceBase<double>> &domainSpace,
    const RCP<Epetra_Operator> &epetraOperator,
    const EOpTransp opTrans = NOTRANS,
    const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
    const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED);

  /** \brief Initialize. */
  void constInitialize(const RCP<const Epetra_Operator> &epetraOperator,
    const EOpTransp opTrans = NOTRANS,
    const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
    const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED);

  void constInitialize(
    const RCP<const VectorSpaceBase<double>> &rangeSpace,
    const RCP<const VectorSpaceBase<double>> &domainSpace,
    const RCP<const Epetra_Operator> &epetraOperator,
    const EOpTransp opTrans = NOTRANS,
    const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
    const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED);

  /** \brief Get embedded non-const Epetra_Operator. */
  RCP<Epetra_Operator>
  getEpetraOperator();

  /** \brief Get embedded const Epetra_Operator. */
  RCP<const Epetra_Operator>
  getConstEpetraOperator() const;

  //@}

  /** \name Overridden from EpetraLinearOpBase */
  //@{

  /** \brief . */
  void getNonconstEpetraOpView(
    const Ptr<RCP<Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    );
  /** \brief . */
  void getEpetraOpView(
    const Ptr<RCP<const Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    ) const;

  //@}

  /** \name Public Overridden functions from LinearOpBase. */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<double>> range() const;

  /** \brief . */
  RCP<const VectorSpaceBase<double>> domain() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;
  
  //@}
  

protected:

  /** \name Protected Overridden functions from LinearOpBase. */
  //@{

  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &X_in,
    const Teuchos::Ptr<MultiVectorBase<double>> &Y_inout,
    const double alpha,
    const double beta
    ) const;

  //@}

  /** \name Protected member functions overridden from ScaledLinearOpBase. */
  //@{

  /** \brief . */
  virtual bool supportsScaleLeftImpl() const;

  /** \brief . */
  virtual bool supportsScaleRightImpl() const;

  /** \brief . */
  virtual void scaleLeftImpl(const VectorBase<double> &row_scaling);

  /** \brief . */
  virtual void scaleRightImpl(const VectorBase<double> &col_scaling);
  
  //@}

  /** \name Protected member functions overridden from RowStatLinearOpBase. */
  //@{

  /** \brief . */
  virtual bool rowStatIsSupportedImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat) const;

  /** \brief . */
  virtual void getRowStatImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat,
    const Ptr<VectorBase<double>> &rowStatVec) const;

  //@}

private:

  RCP<const EpetraVectorSpace>
  rangeSpace_;

  RCP<const EpetraVectorSpace>
  domainSpace_;

  Teuchos::ConstNonconstObjectContainer<Epetra_Operator>
  epetraOperator_;

  EAdjointEpetraOp  adjointSupport_;  // Whether applying the adjoint is supported
  EApplyEpetraOpAs  applyAs_;         // applyImpl may call Apply or ApplyInverse on the stored op.
  EOpTransp         opTrans_;         // applyImpl may transpose the underlying op before calling Apply (or ApplyInverse)

  template<class EpetraOperator_t>
  void initializeImpl(
    const RCP<const VectorSpaceBase<double>> &rangeSpace,
    const RCP<const VectorSpaceBase<double>> &domainSpace,
    const RCP<EpetraOperator_t> &epetraOperator,
    const EOpTransp opTrans,
    const EApplyEpetraOpAs applyAs,
    const EAdjointEpetraOp adjointSupport
  );

};


/** \brief Nonmmeber constructor for EpetraLinearOp.
 *
 * \relates EpetraLinearOp
 */
inline
RCP<EpetraLinearOp>
epetraLinearOp(
  const RCP<Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans = NOTRANS,
  const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED)
{
  const RCP<EpetraLinearOp> op (new EpetraLinearOp());
  op->initialize(epetraOperator, opTrans, applyAs);
  return op;
}

inline
RCP<EpetraLinearOp>
epetraLinearOp(
  const RCP<const VectorSpaceBase<double>> &rangeSpace,
  const RCP<const VectorSpaceBase<double>> &domainSpace,
  const RCP<Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans = NOTRANS,
  const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED)
{
  const RCP<EpetraLinearOp> op (new EpetraLinearOp());
  op->initialize(rangeSpace, domainSpace, epetraOperator, opTrans, applyAs);
  return op;
}


/** \brief Nonmmeber constructor for EpetraLinearOp.
 *
 * \relates EpetraLinearOp
 */
inline
RCP<const EpetraLinearOp>
constEpetraLinearOp(
  const RCP<const Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans = NOTRANS,
  const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED)
{
  const RCP<EpetraLinearOp> op( new EpetraLinearOp() );
  op->constInitialize(epetraOperator, opTrans, applyAs);
  return op;
}

inline
RCP<const EpetraLinearOp>
constEpetraLinearOp(
  const RCP<const VectorSpaceBase<double>> &rangeSpace,
  const RCP<const VectorSpaceBase<double>> &domainSpace,
  const RCP<const Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans = NOTRANS,
  const EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  const EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED)
{
  const RCP<EpetraLinearOp> op( new EpetraLinearOp() );
  op->constInitialize(rangeSpace, domainSpace, epetraOperator, opTrans, applyAs);
  return op;
}

}  // namespace Thyra

#endif	// THYRA_EPETRA_LINEAR_OP_HPP
