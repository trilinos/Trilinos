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

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraVectorSpace.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

namespace Thyra {

// Constructors/initializers


EpetraLinearOp::EpetraLinearOp()
{}


void EpetraLinearOp::initialize(
  const RCP<Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans,
  const EApplyEpetraOpAs applyAs,
  const EAdjointEpetraOp adjointSupport)
{
  Teuchos::RCP<const VectorSpaceBase<double>> rangeSpace, domainSpace;
  rangeSpace  = epetraVectorSpace(Teuchos::rcpFromRef(epetraOperator->OperatorRangeMap()));
  domainSpace = epetraVectorSpace(Teuchos::rcpFromRef(epetraOperator->OperatorDomainMap()));
  
  initializeImpl(rangeSpace, domainSpace, epetraOperator, opTrans, applyAs, adjointSupport);
}

void EpetraLinearOp::initialize(
  const RCP<const VectorSpaceBase<double>> &rangeSpace,
  const RCP<const VectorSpaceBase<double>> &domainSpace,
  const RCP<Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans,
  const EApplyEpetraOpAs applyAs,
  const EAdjointEpetraOp adjointSupport)
{
  initializeImpl(rangeSpace, domainSpace, epetraOperator, opTrans, applyAs, adjointSupport);
}

void EpetraLinearOp::constInitialize(
  const RCP<const Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans,
  const EApplyEpetraOpAs applyAs,
  const EAdjointEpetraOp adjointSupport)
{
  Teuchos::RCP<const VectorSpaceBase<double>> rangeSpace, domainSpace;
  rangeSpace  = epetraVectorSpace(Teuchos::rcpFromRef(epetraOperator->OperatorRangeMap()));
  domainSpace = epetraVectorSpace(Teuchos::rcpFromRef(epetraOperator->OperatorDomainMap()));
  initializeImpl(rangeSpace, domainSpace, epetraOperator, opTrans, applyAs, adjointSupport);
}


void EpetraLinearOp::constInitialize(
  const RCP<const VectorSpaceBase<double>> &rangeSpace,
  const RCP<const VectorSpaceBase<double>> &domainSpace,
  const RCP<const Epetra_Operator> &epetraOperator,
  const EOpTransp opTrans,
  const EApplyEpetraOpAs applyAs,
  const EAdjointEpetraOp adjointSupport)
{
  initializeImpl(rangeSpace, domainSpace, epetraOperator, opTrans, applyAs, adjointSupport);
}


RCP<Epetra_Operator>
EpetraLinearOp::getEpetraOperator()
{
  return epetraOperator_.getNonconstObj();
}


RCP<const Epetra_Operator>
EpetraLinearOp::getConstEpetraOperator() const
{
  return epetraOperator_;
}


// Overridden from EpetraLinearOpBase

void EpetraLinearOp::getNonconstEpetraOpView(
  const Ptr<RCP<Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
  )
{
  *epetraOp = epetraOperator_.getNonconstObj();
  *epetraOpTransp = opTrans_;
  *epetraOpApplyAs = applyAs_;
  *epetraOpAdjointSupport = adjointSupport_;
}


void EpetraLinearOp::getEpetraOpView(
  const Ptr<RCP<const Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
  ) const
{
  *epetraOp = epetraOperator_.getConstObj();
  *epetraOpTransp = opTrans_;
  *epetraOpApplyAs = applyAs_;
  *epetraOpAdjointSupport = adjointSupport_;
}

// Public Overridden functions from LinearOpBase


RCP<const Thyra::VectorSpaceBase<double>>
EpetraLinearOp::range() const
{
  return rangeSpace_;
}


RCP<const Thyra::VectorSpaceBase<double>>
EpetraLinearOp::domain() const
{
  return domainSpace_;
}

// Overridden from Teuchos::Describable


std::string EpetraLinearOp::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{";
  if(Teuchos::nonnull(epetraOperator_.getConstObj())) {
    oss << "op=\'"<<typeName(*epetraOperator_)<<"\'";
    oss << ",rangeDim="<<this->range()->dim();
    oss << ",domainDim="<<this->domain()->dim();
  }
  else {
    oss << "op=NULL";
  }
  oss << "}";
  return oss.str();
}


void EpetraLinearOp::describe(
  FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::includesVerbLevel;
  using Teuchos::as;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Teuchos::describe;
  OSTab tab(out);
  if ( as<int>(verbLevel) == as<int>(Teuchos::VERB_LOW) || epetraOperator_.getConstObj().is_null()) {
    out << this->description() << std::endl;
  }
  else if (includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM)) {
    out
      << Teuchos::Describable::description()
      << "{"
      << "rangeDim=" << this->range()->dim()
      << ",domainDim=" << this->domain()->dim()
      << "}\n";
    OSTab tab2(out);
    if (Teuchos::nonnull(epetraOperator_.getConstObj())) {
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        out << "opTrans="<<toString(opTrans_)<<"\n";
        out << "applyAs="<<toString(applyAs_)<<"\n";
        out << "adjointSupport="<<toString(adjointSupport_)<<"\n";
        out << "op="<<typeName(*epetraOperator_.getConstObj())<<"\n";
      }
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
        OSTab tab3(out);
        RCP<const Epetra_CrsMatrix>
          csr_op = rcp_dynamic_cast<const Epetra_CrsMatrix>(epetraOperator_.getConstObj());
        if (!is_null(csr_op)) {
          csr_op->Print(out);
        }
      }
    } else {
      out << "op=NULL"<<"\n";
    }
  }
}



// Protected Overridden functions from LinearOpBase


bool EpetraLinearOp::opSupportedImpl(
  Thyra::EOpTransp M_trans) const
{
  if (is_null(epetraOperator_)) {
    return false;
  }

  return M_trans==NOTRANS || adjointSupport_==EPETRA_OP_ADJOINT_SUPPORTED;
}


void EpetraLinearOp::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<double> &X_in,
  const Teuchos::Ptr<Thyra::MultiVectorBase<double>> &Y_inout,
  const double alpha,
  const double beta
  ) const
{
  // Short names
  using Teuchos::rcpFromRef;
  using Teuchos::rcpFromPtr;
  using ST = Teuchos::ScalarTraits<double>;
  using ConverterE = EpetraOperatorVectorExtraction;

  // Sanity check
  const EOpTransp real_M_trans = real_trans(M_trans);
  TEUCHOS_TEST_FOR_EXCEPTION(!opSupportedImpl(real_M_trans), std::runtime_error,
                              "Error! Input value for M_trans is not supported.\n")

  // Get Epetra_MultiVector objects for X and Y
  const RCP<const Epetra_MultiVector> X = ConverterE::getConstEpetraMultiVector(rcpFromRef(X_in));
  const RCP<Epetra_MultiVector> Y = ConverterE::getEpetraMultiVector(rcpFromPtr(Y_inout));

  // NOTE: there are two transpose flags: the one passed in (M_trans), and whether this operator
  //       is wrapping the transpose of an Epetra_Operator (opTrans_).
  //       In order to figure out what transpose mode to set on the wrapped operator, we need
  //       to combine the two. After that, check if the resulting mode is already the one
  //       set in the operator: if yes, simply proceed with calling Apply(Inverse); if no,
  //       temporarily switch the UseTranspose flag of the operator, call Apply(Inverse),
  //       and finally restore the original flag in the operator.

  // Store the old setting of UseTranspose, so we can restore it after the call to Apply(...)
  bool use_transpose = epetraOperator_.getConstObj()->UseTranspose();

  // Combine the requested transpose mode with the mode the underlying Epetra op is applied
  // and establishing what needs to be set on the underlying operator
  const EOpTransp overallTrans= real_trans(trans_trans(opTrans_,M_trans));
  bool needs_trans = (overallTrans==TRANS);

  // This is tricky: this method is marked const, but if use_transpose!=overallTrans,
  // we need to call SetUseTranspose on the Epetra_Operator. If this wrapper was created
  // from a const object, the runtime check will not allow that call. Hence, if the operator is const,
  // we do a rcp_const_cast. Yes, it is ugly, but we can't avoid it.
  // If Epetra were designed with Epetra_Operator accepting the transpose flag
  // at the Apply(...) call time instead of requiring one to set it in the
  // operator before Apply(...), then we would be fine.

  Teuchos::RCP<Epetra_Operator> nonConstOp;
  if (use_transpose!=needs_trans) {
    if (epetraOperator_.isConst()) {
      nonConstOp = Teuchos::rcp_const_cast<Epetra_Operator>(epetraOperator_.getConstObj());
    } else {
      nonConstOp = epetraOperator_.getNonconstObj();
    }
    nonConstOp->SetUseTranspose(needs_trans);
  }

  RCP<Epetra_MultiVector> Y_copy;
  if (beta!=ST::zero()) {
    // Epetra does not offer an apply-and-update kind of method. Hence, if
    // beta!=0, we need to store a copy of the original output vector
    Y_copy = Teuchos::rcp( new Epetra_MultiVector(*Y) );
  }

  // NOTE: we could have done this check at the beginning, and avoid all the conversions,
  //       since we could have just called Y_inout->scale(beta). However, this would
  //       allow buggy code to execute; for instance, if the value of M_trans was not supported
  if (alpha==ST::zero()) {
    Y->PutScalar(0.0);
  } else {
    // Apply the operator
    if (applyAs_==EPETRA_OP_APPLY_APPLY) {
      epetraOperator_.getConstObj()->Apply(*X,*Y);
    } else {
      epetraOperator_.getConstObj()->ApplyInverse(*X,*Y);
    }

    // If needed, scale by alpha
    if (alpha!=ST::one()) {
      Y->Scale(alpha);
    }
  }

  // If needed, add beta*Y_old
  if (beta!=ST::zero()) {
    Y->Update(beta,*Y_copy,1.0);
  }

  // Finally, if needed, restore the original transpose flag
  if (use_transpose!=needs_trans) {
    nonConstOp->SetUseTranspose(use_transpose);
  }
}

// Protected member functions overridden from ScaledLinearOpBase


bool EpetraLinearOp::supportsScaleLeftImpl() const
{
  return Teuchos::nonnull(epetraOperator_);
}


bool EpetraLinearOp::supportsScaleRightImpl() const
{
  return Teuchos::nonnull(epetraOperator_);
}


void
EpetraLinearOp::
scaleLeftImpl(const VectorBase<double> &row_scaling_in)
{
  using Teuchos::rcpFromRef;

  const RCP<const Epetra_Vector> row_scaling =
    EpetraOperatorVectorExtraction::getConstEpetraVector(rcpFromRef(row_scaling_in));

  const RCP<Epetra_RowMatrix> rowMatrix =
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetraOperator_.getNonconstObj(),true);

  rowMatrix->LeftScale(*row_scaling);
}


void
EpetraLinearOp::
scaleRightImpl(const VectorBase<double> &col_scaling_in)
{
  using Teuchos::rcpFromRef;

  const RCP<const Epetra_Vector> col_scaling =
    EpetraOperatorVectorExtraction::getConstEpetraVector(rcpFromRef(col_scaling_in));

  const RCP<Epetra_RowMatrix> rowMatrix =
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetraOperator_.getNonconstObj(),true);

  rowMatrix->RightScale(*col_scaling);
}

// Protected member functions overridden from RowStatLinearOpBase


bool EpetraLinearOp::
rowStatIsSupportedImpl(
  const RowStatLinearOpBaseUtils::ERowStat rowStat) const
{
  if (Teuchos::is_null(epetraOperator_)) {
    return false;
  }

  // We can only get row stats if the operator is derived from Epetra_RowMatrix
  const RCP<const Epetra_RowMatrix> matrix = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(epetraOperator_.getConstObj());
  if (matrix.is_null()) {
    return false;
  }

  // NOTE: as of today (7/30/18) these are all the supported options,
  //       so it may be tempting to just return true without a switch.
  //       However, having a switch protects in case of future additions
  //       to the ERowStat enum.

  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:    // Fallthrough
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:        // Fallthrough
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:    // Fallthrough
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      return true;
    default:
      return false;
  }

  TEUCHOS_UNREACHABLE_RETURN(false);
}


void EpetraLinearOp::getRowStatImpl(
  const RowStatLinearOpBaseUtils::ERowStat rowStat,
  const Ptr<VectorBase<double>> &rowStatVec_in
  ) const
{
  TEUCHOS_ASSERT(nonnull(epetraOperator_));
  TEUCHOS_ASSERT(nonnull(rowStatVec_in));
  const RCP<Epetra_Vector> eRowStat = EpetraOperatorVectorExtraction::getEpetraVector(rcpFromPtr(rowStatVec_in));
  const RCP<const Epetra_RowMatrix> matrix = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(epetraOperator_.getConstObj(),true);

  // NOTE: Epetra offers InvRowSums and InvColSums methods. So for both Inv and non-Inv cases,
  //       we call the Inv*Sums. Then we invert element-wise in case non-Inv was asked.
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:  // Fallthrough
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      matrix->InvRowSums(*eRowStat);
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:  // Fallthrough
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      matrix->InvColSums(*eRowStat);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error! Unsupported row stat.\n");
  }

  // If the requests was for non-inverse sums, we need to invert element-wise the vector
  if ( (rowStat == RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM) ||
       (rowStat == RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM) ) {
    eRowStat->Reciprocal(*eRowStat);
  }
}


// private


template<class EpetraOperator_t>
void EpetraLinearOp::initializeImpl(
  const RCP<const VectorSpaceBase<double>> &rangeSpace,
  const RCP<const VectorSpaceBase<double>> &domainSpace,
  const RCP<EpetraOperator_t> &epetraOperator,
  const EOpTransp opTrans,
  const EApplyEpetraOpAs applyAs,
  const EAdjointEpetraOp adjointSupport)
{
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(rangeSpace));
  TEUCHOS_ASSERT(nonnull(domainSpace));
  TEUCHOS_ASSERT(nonnull(epetraOperator));
#endif

  rangeSpace_  = Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(rangeSpace);
  domainSpace_ = Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(domainSpace);

  TEUCHOS_TEST_FOR_EXCEPTION(rangeSpace_.is_null(), std::runtime_error,
                             "Error! Could not cast input range space to EpetraVectorSpace.\n");
  TEUCHOS_TEST_FOR_EXCEPTION(domainSpace_.is_null(), std::runtime_error,
                             "Error! Could not cast input domain space to EpetraVectorSpace.\n");

  epetraOperator_ = epetraOperator;
  applyAs_ = applyAs;
  opTrans_ = opTrans;
  adjointSupport_ = adjointSupport;

  // Now check if this operator supports adjoint
  // Combine TRANS with the mode the underlying Epetra op is applied
  // and establishing what needs to be set on the underlying operator
  const EOpTransp overallTrans= real_trans(trans_trans(opTrans_,TRANS));
  bool is_overallTrans_trans = (overallTrans==TRANS);

  // Store the current value of UseTranspose so we can restore it later
  bool use_transpose = epetraOperator_.getConstObj()->UseTranspose();

  if (use_transpose==is_overallTrans_trans) {
    adjointSupport_ = EPETRA_OP_ADJOINT_SUPPORTED;
  } else {
    // This is tricky, cause the stored operator may be const, in which case
    // we cannot call that method. Hence, if the operator is const, we do a
    // rcp_const_cast. Yes, it is ugly, but we can't avoid it.
    // If Epetra were designed with Epetra_Operator accepting the transpose flag
    // at the Apply(...) call time instead of requiring one to set it in the
    // operator before Apply(...), then we would be fine.
    Teuchos::RCP<Epetra_Operator> nonConstOp;
    if (epetraOperator_.isConst()) {
      nonConstOp = Teuchos::rcp_const_cast<Epetra_Operator>(epetraOperator_.getConstObj());
    } else {
      nonConstOp = epetraOperator_.getNonconstObj();
    }

    int errorCode;
    try {
      errorCode = nonConstOp->SetUseTranspose(is_overallTrans_trans);
    } catch (...) {
      // If we catch an exception, we ASSUME is because the operator does not support !use_transpose,
      // and the developer who implemented it decided it was appropriate to throw in such case.
      // Here we are not really interested in setting the value of UseTranspose, but just in checking
      // whether that functionality is supported.
      errorCode = 1;
    }

    // Immediately restore the original setting on the underlying operator
    nonConstOp->SetUseTranspose(use_transpose);

    if (errorCode==0) {
      // The underlying operator supports the requested transpose flag
      adjointSupport_ = EPETRA_OP_ADJOINT_SUPPORTED;
    } else {
      adjointSupport_ = EPETRA_OP_ADJOINT_UNSUPPORTED;
    }
  }

}

} // namespace Thyra
