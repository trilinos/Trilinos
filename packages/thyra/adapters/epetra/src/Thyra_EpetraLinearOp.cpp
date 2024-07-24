// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h" // Printing and absolute row sums only!


namespace Thyra {


// Constructors / initializers / accessors


EpetraLinearOp::EpetraLinearOp()
  :isFullyInitialized_(false),
   opTrans_(NOTRANS),
   applyAs_(EPETRA_OP_APPLY_APPLY),
   adjointSupport_(EPETRA_OP_ADJOINT_UNSUPPORTED)
{}


void EpetraLinearOp::initialize(
  const RCP<Epetra_Operator> &op,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport,
  const RCP< const VectorSpaceBase<double> > &range_in,
  const RCP< const VectorSpaceBase<double> > &domain_in
  )
{
  
  using Teuchos::rcp_dynamic_cast;
  typedef SpmdVectorSpaceBase<double> SPMDVSB;

  // Validate input, allocate spaces, validate ...
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(op), std::invalid_argument,
    "Thyra::EpetraLinearOp::initialize(...): Error!" );
  // ToDo: Validate spmdRange, spmdDomain against op maps!
#endif

  RCP<const SPMDVSB> l_spmdRange;
  if(!is_null(range_in))
    l_spmdRange = rcp_dynamic_cast<const SPMDVSB>(range_in,true);
  else
    l_spmdRange = ( applyAs==EPETRA_OP_APPLY_APPLY
      ? allocateRange(op,opTrans) : allocateDomain(op,opTrans) );

  RCP<const SPMDVSB> l_spmdDomain;
  if(!is_null(domain_in))
    l_spmdDomain = rcp_dynamic_cast<const SPMDVSB>(domain_in,true);
  else
    l_spmdDomain = ( applyAs==EPETRA_OP_APPLY_APPLY
      ? allocateDomain(op,opTrans) : allocateRange(op,opTrans) );
  
  // Set data (no exceptions should be thrown now)
  isFullyInitialized_ = true;
  op_ = op;
  rowMatrix_ = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op_);
  opTrans_ = opTrans;
  applyAs_ = applyAs;
  adjointSupport_ = adjointSupport;
  range_ = l_spmdRange;
  domain_ = l_spmdDomain;

}


void EpetraLinearOp::partiallyInitialize(
  const RCP<const VectorSpaceBase<double> > &range_in,
  const RCP<const VectorSpaceBase<double> > &domain_in,
  const RCP<Epetra_Operator> &op,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport
  )
{
  
  using Teuchos::rcp_dynamic_cast;
  typedef SpmdVectorSpaceBase<double> SPMDVSB;

  // Validate input, allocate spaces, validate ...
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(range_in), std::invalid_argument,
    "Thyra::EpetraLinearOp::partiallyInitialize(...): Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(domain_in), std::invalid_argument,
    "Thyra::EpetraLinearOp::partiallyInitialize(...): Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(op), std::invalid_argument,
    "Thyra::EpetraLinearOp::partiallyInitialize(...): Error!" );
#endif

  RCP<const SPMDVSB>
    l_spmdRange = rcp_dynamic_cast<const SPMDVSB>(range_in,true);
  RCP<const SPMDVSB>
    l_spmdDomain = rcp_dynamic_cast<const SPMDVSB>(domain_in,true);
  
  // Set data (no exceptions should be thrown now)
  isFullyInitialized_ = false;
  op_ = op;
  rowMatrix_ = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op_);
  opTrans_ = opTrans;
  applyAs_ = applyAs;
  adjointSupport_ = adjointSupport;
  range_ = l_spmdRange;
  domain_ = l_spmdDomain;
  
}


void EpetraLinearOp::setFullyInitialized(bool /* isFullyInitialized */)
{
  // ToDo: Validate that everything matches up!
  isFullyInitialized_ = true;
}


void EpetraLinearOp::uninitialize(
  RCP<Epetra_Operator> *op,
  EOpTransp *opTrans,
  EApplyEpetraOpAs *applyAs,
  EAdjointEpetraOp *adjointSupport,
  RCP<const VectorSpaceBase<double> > *range_out,
  RCP<const VectorSpaceBase<double> > *domain_out
  )
{

  if(op) *op = op_;
  if(opTrans) *opTrans = opTrans_;
  if(applyAs) *applyAs = applyAs_;
  if(adjointSupport) *adjointSupport = adjointSupport_;
  if(range_out) *range_out = range_;
  if(domain_out) *domain_out = domain_;

  isFullyInitialized_ = false;
  op_ = Teuchos::null;
  rowMatrix_ = Teuchos::null;
  opTrans_ = NOTRANS;
  applyAs_ = EPETRA_OP_APPLY_APPLY;
  adjointSupport_ = EPETRA_OP_ADJOINT_SUPPORTED;
  range_ = Teuchos::null;
  domain_ = Teuchos::null;

}


RCP< const SpmdVectorSpaceBase<double> >
EpetraLinearOp::spmdRange() const
{
  return range_;
}


RCP< const SpmdVectorSpaceBase<double> >
EpetraLinearOp::spmdDomain() const
{
  return domain_;
}


RCP<Epetra_Operator>
EpetraLinearOp::epetra_op() 
{
  return op_;
}


RCP<const Epetra_Operator>
EpetraLinearOp::epetra_op() const 
{
  return op_;
}


// Overridden from EpetraLinearOpBase


void EpetraLinearOp::getNonconstEpetraOpView(
  const Ptr<RCP<Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
  )
{
  *epetraOp = op_;
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
  *epetraOp = op_;
  *epetraOpTransp = opTrans_;
  *epetraOpApplyAs = applyAs_;
  *epetraOpAdjointSupport = adjointSupport_;
}


// Overridden from LinearOpBase


RCP<const VectorSpaceBase<double> >
EpetraLinearOp::range() const
{
  return range_;
}


RCP<const VectorSpaceBase<double> >
EpetraLinearOp::domain() const
{
  return domain_;
}


RCP<const LinearOpBase<double> >
EpetraLinearOp::clone() const
{
  assert(0); // ToDo: Implement when needed
  return Teuchos::null;
}


// Overridden from Teuchos::Describable


std::string EpetraLinearOp::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{";
  if(op_.get()) {
    oss << "op=\'"<<typeName(*op_)<<"\'";
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
  if ( as<int>(verbLevel) == as<int>(Teuchos::VERB_LOW) || is_null(op_)) {
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
    if (op_.get()) {
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        out << "opTrans="<<toString(opTrans_)<<"\n";
        out << "applyAs="<<toString(applyAs_)<<"\n";
        out << "adjointSupport="<<toString(adjointSupport_)<<"\n";
        out << "op="<<typeName(*op_)<<"\n";
      }
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
        OSTab tab3(out);
        RCP<const Epetra_CrsMatrix>
          csr_op = rcp_dynamic_cast<const Epetra_CrsMatrix>(op_);
        if (!is_null(csr_op)) {
          csr_op->Print(out);
        }
      }
    }
    else {
      out << "op=NULL"<<"\n";
    }
  }
}


// protected


// Protected member functions overridden from LinearOpBase


bool EpetraLinearOp::opSupportedImpl(EOpTransp M_trans) const
{
  if (!isFullyInitialized_)
    return false;
  return ( M_trans == NOTRANS
    ? true : adjointSupport_==EPETRA_OP_ADJOINT_SUPPORTED );
}


void EpetraLinearOp::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X_in,
  const Ptr<MultiVectorBase<double> > &Y_inout,
  const double alpha,
  const double beta
  ) const
{

  THYRA_FUNC_TIME_MONITOR("Thyra::EpetraLinearOp::euclideanApply");

  const EOpTransp real_M_trans = real_trans(M_trans);

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!isFullyInitialized_);
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "EpetraLinearOp::euclideanApply(...)", *this, M_trans, X_in, Y_inout
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    real_M_trans==TRANS && adjointSupport_==EPETRA_OP_ADJOINT_UNSUPPORTED,
    Exceptions::OpNotSupported,
    "EpetraLinearOp::apply(...): *this was informed that adjoints "
    "are not supported when initialized." 
    );
#endif

  const RCP<const VectorSpaceBase<double> > XY_domain = X_in.domain();
  const int numCols = XY_domain->dim();
 
  //
  // Get Epetra_MultiVector objects for the arguments
  //
  // 2007/08/18: rabartl: Note: After profiling, I found that calling the more
  // general functions get_Epetra_MultiVector(...) was too slow. These
  // functions must ensure that memory is being remembered efficiently and the
  // use of extra data with the RCP and other things is slow.
  //
  RCP<const Epetra_MultiVector> X;
  RCP<Epetra_MultiVector> Y;
  {
    THYRA_FUNC_TIME_MONITOR_DIFF(
      "Thyra::EpetraLinearOp::euclideanApply: Convert MultiVectors", MultiVectors);
    // X
    X = get_Epetra_MultiVector(
      real_M_trans==NOTRANS ? getDomainMap() : getRangeMap(), X_in );
    // Y
    if( beta == 0 ) {
      Y = get_Epetra_MultiVector(
        real_M_trans==NOTRANS ? getRangeMap() : getDomainMap(), *Y_inout );
    }
  }

  //
  // Set the operator mode
  //

  /* We need to save the transpose state here, and then reset it after 
   * application. The reason for this is that if we later apply the 
   * operator outside Thyra (in Aztec, for instance), it will remember
   * the transpose flag set here. */
  bool oldState = op_->UseTranspose();
  op_->SetUseTranspose(
    real_trans(trans_trans(opTrans_,M_trans)) == NOTRANS ? false : true );

  //
  // Perform the apply operation
  //
  {
    THYRA_FUNC_TIME_MONITOR_DIFF(
      "Thyra::EpetraLinearOp::euclideanApply: Apply", Apply);
    if( beta == 0.0 ) {
      // Y = M * X
      if( applyAs_ == EPETRA_OP_APPLY_APPLY ) {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta==0): Apply",
          ApplyApply);
        op_->Apply( *X, *Y );
      }
      else if( applyAs_ == EPETRA_OP_APPLY_APPLY_INVERSE ) {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta==0): ApplyInverse",
          ApplyApplyInverse);
        op_->ApplyInverse( *X, *Y );
      }
      else {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
      }
      // Y = alpha * Y
      if( alpha != 1.0 ) {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta==0): Scale Y",
          Scale);
        Y->Scale(alpha);
      }
    }
    else {  // beta != 0.0
      // Y_inout = beta * Y_inout
      if(beta != 0.0) {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta!=0): Scale Y",
          Scale);
        scale( beta, Y_inout );
      }
      else {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta!=0): Y=0",
          Apply2);
        assign( Y_inout, 0.0 );
      }
      // T = M * X
      Epetra_MultiVector T(op_->OperatorRangeMap(), numCols, false);
      // NOTE: Above, op_->OperatorRange() will be right for either
      // non-transpose or transpose because we have already set the
      // UseTranspose flag correctly.
      if( applyAs_ == EPETRA_OP_APPLY_APPLY ) {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta!=0): Apply",
          Apply2);
        op_->Apply( *X, T );
      }
      else if( applyAs_ == EPETRA_OP_APPLY_APPLY_INVERSE ) {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta!=0): ApplyInverse",
          ApplyInverse);
        op_->ApplyInverse( *X, T );
      }
      else {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
      }
      // Y_inout += alpha * T
      {
        THYRA_FUNC_TIME_MONITOR_DIFF(
          "Thyra::EpetraLinearOp::euclideanApply: Apply(beta!=0): Update Y",
          Update);
        update(
          alpha,
          *create_MultiVector(
            Teuchos::rcp(&Teuchos::getConst(T),false),
            Y_inout->range(),
            XY_domain
            ),
          Y_inout
          );
      }
    }
  }

  // Reset the transpose state
  op_->SetUseTranspose(oldState);

  // 2009/04/14: ToDo: This will not reset the transpose flag correctly if an
  // exception is thrown!

}


// Protected member functions overridden from ScaledLinearOpBase


bool EpetraLinearOp::supportsScaleLeftImpl() const
{
  return nonnull(rowMatrix_);
}


bool EpetraLinearOp::supportsScaleRightImpl() const
{
  return nonnull(rowMatrix_);
}


void EpetraLinearOp::scaleLeftImpl(const VectorBase<double> &row_scaling_in)
{
  using Teuchos::rcpFromRef;
  const RCP<const Epetra_Vector> row_scaling =
    get_Epetra_Vector(getRangeMap(), rcpFromRef(row_scaling_in));
  rowMatrix_->LeftScale(*row_scaling);
}


void EpetraLinearOp::scaleRightImpl(const VectorBase<double> &col_scaling_in)
{
  using Teuchos::rcpFromRef;
  const RCP<const Epetra_Vector> col_scaling =
    get_Epetra_Vector(getDomainMap(), rcpFromRef(col_scaling_in));
  rowMatrix_->RightScale(*col_scaling);
}


// Protected member functions overridden from RowStatLinearOpBase


bool EpetraLinearOp::rowStatIsSupportedImpl(
  const RowStatLinearOpBaseUtils::ERowStat rowStat) const
{
  if (is_null(rowMatrix_)) {
    return false;
  }
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      return true;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  TEUCHOS_UNREACHABLE_RETURN(false);
}


void EpetraLinearOp::getRowStatImpl(
  const RowStatLinearOpBaseUtils::ERowStat rowStat,
  const Ptr<VectorBase<double> > &rowStatVec_in
  ) const
{
  using Teuchos::rcpFromPtr;
  const RCP<Epetra_Vector> rowStatVec =
    get_Epetra_Vector(getRangeMap(), rcpFromPtr(rowStatVec_in));
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
      rowMatrix_->InvRowSums(*rowStatVec);
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      // compute absolute row sum
      computeAbsRowSum(*rowStatVec);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}


// Allocators for domain and range spaces


RCP< const SpmdVectorSpaceBase<double> > 
EpetraLinearOp::allocateDomain(
  const RCP<Epetra_Operator> &op,
  EOpTransp /* op_trans */
  ) const
{
  return Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
    create_VectorSpace(Teuchos::rcp(&op->OperatorDomainMap(),false))
    );
  // ToDo: What about the transpose argument???, test this!!!
}


RCP<const SpmdVectorSpaceBase<double> > 
EpetraLinearOp::allocateRange(
  const RCP<Epetra_Operator> &op,
  EOpTransp /* op_trans */
  ) const
{
  return Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
    create_VectorSpace(Teuchos::rcp(&op->OperatorRangeMap(),false))
    );
  // ToDo: What about the transpose argument???, test this!!!
}

// private


const Epetra_Map& EpetraLinearOp::getRangeMap() const
{
  return ( applyAs_ == EPETRA_OP_APPLY_APPLY
    ? op_->OperatorRangeMap() : op_->OperatorDomainMap() );
  // ToDo: What about the transpose argument???, test this!!!
}


const Epetra_Map& EpetraLinearOp::getDomainMap() const
{
  return ( applyAs_ == EPETRA_OP_APPLY_APPLY
    ? op_->OperatorDomainMap() : op_->OperatorRangeMap() );
  // ToDo: What about the transpose argument???, test this!!!
}

void EpetraLinearOp::computeAbsRowSum(Epetra_Vector & x) const
{
  TEUCHOS_ASSERT(!is_null(rowMatrix_));

  RCP<Epetra_CrsMatrix> crsMatrix 
    = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(rowMatrix_);

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(crsMatrix),
    Exceptions::OpNotSupported,
    "EpetraLinearOp::computeAbsRowSum(...): wrapped matrix must be of type "
    "Epetra_CrsMatrix for this method. Other operator types are not supported."
    );

  //
  // Put inverse of the sum of absolute values of the ith row of A in x[i].
  // (this is a modified copy of Epetra_CrsMatrix::InvRowSums)
  //

  if (crsMatrix->Filled()) {
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(crsMatrix),
      std::invalid_argument,
      "EpetraLinearOp::computeAbsRowSum(...): Epetra_CrsMatrix must be filled"
    );
  } 
  int i, j;
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double * xp = (double*)x.Values();
  if (crsMatrix->Graph().RangeMap().SameAs(x.Map()) && crsMatrix->Exporter() != 0) {
    Epetra_Vector x_tmp(crsMatrix->RowMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for (i=0; i < crsMatrix->NumMyRows(); i++) {
      int      NumEntries = 0;
      double * RowValues  = 0;
      crsMatrix->ExtractMyRowView(i,NumEntries,RowValues);
      for (j=0; j < NumEntries; j++)  x_tmp_p[i] += std::abs(RowValues[j]);
    }
    TEUCHOS_TEST_FOR_EXCEPT(0!=x.Export(x_tmp, *crsMatrix->Exporter(), Add)); //Export partial row sums to x.
  }
  else if (crsMatrix->Graph().RowMap().SameAs(x.Map())) {
    for (i=0; i < crsMatrix->NumMyRows(); i++) {
      int      NumEntries = 0;
      double * RowValues  = 0;
      crsMatrix->ExtractMyRowView(i,NumEntries,RowValues);
      double scale = 0.0;
      for (j=0; j < NumEntries; j++) scale += std::abs(RowValues[j]);
      xp[i] = scale;
    }
  }
  else { // x.Map different than both crsMatrix->Graph().RowMap() and crsMatrix->Graph().RangeMap()
    TEUCHOS_TEST_FOR_EXCEPT(true); // The map of x must be the RowMap or RangeMap of A.
  }
}

}	// end namespace Thyra


// Nonmembers


Teuchos::RCP<Thyra::EpetraLinearOp>
Thyra::nonconstEpetraLinearOp()
{
  return Teuchos::rcp(new EpetraLinearOp());
}


Teuchos::RCP<Thyra::EpetraLinearOp>
Thyra::partialNonconstEpetraLinearOp(
  const RCP<const VectorSpaceBase<double> > &range,
  const RCP<const VectorSpaceBase<double> > &domain,
  const RCP<Epetra_Operator> &op,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport
  )
{
  RCP<EpetraLinearOp> thyraEpetraOp = Teuchos::rcp(new EpetraLinearOp());
  thyraEpetraOp->partiallyInitialize(
    range, domain,op,opTrans, applyAs, adjointSupport
    );
  return thyraEpetraOp;
}


Teuchos::RCP<Thyra::EpetraLinearOp>
Thyra::nonconstEpetraLinearOp(
  const RCP<Epetra_Operator> &op,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport,
  const RCP< const VectorSpaceBase<double> > &range,
  const RCP< const VectorSpaceBase<double> > &domain
  )
{
  RCP<EpetraLinearOp> thyraEpetraOp = Teuchos::rcp(new EpetraLinearOp());
  thyraEpetraOp->initialize(
    op,opTrans, applyAs, adjointSupport, range, domain
    );
  return thyraEpetraOp;
}


Teuchos::RCP<const Thyra::EpetraLinearOp>
Thyra::epetraLinearOp(
  const RCP<const Epetra_Operator> &op,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport,
  const RCP<const VectorSpaceBase<double> > &range,
  const RCP<const VectorSpaceBase<double> > &domain
  )
{
  RCP<EpetraLinearOp> thyraEpetraOp = Teuchos::rcp(new EpetraLinearOp());
  thyraEpetraOp->initialize(
    Teuchos::rcp_const_cast<Epetra_Operator>(op), // Safe cast due to return type!
    opTrans, applyAs, adjointSupport, range, domain
    );
  return thyraEpetraOp;
}


Teuchos::RCP<Thyra::EpetraLinearOp>
Thyra::nonconstEpetraLinearOp(
  const RCP<Epetra_Operator> &op,
  const std::string &label,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport,
  const RCP<const VectorSpaceBase<double> > &range,
  const RCP<const VectorSpaceBase<double> > &domain
  )
{
  RCP<EpetraLinearOp> thyraEpetraOp = Teuchos::rcp(new EpetraLinearOp());
  thyraEpetraOp->initialize(
    op,opTrans, applyAs, adjointSupport, range, domain
    );
  thyraEpetraOp->setObjectLabel(label);
  return thyraEpetraOp;
}


Teuchos::RCP<const Thyra::EpetraLinearOp>
Thyra::epetraLinearOp(
  const RCP<const Epetra_Operator> &op,
  const std::string &label,
  EOpTransp opTrans,
  EApplyEpetraOpAs applyAs,
  EAdjointEpetraOp adjointSupport,
  const RCP< const SpmdVectorSpaceBase<double> > &range,
  const RCP< const SpmdVectorSpaceBase<double> > &domain
  )
{
  RCP<EpetraLinearOp> thyraEpetraOp = Teuchos::rcp(new EpetraLinearOp());
  thyraEpetraOp->initialize(
    Teuchos::rcp_const_cast<Epetra_Operator>(op), // Safe cast due to return type!
    opTrans, applyAs, adjointSupport, range, domain
    );
  thyraEpetraOp->setObjectLabel(label);
  return thyraEpetraOp;
}
