// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DEF_HPP
#define THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DEF_HPP


#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve_decl.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// public


// Constructors/Initializers/Accessors


template<class Scalar>
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::DefaultBlockedTriangularLinearOpWithSolve()
  : blockFillIsActive_(false), numDiagBlocks_(0),
    hasLowerBlocks_(false), hasUpperBlocks_(false)
{}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setNonconstBlocks(
  const RCP<PhysicallyBlockedLinearOpBase<Scalar> > &blocks
  )
{
  assertAndSetBlockStructure(*blocks);
  blocks_.initialize(blocks);
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setBlocks(
  const RCP<const PhysicallyBlockedLinearOpBase<Scalar> > &blocks
  )
{
  assertAndSetBlockStructure(*blocks);
  blocks_.initialize(blocks);
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockTriangular(
  const PhysicallyBlockedLinearOpBase<Scalar> &blocks
  )
{
  bool hasLower, hasUpper;
  getBlockStructure(blocks, hasLower, hasUpper);
}


template<class Scalar>
RCP<PhysicallyBlockedLinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getNonconstBlocks()
{
  return blocks_.getNonconstObj();
}


template<class Scalar>
RCP<const PhysicallyBlockedLinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getBlocks()
{
  return blocks_.getConstObj();
}


// Overridden from PhysicallyBlockedLinearOpWithSolveBase


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::acceptsLOWSBlock(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  return i==j; // Only accept LOWS blocks along the diagonal!
}

template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setNonconstLOWSBlock(
  const int i, const int j,
  const Teuchos::RCP<LinearOpWithSolveBase<Scalar> > &block
  )
{
  setLOWSBlockImpl(i,j,block);
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setLOWSBlock(
  const int i, const int j,
  const Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > &block
  )
{
  setLOWSBlockImpl(i,j,block);
}


// Overridden from PhysicallyBlockedLinearOpBase


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::beginBlockFill()
{
 assertBlockFillIsActive(false);
 TEUCHOS_TEST_FOR_EXCEPT("ToDo: Have not implemented flexible block fill yet!");
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::beginBlockFill(
  const int numRowBlocks, const int numColBlocks
  )
{
  using Teuchos::null;
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT_EQUALITY(numRowBlocks, numColBlocks);
#else
  (void)numColBlocks;
#endif
  assertBlockFillIsActive(false);
  numDiagBlocks_ = numRowBlocks;
  diagonalBlocks_.resize(numDiagBlocks_);
  productRange_ = null;
  productDomain_ = null;
  removeOffDiagonalBlocks();
  blockFillIsActive_ = true;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::beginBlockFill(
  const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productRange_in,
  const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productDomain_in
  )
{
#ifdef THYRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( is_null(productRange_in) );
  TEUCHOS_TEST_FOR_EXCEPT( is_null(productDomain_in) );
  TEUCHOS_TEST_FOR_EXCEPT( productRange_in->numBlocks() != productDomain_in->numBlocks() );
#endif
  assertBlockFillIsActive(false);
  productRange_ = productRange_in;
  productDomain_ = productDomain_in;
  numDiagBlocks_ = productRange_in->numBlocks();
  diagonalBlocks_.resize(numDiagBlocks_);
  removeOffDiagonalBlocks();
  blockFillIsActive_ = true;
}


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::blockFillIsActive() const
{
  return blockFillIsActive_;
}


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::acceptsBlock(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  // The diagonal blocks are set with setLOWSBlock(), and the off-diagonal
  // blocks with setBlocks() after the block fill
  return false;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setNonconstBlock(
  const int /* i */, const int /* j */,
  const Teuchos::RCP<LinearOpBase<Scalar> > &/* block */
  )
{
  assertBlockFillIsActive(true);
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error, DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setNonconstBlock(): "
    "the diagonal blocks are set with setNonconstLOWSBlock(), and the off-diagonal "
    "blocks with setNonconstBlocks() after the block fill!"
    );
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setBlock(
  const int /* i */, const int /* j */,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &/* block */
  )
{
  assertBlockFillIsActive(true);
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error, DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setBlock(): "
    "the diagonal blocks are set with setLOWSBlock(), and the off-diagonal "
    "blocks with setBlocks() after the block fill!"
    );
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::endBlockFill()
{
  assertBlockFillIsActive(true);
  Array<RCP<const VectorSpaceBase<Scalar> > > rangeSpaces;
  Array<RCP<const VectorSpaceBase<Scalar> > > domainSpaces;
  for ( int k = 0; k < numDiagBlocks_; ++k ) {
    const RCP<const LinearOpWithSolveBase<Scalar> > lows_k = 
      diagonalBlocks_[k].getConstObj();
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(lows_k), std::logic_error,
      "Error, the block diagonal k="<<k<<" can not be null when ending block fill!"
      );
    if (is_null(productRange_)) {
      rangeSpaces.push_back(lows_k->range());
      domainSpaces.push_back(lows_k->domain());
    }
  }
  if (is_null(productRange_)) {
    productRange_ = productVectorSpace<Scalar>(rangeSpaces());
    productDomain_ = productVectorSpace<Scalar>(domainSpaces());
  }
  blockFillIsActive_ = false;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::uninitialize()
{
  assertBlockFillIsActive(false);
  productRange_ = Teuchos::null;
  productDomain_ = Teuchos::null;
  numDiagBlocks_ = 0;
  diagonalBlocks_.resize(0);
  removeOffDiagonalBlocks();
}


// Overridden from BlockedLinearOpWithSolveBase


template<class Scalar>
Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getNonconstLOWSBlock(
  const int i, const int j
  )
{
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  if (i!=j)
    return Teuchos::null;
  return diagonalBlocks_[i].getNonconstObj();
} 


template<class Scalar>
Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getLOWSBlock(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(false);
  assertBlockRowCol(i, j);
  if (i != j)
    return Teuchos::null;
  return diagonalBlocks_[i].getConstObj();
}


// Overridden from BlockedLinearOpBase


template<class Scalar>
Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::productRange() const
{
  return productRange_;
}


template<class Scalar>
Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::productDomain() const
{
  return productDomain_;
}


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::blockExists(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  if (i!=j)
    return nonnull(getOffDiagonalBlock(i,j));
  return !is_null(diagonalBlocks_[i].getConstObj());
} 


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::blockIsConst(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  return diagonalBlocks_[i].isConst();
} 


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getNonconstBlock(
  const int i, const int j
  )
{
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  if (i!=j) {
    if (is_null(blocks_.getConstObj()))
      return Teuchos::null;
    // Throws if the blocks were set with setBlocks()
    return blocks_.getNonconstObj()->getNonconstBlock(i,j);
  }
  return this->getNonconstLOWSBlock(i,j);
} 


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getBlock(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  if (i!=j)
    return getOffDiagonalBlock(i,j);
  return this->getLOWSBlock(i,j);
} 


// Overridden from LinearOpBase


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::range() const
{
  return this->productRange();
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::domain() const
{
  return this->productDomain();
}


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null;  // ToDo: Implement clone when needed!
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::description() const
{
  assertBlockFillIsActive(false);
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description() << "{"
    << "numDiagBlocks="<<numDiagBlocks_
    << "}";
  return oss.str();
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  assertBlockFillIsActive(false);
  Teuchos::Describable::describe(out, verbLevel);
  // ToDo: Fill in a better version of this!
}


// protected


// Overridden from LinearOpBase


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::opSupportedImpl(
  EOpTransp M_trans
  ) const
{
  using Thyra::opSupported;
  assertBlockFillIsActive(false);
  for ( int k = 0; k < numDiagBlocks_; ++k ) {
    if ( !opSupported(*diagonalBlocks_[k].getConstObj(),M_trans) )
      return false;
  }
  return offDiagonalBlocksSupport(M_trans);
  // ToDo: To be safe we really should do a collective reduction of all
  // clusters of processes.  However, for the typical use case, every block
  // will return the same info and we should be safe!
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X_in,
  const Ptr<MultiVectorBase<Scalar> > &Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  using Thyra::apply;

#ifdef THYRA_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::apply(...)",
    *this, M_trans, X_in, &*Y_inout
    );
#endif // THYRA_DEBUG  

  //
  // Y = alpha * op(M) * X + beta*Y
  //
  // =>
  //
  // Y[i] = beta*Y[i] + alpha*sum_j op(M)[i,j]*X[j], for i=0...numDiagBlocks-1
  //
  // where op(M)[i,j] is op(M[i,j]) for op = NOTRANS or CONJ, and op(M[j,i])
  // for op = TRANS or CONJTRANS (see getOpOffDiagonalBlock()).
  //

  typedef Teuchos::ScalarTraits<Scalar> ST;

  const ProductMultiVectorBase<Scalar>
    &X = dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in);
  ProductMultiVectorBase<Scalar>
    &Y = dyn_cast<ProductMultiVectorBase<Scalar> >(*Y_inout);
  
  for ( int i = 0; i < numDiagBlocks_; ++ i ) {
    const RCP<MultiVectorBase<Scalar> > Y_i = Y.getNonconstMultiVectorBlock(i);
    Thyra::apply( *diagonalBlocks_[i].getConstObj(), M_trans,
      *X.getMultiVectorBlock(i), Y_i.ptr(), alpha, beta
      );
    for ( int j = 0; j < numDiagBlocks_; ++j ) {
      const RCP<const LinearOpBase<Scalar> >
        opM_ij = getOpOffDiagonalBlock(M_trans,i,j);
      if (nonnull(opM_ij)) {
        Thyra::apply( *opM_ij, M_trans, *X.getMultiVectorBlock(j), Y_i.ptr(),
          alpha, ST::one() );
      }
    }
  }

}


// Overridden from LinearOpWithSolveBase


template<class Scalar>
bool
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::solveSupportsImpl(
  EOpTransp M_trans
  ) const
{
  assertBlockFillIsActive(false);
  for ( int k = 0; k < numDiagBlocks_; ++k ) {
    if ( !Thyra::solveSupports( *diagonalBlocks_[k].getConstObj(), M_trans ) )
      return false;
  }
  // The off-diagonal blocks are applied in the block substitution
  return offDiagonalBlocksSupport(M_trans);
  // ToDo: To be safe we really should do a collective reduction of all
  // clusters of processes.  However, for the typical use case, every block
  // will return the same info and we should be safe!
}


template<class Scalar>
bool
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  using Thyra::solveSupportsSolveMeasureType;
  assertBlockFillIsActive(false);
  for ( int k = 0; k < numDiagBlocks_; ++k ) {
    if (
      !solveSupportsSolveMeasureType(
        *diagonalBlocks_[k].getConstObj(),
        M_trans, solveMeasureType
        )
      )
    {
      return false;
    }
  }
  return true;
}


template<class Scalar>
SolveStatus<Scalar>
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::solveImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &B_in,
  const Ptr<MultiVectorBase<Scalar> > &X_inout,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  ) const
{

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  using Thyra::solve;

#ifdef THYRA_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::apply(...)",
    *this, M_trans, *X_inout, &B_in
    );
  TEUCHOS_TEST_FOR_EXCEPT(!this->solveSupportsImpl(M_trans));
  // TEUCHOS_TEST_FOR_EXCEPTION(
  //   nonnull(solveCriteria) && !solveCriteria->solveMeasureType.useDefault(),
  //   std::logic_error,
  //   "Error, we can't handle any non-default solve criteria yet!"
  //   );
  // ToDo: If solve criteria is to be handled, then we will have to be very
  // carefull how it it interpreted in terms of the individual period solves!
#endif // THYRA_DEBUG  

  //
  // X = inv(op(M)) * B
  //
  // => block substitution.  If op(M) is block lower triangular (M block lower
  // triangular and op = NOTRANS or CONJ, or M block upper triangular and
  // op = TRANS or CONJTRANS):
  //
  // X[i] = inv(op(M)[i,i])*( B[i] - sum_{j<i} op(M)[i,j]*X[j] ),
  //   for i=0...numDiagBlocks-1
  //
  // and if op(M) is block upper triangular:
  //
  // X[i] = inv(op(M)[i,i])*( B[i] - sum_{j>i} op(M)[i,j]*X[j] ),
  //   for i=numDiagBlocks-1...0
  //
  // (op(M)[i,j] as in applyImpl()).  With only diagonal blocks, the blocks
  // are solved independently.  The solve criteria are passed to the solves
  // of the diagonal blocks.
  //

  typedef Teuchos::ScalarTraits<Scalar> ST;

  const ProductMultiVectorBase<Scalar>
    &B = dyn_cast<const ProductMultiVectorBase<Scalar> >(B_in);
  ProductMultiVectorBase<Scalar>
    &X = dyn_cast<ProductMultiVectorBase<Scalar> >(*X_inout);

  const bool opIsUpper =
    ( real_trans(M_trans) == NOTRANS ? hasUpperBlocks_ : hasLowerBlocks_ );

  bool converged = true;
  for ( int k = 0; k < numDiagBlocks_; ++ k ) {
    const int i = ( opIsUpper ? numDiagBlocks_-1-k : k );
    const RCP<const LinearOpWithSolveBase<Scalar> >
      Op_k = diagonalBlocks_[i].getConstObj();
    Op_k->setOStream(this->getOStream());
    Op_k->setVerbLevel(this->getVerbLevel());
    // B[i] minus the off-diagonal blocks times the solved blocks of X
    RCP<const MultiVectorBase<Scalar> > B_i = B.getMultiVectorBlock(i);
    RCP<MultiVectorBase<Scalar> > R_i;
    const int j_begin = ( opIsUpper ? i+1 : 0 );
    const int j_end = ( opIsUpper ? numDiagBlocks_ : i );
    for ( int j = j_begin; j < j_end; ++j ) {
      const RCP<const LinearOpBase<Scalar> >
        opM_ij = getOpOffDiagonalBlock(M_trans,i,j);
      if (nonnull(opM_ij)) {
        if (is_null(R_i))
          R_i = B_i->clone_mv();
        Thyra::apply( *opM_ij, M_trans, *X.getMultiVectorBlock(j), R_i.ptr(),
          -ST::one(), ST::one() );
      }
    }
    if (nonnull(R_i))
      B_i = R_i;
    SolveStatus<Scalar> status =
      Thyra::solve( *Op_k, M_trans, *B_i,
                    X.getNonconstMultiVectorBlock(i).ptr(), solveCriteria );
    if (status.solveStatus != SOLVE_STATUS_CONVERGED)
      converged = false;
  }

   SolveStatus<Scalar> solveStatus;
   solveStatus.solveStatus =
     converged ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNCONVERGED;

  return solveStatus;

}



// private


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockFillIsActive(
  bool blockFillIsActive_in
  ) const
{
#ifdef THYRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!(blockFillIsActive_==blockFillIsActive_in));
#else
  (void)blockFillIsActive_in;
#endif
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockRowCol(
  const int i, const int j
  ) const
{
#ifdef THYRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= i && i < numDiagBlocks_ ), std::logic_error,
    "Error, i="<<i<<" does not fall in the range [0,"<<numDiagBlocks_-1<<"]!"
      );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= j && j < numDiagBlocks_ ), std::logic_error,
    "Error, j="<<j<<" does not fall in the range [0,"<<numDiagBlocks_-1<<"]!"
      );
#else
  (void)i;
  (void)j;
#endif
}


template<class Scalar>
template<class LinearOpWithSolveType>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setLOWSBlockImpl(
  const int i, const int j,
  const Teuchos::RCP<LinearOpWithSolveType> &block
  )
{
  assertBlockFillIsActive(true);
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT_INEQUALITY( i, >=, 0 );
  TEUCHOS_ASSERT_INEQUALITY( j, >=, 0 );
  TEUCHOS_ASSERT_INEQUALITY( i, <, numDiagBlocks_ );
  TEUCHOS_ASSERT_INEQUALITY( j, <, numDiagBlocks_ );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !this->acceptsLOWSBlock(i,j), std::logic_error,
    "Error, this DefaultBlockedTriangularLinearOpWithSolve does not accept\n"
    "LOWSB objects for the block i="<<i<<", j="<<j<<"!"
    );
#else
  (void)j;
#endif
  diagonalBlocks_[i] = block;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(
  const PhysicallyBlockedLinearOpBase<Scalar>& blocks
  )
{
  // Checked in all builds: apply and solve access the blocks (i,j) of blocks
  // for i,j < numDiagBlocks_
  TEUCHOS_TEST_FOR_EXCEPTION(
    blocks.productRange()->numBlocks() != numDiagBlocks_, std::logic_error,
    "Error, DefaultBlockedTriangularLinearOpWithSolve<Scalar>::"
    "assertAndSetBlockStructure(blocks): the blocked operator has "
    << blocks.productRange()->numBlocks() << " block rows, but this operator has "
    << numDiagBlocks_ << " diagonal blocks!"
    );
#ifdef THYRA_DEBUG
  THYRA_ASSERT_VEC_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(blocks)",
    *blocks.range(), *this->range()
    );
  THYRA_ASSERT_VEC_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(blocks)",
    *blocks.domain(), *this->domain()
    );
#endif
  // Sets the flags only if blocks is block triangular, so that a rejected
  // operator leaves this object unchanged
  getBlockStructure(blocks, hasLowerBlocks_, hasUpperBlocks_);
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getBlockStructure(
  const PhysicallyBlockedLinearOpBase<Scalar> &blocks,
  bool &hasLower, bool &hasUpper
  )
{
  const int numBlocks = blocks.productRange()->numBlocks();
  TEUCHOS_TEST_FOR_EXCEPTION(
    blocks.productDomain()->numBlocks() != numBlocks, std::logic_error,
    "Error, DefaultBlockedTriangularLinearOpWithSolve<Scalar>: the blocked "
    "operator has " << numBlocks << " block rows and "
    << blocks.productDomain()->numBlocks() << " block columns!"
    );
  // The off-diagonal blocks must all be above or all be below the diagonal
  bool lower = false, upper = false;
  for ( int i = 0; i < numBlocks; ++i ) {
    for ( int j = 0; j < numBlocks; ++j ) {
      if ( i != j && nonnull(blocks.getBlock(i,j)) ) {
        if (i > j)
          lower = true;
        else
          upper = true;
      }
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION( lower && upper, std::logic_error,
    "Error, DefaultBlockedTriangularLinearOpWithSolve<Scalar>: the blocked "
    "operator has nonzero blocks both below and above the diagonal, it is not "
    "block triangular!"
    );
  hasLower = lower;
  hasUpper = upper;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::removeOffDiagonalBlocks()
{
  blocks_.uninitialize();
  hasLowerBlocks_ = false;
  hasUpperBlocks_ = false;
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getOffDiagonalBlock(
  const int i, const int j
  ) const
{
  const RCP<const PhysicallyBlockedLinearOpBase<Scalar> >
    blocks = blocks_.getConstObj();
  if ( i == j || is_null(blocks) )
    return Teuchos::null;
  return blocks->getBlock(i,j);
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getOpOffDiagonalBlock(
  const EOpTransp M_trans, const int i, const int j
  ) const
{
  return ( real_trans(M_trans) == NOTRANS
    ? getOffDiagonalBlock(i,j) : getOffDiagonalBlock(j,i) );
}


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::offDiagonalBlocksSupport(
  const EOpTransp M_trans
  ) const
{
  for ( int i = 0; i < numDiagBlocks_; ++i ) {
    for ( int j = 0; j < numDiagBlocks_; ++j ) {
      const RCP<const LinearOpBase<Scalar> > M_ij = getOffDiagonalBlock(i,j);
      if ( nonnull(M_ij) && !Thyra::opSupported(*M_ij,M_trans) )
        return false;
    }
  }
  return true;
}


} // namespace Thyra


#endif	// THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DEF_HPP
