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
  : blockFillIsActive_(false), numDiagBlocks_(0)
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
 TEST_FOR_EXCEPT("ToDo: Have not implemented flexible block fill yet!");
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::beginBlockFill(
  const int numRowBlocks, const int numColBlocks
  )
{
  using Teuchos::null;
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT_EQUALITY(numRowBlocks, numColBlocks);
#endif
  assertBlockFillIsActive(false);
  numDiagBlocks_ = numRowBlocks;
  diagonalBlocks_.resize(numDiagBlocks_);
  productRange_ = null;
  productDomain_ = null;
  blockFillIsActive_ = true;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::beginBlockFill(
  const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productRange_in,
  const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productDomain_in
  )
{
#ifdef THYRA_DEBUG
  TEST_FOR_EXCEPT( is_null(productRange_in) );
  TEST_FOR_EXCEPT( is_null(productDomain_in) );
  TEST_FOR_EXCEPT( productRange_in->numBlocks() != productDomain_in->numBlocks() );
#endif
  assertBlockFillIsActive(false);
  productRange_ = productRange_in;
  productDomain_ = productDomain_in;
  numDiagBlocks_ = productRange_in->numBlocks();
  diagonalBlocks_.resize(numDiagBlocks_);
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
  return false; // ToDo: Change this once we accept off-diagonal blocks
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setNonconstBlock(
  const int i, const int j,
  const Teuchos::RCP<LinearOpBase<Scalar> > &block
  )
{
  assertBlockFillIsActive(true);
  TEST_FOR_EXCEPT("Error, we don't support off-diagonal LOB objects yet!");
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::setBlock(
  const int i, const int j,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &block
  )
{
  assertBlockFillIsActive(true);
  TEST_FOR_EXCEPT("Error, we don't support off-diagonal LOB objects yet!");
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
    TEST_FOR_EXCEPTION(is_null(lows_k), std::logic_error,
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
    return false; // ToDo: Update this when off-diagonals are supported!
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
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  if (i!=j)
    return Teuchos::null; // ToDo: Update when off-diagonals are supported!
  return this->getNonconstLOWSBlock(i,j);
} 


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::getBlock(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  if (i!=j)
    return Teuchos::null; // ToDo: Update when off-diagonals are supported!
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
  typedef Teuchos::ScalarTraits<Scalar>  ST;
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
  return true;
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
  typedef Teuchos::ScalarTraits<Scalar> ST;
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
  // Y[i] = beta+Y[i] + alpha*op(Op)[i]*X[i], for i=0...numDiagBlocks-1
  //
  // ToDo: Update to handle off diagonal blocks when needed!
  //

  const ProductMultiVectorBase<Scalar>
    &X = dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in);
  ProductMultiVectorBase<Scalar>
    &Y = dyn_cast<ProductMultiVectorBase<Scalar> >(*Y_inout);
  
  for ( int i = 0; i < numDiagBlocks_; ++ i ) {
    Thyra::apply( *diagonalBlocks_[i].getConstObj(), M_trans,
      *X.getMultiVectorBlock(i), Y.getNonconstMultiVectorBlock(i).ptr(),
      alpha, beta
      );
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
  return true;
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
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::solve;

#ifdef THYRA_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::apply(...)",
    *this, M_trans, *X_inout, &B_in
    );
  TEST_FOR_EXCEPT(!this->solveSupportsImpl(M_trans));
  TEST_FOR_EXCEPTION(
    nonnull(solveCriteria) && !solveCriteria->solveMeasureType.useDefault(),
    std::logic_error,
    "Error, we can't handle any non-default solve criteria yet!"
    );
  // ToDo: If solve criteria is to be handled, then we will have to be very
  // carefull how it it interpreted in terms of the individual period solves!
#endif // THYRA_DEBUG  

  //
  // Y = alpha * inv(op(M)) * X + beta*Y
  //
  // =>
  //
  // X[i] = inv(op(Op[i]))*B[i], for i=0...numDiagBlocks-1
  //
  // ToDo: Update to handle off diagonal blocks when needed!
  //

  const ProductMultiVectorBase<Scalar>
    &B = dyn_cast<const ProductMultiVectorBase<Scalar> >(B_in);
  ProductMultiVectorBase<Scalar>
    &X = dyn_cast<ProductMultiVectorBase<Scalar> >(*X_inout);
  
  for ( int i = 0; i < numDiagBlocks_; ++ i ) {
    const RCP<const LinearOpWithSolveBase<Scalar> >
      Op_k = diagonalBlocks_[i].getConstObj();
    Op_k->setOStream(this->getOStream());
    Op_k->setVerbLevel(this->getVerbLevel());
    Thyra::solve( *Op_k, M_trans, *B.getMultiVectorBlock(i),
      X.getNonconstMultiVectorBlock(i).ptr() );
    // ToDo: Pass in solve criteria when needed!
  }

  return SolveStatus<Scalar>();

}



// private


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockFillIsActive(
  bool blockFillIsActive_in
  ) const
{
#ifdef THYRA_DEBUG
  TEST_FOR_EXCEPT(!(blockFillIsActive_==blockFillIsActive_in));
#endif
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockRowCol(
  const int i, const int j
  ) const
{
#ifdef THYRA_DEBUG
  TEST_FOR_EXCEPTION(
    !( 0 <= i && i < numDiagBlocks_ ), std::logic_error,
    "Error, i="<<i<<" does not fall in the range [0,"<<numDiagBlocks_-1<<"]!"
      );
  TEST_FOR_EXCEPTION(
    !( 0 <= j && j < numDiagBlocks_ ), std::logic_error,
    "Error, j="<<j<<" does not fall in the range [0,"<<numDiagBlocks_-1<<"]!"
      );
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
  TEST_FOR_EXCEPTION(
    !this->acceptsLOWSBlock(i,j), std::logic_error,
    "Error, this DefaultBlockedTriangularLinearOpWithSolve does not accept\n"
    "LOWSB objects for the block i="<<i<<", j="<<j<<"!"
    );
#endif
  diagonalBlocks_[i] = block;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(
  const PhysicallyBlockedLinearOpBase<Scalar>& blocks
  )
{
#ifdef THYRA_DEBUG
  THYRA_ASSERT_VEC_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(blocks)",
    *blocks.range(), *this->range()
    );
  THYRA_ASSERT_VEC_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(blocks)",
    *blocks.domain(), *this->domain()
    );
  // ToDo: Make sure that all of the blocks are above or below the diagonal
  // but not both!
#endif
  // ToDo: Set if this is an upper or lower triangular block operator.
}


} // namespace Thyra


#endif	// THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DEF_HPP
