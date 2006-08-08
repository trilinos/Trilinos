// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP
#define THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP

#include "Thyra_DefaultBlockedLinearOpDecl.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"

namespace Thyra {

// Constructors

template<class Scalar>
DefaultBlockedLinearOp<Scalar>::DefaultBlockedLinearOp()
  :blockFillIsActive_(false)
{}

// Overridden from PhysicallyBlockedLinearOpBase

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill()
{
  assertBlockFillIsActive(false);
  resetStorage(0,0);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill(
  const int numRowBlocks, const int numColBlocks
  )
{
  assertBlockFillIsActive(false);
  resetStorage(numRowBlocks,numColBlocks);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill(
  const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >  &productRange
  ,const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > &productDomain
  )
{
  assertBlockFillIsActive(false);
  productRange_ = productRange.assert_not_null();
  productDomain_ = productDomain.assert_not_null();
  // Note: the above spaces must be set before calling the next function!
  resetStorage(productRange_->numBlocks(),productDomain_->numBlocks());
}

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockFillIsActive() const
{
  return blockFillIsActive_;
}

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::acceptsBlock(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  return true;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::setNonconstBlock(
  const int i, const int j
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &block
  )
{
  setBlockImpl(i,j,block);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::setBlock(
  const int i, const int j
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &block
  )
{
  setBlockImpl(i,j,block);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::endBlockFill()
{
  using Teuchos::rcp;
  using Teuchos::arrayArg;
  assertBlockFillIsActive(true);
  // Create the product range and domain spaces if these are not already set.
  if(!productRange_.get()) {
#ifdef TEUCHOS_DEBUG
    for(int i = 0; i < numRowBlocks_; ++i) {
      TEST_FOR_EXCEPTION(
        !rangeBlocks_[i].get(), std::logic_error
        ,"DefaultBlockedLinearOp<Scalar>::endBlockFill():"
        " Error, no linear operator block for the i="<<i<<" block row was added"
        " and we can not complete the block fill!"
        );
    }
    for(int j = 0; j < numColBlocks_; ++j) {
      TEST_FOR_EXCEPTION(
        !domainBlocks_[j].get(), std::logic_error
        ,"DefaultBlockedLinearOp<Scalar>::endBlockFill():"
        " Error, no linear operator block for the j="<<j<<" block column was added"
        " and we can not complete the block fill!"
        );
    }
#endif
    productRange_
      = rcp(
        new DefaultProductVectorSpace<Scalar>(
          numRowBlocks_,&rangeBlocks_[0]
          )
        );
    productDomain_
      = rcp(
        new DefaultProductVectorSpace<Scalar>(
          numColBlocks_,&domainBlocks_[0]
          )
        );
  }
  numRowBlocks_ = productRange_->numBlocks();
  numColBlocks_ = productDomain_->numBlocks();
  rangeBlocks_.resize(0);
  domainBlocks_.resize(0);
  // Insert the block LOB objects if doing a flexible fill.
  if(Ops_stack_.size()) {
    Ops_.resize(numRowBlocks_*numColBlocks_);
    for( int k = 0; k < static_cast<int>(Ops_stack_.size()); ++k ) {
      const BlockEntry<Scalar> &block_i_j = Ops_stack_[k];
      Ops_[numRowBlocks_*block_i_j.j+block_i_j.i] = block_i_j.block;
    }
    Ops_stack_.resize(0);
  }
  blockFillIsActive_ = false;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::uninitialize()
{
  productRange_ = Teuchos::null;
  productDomain_ = Teuchos::null;
  numRowBlocks_ = 0;
  numColBlocks_ = 0;
  Ops_.resize(0);
  Ops_stack_.resize(0);
  rangeBlocks_.resize(0);
  domainBlocks_.resize(0);
  blockFillIsActive_ = false;
}

// Overridden from BlockedLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::productRange() const
{
  return productRange_;
}

template<class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::productDomain() const
{
  return productDomain_;
}

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockExists(
  const int i, const int j
  ) const
{
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  return true;
} 

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockIsConst(
  const int i, const int j
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!blockExists(i,j));
#endif
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  return Ops_[numRowBlocks_*j+i].isConst();
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::getNonconstBlock(const int i, const int j)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!blockExists(i,j));
#endif
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  return Ops_[numRowBlocks_*j+i].getNonconstObj();
} 

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::getBlock(const int i, const int j) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!blockExists(i,j));
#endif
  assertBlockFillIsActive(false);
  assertBlockRowCol(i,j);
  return Ops_[numRowBlocks_*j+i].getConstObj();
} 

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::range() const
{
  return productRange_;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::domain() const
{
  return productDomain_;
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // ToDo: Implement this when needed!
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultBlockedLinearOp<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  assertBlockFillIsActive(false);
  std::ostringstream oss;
  oss
    << "Thyra::DefaultBlockedLinearOp<" << ST::name() << ">"
    << "{"
    << "numRowBlocks="<<numRowBlocks_
    << ",numColBlocks="<<numColBlocks_
    << "}";
  return oss.str();
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  assertBlockFillIsActive(false);
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << "Thyra::DefaultBlockedLinearOp<" << ST::name() << ">{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim=" << this->domain()->dim()
        << ",numRowBlocks=" << numRowBlocks_
        << ",numColBlocks=" << numColBlocks_
        << "}\n";
      OSTab tab(out);
      *out
        <<  "Constituent LinearOpBase objects for M = [ Op[0,0] ... ; ... ; ... Op[numRowBlocks-1,numColBlocks-1]:\n";
      tab.incrTab();
      for( int i = 0; i < numRowBlocks_; ++i ) {
        for( int j = 0; j < numColBlocks_; ++j ) {
          *out << "Op["<<i<<","<<j<<"] = ";
          Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
            block_i_j = getBlock(i,j);
          if(block_i_j.get())
            *out << "\n" << Teuchos::describe(*getBlock(i,j),verbLevel);
          else
            *out << "NULL\n";
        }
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::opSupported(
  ETransp M_trans
  ) const
{
  bool opSupported = true;
  for( int i = 0; i < numRowBlocks_; ++i ) {
    for( int j = 0; j < numColBlocks_; ++j ) {
      Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
        block_i_j = getBlock(i,j);
      if( block_i_j.get() && !Thyra::opSupported(*block_i_j,M_trans) )
        opSupported = false;
    }
  }
  return opSupported;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X_in
  ,MultiVectorBase<Scalar>          *Y_inout
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  using Teuchos::RefCountPtr;
  using Teuchos::dyn_cast;
  typedef Teuchos::ScalarTraits<Scalar>                ST;
  typedef RefCountPtr<MultiVectorBase<Scalar> >        MultiVectorPtr;
  typedef RefCountPtr<const MultiVectorBase<Scalar> >  ConstMultiVectorPtr;
  typedef RefCountPtr<const LinearOpBase<Scalar> >     ConstLinearOpPtr;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultBlockedLinearOp<Scalar>::apply(...)",*this,M_trans,X_in,Y_inout
    );
#endif // TEUCHOS_DEBUG  
  const bool
    struct_transp = (real_trans(M_trans)!=NOTRANS); // Structural transpose?
  const int
    opNumRowBlocks = ( !struct_transp ? numRowBlocks_ : numColBlocks_ ), 
    opNumColBlocks = ( !struct_transp ? numColBlocks_ : numRowBlocks_ ); 
  //
  // Y = alpha * op(M) * X + beta*Y
  //
  // =>
  //
  // Y[i] = beta+Y[i] + sum(alpha*op(Op)[i,j]*X[j],j=0...opNumColBlocks-1)
  //
  //   , for i=0...opNumRowBlocks-1
  //
  const ProductMultiVectorBase<Scalar>
    &X = dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in);
  ProductMultiVectorBase<Scalar>
    &Y = dyn_cast<ProductMultiVectorBase<Scalar> >(*Y_inout);
  for( int i = 0; i < opNumRowBlocks; ++i ) {
    MultiVectorPtr Y_i = Y.getNonconstMultiVectorBlock(i);
    for( int j = 0; j < opNumColBlocks; ++j ) {
      ConstLinearOpPtr     Op_i_j = ( !struct_transp ? getBlock(i,j) : getBlock(j,i) );
      ConstMultiVectorPtr  X_j    = X.getMultiVectorBlock(j);
      if(j==0) {
        if(Op_i_j.get())  Thyra::apply(*Op_i_j,M_trans,*X_j,&*Y_i,alpha,beta);
        else              scale(beta,&*Y_i);
      }
      else {
        if(Op_i_j.get())  Thyra::apply(*Op_i_j,M_trans,*X_j,&*Y_i,alpha,ST::one());
      }
    }
  }
}

// private

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::resetStorage(
  const int numRowBlocks, const int numColBlocks
  )
{
  uninitialize();
  numRowBlocks_ = numRowBlocks;
  numColBlocks_ = numColBlocks;
  Ops_.resize(numRowBlocks_*numColBlocks_);
  if(!productRange_.get()) {
    rangeBlocks_.resize(numRowBlocks);
    domainBlocks_.resize(numColBlocks);
  }
  blockFillIsActive_ = true;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::assertBlockFillIsActive(
  bool blockFillIsActive
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!(blockFillIsActive_==blockFillIsActive));
#endif
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::assertBlockRowCol(
  const int i, const int j
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !( 0 <= i ), std::logic_error
    ,"Error, i="<<i<<" is invalid!"
    );
  TEST_FOR_EXCEPTION(
    !( 0 <= j ), std::logic_error
    ,"Error, j="<<j<<" is invalid!"
    );
  // Only validate upper range if the number of row and column blocks is
  // fixed!
  if(Ops_.size()) {
    TEST_FOR_EXCEPTION(
      !( 0 <= i && i < numRowBlocks_ ), std::logic_error
      ,"Error, i="<<i<<" does not fall in the range [0,"<<numRowBlocks_-1<<"]!"
      );
    TEST_FOR_EXCEPTION(
      !( 0 <= j && j < numColBlocks_ ), std::logic_error
      ,"Error, j="<<j<<" does not fall in the range [0,"<<numColBlocks_-1<<"]!"
      );
  }
#endif
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::setBlockSpaces(
  const int i, const int j, const LinearOpBase<Scalar> &block
  )
{
  typedef std::string s;
  using Teuchos::toString;
  assertBlockFillIsActive(true);
  assertBlockRowCol(i,j);
  // Validate that if the vector space block is already set that it is
  // compatible with the block that is being set.
  if( i < numRowBlocks_ && j < numColBlocks_ ) {
#ifdef TEUCHOS_DEBUG
    Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
      rangeBlock = (
        productRange_.get()
        ? productRange_->getBlock(i)
        : rangeBlocks_[i]
        ),
      domainBlock = (
        productDomain_.get()
        ? productDomain_->getBlock(j)
        : domainBlocks_[j]
        );
    if(rangeBlock.get()) {
      THYRA_ASSERT_VEC_SPACES_NAMES(
        "DefaultBlockedLinearOp<Scalar>::setBlockSpaces(i,j,block)"
        ,*rangeBlock,("(*productRange->getBlock("+toString(i)+"))")
        ,*block.range(),("(*block["+toString(i)+","+toString(j)+"].range())")
        );
    }
    if(domainBlock.get()) {
      THYRA_ASSERT_VEC_SPACES_NAMES(
        "DefaultBlockedLinearOp<Scalar>::setBlockSpaces(i,j,block)"
        ,*domainBlock,("(*productDomain->getBlock("+toString(j)+"))")
        ,*block.domain(),("(*block["+toString(i)+","+toString(j)+"].domain())")
        );
    }
#endif // TEUCHOS_DEBUG
  }
  // Add spaces missing range and domain space blocks if we are doing a
  // flexible fill (otherwise these loops will not be executed)
  for( int k = numRowBlocks_; k <= i; ++k )
    rangeBlocks_.push_back(Teuchos::null);
  for( int k = numColBlocks_; k <= j; ++k )
    domainBlocks_.push_back(Teuchos::null);
  // Set the incoming range and domain blocks if not already set
  if(!productRange_.get()) {
    if(!rangeBlocks_[i].get())
      rangeBlocks_[i] = block.range();
    if(!domainBlocks_[j].get()) {
      domainBlocks_[j] = block.domain();
    }
  }
  // Update the current number of row and columns blocks if doing a flexible
  // fill.
  if(!Ops_.size()) {
    numRowBlocks_ = rangeBlocks_.size();
    numColBlocks_ = domainBlocks_.size();
  }
}

template<class Scalar>
template<class LinearOpType>
void DefaultBlockedLinearOp<Scalar>::setBlockImpl(
  const int i, const int j
  ,const Teuchos::RefCountPtr<LinearOpType> &block
  )
{
  setBlockSpaces(i,j,*block);
  if(Ops_.size()) {
    // We are doing a fill with a fixed number of row and column blocks so we
    // can just set this.
    Ops_[numRowBlocks_*j+i] = block;
  }
  else {
    // We are doing a flexible fill so add the block to the stack of blocks or
    // replace a block that already exists.
    bool foundBlock = false;
    for( unsigned int k = 0; k < Ops_stack_.size(); ++k ) {
      BlockEntry<Scalar> &block_i_j = Ops_stack_[k];
      if( block_i_j.i == i && block_i_j.j == j ) {
        block_i_j.block = block;
        foundBlock = true;
        break;
      }
    }
    if(!foundBlock)
      Ops_stack_.push_back(BlockEntry<Scalar>(i,j,block));
  }
}

} // namespace Thyra

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::block2x2(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &A00
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &A01
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &A10
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &A11
  )
{
  Teuchos::RefCountPtr<PhysicallyBlockedLinearOpBase<Scalar> >
    M = Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
  M->beginBlockFill(2,2);
  if(A00.get()) M->setBlock(0,0,A00);
  if(A01.get()) M->setBlock(0,1,A01);
  if(A10.get()) M->setBlock(1,0,A10);
  if(A11.get()) M->setBlock(1,1,A11);
  M->endBlockFill();
  return M;
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::block2x1(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &A00
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &A10
  )
{
  Teuchos::RefCountPtr<PhysicallyBlockedLinearOpBase<Scalar> >
    M = Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
  M->beginBlockFill(2,1);
  if(A00.get()) M->setBlock(0,0,A00);
  if(A10.get()) M->setBlock(1,0,A10);
  M->endBlockFill();
  return M;
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::block1x2(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &A00
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &A01
  )
{
  Teuchos::RefCountPtr<PhysicallyBlockedLinearOpBase<Scalar> >
    M = Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
  M->beginBlockFill(1,2);
  if(A00.get()) M->setBlock(0,0,A00);
  if(A01.get()) M->setBlock(0,1,A01);
  M->endBlockFill();
  return M;
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock2x2(
  const Teuchos::RefCountPtr<LinearOpBase<Scalar> >    &A00
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> >   &A01
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> >   &A10
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> >   &A11
  )
{
  Teuchos::RefCountPtr<PhysicallyBlockedLinearOpBase<Scalar> >
    M = Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
  M->beginBlockFill(2,2);
  if(A00.get()) M->setNonconstBlock(0,0,A00);
  if(A01.get()) M->setNonconstBlock(0,1,A01);
  if(A10.get()) M->setNonconstBlock(1,0,A10);
  if(A11.get()) M->setNonconstBlock(1,1,A11);
  M->endBlockFill();
  return M;
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock2x1(
  const Teuchos::RefCountPtr<LinearOpBase<Scalar> >    &A00
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> >   &A10
  )
{
  Teuchos::RefCountPtr<PhysicallyBlockedLinearOpBase<Scalar> >
    M = Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
  M->beginBlockFill(2,1);
  if(A00.get()) M->setNonconstBlock(0,0,A00);
  if(A10.get()) M->setNonconstBlock(1,0,A10);
  M->endBlockFill();
  return M;
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock1x2(
  const Teuchos::RefCountPtr<LinearOpBase<Scalar> >    &A00
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> >   &A01
  )
{
  Teuchos::RefCountPtr<PhysicallyBlockedLinearOpBase<Scalar> >
    M = Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
  M->beginBlockFill(1,2);
  if(A00.get()) M->setNonconstBlock(0,0,A00);
  if(A01.get()) M->setNonconstBlock(0,1,A01);
  M->endBlockFill();
  return M;
}

#endif	// THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP
