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
void DefaultBlockedLinearOp<Scalar>::beginBlockFill(
  const int numRowBlocks, const int numColBlocks
  )
{
  assertBlockFillIsActive(false);
  resetStorage(numRowBlocks,numColBlocks);
  rangeBlocks_.resize(numRowBlocks);
  domainBlocks_.resize(numColBlocks);
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
  setBlockSpaces(i,j,*block);
  Ops_[numRowBlocks_*j+i] = block;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::setBlock(
  const int i, const int j
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &block
  )
{
  setBlockSpaces(i,j,*block);
  Ops_[numRowBlocks_*j+i] = block;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::endBlockFill()
{
  using Teuchos::rcp;
  using Teuchos::arrayArg;
  assertBlockFillIsActive(true);
  if(!productRange_.get()) {
#ifdef TEUCHOS_DEBUG
    for(int i = 0; i < numRowBlocks_; ++i) {
      TEST_FOR_EXCEPTION(
        !rangeBlocks_[i].get(), std::logic_error
        ,"DefaultBlockedLinearOp<Scalar>::endBlockFill():"
        " Error, the block matrix for the i="<<i<<" block row space is missing"
        " and we can not complete the block fill!"
        );
    }
    for(int j = 0; j < numRowBlocks_; ++j) {
      TEST_FOR_EXCEPTION(
        !domainBlocks_[j].get(), std::logic_error
        ,"DefaultBlockedLinearOp<Scalar>::endBlockFill():"
        " Error, the block matrix for the j="<<j<<" block column space is missing"
        " and we can not complete the block fill!"
        );
    }
#endif
  }
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
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  TEST_FOR_EXCEPT(true);
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
  if(productRange_.get()) {
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
    !( 0 <= i && i < numRowBlocks_ ), std::logic_error
    ,"Error, i="<<i<<" does not fall in the range [0,"<<numRowBlocks_-1<<"]!"
    );
  TEST_FOR_EXCEPTION(
    !( 0 <= j && j < numColBlocks_ ), std::logic_error
    ,"Error, j="<<j<<" does not fall in the range [0,"<<numColBlocks_-1<<"]!"
    );
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
      ,*rangeBlock,("productRange->getBlock("+toString(i)+")")
      ,*block.range(),("block["+toString(i)+","+toString(j)+"].range()")
      );
  }
  if(domainBlock.get()) {
    THYRA_ASSERT_VEC_SPACES_NAMES(
      "DefaultBlockedLinearOp<Scalar>::setBlockSpaces(i,j,block)"
      ,*domainBlock,("productDomain->getBlock("+toString(j)+")")
      ,*block.domain(),("block["+toString(i)+","+toString(j)+"].domain()")
      );
  }
#endif
  if(!productRange_.get()) {
    if(!rangeBlocks_[i].get())
      rangeBlocks_[i] = block.range();
    if(!domainBlocks_[j].get()) {
      domainBlocks_[j] = block.domain();
    }
  }
}

} // namespace Thyra

#endif	// THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP
