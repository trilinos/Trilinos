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

#ifndef THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DEF_HPP
#define THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DEF_HPP


#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve_decl.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
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
  assertBlockFillIsActive(false);
  TEST_FOR_EXCEPT("ToDo: Have not implemented block fill with just numBlocks yet!!");
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::beginBlockFill(
  const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >  &productRange,
  const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productDomain
  )
{
  assertBlockFillIsActive(false);
  TEST_FOR_EXCEPT( is_null(productRange) );
  TEST_FOR_EXCEPT( is_null(productDomain) );
  TEST_FOR_EXCEPT( productRange->numBlocks() != productDomain->numBlocks() );
  productRange_ = productRange.assert_not_null();
  productDomain_ = productDomain.assert_not_null();
  numDiagBlocks_ = productRange->numBlocks();
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
  for ( int k = 0; k < numDiagBlocks_; ++k ) {
    TEST_FOR_EXCEPTION(
      is_null(diagonalBlocks_[k].getConstObj()), std::logic_error,
      "Error, the block diagonal k="<<k<<" can not be null when ending block fill!"
      );
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
  assertBlockRowCol(i,j);
  if (i!=j)
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
  Teuchos::Describable::describe(out,verbLevel);
  // ToDo: Fill in a better version of this!
}


// protected


// Overridden from SingleScalarLinearOpWithSolveBase


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::solveSupportsTrans(
  EOpTransp M_trans
  ) const
{
  assertBlockFillIsActive(false);
  for ( int k = 0; k < numDiagBlocks_; ++k ) {
    if ( !solveSupports( *diagonalBlocks_[k].getConstObj(), M_trans ) )
      return false;
  }
  return true;
  // ToDo: To be safe we really should do a collective reduction of all
  // clusters of processes.  However, for the typical use case, every block
  // will return the same info and we should be safe!
}


template<class Scalar>
bool
DefaultBlockedTriangularLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureType(
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
  // ToDo: To be safe we really should do a collective reduction of all
  // clusters of processes.  However, for the typical use case, every block
  // will return the same info and we should be safe!
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::solve(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &B_in,
  MultiVectorBase<Scalar> *X_inout,
  const int numBlocks,
  const BlockSolveCriteria<Scalar> blockSolveCriteria[],
  SolveStatus<Scalar> blockSolveStatus[]
  ) const
{

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::solve;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==X_inout);
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::apply(...)",*this,M_trans,*X_inout,&B_in
    );
  TEST_FOR_EXCEPT(!this->solveSupportsTrans(M_trans));
  TEST_FOR_EXCEPT( numBlocks > 1 );
  TEST_FOR_EXCEPTION(
    blockSolveCriteria && !blockSolveCriteria[0].solveCriteria.solveMeasureType.useDefault(),
    std::logic_error,
    "Error, we can't handle any non-default solve criteria yet!"
    );
  // ToDo: If solve criteria is to be handled, then we will have to be very
  // carefull how it it interpreted in terms of the individual period solves!
#endif // TEUCHOS_DEBUG  

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
      &*X.getNonconstMultiVectorBlock(i) );
    // ToDo: Pass in solve criteria when needed!
  }
  
  // ToDo: We really need to collect the SolveStatus objects across clusters
  // in order to really implement this interface correctly!  If a solve fails
  // on some cluster but not another, then different solve status information
  // will be returned.  This may result in different decisions being taken in
  // different clusters with bad consequences.  To really handle multiple
  // clusters correctly, we need a strategy object that knows how to do these
  // reductions correctly.  Actually, we could use Teuchos:Comm along with a
  // user-defined reduction object to do these reductions correctly!

}


// Overridden from SingleScalarLinearOpBase


template<class Scalar>
bool DefaultBlockedTriangularLinearOpWithSolve<Scalar>::opSupported(
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
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::apply(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X_in,
  MultiVectorBase<Scalar> *Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::apply;

#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultBlockedTriangularLinearOpWithSolve<Scalar>::apply(...)",*this,M_trans,X_in,Y_inout
    );
#endif // TEUCHOS_DEBUG  

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
      *X.getMultiVectorBlock(i),
      &*Y.getNonconstMultiVectorBlock(i)
      );
  }

}


// private


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockFillIsActive(
  bool blockFillIsActive
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!(blockFillIsActive_==blockFillIsActive));
#endif
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertBlockRowCol(
  const int i, const int j
  ) const
{
#ifdef TEUCHOS_DEBUG
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
  TEST_FOR_EXCEPTION(
    !this->acceptsLOWSBlock(i,j), std::logic_error,
    "Error, this DefaultBlockedTriangularLinearOpWithSolve does not accept\n"
    "LOWSB objects for the block i="<<i<<", j="<<j<<"!"
    );
  diagonalBlocks_[i] = block;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolve<Scalar>::assertAndSetBlockStructure(
  const PhysicallyBlockedLinearOpBase<Scalar>& blocks
  )
{
#ifdef TEUCHOS_DEBUG
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
