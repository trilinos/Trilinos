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

#ifndef THYRA_REAL_COMPLEX_FFT_LINEAR_OP_HPP
#define THYRA_REAL_COMPLEX_FFT_LINEAR_OP_HPP

#include "ComplexFFTLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DefaultSerialVectorSpaceConverter.hpp"
#include "serial_1D_FFT.hpp"

/** \brief Simple example concrete subclass for a serial real-to-complex FFT.
 *
 * This is a very bad but straightforward implementation of a real-to-complex
 * FFT operator that simply uses the implemention <tt>ComplexFFTLinearOp</tt>.
 *
 * \ingroup Thyra_Op_Vec_examples_fft_grp
 */
template<class RealScalar>
class RealComplexFFTLinearOp : virtual public Thyra::LinearOpWithSolveBase< std::complex<RealScalar>, RealScalar >
{
public:

  /** \brief . */
  typedef std::complex<RealScalar>                                                 RangeScalar;
  /** \brief . */
  typedef RealScalar                                                               DomainScalar;
  /** \brief .*/
  typedef typename Teuchos::PromotionTraits<RangeScalar,DomainScalar>::promote     Scalar;

  /** \brief . */
  RealComplexFFTLinearOp( const int N );

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<RangeScalar> > range() const;
  /** \brief . */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<DomainScalar> > domain() const;
  /** \brief . */
  bool applySupports( const Thyra::EConj conj ) const;
  /** \brief . */
  bool applyTransposeSupports( const Thyra::EConj conj ) const;
  /** \brief . */
  void apply(
    const Thyra::EConj                            conj
    ,const Thyra::MultiVectorBase<DomainScalar>   &X
    ,Thyra::MultiVectorBase<RangeScalar>          *Y
    ,const RangeScalar                            alpha
    ,const RangeScalar                           beta
    ) const;
  /** \brief . */
  void applyTranspose(
    const Thyra::EConj                            conj
    ,const Thyra::MultiVectorBase<RangeScalar>    &X
    ,Thyra::MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                           alpha
    ,const DomainScalar                           beta
    ) const;

  //@}

  /** @name Overridden from LinearOpWithSolveBase */
  //@{

  /** \brief . */
  bool solveSupportsConj(Thyra::EConj conj) const;
  /** \brief . */
  bool solveTransposeSupportsConj(Thyra::EConj conj) const;
  /** \brief . */
  bool solveSupportsSolveMeasureType(Thyra::EConj conj, const Thyra::SolveMeasureType solveMeasureType) const;
  /** \brief . */
  bool solveTransposeSupportsSolveMeasureType(Thyra::EConj conj, const Thyra::SolveMeasureType solveMeasureType) const;
  /** \brief . */
  void solve(
    const Thyra::EConj                           conj
    ,const Thyra::MultiVectorBase<RangeScalar>   &B
    ,Thyra::MultiVectorBase<DomainScalar>        *X
    ,const int                                   numBlocks
    ,const Thyra::BlockSolveCriteria<Scalar>     blockSolveCriteria[]
    ,Thyra::SolveStatus<Scalar>                  blockSolveStatus[]
    ) const;
  /** \brief . */
  void solveTranspose(
    const Thyra::EConj                           conj
    ,const Thyra::MultiVectorBase<DomainScalar>  &B
    ,Thyra::MultiVectorBase<RangeScalar>         *X
    ,const int                                   numBlocks
    ,const Thyra::BlockSolveCriteria<Scalar>     blockSolveCriteria[]
    ,Thyra::SolveStatus<Scalar>                  blockSolveStatus[]
    ) const;

  //@}

private:

  Thyra::DefaultSerialVectorSpaceConverter<DomainScalar,RangeScalar>      realToComplexConverter_;
  ComplexFFTLinearOp<RealScalar>                                          complexFFTOp_;
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase< RealScalar > >      domain_;
  
  /// V = beta*V + real(V_complex)
  static void updateReal(
    const Thyra::MultiVectorBase<RangeScalar>      &V_complex
    ,const DomainScalar                            &beta
    ,Thyra::MultiVectorBase<DomainScalar>          *V
    );

};

// /////////////////////////
// Definitions

template<class RealScalar>
RealComplexFFTLinearOp<RealScalar>::RealComplexFFTLinearOp( const int N )
  :complexFFTOp_(N), domain_(realToComplexConverter_.createVectorSpaceFrom(*complexFFTOp_.domain()))
{}

// Overridden from LinearOpBase

template<class RealScalar>
Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<typename RealComplexFFTLinearOp<RealScalar>::RangeScalar > >
RealComplexFFTLinearOp<RealScalar>::range() const
{
  return complexFFTOp_.range();
}

template<class RealScalar>
Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<typename RealComplexFFTLinearOp<RealScalar>::DomainScalar > >
RealComplexFFTLinearOp<RealScalar>::domain() const
{
  return domain_;
}

template<class RealScalar>
bool RealComplexFFTLinearOp<RealScalar>::applySupports( const Thyra::EConj conj ) const
{
  return (conj == Thyra::NONCONJ_ELE);
}

template<class RealScalar>
bool RealComplexFFTLinearOp<RealScalar>::applyTransposeSupports( const Thyra::EConj conj ) const
{
  return (conj == Thyra::CONJ_ELE);
}

template<class RealScalar>
void RealComplexFFTLinearOp<RealScalar>::apply(
  const Thyra::EConj                            conj
  ,const Thyra::MultiVectorBase<DomainScalar>   &X
  ,Thyra::MultiVectorBase<RangeScalar>          *Y
  ,const RangeScalar                            alpha
  ,const RangeScalar                            beta
  ) const
{
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<RangeScalar> >
    X_complex = Thyra::createMembers(complexFFTOp_.domain(),X.domain()->dim());
  realToComplexConverter_.convert( X, &*X_complex );
  Thyra::apply(complexFFTOp_,conj,*X_complex,Y,alpha,beta);
}

template<class RealScalar>
void RealComplexFFTLinearOp<RealScalar>::applyTranspose(
  const Thyra::EConj                            conj
  ,const Thyra::MultiVectorBase<RangeScalar>    &X
  ,Thyra::MultiVectorBase<DomainScalar>         *Y
  ,const DomainScalar                           alpha
  ,const DomainScalar                           beta
  ) const
{
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<RangeScalar> >
    Y_complex = Thyra::createMembers(complexFFTOp_.range(),X.domain()->dim());
  Thyra::applyTranspose(complexFFTOp_,conj,X,&*Y_complex,RangeScalar(alpha));
  updateReal( *Y_complex, beta, Y );
}

// Overridden from LinearOpWithSolveBase

template<class RealScalar>
bool RealComplexFFTLinearOp<RealScalar>::solveSupportsConj(Thyra::EConj conj) const
{
  return (conj == Thyra::NONCONJ_ELE);
}

template<class RealScalar>
bool RealComplexFFTLinearOp<RealScalar>::solveTransposeSupportsConj(Thyra::EConj conj) const
{
  return (conj == Thyra::CONJ_ELE);
}

template<class RealScalar>
bool RealComplexFFTLinearOp<RealScalar>::solveSupportsSolveMeasureType(Thyra::EConj conj, const Thyra::SolveMeasureType solveMeasureType) const
{
  return (conj == Thyra::NONCONJ_ELE);
}

template<class RealScalar>
bool RealComplexFFTLinearOp<RealScalar>::solveTransposeSupportsSolveMeasureType(Thyra::EConj conj, const Thyra::SolveMeasureType solveMeasureType) const
{
  return (conj == Thyra::CONJ_ELE);
}

template<class RealScalar>
void RealComplexFFTLinearOp<RealScalar>::solve(
  const Thyra::EConj                           conj
  ,const Thyra::MultiVectorBase<RangeScalar>   &B
  ,Thyra::MultiVectorBase<DomainScalar>        *X
  ,const int                                   numBlocks
  ,const Thyra::BlockSolveCriteria<Scalar>     blockSolveCriteria[]
  ,Thyra::SolveStatus<Scalar>                  blockSolveStatus[]
  ) const
{
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<RangeScalar> >
    X_complex = Thyra::createMembers(complexFFTOp_.domain(),B.domain()->dim());
  Thyra::solve(complexFFTOp_,conj,B,&*X_complex,numBlocks,blockSolveCriteria,blockSolveStatus);
  updateReal( *X_complex, Teuchos::ScalarTraits<DomainScalar>::zero(), X );
}

template<class RealScalar>
void RealComplexFFTLinearOp<RealScalar>::solveTranspose(
  const Thyra::EConj                           conj
  ,const Thyra::MultiVectorBase<DomainScalar>  &B
  ,Thyra::MultiVectorBase<RangeScalar>         *X
  ,const int                                   numBlocks
  ,const Thyra::BlockSolveCriteria<Scalar>     blockSolveCriteria[]
  ,Thyra::SolveStatus<Scalar>                  blockSolveStatus[]
  ) const
{
  Teuchos::RefCountPtr< Thyra::MultiVectorBase<RangeScalar> >
    B_complex = Thyra::createMembers(complexFFTOp_.domain(),B.domain()->dim());
  realToComplexConverter_.convert( B, &*B_complex );
  Thyra::solveTranspose(complexFFTOp_,conj,*B_complex,X,numBlocks,blockSolveCriteria,blockSolveStatus);
}

// private

template<class RealScalar>
void RealComplexFFTLinearOp<RealScalar>::updateReal(
  const Thyra::MultiVectorBase<RangeScalar>      &V_complex
  ,const DomainScalar                            &beta
  ,Thyra::MultiVectorBase<DomainScalar>          *V
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( V==NULL );
#endif  
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Thyra::ConstDetachedMultiVectorView<RangeScalar>           ev_V_complex(V_complex);
  Thyra::DetachedMultiVectorView<DomainScalar>   ev_V(*V);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( ev_V_complex.subDim() != ev_V.subDim() || ev_V_complex.numSubCols() != ev_V.numSubCols() );
#endif  
  if( beta == ST::zero() ) {
    for( Thyra::Index j = 0; j < ev_V.numSubCols(); ++j ) {
      for( Thyra::Index i = 0; i < ev_V.subDim(); ++i ) {
        ev_V(i,j) = ev_V_complex(i,j).real();
      }
    }
  }
  else if( beta == ST::one() ) {
    for( Thyra::Index j = 0; j < ev_V.numSubCols(); ++j ) {
      for( Thyra::Index i = 0; i < ev_V.subDim(); ++i ) {
        ev_V(i,j) += ev_V_complex(i,j).real();
      }
    }
  }
  else {
    for( Thyra::Index j = 0; j < ev_V.numSubCols(); ++j ) {
      for( Thyra::Index i = 0; i < ev_V.subDim(); ++i ) {
        ev_V(i,j) = beta*ev_V(i,j) + ev_V_complex(i,j).real();
      }
    }
  }
}

#endif	// THYRA_REAL_COMPLEX_FFT_LINEAR_OP_HPP
