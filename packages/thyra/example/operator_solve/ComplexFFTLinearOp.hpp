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

#ifndef THYRA_COMPLEX_FFT_LINEAR_OP_HPP
#define THYRA_COMPLEX_FFT_LINEAR_OP_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "serial_1D_FFT.hpp"

/** \brief Simple concrete subclass for a serial complex-to-complex FFT.
 *
 * This implementation uses orthonormal columns and rows and therefore the
 * adjoint is the same as the inverse.
 *
 * \ingroup Thyra_Op_Vec_examples_fft_grp
 */
template<class RealScalar>
class ComplexFFTLinearOp
  : virtual public Thyra::LinearOpWithSolveBase< std::complex<RealScalar> >
  , virtual protected Thyra::SingleRhsLinearOpWithSolveBase< std::complex<RealScalar> >
{
public:

  /** \brief . */
  typedef std::complex<RealScalar> Scalar;
  /** \brief . */
  ComplexFFTLinearOp( const int N );

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase< std::complex<RealScalar> > > range() const;
  /** \brief . */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase< std::complex<RealScalar> > > domain() const;
  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{
  /** \brief . */
  bool opSupported(Thyra::ETransp M_trans) const;
  //@}

  /** @name Overridden from SingleRhsLinearOpBase */
  //@{
  /** \brief . */
  void apply(
    const Thyra::ETransp                                    M_trans
    ,const Thyra::VectorBase< std::complex<RealScalar> >    &x
    ,Thyra::VectorBase< std::complex<RealScalar> >          *y
    ,const std::complex<RealScalar>                         alpha
    ,const std::complex<RealScalar>                         beta
    ) const;
  //@}

  /** @name Overridden from SingleScalarLinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsTrans(Thyra::ETransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureType(
    Thyra::ETransp M_trans, const Thyra::SolveMeasureType& solveMeasureType
    ) const;
  //@}

  /** @name Overridden from SingleRhsLinearOpWithSolveBase */
  //@{
  /** \brief . */
  Thyra::SolveStatus< std::complex<RealScalar> > solve(
    const Thyra::ETransp                                             M_trans
    ,const Thyra::VectorBase< std::complex<RealScalar> >             &b
    ,Thyra::VectorBase< std::complex<RealScalar> >                   *x
    ,const Thyra::SolveCriteria< std::complex<RealScalar> >          *solveCriteria
    ) const;
  //@}

private:

  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase< std::complex<RealScalar> > > space_;

};

// /////////////////////////
// Definitions

template<class RealScalar>
ComplexFFTLinearOp<RealScalar>::ComplexFFTLinearOp( const int N )
{
  space_ = Teuchos::rcp(
    new Thyra::DefaultSpmdVectorSpace< std::complex<RealScalar> >(
      int(std::pow(double(2),N))
      )
    );
}

// Overridden from LinearOpBase

template<class RealScalar>
Teuchos::RefCountPtr< const Thyra::VectorSpaceBase< std::complex<RealScalar> > >
ComplexFFTLinearOp<RealScalar>::range() const
{
  return space_;
}

template<class RealScalar>
Teuchos::RefCountPtr< const Thyra::VectorSpaceBase< std::complex<RealScalar> > >
ComplexFFTLinearOp<RealScalar>::domain() const
{
  return space_;
}

// Overridden from SingleScalarLinearOpBase

template<class RealScalar>
bool ComplexFFTLinearOp<RealScalar>::opSupported(Thyra::ETransp M_trans) const
{
  return ( M_trans == Thyra::NOTRANS || M_trans == Thyra::CONJTRANS );
}

// Overridden from SingleRhsLinearOpBase

template<class RealScalar>
void ComplexFFTLinearOp<RealScalar>::apply(
  const Thyra::ETransp                                    M_trans
  ,const Thyra::VectorBase< std::complex<RealScalar> >    &x
  ,Thyra::VectorBase< std::complex<RealScalar> >          *y
  ,const std::complex<RealScalar>                         alpha
  ,const std::complex<RealScalar>                         beta
  ) const
{
  typedef Teuchos::ScalarTraits< std::complex<RealScalar> > ST;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( M_trans == Thyra::NOTRANS || M_trans == Thyra::CONJTRANS ) );
#endif
  // Update y first
  Thyra::Vt_S( y, beta );
  // Translate from input x into one long array with data[] that will be
  // passed to and from the FFT routine
  const Thyra::ConstDetachedVectorView<Scalar>  x_ev(x);
  std::vector<RealScalar> data(2*x_ev.subDim());
  for( int k = 0; k < x_ev.subDim(); ++k ) {
    data[2*k]   = x_ev[k].real();
    data[2*k+1] = x_ev[k].imag();
  }
  // Call the FFT rountine
  serial_1D_FFT(
    &data[0]-1                                            // This function is 1-based of all things!
    ,x_ev.subDim()                                        // 1/2 length of data[]
    ,Thyra::real_trans(M_trans)==Thyra::NOTRANS ? +1 : -1 // +1 = fwd, -1 = adjoint
    );
  // Add the scaled result into y
  const Thyra::DetachedVectorView<Scalar>  y_ev(*y);
  const Scalar scalar = alpha * Scalar(1)/ST::squareroot(x_ev.subDim()); // needed to make adjoint == inverse!
  for( int k = 0; k < y_ev.subDim(); ++k ) {
    y_ev[k] += ( scalar * Scalar(data[2*k],data[2*k+1]) );
  }
}

// Overridden from SingleScalarLinearOpWithSolveBase

template<class RealScalar>
bool ComplexFFTLinearOp<RealScalar>::solveSupportsTrans(Thyra::ETransp M_trans) const
{
  return ( M_trans == Thyra::NOTRANS || M_trans == Thyra::CONJTRANS );
}

template<class RealScalar>
bool ComplexFFTLinearOp<RealScalar>::solveSupportsSolveMeasureType(
  Thyra::ETransp M_trans, const Thyra::SolveMeasureType& solveMeasureType
  ) const
{
  return ( M_trans == Thyra::NOTRANS || M_trans == Thyra::CONJTRANS );
}

// Overridden from SingleRhsLinearOpWithSolveBase

template<class RealScalar>
Thyra::SolveStatus< std::complex<RealScalar> >
ComplexFFTLinearOp<RealScalar>::solve(
  const Thyra::ETransp                                             M_trans
  ,const Thyra::VectorBase< std::complex<RealScalar> >             &b
  ,Thyra::VectorBase< std::complex<RealScalar> >                   *x
  ,const Thyra::SolveCriteria< std::complex<RealScalar> >          *solveCriteria
  ) const
{
  typedef Teuchos::ScalarTraits< std::complex<RealScalar> > ST;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( M_trans == Thyra::NOTRANS || M_trans == Thyra::CONJTRANS ) );
#endif
  Thyra::apply( *this, M_trans==Thyra::NOTRANS?Thyra::CONJTRANS:Thyra::NOTRANS, b, x );
  typedef Thyra::SolveCriteria< std::complex<RealScalar> >  SC;
  typedef Thyra::SolveStatus< std::complex<RealScalar> >    SS;
  SS solveStatus;
  solveStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
  solveStatus.achievedTol = SS::unknownTolerance();
  return solveStatus;
}

#endif	// THYRA_COMPLEX_FFT_LINEAR_OP_HPP
