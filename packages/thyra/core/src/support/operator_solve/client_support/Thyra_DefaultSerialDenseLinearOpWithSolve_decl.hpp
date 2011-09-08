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


#ifndef THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_DECL_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "RTOpPack_LapackWrappers.hpp"


namespace Thyra {


/* \brief . */
inline RTOpPack::ETransp convertToRTOpPackETransp( const EOpTransp transp )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(transp == CONJ);
#endif
  switch(transp) {
    case NOTRANS:
      return RTOpPack::NOTRANS;
    case TRANS:
      return RTOpPack::TRANS;
    case CONJTRANS:
      return RTOpPack::CONJTRANS;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return RTOpPack::NOTRANS; // will never be called!
}
// ToDo: Move the above function into Thyra_OperatorVectorTypes.hpp


/** \brief Simple concreate subclass of <tt>LinearOpWithSolveBase</tt> for
 * serial dense matrices implemented using LAPACK.
 *
 * This class uses the helper class <tt>DetachedMultiVectorView</tt> to
 * extract an explicit view of the matrix elements and then uses
 * <tt>Teuchos::LAPACK</tt> to factor <tt>M = L * U</tt> and then do
 * back-solves with the factors <tt>L</tt> and <tt>U</tt>.
 *
 * Even through this class accesses explicit matrix entries and is called
 * <tt>SerialDense</tt>, it is still considered an ANA subclass since it does
 * not have any direct dependance on a specific computing environment or
 * concreate operator/vector/vectorspace implementation.
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class DefaultSerialDenseLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<Scalar>
{
public:
  
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  DefaultSerialDenseLinearOpWithSolve();

  /** \brief . */
  void initialize( const RCP<const MultiVectorBase<Scalar> > &M );

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > getFwdOp() const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > domain() const;

  //@}

protected:
  
  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
 void applyImpl(
   const EOpTransp M_trans,
   const MultiVectorBase<Scalar> &X,
   const Ptr<MultiVectorBase<Scalar> > &Y,
   const Scalar alpha,
   const Scalar beta
   ) const;

  //@}

  /** @name Overridden from LinearOpWithSolveBase */
  //@{

  /** \brief . */
  bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const;
  /** \brief . */
  SolveStatus<Scalar> solveImpl(
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;

 //@}

private:

  // /////////////////////////
  // Private data members

  RCP<const MultiVectorBase<Scalar> > M_;
  RTOpPack::ConstSubMultiVectorView<Scalar> LU_;
  Array<int> ipiv_;

  // /////////////////////////
  // Private member functions

  static void factorize(
    const MultiVectorBase<Scalar> &M,
    const Ptr<RTOpPack::ConstSubMultiVectorView<Scalar> > &LU,
    const Ptr<Array<int> > &ipiv
    );

  static void backsolve(
    const RTOpPack::ConstSubMultiVectorView<Scalar> &LU,
    const ArrayView<const int> ipiv,
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X
    );

  // Not defined and not to be called
  DefaultSerialDenseLinearOpWithSolve(const DefaultSerialDenseLinearOpWithSolve&);
  DefaultSerialDenseLinearOpWithSolve& operator=(const DefaultSerialDenseLinearOpWithSolve&);

};


/** \brief Nonmember constructor.
 *
 * \relates DefaultSerialDenseLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultSerialDenseLinearOpWithSolve<Scalar> >
defaultSerialDenseLinearOpWithSolve()
{
  return Teuchos::rcp(new DefaultSerialDenseLinearOpWithSolve<Scalar>);
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultSerialDenseLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultSerialDenseLinearOpWithSolve<Scalar> >
defaultSerialDenseLinearOpWithSolve( const RCP<const MultiVectorBase<Scalar> > &M )
{
  RCP<DefaultSerialDenseLinearOpWithSolve<Scalar> >
    M_lows = Teuchos::rcp(new DefaultSerialDenseLinearOpWithSolve<Scalar>());
  M_lows->initialize(M);  // With throw if singular
  return M_lows;
}


}	// end namespace Thyra


#endif	// THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_DECL_HPP
