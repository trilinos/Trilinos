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

#ifndef THYRA_LINEAROPERATOR_DECL_HPP
#define THYRA_LINEAROPERATOR_DECL_HPP

#include "Teuchos_Handle.hpp"
#include "Thyra_VectorSpaceDecl.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_LinearOpBaseDecl.hpp"

namespace Thyra
{

  /** Read-only handle class for <tt>Thyra::LinearOpBase</tt> objects which
   * supports operator-overloading implicit linear operator construction.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
 template <class RangeScalar, class DomainScalar=RangeScalar>
  class ConstLinearOperator 
    : public virtual Teuchos::ConstHandle<LinearOpBase<RangeScalar, DomainScalar> >
  {
  public:

    typedef RangeScalar Scalar;

    /** \brief . */
    ConstLinearOperator( const Teuchos::ENull _null = Teuchos::null )
      : Teuchos::ConstHandle<LinearOpBase<RangeScalar, DomainScalar> >(){;}

    /** \brief Construct from a raw pointer */
    ConstLinearOperator(Teuchos::ConstHandleable<LinearOpBase<RangeScalar, DomainScalar> >* rawPtr) 
      : Teuchos::ConstHandle<LinearOpBase<RangeScalar, DomainScalar> >(rawPtr){;}

    /** \brief Construct from a smart pointer */
    ConstLinearOperator(const Teuchos::RCP<const LinearOpBase<RangeScalar, DomainScalar> >& smartPtr) 
      : Teuchos::ConstHandle<LinearOpBase<RangeScalar, DomainScalar> >(smartPtr){;}

    /** \brief Return the domain space */
    const VectorSpace<DomainScalar> domain() const ;

    /** \brief Return the range space */
    const VectorSpace<RangeScalar> range() const ;

    /** \brief Apply the linear operator
     *
     * Compute
     * \code
     * out = beta*out + alpha*op*in;
     * \endcode
     **/
    void apply(const ConstVector<DomainScalar>& in,
               Vector<RangeScalar>& out,
               const RangeScalar& alpha = 1.0,
               const RangeScalar& beta = 0.0) const ;

    /** \brief Apply the transpose of the linear operator
     *
     * Compute
     * \code
     * out = beta*out + alpha*op^T*in;
     * \endcode
     **/
    void applyTranspose(const ConstVector<RangeScalar>& in,
                        Vector<DomainScalar>& out,
                        const DomainScalar& alpha = 1.0,
                        const DomainScalar& beta = 0.0) const ;


    /** \brief Return the number of block rows */
    int numBlockRows() const ;

    /** \brief Return the number of block columns */
    int numBlockCols() const ;

    /** \brief Return the (blockRow, blockCol)-th subblock */
    ConstLinearOperator<RangeScalar, DomainScalar> getBlock(int blockRow, int blockCol) const ;

  };

  /** Handle class for <tt>Thyra::LinearOpBase</tt> objects which supports
   * operator-overloading implicit linear operator construction.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
  template <class RangeScalar, class DomainScalar=RangeScalar>
  class LinearOperator 
    : public Teuchos::Handle<LinearOpBase<RangeScalar, DomainScalar> >,
      public ConstLinearOperator<RangeScalar, DomainScalar>
  {
  public:

    typedef RangeScalar Scalar;

    /** \brief .  */
    LinearOperator( const Teuchos::ENull _null = Teuchos::null )
      : Teuchos::Handle<LinearOpBase<RangeScalar, DomainScalar> >(){;}

    /** \brief Construct from a raw pointer */
    LinearOperator(Teuchos::Handleable<LinearOpBase<RangeScalar, DomainScalar> >* rawPtr) 
      : Teuchos::Handle<LinearOpBase<RangeScalar, DomainScalar> >(rawPtr){;}

    /** \brief Construct from a smart pointer */
    LinearOperator(const Teuchos::RCP<LinearOpBase<RangeScalar, DomainScalar> >& smartPtr) 
      : Teuchos::Handle<LinearOpBase<RangeScalar, DomainScalar> >(smartPtr){;}

    /** \brief Return the (blockRow, blockCol)-th subblock */
    LinearOperator<RangeScalar, DomainScalar> getBlock(int blockRow, int blockCol) ;

  };

  /** \brief Implicitly scale a linear operator.
   *
   * \relates ConstLinearOperator
   */
  template <class Scalar>
  ConstLinearOperator<Scalar>
  operator*(const Scalar& a, 
            const ConstLinearOperator<Scalar>& A);

  /** \brief Implicitly scale a linear operator.
   *
   * \relates LinearOperator
   */
  template <class Scalar>
  LinearOperator<Scalar>
  operator*(const Scalar& a, 
            const LinearOperator<Scalar>& A);

  /** \brief Implicitly scale a linear operator.
   *
   * \relates ConstLinearOperator
   */
  template <class Scalar>
  ConstLinearOperator<Scalar>
  operator*(const ConstLinearOperator<Scalar>& A,
            const Scalar& a);

  /** \brief Implicitly scale a linear operator.
   *
   * \relates LinearOperator
   */
  template <class Scalar>
  LinearOperator<Scalar>
  operator*(const LinearOperator<Scalar>& A,
            const Scalar& a);

  /** \brief Implicitly multiply two linear operators.
   *
   * \relates ConstLinearOperator
   */
  template <class Scalar>
  ConstLinearOperator<Scalar>
  operator*(const ConstLinearOperator<Scalar>& A,
            const ConstLinearOperator<Scalar>& B);

  /** \brief Implicitly multiply two linear operators.
   *
   * \relates LinearOperator
   */
  template <class Scalar>
  LinearOperator<Scalar>
  operator*(const LinearOperator<Scalar>& A,
            const LinearOperator<Scalar>& B);

  /** \brief Implicitly add two linear operators.
   *
   * \relates ConstLinearOperator
   */
  template <class Scalar>
  ConstLinearOperator<Scalar>
  operator+(const ConstLinearOperator<Scalar>& A,
            const ConstLinearOperator<Scalar>& B);

  /** \brief Implicitly add two linear operators.
   *
   * \relates LinearOperator
   */
  template <class Scalar>
  LinearOperator<Scalar>
  operator+(const LinearOperator<Scalar>& A,
            const LinearOperator<Scalar>& B);
    
  /** \brief Form an implicit block 2x2 linear operator <tt>[ A00, A01; A10,
   * A11 ]</tt>.
   *
   * \relates ConstLinearOperator
   */
  template<class Scalar>
  ConstLinearOperator<Scalar>
  block2x2(
           const ConstLinearOperator<Scalar>&   A00,
           const ConstLinearOperator<Scalar>&   A01,
           const ConstLinearOperator<Scalar>&   A10,
           const ConstLinearOperator<Scalar>&   A11
           );

  /** \brief Form an implicit block 2x1 linear operator <tt>[ A00; A10 ]</tt>.
   *
   * \relates ConstLinearOperator
   */
  template<class Scalar>
  ConstLinearOperator<Scalar>
  block2x1(
           const ConstLinearOperator<Scalar>&    A00,
           const ConstLinearOperator<Scalar>&   A10
           );

  /** \brief Form an implicit block 1x2 linear operator <tt>[ A00, A01 ]</tt>.
   *
   * \relates ConstLinearOperator
   */
  template<class Scalar>
  ConstLinearOperator<Scalar>
  block1x2(
           const ConstLinearOperator<Scalar>&    A00,
           const ConstLinearOperator<Scalar>&   A01
           );
  
  /** \brief Form an implicit block 2x2 linear operator <tt>[ A00, A01; A10,
   * A11 ]</tt>.
   *
   * \relates LinearOperator
   */
  template<class Scalar>
  LinearOperator<Scalar>
  block2x2(
           const LinearOperator<Scalar>&    A00,
           const LinearOperator<Scalar>&   A01,
           const LinearOperator<Scalar>&   A10,
           const LinearOperator<Scalar>&   A11
           );

  /** \brief Form an implicit block 2x1 linear operator <tt>[ A00; A10 ]</tt>.
   *
   * \relates LinearOperator
   */
  template<class Scalar>
  LinearOperator<Scalar>
  block2x1(
           const LinearOperator<Scalar>&    A00,
           const LinearOperator<Scalar>&   A10
           );

  /** \brief Form an implicit block 1x2 linear operator <tt>[ A00, A01 ]</tt>.
   *
   * \relates LinearOperator
   */
  template<class Scalar>
  LinearOperator<Scalar>
  block1x2(
           const LinearOperator<Scalar>&    A00,
           const LinearOperator<Scalar>&   A01
           );

  /** \brief Create an identity operator. */
  template<class Scalar>
  ConstLinearOperator<Scalar>
  identity( const VectorSpace<Scalar> &space );

  /** \brief Create an identity operator. */
  template<class Scalar>
  ConstLinearOperator<Scalar>
  zero( const VectorSpace<Scalar> &range, const VectorSpace<Scalar> &domain );

} // namespace Thyra

#endif
