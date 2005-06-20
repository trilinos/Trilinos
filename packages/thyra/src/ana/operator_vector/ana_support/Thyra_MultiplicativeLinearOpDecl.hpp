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

#ifndef THYRA_MULTIPLICATIVE_LINEAR_OP_DECL_HPP
#define THYRA_MULTIPLICATIVE_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Teuchos_arrayArg.hpp"

namespace Thyra {

/** \brief Concrete composite <tt>LinearOpBase</tt> subclass that creates a
 * multiplicative linear operator out of one or more constituent
 * <tt>LinearOpBase</tt> objects.
 *
 * This class represents a multiplicative linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = gamma * Op[0] * Op[1] * ... * Op[numOps-1]
 \endverbatim
 *
 * where <tt>Op[]</tt> is an array of <tt>numOps</tt>
 * <tt>LinearOp</tt> objects and <tt>gamma</tt> is a scalar.
 * Of course the operator <tt>M</tt> is not constructed explicitly but
 * instead just applies the constituent linear operators accordingly
 * using temporaries.
 *
 * In other words, this class defines <tt>apply()</tt> as:
 *
 \verbatim

 y = alpha*M*x + beta*y
   = (alpha*gamma) * ( Op[0] * ( Op[1] * ( .... ( Op[numOps-1] * x ) ... ) ) ) + beta * y
 \endverbatim
 *
 * for the case where <tt>M_trans==NOTRANS</tt> and as:
 *
 \verbatim

 y = alpha*M'*x + beta*y
   = (alpha*gamma) * ( Op[numOps-1]' * ( .... ( Op[1]' * ( Op[0]' * x ) ) ... ) ) + beta * y
 \endverbatim
 *
 * for the case where <tt>M_trans!=NOTRANS</tt> (where the transpose
 * <tt>'</tt> either defines <tt>TRANS</tt> or <tt>CONJTRANS</tt>).
 *
 * Constructing a multiplicative operator is easy.  For example, suppose one
 * wants to construct the multiplicative operator <tt>D = gamma * A * B' * C</tt>.
 * To do so one would do:

 \code
 template<class Scalar>
 void constructD(
    const Scalar                                                   &gamma
    ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >  &A
    ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >  &B
    ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >  &C
    ,Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >        *D
    )
 {
   typedef Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> > LOB;
   *D = Teuchos::rcp(
     new Thyra::MultiplicativeLinearOp<Scalar>(
       3, Teuchos::arrayArg<LOB>(A,adjoin(B),C)(), gamma
       )
     );
 }
 \endcode
 *
 * Rather than calling the constructor directly, consider using the non-member helper
 * functions described \ref Thyra_Op_Vec_MultipliciateLinearOp_helpers_grp "here".
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class MultiplicativeLinearOp : virtual public LinearOpBase<Scalar> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->numOps()==0</tt>
   * </ul>
   */
  MultiplicativeLinearOp();

  /** Calls <tt>initialize()</tt>.
   *
   * Rather than calling this constructor directly, consider using the non-member helper
   * functions described \ref Thyra_Op_Vec_MultipliciateLinearOp_helpers_grp "here".
   */
  MultiplicativeLinearOp(
    const int                                                   numOps
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    Ops[]
    ,const Scalar                                               &gamma = Teuchos::ScalarTraits<Scalar>::one()
    );

  /** \brief Initialize given a list of linear operators.
   *
   * @param  numOps  [in] Number of constituent operators.
   * @param  Ops     [in] Array (length <tt>numOps</tt>) of
   *                 constituent linear operators and their
   *                 aggregated default definitions of the
   *                 non-transposed operator.
   * @param  gamma   [in] Scalar multiplier
   *
   * Preconditions:<ul>
   * <li><tt>numOps > 0</tt>
   * <li><tt>Ops != NULL</tt>
   * <li><tt>Ops[k].op().get()!=NULL</tt>, for <tt>k=0...numOps-1</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->numOps()==numOps</tt>
   * <li><tt>this->getOp(k).op().get()==Ops[k].op().get()</tt>, for <tt>k=0...numOps-1</tt>
   * <li><tt>this->gamma()==gamma</tt>
   * </ul>
   */
  void initialize(
    const int                                                   numOps
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    Ops[]
    ,const Scalar                                               &gamma = Teuchos::ScalarTraits<Scalar>::one()
    );

  /** \brief Returns the current number of constituent operators.
   *
   * A return value of <tt>0</tt> indicates that <tt>this</tt> is not
   * fully initialized.
   */
  int numOps() const;

  /** \brief Return the <tt>k</tt>th constituent operator.
   *
   * @param  k  [in] The zero-based index of the constituent operator to return.
   *
   * Preconditions:<ul>
   * <li><tt> 0 <= k < this->numOps()</tt>
   * </ul>
   */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > getOp(const int k) const;

  /** \brief Set to uninitialized.
   *
   * @param  numOps  [in] Number of operators (must be <tt>numOps==this->numOps()</tt>).
   * @param  Ops     [out] Array (length <tt>numOps</tt>) that if <tt>Ops!=NULL</tt>
   *                 then <tt>Ops[k]</tt> will be set to <tt>this->getOp(k)</tt>, for
   *                 <tt>k=0...numOps-1</tt>.
   * @param  gamma [out] Optional pointer to scalar <tt>gamma</tt>.
   *               If <tt>gamma!=NULL</tt> then on output <tt>*gamma</tt>
   *               is set to <tt>this->gamma()</tt> (before call).
   *
   * Preconditions:<ul>
   * <li>[<tt>Ops!=NULL</tt>] <tt>numOps==this->numOps()</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->numOps()==0</tt>
   * </ul>
   */
  void uninitialize(
    const int                                             numOps   = 0
    ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    Ops[]    = NULL
    ,Scalar                                               *gamma   = NULL
    );

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Prints just the name <tt>MultiplicativeLinearOp</tt> along with
   * the overall dimensions and the number of constituent operators.
   */
  std::string description() const;

  /** \brief Prints the details about the constituent linear operators.
   *
   * This function outputs different levels of detail based on the value passed in
   * for <tt>verbLevel</tt>:
   *
   * ToDo: Finish documentation!
   */
  std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const;

  //@}

  /** @name Overridden from OpBase */
  //@{
  /** \brief Returns <tt>this->getOp(0).range()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>this->numOps()==0</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  /** \brief Returns <tt>this->getOp(this->numOps()-1).domain()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>this->numOps()==0</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief Returns <tt>true</tt> only if all constituent operators support
   * <tt>M_trans</tt>.
   */
  bool opSupported(ETransp M_trans) const;
  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const;
  /** \brief . */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;

  //@}

private:

  std::vector<Teuchos::RefCountPtr<const LinearOpBase<Scalar> > >   Ops_;
  Scalar                                                            gamma_;

  void assertInitialized() const;

};

/** \defgroup Thyra_Op_Vec_MultipliciateLinearOp_helpers_grp  Helper functions for creating MultiplicativeLinearOp objects.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */

/** \brief Form an implicit multiplication of two linear operators: <tt>return = A * B</tt>.
 *
 * This function simply creates a <tt>MultiplicativeLinearOp</tt> given
 * <tt>A</tt> and <tt>B</tt>.
 *
 * \ingroup Thyra_Op_Vec_MultipliciateLinearOp_helpers_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
multiply(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &A
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &B
  ,const Scalar                                              &gamma
#ifndef __sun
                                                                    = Teuchos::ScalarTraits<Scalar>::one()
#endif
  );

#ifdef __sun
template<class Scalar>
inline Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
multiply(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &A
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &B
  )
{
  return multiply(A,B,Teuchos::ScalarTraits<Scalar>::one());
}
#endif

// /////////////////////////////////
// Inline members

template<class Scalar>
inline
int MultiplicativeLinearOp<Scalar>::numOps() const
{
  return Ops_.size();
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const LinearOpBase<Scalar> > MultiplicativeLinearOp<Scalar>::getOp(const int k) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k];
}

template<class Scalar>
inline
void MultiplicativeLinearOp<Scalar>::assertInitialized() const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( !( numOps() > 0 ) );
#endif
}

}	// end namespace Thyra

// ///////////////////////
// Inline non-members

template<class Scalar>
inline
Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
Thyra::multiply(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &A
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &B
  ,const Scalar                                              &gamma
  )
{
  return Teuchos::rcp(
    new MultiplicativeLinearOp<Scalar>(
      2
      ,Teuchos::arrayArg<Teuchos::RefCountPtr<const LinearOpBase<Scalar> > >(A,B)()
      ,gamma
      )
    );
}

#endif	// THYRA_MULTIPLICATIVE_LINEAR_OP_DECL_HPP
