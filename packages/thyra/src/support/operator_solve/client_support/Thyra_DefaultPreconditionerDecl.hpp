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

#ifndef THYRA_DEFUALT_PRECONDITIONER_DECL_HPP
#define THYRA_DEFUALT_PRECONDITIONER_DECL_HPP

#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief Default implementation of a <tt>PreconditionerBase</tt> that just
 * accepts precreated preconditioner linear operators.
 *
 * Here is how to construct a precondtioner for the four different types of preconditioners:
 * <ul>
 * <li>Single preconditioner linear operator designed or targeted to be applied on the left:
 *   <ul>
 *   <li><tt>DefaultPreconditioner(leftPrecOp,Teuchos::null);
 *   </ul>
 * <li>Single preconditioner linear operator designed or targeted to be applied on the right:
 *   <ul>
 *   <li><tt>DefaultPreconditioner(Teuchos::null,righPrecOp);
 *   </ul>
 * <li>Split two-sided preconditioner with linear operators designed or
 *   targeted to be applied on the left and the right:
 *   <ul>
 *   <li><tt>DefaultPreconditioner(leftPrecOp,rightPrecOp);
 *   </ul>
 * <li>Single preconditioner linear operator not designed or targeted to be
 *   applied on the left or the right:
 *   <ul>
 *   <li><tt>DefaultPreconditioner(precOp);
 *   </ul>
 * <\ul>
 *
 * ToDo: Finish documentation!
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class DefaultPreconditioner : virtual public PreconditionerBase<RangeScalar,DomainScalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized.
   */
  DefaultPreconditioner();

  /** \brief Construct a left-only, or right-only, or split left/right
   * preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
    ,const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
    );

  /** \brief Construct a const-only left-only, or right-only, or split
   * left/right preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
    );

  /** \brief Construct a single unspecified preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
    );

  /** \brief Construct a const-only single unspecified preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
    );

  /** \brief Initialize a left preconditioner.
   */
  void initializeLeft(
    const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
    );

  /** \brief Initialize a const-only left preconditioner.
   */
  void initializeLeft(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
    );

  /** \brief Initialize a right preconditioner.
   */
  void initializeRight(
    const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &rightPrecOp
    );

  /** \brief Initialize a const-only right preconditioner.
   */
  void initializeRight(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &rightPrecOp
    );

  /** \brief Initialize a split left/right preconditioner.
   */
  void initializeLeftRight(
    const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
    ,const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
    );

  /** \brief Initialize a const-only split left/right preconditioner.
   */
  void initializeLeftRight(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
    );

  /** \brief Initialize a single unspecified preconditioner
   * operator.
   */
  void initializeUnspecified(
    const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
    );

  /** \brief Initialize a const-only single unspecified preconditioner
   * operator.
   */
  void initializeUnspecified(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
    );

  /** \brief Uninitialize.
   *
   * Note: If the client wants to access the underlying preconditioner
   * operators, then it had better grab them with the below access functions
   * before calling this function.
   */
  void uninitialize();

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;
  //@}

  /** @name Overridden from PreconditionerBase */
  //@{
  /** \brief . */
  bool isLeftPrecOpConst() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> > getNonconstLeftPrecOp();
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> > getLeftPrecOp() const;
  /** \brief . */
  bool isRightPrecOpConst() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> > getNonconstRightPrecOp();
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> > getRightPrecOp() const;
  /** \brief . */
  bool isUnspecifiedPrecOpConst() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> > getNonconstUnspecifiedPrecOp();
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> > getUnspecifiedPrecOp() const;
  //@}
  
private:
  
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<RangeScalar,DomainScalar> >  leftPrecOp_;
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<RangeScalar,DomainScalar> >  rightPrecOp_;
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<RangeScalar,DomainScalar> >  unspecifiedPrecOp_;
  
};

// ///////////////////////
// Related functions

/** \brief Create a precondioner from a single linear operator not targeted to
 * be used on the left or the right.
 *
 * \relates DefaultPreconditioner
 */
template <class Scalar>
Teuchos::RefCountPtr<const DefaultPreconditioner<Scalar> >
unspecifiedPrec(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &unspecifiedPrecOp
  )
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>(unspecifiedPrecOp));
}

/** \brief Create a precondioner from a single linear operator targeted to be
 * used on the left.
 *
 * \relates DefaultPreconditioner
 */
template <class Scalar>
Teuchos::RefCountPtr<const DefaultPreconditioner<Scalar> >
leftPrec(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &leftPrecOp
  )
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>(leftPrecOp,Teuchos::null));
}

/** \brief Create a precondioner from a single linear operator targeted to be
 * used on the right.
 *
 * \relates DefaultPreconditioner
 */
template <class Scalar>
Teuchos::RefCountPtr<const DefaultPreconditioner<Scalar> >
rightPrec(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &rightPrecOp
  )
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>(Teuchos::null,rightPrecOp));
}

/** \brief Create a split precondioner from two linear operators, one to be
 * applied on the left and one to be applied on the right.
 *
 * \relates DefaultPreconditioner
 */
template <class Scalar>
Teuchos::RefCountPtr<const DefaultPreconditioner<Scalar> >
splitPrec(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >     &leftPrecOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &rightPrecOp
  )
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>(leftPrecOp,rightPrecOp));
}

} // namespace Thyra

#endif // THYRA_DEFUALT_PRECONDITIONER_DECL_HPP
