// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFUALT_PRECONDITIONER_DECL_HPP
#define THYRA_DEFUALT_PRECONDITIONER_DECL_HPP

#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Default implementation of a <tt>PreconditionerBase</tt> that just
 * accepts precreated preconditioner linear operators.
 *
 * Here is how to construct a preconditioner for the four different types of preconditioners:
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
template <class Scalar>
class DefaultPreconditioner : virtual public PreconditionerBase<Scalar>
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
    const Teuchos::RCP<LinearOpBase<Scalar> > &leftPrecOp,
    const Teuchos::RCP<LinearOpBase<Scalar> > &rightPrecOp
    );

  /** \brief Construct a const-only left-only, or right-only, or split
   * left/right preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp,
    const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
    );

  /** \brief Construct a single unspecified preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RCP<LinearOpBase<Scalar> > &unspecifiedPrecOp
    );

  /** \brief Construct a const-only single unspecified preconditioner.
   */
  DefaultPreconditioner(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &unspecifiedPrecOp
    );

  /** \brief Initialize a left preconditioner.
   */
  void initializeLeft(
    const Teuchos::RCP<LinearOpBase<Scalar> > &leftPrecOp
    );

  /** \brief Initialize a const-only left preconditioner.
   */
  void initializeLeft(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp
    );

  /** \brief Initialize a right preconditioner.
   */
  void initializeRight(
    const Teuchos::RCP<LinearOpBase<Scalar> > &rightPrecOp
    );

  /** \brief Initialize a const-only right preconditioner.
   */
  void initializeRight(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
    );

  /** \brief Initialize a split left/right preconditioner.
   */
  void initializeLeftRight(
    const Teuchos::RCP<LinearOpBase<Scalar> > &leftPrecOp
    ,const Teuchos::RCP<LinearOpBase<Scalar> > &rightPrecOp
    );

  /** \brief Initialize a const-only split left/right preconditioner.
   */
  void initializeLeftRight(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp
    ,const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
    );

  /** \brief Initialize a single unspecified preconditioner
   * operator.
   */
  void initializeUnspecified(
    const Teuchos::RCP<LinearOpBase<Scalar> > &unspecifiedPrecOp
    );

  /** \brief Initialize a const-only single unspecified preconditioner
   * operator.
   */
  void initializeUnspecified(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &unspecifiedPrecOp
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
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstLeftPrecOp();
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getLeftPrecOp() const;
  /** \brief . */
  bool isRightPrecOpConst() const;
  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstRightPrecOp();
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getRightPrecOp() const;
  /** \brief . */
  bool isUnspecifiedPrecOpConst() const;
  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstUnspecifiedPrecOp();
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getUnspecifiedPrecOp() const;
  //@}
  
private:
  
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > leftPrecOp_;
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > rightPrecOp_;
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > unspecifiedPrecOp_;
  
};

// ///////////////////////
// Related functions


/** \brief Create a precondioner from a single linear operator not targeted to
 * be used on the left or the right.
 *
 * \relates DefaultPreconditioner
 */
template <class Scalar>
Teuchos::RCP<const DefaultPreconditioner<Scalar> >
unspecifiedPrec(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &unspecifiedPrecOp
  )
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>(unspecifiedPrecOp));
}


/** \brief Create a precondioner from a single linear operator not targeted to
 * be used on the left or the right.
 *
 * \relates DefaultPreconditioner
 */
template <class Scalar>
Teuchos::RCP<DefaultPreconditioner<Scalar> >
nonconstUnspecifiedPrec(
  const Teuchos::RCP<LinearOpBase<Scalar> > &unspecifiedPrecOp
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
Teuchos::RCP<const DefaultPreconditioner<Scalar> >
leftPrec(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp
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
Teuchos::RCP<const DefaultPreconditioner<Scalar> >
rightPrec(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
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
Teuchos::RCP<const DefaultPreconditioner<Scalar> >
splitPrec(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp
  ,const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
  )
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>(leftPrecOp,rightPrecOp));
}


} // namespace Thyra


#endif // THYRA_DEFUALT_PRECONDITIONER_DECL_HPP
