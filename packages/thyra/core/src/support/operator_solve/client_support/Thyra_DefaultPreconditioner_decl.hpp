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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
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
