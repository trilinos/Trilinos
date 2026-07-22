// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFUALT_PRECONDITIONER_DEF_HPP
#define THYRA_DEFUALT_PRECONDITIONER_DEF_HPP

#include "Thyra_DefaultPreconditioner_decl.hpp"
#include "Thyra_LinearOpBase.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultPreconditioner<Scalar>::DefaultPreconditioner()
{}


template <class Scalar>
DefaultPreconditioner<Scalar>::DefaultPreconditioner(
  const Teuchos::RCP<LinearOpBase<Scalar> > &leftPrecOp,
  const Teuchos::RCP<LinearOpBase<Scalar> > &rightPrecOp
 )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( is_null(leftPrecOp) && is_null(rightPrecOp) );
#endif
  if (!is_null(leftPrecOp))
    leftPrecOp_.initialize(leftPrecOp);
  if (!is_null(rightPrecOp))
    rightPrecOp_.initialize(rightPrecOp);
}


template <class Scalar>
DefaultPreconditioner<Scalar>::DefaultPreconditioner(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( is_null(leftPrecOp) && is_null(rightPrecOp) );
#endif
  if (!is_null(leftPrecOp))
    leftPrecOp_.initialize(leftPrecOp);
  if (!is_null(rightPrecOp))
    rightPrecOp_.initialize(rightPrecOp);
}


template <class Scalar>
DefaultPreconditioner<Scalar>::DefaultPreconditioner(
  const Teuchos::RCP<LinearOpBase<Scalar> > &unspecifiedPrecOp
  )
{
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}


template <class Scalar>
DefaultPreconditioner<Scalar>::DefaultPreconditioner(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &unspecifiedPrecOp
  )
{
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeLeft(
  const Teuchos::RCP<LinearOpBase<Scalar> > &leftPrecOp
  )
{
  uninitialize();
  leftPrecOp_.initialize(leftPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeLeft(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp
  )
{
  uninitialize();
  leftPrecOp_.initialize(leftPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeRight(
  const Teuchos::RCP<LinearOpBase<Scalar> > &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeRight(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeLeftRight(
  const Teuchos::RCP<LinearOpBase<Scalar> > &leftPrecOp,
  const Teuchos::RCP<LinearOpBase<Scalar> > &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
  leftPrecOp_.initialize(leftPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeLeftRight(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &leftPrecOp,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
  leftPrecOp_.initialize(leftPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeUnspecified(
  const Teuchos::RCP<LinearOpBase<Scalar> > &unspecifiedPrecOp
  )
{
  uninitialize();
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::initializeUnspecified(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &unspecifiedPrecOp
  )
{
  uninitialize();
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::uninitialize()
{
  leftPrecOp_.uninitialize();
  rightPrecOp_.uninitialize();
  unspecifiedPrecOp_.uninitialize();
}


// Overridden from PreconditionerBase


template <class Scalar>
bool DefaultPreconditioner<Scalar>::isLeftPrecOpConst() const
{
  return leftPrecOp_.isConst();
}


template <class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultPreconditioner<Scalar>::getNonconstLeftPrecOp()
{
  return leftPrecOp_.getNonconstObj();
}


template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultPreconditioner<Scalar>::getLeftPrecOp() const
{
  return leftPrecOp_.getConstObj();
}


template <class Scalar>
bool DefaultPreconditioner<Scalar>::isRightPrecOpConst() const
{
  return rightPrecOp_.isConst();
}


template <class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultPreconditioner<Scalar>::getNonconstRightPrecOp()
{
  return rightPrecOp_.getNonconstObj();
}


template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultPreconditioner<Scalar>::getRightPrecOp() const
{
  return rightPrecOp_.getConstObj();
}


template <class Scalar>
bool DefaultPreconditioner<Scalar>::isUnspecifiedPrecOpConst() const
{
  return unspecifiedPrecOp_.isConst();
}


template <class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultPreconditioner<Scalar>::getNonconstUnspecifiedPrecOp()
{
  return unspecifiedPrecOp_.getNonconstObj();
}


template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultPreconditioner<Scalar>::getUnspecifiedPrecOp() const
{
  return unspecifiedPrecOp_.getConstObj();
}


// Overridden from Teuchos::Describable

                                                
template <class Scalar>
std::string DefaultPreconditioner<Scalar>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{";
  bool wroteOne = false;
  if(!is_null(leftPrecOp_.getConstObj())) {
    if(wroteOne) oss << ",";
    oss << "leftPrecOp=" << leftPrecOp_.getConstObj()->description();
    wroteOne = true;
  }
  if(!is_null(rightPrecOp_.getConstObj())) {
    if(wroteOne) oss << ",";
    oss << "rightPrecOp=" << rightPrecOp_.getConstObj()->description();
    wroteOne = true;
  }
  if(!is_null(unspecifiedPrecOp_.getConstObj())) {
    if(wroteOne) oss << ",";
    oss << "unspecifiedPrecOp=" << unspecifiedPrecOp_.getConstObj()->description();
    wroteOne = true;
  }
  oss << "}";
  return oss.str();
}


template <class Scalar>
void DefaultPreconditioner<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << Teuchos::Describable::description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      out
        << Teuchos::typeName(*this) << "\n";
      OSTab tab2(out);
      if(!is_null(leftPrecOp_.getConstObj()))
        out << "leftPrecOp=" << describe(*leftPrecOp_.getConstObj(),verbLevel);
      if(!is_null(rightPrecOp_.getConstObj()))
        out << "rig htPrecOp=" << describe(*rightPrecOp_.getConstObj(),verbLevel);
      if(!is_null(unspecifiedPrecOp_.getConstObj()))
        out << "unspecifiedPrecOp=" << describe(*unspecifiedPrecOp_.getConstObj(),verbLevel);
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


} // namespace Thyra


#endif // THYRA_DEFUALT_PRECONDITIONER_DEF_HPP
