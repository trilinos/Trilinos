// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFUALT_LINEAR_OP_SOURCE_DECL_HPP
#define THYRA_DEFUALT_LINEAR_OP_SOURCE_DECL_HPP

#include "Thyra_LinearOpSourceBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief Default implementation of a <tt>LinearOpSourceBase</tt> that just
 * accepts and gives up a single linear operator object.
 */
template<class Scalar>
class DefaultLinearOpSource : virtual public LinearOpSourceBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized.
   */
  DefaultLinearOpSource();

  /** \brief Construct with a non-const linear operator.
   */
  DefaultLinearOpSource(
    const Teuchos::RCP<LinearOpBase<Scalar> >    &op
    );

  /** \brief Construct with a const linear operator.
   */
  DefaultLinearOpSource(
    const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
    );

  /** \brief Initialize with a non-const linear operator.
   */
  void initialize(
    const Teuchos::RCP<LinearOpBase<Scalar> >    &op
    );

  /** \brief Initialize with a const linear operator.
   */
  void initialize(
    const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
    );

  /** \brief Uninitialize.
   *
   * Note: If the client wants to access the underlying linear operator, then
   * it had better grab them with the below access functions before calling
   * this function.
   */
  void uninitialize();

  //@}

  // ToDo: Override functions from Describable!

  /** @name Overridden from LinearOpSourceBase */
  //@{
  /** \brief . */
  bool isOpConst() const;
  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstOp();
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getOp() const;
  //@}
  
private:
  
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> >  op_;
  
};

// //////////////////////////////
// Related functions

/** \brief Create a <tt>DefaultLinearOpSource</tt> object out of a
 * <tt>LinearOpBase</tt> object.
 *
 * \relates DefaultLinearOpSource
 */
template <class Scalar>
Teuchos::RCP<const DefaultLinearOpSource<Scalar> >
defaultLinearOpSource(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
  )
{
  return Teuchos::rcp(new DefaultLinearOpSource<Scalar>(op));
}

} // namespace Thyra

#endif // THYRA_DEFUALT_LINEAR_OP_SOURCE_DECL_HPP
