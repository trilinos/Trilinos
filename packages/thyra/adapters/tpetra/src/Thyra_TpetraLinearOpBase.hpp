
#ifndef THYRA_TPETRA_LINEAR_OP_BASE_HPP
#define THYRA_TPETRA_LINEAR_OP_BASE_HPP

#include "Thyra_TpetraTypes.hpp"
#include "Thyra_SingleScalarEuclideanLinearOpBase.hpp"

namespace Thyra {

/** \brief Abstract base class for all <tt>LinearOpBase</tt> objects that can
 * return an <tt>Tpetra::Operator</tt> view of themselves.
 *
 * This interface defines a key interoperability interface that allows a
 * general client to extract an <tt>Tpetra::Operator</tt> view of a linear
 * operator object.  In some cases, <tt>*this</tt> linear operator can be
 * modifed through the <tt>Tpetra::Operator</tt> view and in some cases, it
 * can not.
 *
 * ToDo: Finish documentatation!
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template<class Ordinal, class Scalar>
class TpetraLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** \name Pure virtual functions that must be overridden in subclasses. */
  //@{

  /** \brief Return a smart pointer to a non-<tt>const</tt>
   * <tt>Tpetra::Operator</tt> view of this object.
   *
   * \param  tpetraOp  [out] The non-<tt>const</tt> tpetra operator view of <tt>*this</tt>.
   * \param  tpetraOpAdjointSupport
   *                  [out] Determines if the operator supports transposes or not.
   *
   * <b>Preconditions:</b></ul>
   * <li><tt>tpetraOp!=NULL</tt>
   * <li><tt>tpetraOpAdjointSupport!=NULL</tt>
   * </ul>
   *
   * <b>Posconditions:</b></ul>
   * <li><tt>tpetraOp->get() != NULL</tt>
   * </ul>
   *
   * The object accessed from <tt>*tpetraOp</tt> is only guaranteed to be
   * valid while the returned <tt>Teuchos::RCP</tt> object exits.
   * This allows for some very specialized implementations where a
   * <tt>Tpetra::Operator</tt> view of <tt>*this</tt> can be acquired and
   * released according to the lifetime of the returned
   * <tt>Teuchos::RCP</tt> object.
   *
   * The <tt>Tpetra::Operator</tt> object may be dynamic casted to more
   * specialized interfaces and therefore modified.  Then, when the last
   * <tt>RCP</tt> object ancestor returned from this function goes
   * away, then <tt>*this</tt> will be updated to relect the change.
   *
   * <b>Warning!</b> The client can not assume that the view in
   * <tt>*(*tpetraOp)</tt> will be valid past the lifetime of <tt>*this</tt>
   * object which is providing the view!  The client must take special care in
   * the case!
   */
  virtual void getTpetraOpView(
    Teuchos::RCP<Tpetra::Operator<Ordinal,Scalar> >   *tpetraOp
    ,EAdjointTpetraOp                                         *tpetraOpAdjointSupport
    ) = 0;

  /** \brief Return a smart pointer to a <tt>const</tt>
   * <tt>Tpetra::Operator</tt> view of this object and how the object is
   * applied to implement <tt>*this</tt> linear operator.
   *
   * \param  tpetraOp  [out] The <tt>const</tt> tpetra operator view of <tt>*this</tt>.
   * \param  tpetraOpAdjointSupport
   *                  [out] Determines if the operator supports transposes or not.
   *
   * <b>Preconditions:</b></ul>
   * <li><tt>tpetraOp!=NULL</tt>
   * <li><tt>tpetraOpAdjointSupport!=NULL</tt>
   * </ul>
   *
   * <b>Posconditions:</b></ul>
   * <li><tt>tpetraOp->get() != NULL</tt>
   * </ul>
   *
   * The object accessed from <tt>*return</tt> is only guaranteed to be valid
   * while the returned <tt>Teuchos::RCP</tt> object exits.  This
   * allows for some very specialized implementations where a
   * <tt>Tpetra::Operator</tt> view of <tt>*this</tt> can be acquired and
   * released according to the lifetime of the returned
   * <tt>Teuchos::RCP</tt> object.
   *
   * Note that if the client tries to constant cast the returned object and
   * modify it that this returned view is not guaranteed to update
   * <tt>*this</tt>.  If the goal is to modify <tt>*this</tt> then the client
   * should call the non-<tt>const</tt> version of this function.
   *
   * <b>Warning!</b> The client can not assume that the view in
   * <tt>*(*tpetraOp)</tt> will be valid past the lifetime of <tt>*this</tt>
   * object which is providing the view!  The client must take special care in
   * the case!
   */
  virtual void getTpetraOpView(
    Teuchos::RCP<const Tpetra::Operator<Ordinal,Scalar> >   *tpetraOp
    ,EAdjointTpetraOp                                               *tpetraOpAdjointSupport
    ) const = 0;

  //@}

};

}	// end namespace Thyra

#endif	// THYRA_TPETRA_LINEAR_OP_BASE_HPP
