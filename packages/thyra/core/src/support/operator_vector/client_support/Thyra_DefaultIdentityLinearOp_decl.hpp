// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_IDENTITY_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_IDENTITY_LINEAR_OP_DECL_HPP

#include "Thyra_IdentityLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief Represents a identity linear operator <tt>M = I</tt>.
 *
 * This class implements:

 \verbatim

 y = alpha*op(M)*x + beta*y

 =>

 y = alpha*x + beta*y

 \endverbatim

 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultIdentityLinearOp : virtual public IdentityLinearOpBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  DefaultIdentityLinearOp();

  /** Calls <tt>initialize()</tt>.
   */
  DefaultIdentityLinearOp(const RCP<const VectorSpaceBase<Scalar> > &space);

  /** \brief Initialize given a list of non-const linear operators.
   *
   * \param range [in] Range vector space.
   *
   * \param range [in] Domain vector space.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.get()!=NULL</tt>
   * <li><tt>domain.get()!=NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range().get()==range.get()</tt>
   * <li><tt>this->domain().get()==domain.get()</tt>
   * </ul>
   */
  void initialize(const RCP<const VectorSpaceBase<Scalar> > &space);

  /** \brief Set to uninitialized.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  void uninitialize();

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  
  /** \brief Returns <tt>Teuchos::null</tt> if uninitialized. */
  RCP< const VectorSpaceBase<Scalar> > range() const;
  
  /** \brief Returns <tt>Teuchos::null</tt> if uninitialized. */
  RCP< const VectorSpaceBase<Scalar> > domain() const;
  
  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;
  
  //@}
  
  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Prints just the name <tt>DefaultIdentityLinearOp</tt> along with the
   * overall dimensions.
   */
  std::string description() const;

  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> . */
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

private:

  RCP<const VectorSpaceBase<Scalar> >  space_;

  // Not defined and not to be called
  DefaultIdentityLinearOp(const DefaultIdentityLinearOp&);
  DefaultIdentityLinearOp& operator=(const DefaultIdentityLinearOp&);

};


/** \brief Create an identity linear operator with given a vector space.
 *
 * \relates DefaultIdentityLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
identity(
  const RCP<const VectorSpaceBase<Scalar> > &space,
  const std::string &label = ""
  );


}	// end namespace Thyra


#endif	// THYRA_DEFAULT_IDENTITY_LINEAR_OP_DECL_HPP
