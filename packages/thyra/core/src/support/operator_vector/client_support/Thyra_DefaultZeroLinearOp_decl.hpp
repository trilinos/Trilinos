// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP

#include "Thyra_ZeroLinearOpBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Represents a zero linear operator <tt>M = 0</tt>.
 *
 * This class implements:

 \verbatim

 y = alpha*op(M)*x + beta*y

 =>

 y = beta*y

 \endverbatim

 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultZeroLinearOp 
  : virtual public ZeroLinearOpBase<Scalar>
  , virtual public RowStatLinearOpBase<Scalar>
  , virtual public ScaledLinearOpBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  DefaultZeroLinearOp();

  /** Calls <tt>initialize()</tt>.
   */
  DefaultZeroLinearOp(
    const RCP<const VectorSpaceBase<Scalar> > &range,
    const RCP<const VectorSpaceBase<Scalar> > &domain
    );

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
  void initialize(
    const RCP<const VectorSpaceBase<Scalar> > &range,
    const RCP<const VectorSpaceBase<Scalar> > &domain
    );

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
                                                
  /** \brief Prints just the name <tt>DefaultZeroLinearOp</tt> along with the
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

  /** @name Overridden from RowStatLinearOpBase */
  //@{

  /** \brief . */
  virtual bool rowStatIsSupportedImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat) const;
 	
  /** \brief . */
  virtual void getRowStatImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat, 
    const Teuchos::Ptr<VectorBase< Scalar> > &rowStatVec) const;

  //@}

  /** @name Overridden from ScaledLinearOpBase */
  //@{

  /** \brief . */
  virtual bool supportsScaleLeftImpl() const
  { return true; }

  /** \brief . */
  virtual bool supportsScaleRightImpl() const
  { return true; }
 
  /** \brief . */ // Meaningless operation
  virtual void scaleLeftImpl(const VectorBase< Scalar > &/* row_scaling */)
  { }

  /** \brief . */ // Meaningless operation
  virtual void scaleRightImpl(const VectorBase< Scalar > &/* col_scaling */)
  { }

  //@}

private:

  RCP<const VectorSpaceBase<Scalar> > range_;
  RCP<const VectorSpaceBase<Scalar> > domain_;

  // Not defined and not to be called
  DefaultZeroLinearOp(const DefaultZeroLinearOp&);
  DefaultZeroLinearOp& operator=(const DefaultZeroLinearOp&);

};


/** \brief Create a zero linear operator with given range and domain spaces.
 *
 * \relates DefaultZeroLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
zero(
  const RCP<const VectorSpaceBase<Scalar> > &range,
  const RCP<const VectorSpaceBase<Scalar> > &domain
  );

/** \brief Create a nonconst zero linear operator with given range and domain spaces.
 *
 * This is to enable support for using the ScaledLinearOp interface. Which does nothing
 * yet still requires nonconstant operators.
 *
 * \relates DefaultZeroLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
nonconstZero(
  const RCP<const VectorSpaceBase<Scalar> > &range,
  const RCP<const VectorSpaceBase<Scalar> > &domain
  );


}	// end namespace Thyra


#endif	// THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP
