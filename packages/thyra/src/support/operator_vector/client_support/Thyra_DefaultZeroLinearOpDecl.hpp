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

#ifndef THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP

#include "Thyra_ZeroLinearOpBase.hpp"
#include "Thyra_SingleScalarLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_arrayArg.hpp"

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
  : virtual public ZeroLinearOpBase<Scalar>              // Public interface
  , virtual protected SingleScalarLinearOpBase<Scalar>   // Implementation detail
{
public:

  /** \brief . */
  using SingleScalarLinearOpBase<Scalar>::apply;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  DefaultZeroLinearOp();

  /** Calls <tt>initialize()</tt>.
   */
  DefaultZeroLinearOp(
    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &domain
    );

  /** \brief Initialize given a list of non-const linear operators.
   *
   * \param  range  [in] Range vector space.
   * \param  range  [in] Domain vector space.
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
    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &domain
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
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  
  /** \brief Returns <tt>Teuchos::null</tt> if uninitialized. */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;
  
  //@}
  
  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Prints just the name <tt>DefaultZeroLinearOp</tt> along with the
   * overall dimensions.
   */
  std::string description() const;

  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> . */
  bool opSupported(ETransp M_trans) const;

  /** \brief . */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

private:

  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  range_;
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  domain_;

  // Not defined and not to be called
  DefaultZeroLinearOp(const DefaultZeroLinearOp&);
  DefaultZeroLinearOp& operator=(const DefaultZeroLinearOp&);

};

/** \brief Create a zero linear operator with given range and domain spaces.
 *
 * \relates DefaultZeroLinearOp
 */
template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
zero(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &domain
  );

}	// end namespace Thyra

#endif	// THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP
