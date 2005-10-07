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

#ifndef THYRA_SERIAL_MULTI_VECTOR_STD_DECL_HPP
#define THYRA_SERIAL_MULTI_VECTOR_STD_DECL_HPP

#include "Thyra_SerialMultiVectorBaseDecl.hpp"

namespace Thyra {

/** \brief General, yet efficient, concrete <tt>MultiVectorBase</tt>
 * implementation subclass for serial shared-memory multi-vectors.
 *
 * Objects of this type generally should not be constructed directly by a
 * client but instead by using the concrete vector space subclass
 * <tt>SerialVectorSpaceStd</tt> using the function
 * <tt>Thyra::SerialVectorSpaceStd::createMembers()</tt>.
 *
 * The storage type can be anything since a <tt>Teuchos::RefCountPtr</tt> is
 * used to pass in the values pointer into the constructor and
 * <tt>initialize()</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_serial_concrete_std_grp
 */
template<class Scalar>
class SerialMultiVectorStd : virtual public SerialMultiVectorBase<Scalar> {
public:

  /** \brief . */
  using SerialMultiVectorBase<Scalar>::subView;
  /** \brief . */
  using SerialMultiVectorBase<Scalar>::col;

  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct to uninitialized
  SerialMultiVectorStd();

  /// Calls <tt>initialize()</tt>
  SerialMultiVectorStd(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    );

  /// Calls <tt>initialize()</tt>
  SerialMultiVectorStd(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    ,const Teuchos::RefCountPtr<Scalar>                                    &values
    ,const Index                                                           leadingDim
    );

  /** \brief Initialize.
   *
   * @param  range     [in] Smart pointer to the vector space object
   *                   that defines the range.
   * @param  domain    [in] Smart pointer to the vector space object
   *                   that defines the domain.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.get()!=NULL</tt>
   * <li><tt>domain.get()!=NULL</tt>
   * </ul>
   *
   * This function simply calls <tt>initialize(range,domain,...)</tt>
   * and passes in dynamically allocated data for the values.
   */
  void initialize(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    );

  /** \brief Initialize.
   *
   * @param  range     [in] Smart pointer to the vector space object
   *                   that defines the range.
   * @param  domain    [in] Smart pointer to the vector space object
   *                   that defines the domain.
   * @param  values    [in] Smart pointer to beginning of Fortran-style column-major
   *                   array that defines the local values in the multi-vector.
   *                   This array must be at least of dimension <tt>leadingDim*domain->dim()</tt>
   *                   and <tt>(&*values)[ (i-1) + (j-1)*leadingDim ]</tt> gives the local value
   *                   of the one-based <tt>(i,j)</tt> entry where <tt>i=1...range()->dim()</tt>
   *                   and <tt>j=1...domain->dim()</tt>.
   * @param  leadingDim
   *                   [in] The leading dimension of the multi-vector.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.get()!=NULL</tt>
   * <li><tt>domain.get()!=NULL</tt>
   * <li><tt>values.get()!=NULL</tt>
   * <li><tt>leadingDim >= range->dim()</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->rangeScalarProdVecSpc().get() == range.get()</tt>
   * <li><tt>this->domainScalarProdVecSpc().get() == domain.get()</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    ,const Teuchos::RefCountPtr<Scalar>                                    &values
    ,const Index                                                           leadingDim
    );

  /** \brief Set to an uninitialized state.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->rangeScalarProdVecSpc().get() == NULL</tt>
   * <li><tt>this->domainScalarProdVecSpc().get() == NULL</tt>
   * </ul>
   */
  void uninitialize(
    Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >         *range         = NULL
    ,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >        *domain        = NULL
    ,Teuchos::RefCountPtr<Scalar>                                          *values        = NULL
    ,Index                                                                 *leadingDim    = NULL
    );

  //@}

  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Overridden from EuclideanLinearOpBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;
  /** \brief . */
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
  //@}

  /** @name Overridden from MultiVectorBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr<VectorBase<Scalar> > col(Index j);
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const Range1D& col_rng );
  //@}

  /** @name Overridden from SerialMultiVectorBase */
  //@{

  /** \brief . */
  void getData( const Scalar **values, Index *leadingDim ) const;
  /** \brief . */
  void freeData( const Scalar *values ) const;
  /** \brief . */
  void getData( Scalar **values, Index *leadingDim );
  /** \brief . */
  void commitData( Scalar *values );
  //@}
  
private:
  
  // ///////////////////////////////////////
  // Private data members

  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >     range_;
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >     domain_;
  Teuchos::RefCountPtr<Scalar>                                       values_;
  Index                                                              leadingDim_;
  Index                                                              numRows_; // Cached
  Index                                                              numCols_; // Cached
  
}; // end class SerialMultiVectorStd

} // end namespace Thyra

#endif // THYRA_SERIAL_MULTI_VECTOR_STD_DECL_HPP
