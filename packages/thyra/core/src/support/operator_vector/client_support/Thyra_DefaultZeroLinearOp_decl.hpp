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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP

#include "Thyra_ZeroLinearOpBase.hpp"
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
class DefaultZeroLinearOp : virtual public ZeroLinearOpBase<Scalar>
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


}	// end namespace Thyra


#endif	// THYRA_DEFAULT_ZERO_LINEAR_OP_DECL_HPP
