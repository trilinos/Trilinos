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

#ifndef THYRA_PRECONDITIONER_BASE_HPP
#define THYRA_PRECONDITIONER_BASE_HPP

#include "Thyra_OperatorSolveTypes.hpp"
#include "Teuchos_Describable.hpp"


namespace Thyra {


/** \brief Simple interface class to access a precreated preconditioner as one
 * or more linear operators objects and information on how they are meant to
 * be applied.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 *
 * This class supports four different kinds of preconditioners:
 * <ul>
 * <li>Single preconditioner linear operator designed or targeted to be applied on the left:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()!=NULL</tt>
 *   <li><tt>getRightPrecOp().get()==NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()==NULL</tt>
 *   </ul>
 * <li>Single preconditioner linear operator designed or targeted to be applied on the right:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()==NULL</tt>
 *   <li><tt>getRightPrecOp().get()!=NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()==NULL</tt>
 *   </ul>
 * <li>Split two-sided preconditioner with linear operators designed or
 *   targeted to be applied on the left and the right:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()!=NULL</tt>
 *   <li><tt>getRightPrecOp().get()!=NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()==NULL</tt>
 *   </ul>
 * <li>Single preconditioner linear operator not designed or targeted to be
 *   applied on the left or the right:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()==NULL</tt>
 *   <li><tt>getRightPrecOp().get()==NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()!=NULL</tt>
 *   </ul>
 * </ul>
 */
template<class Scalar>
class PreconditionerBase : virtual public Teuchos::Describable
{
public:

  /** @name Pure virtual public functions that must be overridden in subclasses */
  //@{

  /** \brief Return if the underlying left preconditioner operator is
   * const-only or allows non-const access.
   */
  virtual bool isLeftPrecOpConst() const = 0;

  /** \brief Return a non-const left preconditioner linear operator if one is
   * designed or targeted to be applied on the left.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isLeftPrecOpConst()==true</tt>] <tt>getLeftPrecOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstLeftPrecOp() = 0;

  /** \brief Return a const left preconditioner linear operator if one is
   * designed or targeted to be applied on the left.
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getLeftPrecOp()const = 0;

  /** \brief Return if the underlying right preconditioner operator is
   * const-only or allows non-const access.
   */
  virtual bool isRightPrecOpConst() const = 0;

  /** \brief Return a non-const right preconditioner linear operator if one is
   * designed or targeted to be applied on the right.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isRightPrecOpConst()==true</tt>] <tt>getRightPrecOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstRightPrecOp() = 0;

  /** \brief Return a const right preconditioner linear operator if one is
   * designed or targeted to be applied on the right.
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getRightPrecOp() const = 0;

  /** \brief Return if the underlying unspecified preconditioner operator is
   * const-only or allows non-const access.
   */
  virtual bool isUnspecifiedPrecOpConst() const = 0;

  /** \brief Return a non-const generic preconditioner linear operator that is
   * not designed or targeted to be applied on the left or on the right.
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstUnspecifiedPrecOp() = 0;

  /** \brief Return a const generic preconditioner linear operator that is not
   * designed or targeted to be applied on the left or on the right.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isUnspecifiedPrecOpConst()==true</tt>] <tt>getUnspecifiedPrecOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getUnspecifiedPrecOp() const = 0;
  
  //@}

private:
  
  // Not defined and not to be called
  PreconditionerBase<Scalar>&
  operator=(const PreconditionerBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_PRECONDITIONER_BASE_HPP
