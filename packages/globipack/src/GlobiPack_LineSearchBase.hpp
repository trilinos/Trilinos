/*
// @HEADER
// ***********************************************************************
// 
//    GlobiPack: Collection of Scalar 1D globalizaton utilities
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef GLOBIPACK_LINE_SEARCH_BASE_HPP
#define GLOBIPACK_LINE_SEARCH_BASE_HPP


#include "GlobiPack_MeritFunc1DBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"


namespace GlobiPack {


/** \brief Base class for 1D linearsearch algorithms.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class LineSearchBase
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<LineSearchBase<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief Determines if the linesearch algorithm requires the base
   * derivative at <tt>Dphi(0)</tt> or not.
   */
  virtual bool requiresBaseDeriv() const = 0;

  /** \brief Determines if the linesearch algorithm requires that
   * <tt>Dphi(alpha)</tt> can be computed or not.
   */
  virtual bool requiresDerivEvals() const = 0;

  /** \brief Called to perform a linesearch.
   *
   * \param phi [in] The merit function object that will compute the initial
   * descent direction <tt>Dphi</tt>, and the value <tt>phi(alpha)</tt> and/or
   * derivative <tt>Dphi(alpha) at different points <tt>alpha</tt>.  The last
   * call to <tt>phi.eval(...)</tt> will always be with the value of
   * <tt>alpha_k</tt> returned.
   *
   * \param phi_k [in] The value of <tt>phi(0)</tt>.  This must be computed by
   * the client externally and passed in.
   *
   * \param alpha_k [in/out] On input, <tt>alpha_k</tt> is the initial value
   * to try out (usually 1.0 for most Newton-based algorithms).  On output,
   * <tt>alpha_k</tt> is the accepted value for a successful line search, or
   * it will be the <tt>alpha</tt> for the minimum <tt>phi(alpha)</tt> found
   * during a failed line search algorithm.  The initial value of
   * <tt>alpha_k</tt> on input will typically not be exceeded internally. but
   * that is not guaranteed.
   *
   * \param phi_kp1 [in/out] Merit function at new point
   * <tt>phi(alpha_k)</tt>.  On input it must be the value of
   * <tt>phi(alpha_k)</tt> for <tt>alpha_k</tt> on input.  On output it is set
   * to <tt>phi(alpha_k)</tt>.
   *
   * \param Dphi_kp1 [in/out] Derivative Merit function at new point
   * <tt>phi(alpha_k)</tt>.  On input it can be set the the value of
   * <tt>Dphi(alpha_k)</tt> for <tt>alpha_k</tt> on input.  On output, it is
   * set to <tt>phi.value(alpha_k)</tt>.  Must be <tt>null</tt> if
   * <tt>this->requiresDerivEvals()==false</tt>.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li> [<tt>this->requiresBaseDeriv()==true</tt>]
   * <tt>phi.supportsBaseDeriv()==true</tt>
   *
   * <li> [<tt>this->requiresDerivEvals()==true</tt>]
   * <tt>phi.supportsDerivEvals()==true</tt>
   *
   * <li> [<tt>this->requiresBaseDeriv()</tt>==true] <tt>phi.baseDeriv() <
   * 0.0</tt> (throw <tt>Exceptions::NotDescentDirection</tt>)
   *
   * <li> <tt>!is_null(alpha_k)</tt>
   *
   * <li> <tt>!is_null(phi_kp1)</tt>
   *
   * <li> [<tt>this->requiresDerivEvals()==true</tt>]
   * <tt>!is_null(Dphi_kp1)==true</tt>
   *
   * <li> [<tt>this->requiresDerivEvals()==false</tt>]
   * <tt>is_null(Dphi_kp1)==true</tt>
   *
   * </ul>
   *
   * \returns <tt>true</tt> for successful line search or <tt>false</tt> for a
   * line search failure.
   *
   * This function computes the approximate minimum to 1D merit function
   * <tt>phi(alpha)</tt>.  More specifically the following problem is
   * approximately solved:

   \verbatim
     min  phi(alpha)  s.t. alpha = [0, alpha_upper]<br>
   \endverbatim

   * Actually, if the initial <tt>alpha_k</tt> satisfies an internally defied
   * descent requirement, then it will be choosen over smaller values of
   * <tt>alpha_k</tt> that may result in a greater reduction in the given
   * merit function..
   *
   * If the maximum number of iterations is exceeded without finding an
   * acceptable point, then the subclass object will return <tt>false</tt> and
   * will return values of <tt>alpha_k</tt> and <tt>phi_kp1</tt> will be for
   * the lowest value of <tt>phi_kp1 = phi(alpha_k)</tt> found.  In this case,
   * the last call to <tt>phi(alpha_k)</tt> will be this best value of
   * <tt>phi_kp1</tt>.
   */
  virtual bool doLineSearch(
    const MeritFunc1DBase<Scalar> &phi,
    const Scalar &phi_k,
    const Ptr<Scalar> &alpha_k,
    const Ptr<Scalar> &phi_kp1,
    const Ptr<Scalar> &Dphi_kp1
    ) const = 0;

};


} // namespace GlobiPack


#endif // GLOBIPACK_LINE_SEARCH_BASE_HPP
