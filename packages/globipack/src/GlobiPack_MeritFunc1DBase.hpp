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

#ifndef GLOBIPACK_MERIT_FUNC_1D_BASE_HPP
#define GLOBIPACK_MERIT_FUNC_1D_BASE_HPP


#include "GlobiPack_Types.hpp"
#include "Teuchos_Describable.hpp"


namespace GlobiPack {


/** \brief Base class for 1D merit fucntions used in globalization methods.
 *
 * NOTE: The Scalar type must be a real type since comparions are performed
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class MeritFunc1DBase : virtual public Teuchos::Describable
{
public:

  /** \brief Determine if derivative evaluations are supported or not. */
  virtual bool supportsBaseDeriv() const = 0;

  /** \brief Determine if derivative evaluations are supported or not. */
  virtual bool supportsDerivEvals() const = 0;

  /** \brief Return the base derivative at <tt>alpha=0</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>supportsBaseDeriv() == true</tt>
   * </ul>
   */
  virtual Scalar baseDeriv() const = 0;

  /** \brief Evaluate the merit function at <tt>alpha</tt>.
   *
   * \param alpha [in] The value of the independent variable determining the
   * step length.  Typically <tt>alpha > 0</tt>.
   *
   * \param phi [out] The value of the merit function evaluated at
   * <tt>alpha</tt>.
   *
   * \param Dphi [out] The value of the derivative of the merit function
   * evaluated at <tt>alpha</tt>.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li> If <tt>!is_null(Dphi)</tt> then <tt>supportsDerivEvals()</tt> must
   * be <tt>true</tt>.
   *
   * </ul>
   */
  virtual void eval( const Scalar &alpha, const Ptr<Scalar> &phi,
    const Ptr<Scalar> &Dphi ) const = 0;

};


/** \brief Compute the value of the merit function <tt>phi(alpha)</tt>.
 *
 * \relates MeritFunc1DBase
 */
template<typename Scalar>
typename ScalarTraits<Scalar>::magnitudeType
computeValue(const MeritFunc1DBase<Scalar> &phi, const Scalar &alpha)
{
  Scalar phi_val = ScalarTraits<Scalar>::zero();
  phi.eval(alpha, Teuchos::outArg(phi_val), Teuchos::null);
  return phi_val;
}


/** \brief Compute a point as an object.
 *
 * \relates MeritFunc1DBase
 */
template<typename Scalar>
PointEval1D<Scalar>
computePoint(const MeritFunc1DBase<Scalar> &phi, const Scalar &alpha,
  const bool compute_phi = true, const bool compute_Dphi = false)
{
  using Teuchos::null;
  using Teuchos::outArg;
  PointEval1D<Scalar> p;
  p.alpha = alpha;
  phi.eval( alpha, compute_phi ? outArg(p.phi) : null ,
    compute_Dphi ? outArg(p.Dphi) : null );
  return p;
}


} // namespace GlobiPack


#endif // GLOBIPACK_MERIT_FUNC_1D_BASE_HPP
