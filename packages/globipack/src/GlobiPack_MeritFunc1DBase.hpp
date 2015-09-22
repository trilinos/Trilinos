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
  virtual bool supportsDerivEvals() const = 0;

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
