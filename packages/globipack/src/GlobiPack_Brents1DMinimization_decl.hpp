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

#ifndef GLOBIPACK_BRENTS_1D_MINIMIZATION_DECL_HPP
#define GLOBIPACK_BRENTS_1D_MINIMIZATION_DECL_HPP


#include "GlobiPack_MeritFunc1DBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace GlobiPack {


/** \brief Simple concrete class that implements a 1D algorithm to mimimize a
 * 1D function.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class Brents1DMinimization
  : public Teuchos::Describable,
    public Teuchos::VerboseObject<Brents1DMinimization<Scalar> >,
    public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \name Constructor/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters. */
  Brents1DMinimization();

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Bracket. */
  //@{

  /** \brief Approximatly mimimize a 1D function.
   *
   * \param phi [in] The evaluator object for the merit function which
   * evaluates <tt>phi(alpha)</tt>).
   *
   * \param pointLower [in] Gives the lower bound for the point.  This lower
   * bound will be respected and will never be violated.  The derivative field
   * <tt>pointLower->Dphi</tt> is not accessed.
   *
   * \param pointMiddle [in/out] In input, <tt>*pointMiddle</tt> give the
   * initial guess for the point.  On output, <tt>*pointMiddle<tt> gives the
   * approximate minumum solution.  The derivative field
   * <tt>pointUpper->Dphi</tt> is not accessed.
   *
   * \param pointUpper [in] Gives the upper bound for the point.  The
   * derivative field <tt>pointUpper->Dphi</tt> is not accessed.
   *
   * \param numIters [out] If not null, on output, <tt>numIters<tt> gives the
   * number of iterations used.
   *
   * \return Returns <tt>true</tt> if an approximate local minimum has been
   * found, <tt>false</tt> otherwise.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>pointLower.alpha < pointMiddle->alpha && pointMiddle->alpha <
   * pointUpper.alpha</tt>
   *
   * <li><tt>pointLower.phi > pointMiddle->phi && pointMiddle->phi <
   * pointUpper.phi</tt>
   *
   * </ul>
   */
  bool approxMinimize(
    const MeritFunc1DBase<Scalar> &phi,
    const PointEval1D<Scalar> &pointLower,
    const Ptr<PointEval1D<Scalar> > &pointMiddle,
    const PointEval1D<Scalar> &pointUpper,
    const Ptr<int> &numIters = Teuchos::null
    ) const;

  //@}

private:

  // //////////////////////
  // Private data members

  Scalar rel_tol_;
  Scalar bracket_tol_;
  int max_iters_;

};


/** \brief Nonmember constructor.
 *
 * \relates Brents1DMinimization
 */
template<typename Scalar>
inline
const RCP<Brents1DMinimization<Scalar> > brents1DMinimization()
{
  return Teuchos::rcp(new Brents1DMinimization<Scalar>());
}


// Default values are exposed here for unit testing purposes


namespace Brents1DMinimizationUtils {


const std::string rel_tol_name = "Relative Tol";
const double rel_tol_default = 1.0e-5;

const std::string bracket_tol_name = "Bracket Tol";
const double bracket_tol_default = 1.0e-5;

const std::string max_iters_name = "Max Iterations";
const int max_iters_default = 10;


} // namespace Brents1DMinimizationUtils



} // namespace GlobiPack


#endif // GLOBIPACK_BRENTS_1D_MINIMIZATION_DECL_HPP
