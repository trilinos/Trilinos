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
