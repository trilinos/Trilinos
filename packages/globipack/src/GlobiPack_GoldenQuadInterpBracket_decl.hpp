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

#ifndef GLOBIPACK_GOLDEN_BRACKET_QUAD_INTERP_DECL_HPP
#define GLOBIPACK_GOLDEN_BRACKET_QUAD_INTERP_DECL_HPP


#include "GlobiPack_MeritFunc1DBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace GlobiPack {


/** \brief Simple concrete class that implements a 1D algorithm to bracket the
 * minimum of a 1D merit function.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class GoldenQuadInterpBracket
  : public Teuchos::Describable,
    public Teuchos::VerboseObject<GoldenQuadInterpBracket<Scalar> >,
    public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \name Constructor/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters. */
  GoldenQuadInterpBracket();

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

  /** \brief Bracket the minimum of a 1D function.
   *
   * \param phi [in] The evaluator object for the merit function which
   * evaluates <tt>phi(alpha)</tt>).
   *
   * \param pointLower [in/out] In input, <tt>*pointLower</tt> give the
   * initial guess for the lower bound for the point.  This lower bound will
   * be respected and will never be violated.  On output, <tt>*pointLower<tt>
   * gives the lower bound for the bracket of the minimum.  The derivative
   * field <tt>pointLower->Dphi</tt> is not accessed.
   *
   * \param pointMiddle [in/out] In input, <tt>*pointMiddle</tt> give the
   * initial guess for the point.  On output, <tt>*pointMiddle<tt> gives the
   * bracketed minimum.  The derivative field <tt>pointUpper->Dphi</tt> is not
   * accessed.
   *
   * \param pointUpper [out] On output, <tt>*pointUpper<tt> gives the upper
   * bound for the bracket of the minimum.  The derivative field
   * <tt>pointUpper->Dphi</tt> is not accessed.
   *
   * \param numIters [out] If not null, on output, <tt>numIters<tt> gives the
   * number of iterations used.
   *
   * \return Returns <tt>true</tt> if a bracket has been found, <tt>false</tt>
   * otherwise.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>pointLower->alpha < pointMiddle->alpha</tt>
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>[<tt>returnVal==true</tt>] <tt>pointLower->alpha < pointMiddle->alpha
   * && pointMiddle->alpha < pointUpper->alpha</tt>
   *
   * <li>[<tt>returnVal==true</tt>] <tt>pointLower->phi > pointMiddle->phi
   * && pointMiddle->phi < pointUpper->phi</tt>
   *
   * </ul>
   */
  bool bracketMinimum(
    const MeritFunc1DBase<Scalar> &phi,
    const Ptr<PointEval1D<Scalar> > &pointLower,
    const Ptr<PointEval1D<Scalar> > &pointMiddle,
    const Ptr<PointEval1D<Scalar> > &pointUpper,
    const Ptr<int> &numIters = Teuchos::null
    ) const;

  //@}

private:

  // //////////////////////
  // Private data members

};


/** \brief Nonmember constructor.
 *
 * \relates GoldenQuadInterpBracket
 */
template<typename Scalar>
inline
const RCP<GoldenQuadInterpBracket<Scalar> > goldenQuadInterpBracket()
{
  return Teuchos::rcp(new GoldenQuadInterpBracket<Scalar>());
}


// Default values are exposed here for unit testing purposes


namespace GoldenQuadInterpBracketUtils {


const std::string eta_name = "Armijo Slope Fraction";
const double eta_default = 1.0e-4;


} // namespace GoldenQuadInterpBracketUtils



} // namespace GlobiPack


#endif // GLOBIPACK_GOLDEN_BRACKET_QUAD_INTERP_DECL_HPP
