/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
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

#ifndef OPTIPACK_NONLINEAR_CG_DECL_HPP
#define OPTIPACK_NONLINEAR_CG_DECL_HPP


#include "OptiPack_Types.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "GlobiPack_LineSearchBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace OptiPack {


namespace NonlinearCGUtils {


/** \brief . */
enum ESolveReturn {
  SOLVE_SOLUTION_FOUND,      ///< .
  SOLVE_LINSEARCH_FAILURE,   ///< .
  SOLVE_MAX_ITERS_EXCEEDED   ///< .
};


} // namespace NonlinearCGUtils


/** \brief Concrete class implementing several nonlinear CG algorithms.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class NonlinearCG
  : public Teuchos::Describable,
    public Teuchos::VerboseObject<NonlinearCG<Scalar> >,
    protected Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructor/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters. */
  NonlinearCG();

  /** \brief Initialize. */
  void initialize(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const int paramIndex,
    const int responseIndex,
    const RCP<GlobiPack::LineSearchBase<Scalar> > &linesearch
    );

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Solve. */
  //@{

  /** \brief Perform a solve.
   *
   * \param p [in/out] On input <tt>p</tt> is the initial guess for the
   * solution.  On output, will be the final estimate for the solution.
   *
   * \param g_opt [out] On output, <tt>*g_opt</tt> will be set to the final
   * value of the objective function.
   *
   * \param tol [in] If <tt>!is_null(tol)</tt>, then <tt>*tol</tt> will be the
   * tolerance used to determine the convergence of the algorithm by comparing
   * to <tt>norm(g_grad)</tt> (where <tt>norm(...)</tt> is the natural norm
   * defined by the vector spaces scalar product).  If <tt>is_null(tol)</tt>,
   * then the tolerance will be determined in some other way.
   *
   * \param alpha_init [in] If <tt>!is_null(alpha_init)</tt>, then
   * <tt>*alpha_init</tt> will be the initial line search step length on the
   * very first nonlinear CG iteration.  Of <tt>is_null(alpha_init)</tt>, the
   * initial step length will be determined automatically.
   *
   * \return Returns <tt>true</tt> if the solution was found.  Returns
   * <tt>false</tt> if a line search failure is encountered.
   *
   */
  NonlinearCGUtils::ESolveReturn
  doSolve(
    const Ptr<Thyra::VectorBase<Scalar> > &p,
    const Ptr<ScalarMag> &g_opt,
    const Ptr<const ScalarMag> &tol = Teuchos::null,
    const Ptr<const ScalarMag> &alpha_init = Teuchos::null
    );

  //@}

private:

  // //////////////////////
  // Private data members

  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  int paramIndex_;
  int responseIndex_;
  RCP<GlobiPack::LineSearchBase<Scalar> > linesearch_;

  mutable int numIters_;

};


/** \brief Nonmember constructor.
 *
 * \relates NonlinearCG
 */
template<typename Scalar>
const RCP<NonlinearCG<Scalar> >
nonlinearCG(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const int paramIndex,
  const int responseIndex,
  const RCP<GlobiPack::LineSearchBase<Scalar> > &linesearch
  )
{
  const RCP<NonlinearCG<Scalar> > solver = 
    Teuchos::rcp(new NonlinearCG<Scalar>);
  solver->initialize(model, paramIndex, responseIndex, linesearch);
  return solver;
}


// Default values are exposed here for unit testing purposes


namespace NonlinearCGUtils {


const std::string eta_name = "Armijo Slope Fraction";
const double eta_default = 1.0e-4;


} // namespace NonlinearCGUtils



} // namespace OptiPack


#endif // OPTIPACK_NONLINEAR_CG_DECL_HPP
