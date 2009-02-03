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

#ifndef GLOBIPACK_POLY_INTERP_LINE_SEARCH_DECL_HPP
#define GLOBIPACK_POLY_INTERP_LINE_SEARCH_DECL_HPP


#include "GlobiPack_LineSearchBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace GlobiPack {


/** \brief Linesearch subclass implementing a bracketing line search using an
 * Armijo cord test condition and a quadratic interploation.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class ArmijoPolyInterpLineSearch
  : public LineSearchBase<Scalar>,
    protected Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \name Constructor/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters. */
  ArmijoPolyInterpLineSearch();

  /** \brief Return the number of iterations after a line search. */
  int numIterations() const;

  /** \brief . */
  Scalar eta() const;
  /** \brief . */
  Scalar minFrac() const;
  /** \brief . */
  Scalar maxFrac() const;
  /** \brief . */
  int	minIters() const;
  /** \brief . */
  int	maxIters() const;
  /** \brief . */
  bool doMaxIters() const;

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overrridden from LineSearchBase. */
  //@{

  /** \brief Returns true. */
  virtual bool requiresBaseDeriv() const;

  /** \brief Returns false. */
  virtual bool requiresDerivEvals() const;

  /** \brief . */
  virtual bool doLineSearch(
    const MeritFunc1DBase<Scalar> &phi,
    const Scalar &phi_k,
    const Ptr<Scalar> &alpha_k,
    const Ptr<Scalar> &phi_kp1,
    const Ptr<Scalar> &Dphi_kp1
    ) const;

  //@}

private:

  // //////////////////////
  // Private data members

  Scalar eta_;
  Scalar minFrac_;
  Scalar maxFrac_;
  int	minIters_;
  int	maxIters_;
  bool doMaxIters_;

  mutable int numIters_;

};


/** \brief Nonmember constructor.
 *
 * \relates ArmijoPolyInterpLineSearch
 */
template<typename Scalar>
const RCP<ArmijoPolyInterpLineSearch<Scalar> > armijoQuadraticLineSearch()
{
  return Teuchos::rcp(new ArmijoPolyInterpLineSearch<Scalar>());
}


// Default values are exposed here for unit testing purposes


namespace ArmijoPolyInterpLineSearchUtils {


const std::string eta_name = "Armijo Slope Fraction";
const double eta_default = 1.0e-4;

const std::string minFrac_name = "Min Backtrack Fraction";
const double minFrac_default = 0.1;

const std::string maxFrac_name = "Max Backtrack Fraction";
const double maxFrac_default = 0.5;

const std::string minIters_name = "Min Num Iterations";
const int minIters_default = 0;

const std::string maxIters_name = "Max Num Iterations";
const int maxIters_default = 20;

const std::string doMaxIters_name = "Do Max Iterations";
const bool doMaxIters_default = false;


} // namespace ArmijoPolyInterpLineSearchUtils



} // namespace GlobiPack


#endif // GLOBIPACK_POLY_INTERP_LINE_SEARCH_DECL_HPP
