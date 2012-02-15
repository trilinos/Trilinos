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

#ifndef GLOBIPACK_POLY_INTERP_LINE_SEARCH_DECL_HPP
#define GLOBIPACK_POLY_INTERP_LINE_SEARCH_DECL_HPP


#include "GlobiPack_LineSearchBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace GlobiPack {


/** \brief Linesearch subclass implementing a backtracking-only line search
 * using an Armijo cord test condition and a quadratic interploation.
 *
 * This linesearch class is really designed for (quasi) Newton methods where a
 * backtracking only linesearch is the only thing the makes sense.
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
    const PointEval1D<Scalar> &point_k,
    const Ptr<PointEval1D<Scalar> > &point_kp1,
    const Ptr<int> &numIters
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
