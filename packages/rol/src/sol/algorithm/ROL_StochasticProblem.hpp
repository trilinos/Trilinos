// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_STOCHASTICPROBLEM_HPP
#define ROL_STOCHASTICPROBLEM_HPP

/** @ingroup stochastic_group
    \class ROL::StochasticProblem
    \brief Provides an interface for stochastic optimization problems.

    Let \f$\mathcal{X}_0\f$ denote a space of random variables defined on a
    probability space \f$(\Omega,\mathcal{F},\mathbb{P})\f$.
    ROL::StochasticProblem provides the interface for problems of the form
    \f[
       \begin{array}{l l}
         \text{minimize}   & \mathcal{R}_0(f_0(z))                   \\
         \text{subject to} & \ell \le z \le u                        \\
                           & g(z) \{\le,\, =\} 0                     \\
                           & h(z,\cdot) \{\le,\, =\} 0 \text{w.p.~1} \\
       \end{array}
    \f]
    where \f$\mathcal{R}_0:\mathcal{X}_0\to\mathbb{R}\f$ ``scalarizes'' the
    random variable objective function \f$f_0:Z \to \mathcal{X}_0\f$,
    \f$\ell\le u\f$ are vectors in \f$Z\f$, \f$g:Z \to C_g\f$ is a
    deterministic constraint and \f$h:Z\times\Omega C_h\f$ is a stochastic
    constraint required to hold with probability 1.

    To Do: Incorporate constraints of the form
    \f[
       \mathcal{R}_i(f_i(z)) \{\le,\, =\} 0
    \f]
    where \f$\mathcal{R}_i:\mathcal{X}_i\to\mathbb{R}\f$ ``scalarizes'' the
    random variable ``costs'' \f$f_i:Z \to\mathcal{X}_i\f$ for
    \f$i=1,\ldots,N\f$.  Here, \f$\mathcal{X}_i\f$ denote spaces of random
    variables defined on \f$(\Omega,\mathcal{F},\mathbb{P})\f$.
*/

#include "ROL_OptimizationProblem.hpp"
#include "ROL_SampleGenerator.hpp"

// Risk-Neutral Includes
#include "ROL_RiskNeutralObjective.hpp"

// Risk-Averse Includes
#include "ROL_RiskAverseObjective.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_RiskBoundConstraint.hpp"

// Risk-Less Includes
#include "ROL_RiskLessConstraint.hpp"

// Almost Sure Includes
#include "ROL_AlmostSureConstraint.hpp"
#include "ROL_SimulatedVector.hpp"

// BPOE Includes
#include "ROL_BPOEObjective.hpp"

#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class StochasticProblem : public OptimizationProblem<Real> {
private:
  Teuchos::RCP<Objective<Real> >                     INPUT_obj_;
  Teuchos::RCP<Vector<Real> >                        INPUT_sol_;
  Teuchos::RCP<BoundConstraint<Real> >               INPUT_bnd_;
  std::vector<Teuchos::RCP<Constraint<Real> > >      INPUT_econ_;
  std::vector<Teuchos::RCP<Vector<Real> > >          INPUT_emul_;
  std::vector<Teuchos::RCP<Constraint<Real> > >      INPUT_icon_;
  std::vector<Teuchos::RCP<Vector<Real> > >          INPUT_imul_;
  std::vector<Teuchos::RCP<BoundConstraint<Real> > > INPUT_ibnd_;

  std::vector<Teuchos::RCP<SampleGenerator<Real> > > esampler_;
  std::vector<Teuchos::RCP<SampleGenerator<Real> > > isampler_;
  Teuchos::RCP<SampleGenerator<Real> >               vsampler_;
  Teuchos::RCP<SampleGenerator<Real> >               gsampler_;
  Teuchos::RCP<SampleGenerator<Real> >               hsampler_;

  bool setVector_;
  bool setObjective_;
  bool setBound_;
  bool setMultiplier_;
  bool isCon_;

  std::vector<Real> computeSampleMean(Teuchos::RCP<SampleGenerator<Real> > &sampler) {
    // Compute mean value of inputs and set parameter in objective
    int dim = sampler->getMyPoint(0).size(), nsamp = sampler->numMySamples();
    std::vector<Real> loc(dim), mean(dim), pt(dim);
    Real wt(0);
    for (int i = 0; i < nsamp; i++) {
      pt = sampler->getMyPoint(i);
      wt = sampler->getMyWeight(i);
      for (int j = 0; j < dim; j++) {
        loc[j] += wt*pt[j];
      }
    }
    sampler->sumAll(&loc[0],&mean[0],dim);
    return mean;
  }

public:
  /** \brief Default constructor.
  */
  StochasticProblem(void)
    : OptimizationProblem<Real>(),
      setVector_(false), setObjective_(false), setBound_(false),
      setMultiplier_(false), isCon_(false) {}

  /***************************************************************************/
  /********* FULL OPTION LIST ************************************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(civec,isampler,parlist);
    setInequalityMultiplier(livec,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(econ,esampler,parlist);
    setEqualityMultiplier(emul,esampler,parlist);
    setInequalityConstraint(civec,isampler,parlist);
    setInequalityMultiplier(livec,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(icon,isampler,parlist);
    setInequalityMultiplier(imul,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(emul,esampler,parlist);
    setEqualityMultiplier(econ,esampler,parlist);
    setInequalityConstraint(icon,isampler,parlist);
    setInequalityMultiplier(imul,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > esampler(cevec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(civec,isampler,parlist);
    setInequalityMultiplier(livec,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,econ,emul,Teuchos::null,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > esampler(cevec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(icon,isampler,parlist);
    setInequalityMultiplier(imul,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,econ,emul,Teuchos::null,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > isampler(civec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(civec,isampler,parlist);
    setInequalityMultiplier(livec,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > isampler(civec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(icon,isampler,parlist);
    setInequalityMultiplier(imul,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,econ,emul,esampler,civec,livec,bivec,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,econ,emul,esampler,icon,imul,ibnd,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > esampler(cevec.size(),Teuchos::null);
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > isampler(civec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(cevec,esampler,parlist);
    setEqualityMultiplier(levec,esampler,parlist);
    setInequalityConstraint(civec,isampler,parlist);
    setInequalityMultiplier(livec,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > isampler(civec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(econ,Teuchos::null,parlist);
    setEqualityMultiplier(emul,Teuchos::null,parlist);
    setInequalityConstraint(icon,isampler,parlist);
    setInequalityMultiplier(imul,isampler,parlist);
    setInequalityBound(bivec,isampler,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>() {
    std::vector<Teuchos::RCP<SampleGenerator<Real> > > esampler(cevec.size(),Teuchos::null);
    setObjectiveSampleGenerators(vsampler,gsampler,hsampler);
    setObjective(obj,parlist);
    setSolutionVector(x,parlist);
    setBoundConstraint(bnd,parlist);
    setEqualityConstraint(econ,esampler,parlist);
    setEqualityMultiplier(emul,esampler,parlist);
    setInequalityConstraint(icon,Teuchos::null,parlist);
    setInequalityMultiplier(imul,Teuchos::null,parlist);
    setInequalityBound(bivec,Teuchos::null,parlist);
  }

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,econ,emul,Teuchos::null,icon,imul,ibnd,Teuchos::null,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO BOUNDS *******************************************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,esampler,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,esampler,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,esampler,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,esampler,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,Teuchos::null,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    isampler  is a sampler for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,Teuchos::null,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,esampler,civec,livec,bivec,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,esampler,civec,livec,bivec,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,esampler,icon,imul,ibnd,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,esampler,icon,imul,ibnd,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,esampler,civec,livec,bivec,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,esampler,civec,livec,bivec,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraints
      @param[in]    ibnd      is a bound constraint for the inequality constraints
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,icon,imul,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    econ      is an equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint
      @param[in]    icon      is an inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<Constraint<Real > >                   &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,Teuchos::null,icon,imul,ibnd,Teuchos::null,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO INEQUALITIES *************************************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic equality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler = std::vector<Teuchos::RCP<SampleGenerator<Real> > >(),
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,cevec,levec,csampler,Teuchos::null,Teuchos::null,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic equality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is a equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,econ,emul,esampler,Teuchos::null,Teuchos::null,Teuchos::null,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO EQUALITIES ***************************************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic inequality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler = std::vector<Teuchos::RCP<SampleGenerator<Real> > >(),
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,Teuchos::null,Teuchos::null,Teuchos::null,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with bound,
             detereministic and almost sure stochastic inequality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    icon      is a inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    isampler  is a sampler for the inequality constraint (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<Constraint<Real> >                    &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,Teuchos::null,Teuchos::null,Teuchos::null,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO INEQUALITIES OR BOUNDS ***************************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic equality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    cevec     is a vector of equality constraints
      @param[in]    levec     is a vector of multipliers for the equality constraints
      @param[in]    esampler  is a vector of samplers for the equality constraints (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &cevec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &levec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &esampler = std::vector<Teuchos::RCP<SampleGenerator<Real> > >(),
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,cevec,levec,csampler,Teuchos::null,Teuchos::null,Teuchos::null,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic equality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    econ      is a equality constraint
      @param[in]    emul      is a multiplier for the equality constraint
      @param[in]    esampler  is a sampler for the equality constraint (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &econ,
                     const Teuchos::RCP<Vector<Real> >                        &emul,
                     const Teuchos::RCP<SampleGenerator<Real> >               &esampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,econ,emul,esampler,Teuchos::null,Teuchos::null,Teuchos::null,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO EQUALITIES OR BOUNDS *****************************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic inequality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    civec     is a vector of inequality constraints
      @param[in]    livec     is a vector of multipliers for the inequality constraints
      @param[in]    bivec     is a vector of bound constraints for the inequality constraints
      @param[in]    isampler  is a vector of samplers for the inequality constraints (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const std::vector<Teuchos::RCP<Constraint<Real> > >      &civec,
                     const std::vector<Teuchos::RCP<Vector<Real> > >          &livec,
                     const std::vector<Teuchos::RCP<BoundConstraint<Real> > > &bivec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &isampler = std::vector<Teuchos::RCP<SampleGenerator<Real> > >(),
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,civec,livec,bivec,isampler,gsampler,hsampler) {}

  /** \brief Construct stochastic optimization problem with
             detereministic and almost sure stochastic inequality
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    icon      is a inequality constraint
      @param[in]    imul      is a multiplier for the inequality constraint
      @param[in]    ibnd      is a bound constraint for the inequality constraint
      @param[in]    isampler  is a sampler for the inequality constraint (optional)
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<Constraint<Real> >                    &icon,
                     const Teuchos::RCP<Vector<Real> >                        &imul,
                     const Teuchos::RCP<BoundConstraint<Real> >               &ibnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &isampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,icon,imul,ibnd,isampler,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO INEQUALITIES OR EQUALITIES ***********************************/
  /***************************************************************************/

  /** \brief Construct stochastic optimization problem with bound
             constraints.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    bnd       is the bound constraint
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<BoundConstraint<Real> >               &bnd,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,bnd,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,gsampler,hsampler) {}

  /***************************************************************************/
  /********* NO INEQUALITIES, EQUALITIES OR BOUNDS ***************************/
  /***************************************************************************/

  /** \brief Construct unconstrained stochastic optimization problem.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
      @param[in]    x         is the solution vector
      @param[in]    obj       is the random variable objective function
      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    gsampler  is a sampler for the gradient of the objective function (optional)
      @param[in]    hsampler  is a sampler for the hessian of the objective function (optional)
  */
  StochasticProblem( Teuchos::ParameterList                                   &parlist,
                     const Teuchos::RCP<Vector<Real> >                        &x,
                     const Teuchos::RCP<Objective<Real> >                     &obj,
                     const Teuchos::RCP<SampleGenerator<Real> >               &vsampler,
                     const Teuchos::RCP<SampleGenerator<Real> >               &gsampler = Teuchos::null,
                     const Teuchos::RCP<SampleGenerator<Real> >               &hsampler = Teuchos::null)
    : StochasticProblem<Real>(parlist,x,obj,vsampler,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,gsampler,hsampler) {}


  /** \brief Set objective value sample generator.  Must reset objective after
             calling this function.

      @param[in]    vsampler  is a sampler for the value of the objective function
  */
  void setValueSampleGenerator(const Teuchos::RCP<SampleGenerator<Real> > &vsampler) {
    vsampler_ = vsampler;
    if ( gsampler_ == Teuchos::null ) {
      gsampler_ = vsampler_;
    }
    if ( hsampler_ == Teuchos::null ) {
      hsampler_ = gsampler_;
    }
    setObjective_ = false;
  }

  /** \brief Set objective gradient generator.  Must reset objective after
             calling this function.

      @param[in]    gsampler  is a sampler for the gradient of the objective function
  */
  void setGradientSampleGenerator(const Teuchos::RCP<SampleGenerator<Real> > &gsampler) {
    gsampler_ = gsampler;
    if ( hsampler_ == Teuchos::null ) {
      hsampler_ = gsampler_;
    }
    setObjective_ = false;
  }

  /** \brief Set objective hessian sample generator.  Must reset objective after
             calling this function.

      @param[in]    hsampler  is a sampler for the hessian of the objective function
  */
  void setHessianSampleGenerator(const Teuchos::RCP<SampleGenerator<Real> > &hsampler) {
    hsampler_ = hsampler;
    setObjective_ = false;
  }

  /** \brief Set objective sample generators.  Must reset objective after
             calling this function.

      @param[in]    vsampler  is a sampler for the value of the objective function
      @param[in]    gsampler  is a sampler for the gradient of the objective function
      @param[in]    hsampler  is a sampler for the hessian of the objective function
  */
  void setObjectiveSampleGenerators(const Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                                    const Teuchos::RCP<SampleGenerator<Real> > &gsampler = Teuchos::null,
                                    const Teuchos::RCP<SampleGenerator<Real> > &hsampler = Teuchos::null) {
    setHessianSampleGenerator(hsampler);
    setGradientSampleGenerator(gsampler);
    setValueSampleGenerator(vsampler);
  }

  /** \brief Set objective function.  Must reset solution vector and the
 *           bound constraint (if using) after calling this function.

      @param[in]    obj       is the random variable objective function
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setObjective(const Teuchos::RCP<Objective<Real> > &obj,
                    Teuchos::ParameterList               &parlist) {
    INPUT_obj_ = obj;
    if ( vsampler_ == Teuchos::null ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL::StochasticProblem): value sampler not set!");
    }
    else {
      // Determine Stochastic Optimization Type
      std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      Teuchos::RCP<Objective<Real> > obj0;
      if ( type == "Risk Neutral" ) {
        bool storage = parlist.sublist("SOL").get("Store Sampled Value and Gradient",true);
        obj0 = Teuchos::rcp(new RiskNeutralObjective<Real>(obj,vsampler_,gsampler_,hsampler_,storage));
      }
      else if ( type == "Risk Averse" ) {
        obj0 = Teuchos::rcp(new RiskAverseObjective<Real>(obj,parlist,vsampler_,gsampler_,hsampler_));
      }
      else if ( type == "BPOE" ) {
        Real order     = parlist.sublist("SOL").sublist("BPOE").get("Moment Order",1.);
        Real threshold = parlist.sublist("SOL").sublist("BPOE").get("Threshold",0.);
        obj0 = Teuchos::rcp(new BPOEObjective<Real>(obj,order,threshold,vsampler_,gsampler_,hsampler_)); 
      }
      else if ( type == "Mean Value" ) {
        std::vector<Real> mean = computeSampleMean(vsampler_);
        obj->setParameter(mean);
        obj0 = obj;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                   "Invalid stochastic optimization type" << type);
      }
      // Set OptimizationProblem data
      OptimizationProblem<Real>::setObjective(obj0);
    }
    setObjective_ = true;
    setVector_    = false;
    setBound_     = false;
  }

  /** \brief Set solution vector.

      @param[in]    x         is the solution vector
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setSolutionVector(const Teuchos::RCP<Vector<Real> > &sol,
                         Teuchos::ParameterList            &parlist) {
    INPUT_sol_ = sol;
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    Teuchos::RCP<Vector<Real> > sol0;
    if ( type == "Risk Neutral" || type == "Mean Value" ) {
      sol0 = sol;
    }
    else if ( type == "Risk Averse" ) {
      sol0 = Teuchos::rcp(new RiskVector<Real>(parlist,sol));
    }
    else if ( type == "BPOE" ) {
      std::vector<Real> stat(1,1);
      sol0 = Teuchos::rcp(new RiskVector<Real>(sol,stat,true));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setSolutionVector(sol0);
    setVector_ = true;
  }

  /** \brief Set solution statistic.  For risk-averse or bPOE minimization.

      @param[in]    stat      is the solution statistic guess
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setSolutionStatistic(const Real stat,
                            Teuchos::ParameterList &parlist) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Averse" || type == "BPOE" ) {
      Teuchos::dyn_cast<RiskVector<Real> >(*(OptimizationProblem<Real>::getSolutionVector())).setStatistic(stat);
    }
  }

  /** \brief Set bound constraint.

      @param[in]    bnd       is the bound constraint
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setBoundConstraint(const Teuchos::RCP<BoundConstraint<Real> > &bnd,
                          Teuchos::ParameterList                     &parlist) {
    INPUT_bnd_ = bnd;
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    Teuchos::RCP<BoundConstraint<Real> > bnd0;
    if ( type == "Risk Neutral" || type == "Mean Value" ) {
      bnd0 = bnd;
    }
    else if ( type == "Risk Averse" || type == "BPOE" ) {
      bnd0 = Teuchos::rcp(new RiskBoundConstraint<Real>(parlist,bnd));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setBoundConstraint(bnd0);
    setBound_ = true;
  }

  /** \brief Set constraint.  Must reset multiplier vector after calling.

      @param[in]    con       is the vector of constraints
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setConstraint(const std::vector<Teuchos::RCP<Constraint<Real> > > &cvec,
                     Teuchos::ParameterList                              &parlist) {
    isCon_ = false;
    INPUT_con_ = cvec; csampler_.clear();
    int nc = static_cast<int>(cvec.size());
    std::vector<Teuchos::RCP<Constraint<Real> > > cvec0(nc);
    for (int i = 0; i < nc; ++i) {
      if ( cvec[i] != Teuchos::null ) {
        // Determine Stochastic Optimization Type
        std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
        if ( type == "Risk Neutral" || type == "Mean Value" ) {
          cvec0[i] = cvec[i];
        }
        else if ( type == "Risk Averse" || type == "BPOE" ) {
          cvec0[i] = Teuchos::rcp(new RiskLessConstraint<Real>(cvec[i]));
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                     "Invalid stochastic optimization type" << type);
        }
        isCon_ = true;
      }
      else {
        cvec0[i] = cvec[i];
      }
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setConstraint(cvec0);
    setMultiplier_ = false;
  }

  /** \brief Set constraint.  Must reset multiplier vector after calling.

      @param[in]    con       is the vector of constraints
      @param[in]    csampler  is the vector of samplers for the constraints
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setConstraint(const std::vector<Teuchos::RCP<Constraint<Real> > >      &cvec,
                     const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &csampler,
                     Teuchos::ParameterList                                   &parlist) {
    isCon_ = false;
    int nc = static_cast<int>(cvec.size());
    if ( nc != static_cast<int>(csampler.size()) ) {
      throw Exception::NotImplemented(">>> ROL::StochasticProblem::setConstraint: Constraint vector and SampleGenerator vector are not the same size!");
    }
    INPUT_con_ = cvec; csampler_ = csampler;
    std::vector<Teuchos::RCP<Constraint<Real> > > cvec0(nc);
    for (int i = 0; i < nc; ++i) {
      if ( cvec[i] != Teuchos::null ) {
        if ( csampler[i] != Teuchos::null ) {
          // Almost Sure Constraint
          cvec0[i] = Teuchos::rcp(new AlmostSureConstraint<Real>(csampler[i],cvec[i]));
        }
        else {
          // Deterministic Constraint
          // Determine Stochastic Optimization Type
          std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
          if ( type == "Risk Neutral" || type == "Mean Value" ) {
            cvec0[i] = cvec[i];
          }
          else if ( type == "Risk Averse" || type == "BPOE" ) {
            cvec0[i] = Teuchos::rcp(new RiskLessConstraint<Real>(cvec[i]));
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                       "Invalid stochastic optimization type" << type);
          }
        }
        isCon_ = true;
      }
      else {
        cvec0[i] = cvec[i];
      }
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setConstraint(cvec0);
    setMultiplier_ = false;
  }

  /** \brief Set constraint.  Must reset multiplier vector after calling.

      @param[in]    con       is the constraint
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setConstraint(const Teuchos::RCP<Constraint<Real> > &con,
                     Teuchos::ParameterList                &parlist) {
    isCon_ = false;
    INPUT_con_.clear(); INPUT_con_.push_back(con);
    Teuchos::RCP<Constraint<Real> > con0 = Teuchos::null;
    if ( con != Teuchos::null ) {
      // Determine Stochastic Optimization Type
      std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      if ( type == "Risk Neutral" || type == "Mean Value" ) {
        con0 = con;
      }
      else if ( type == "Risk Averse" || type == "BPOE" ) {
        con0 = Teuchos::rcp(new RiskLessConstraint<Real>(con));
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                   "Invalid stochastic optimization type" << type);
      }
      isCon_ = true;
    }
    else {
      con0 = con;
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setConstraint(con0);
    setMultiplier_ = false;
  }

  /** \brief Set constraint.  Must reset multiplier vector after calling.

      @param[in]    con       is the constraint
      @param[in]    csampler  is the sampler for the constraint
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setConstraint(const Teuchos::RCP<Constraint<Real> >      &con,
                     const Teuchos::RCP<SampleGenerator<Real> > &csampler,
                     Teuchos::ParameterList                     &parlist) {
    isCon_ = false;
    INPUT_con_.clear(); INPUT_con_.push_back(con);
    csampler_.clear(); csampler_.push_back(csampler);
    Teuchos::RCP<Constraint<Real> > con0 = Teuchos::null;
    if ( con != Teuchos::null ) {
      if ( csampler != Teuchos::null ) {
        // Almost Sure Constraint
        con0 = Teuchos::rcp(new AlmostSureConstraint(csampler,con));
      }
      else {
        // Deterministic Constraint
        // Determine Stochastic Optimization Type
        std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
        if ( type == "Risk Neutral" || type == "Mean Value" ) {
          con0 = con;
        }
        else if ( type == "Risk Averse" || type == "BPOE" ) {
          con0 = Teuchos::rcp(new RiskLessConstraint<Real>(con));
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                     "Invalid stochastic optimization type" << type);
        }
      }
      isCon_ = true;
    }
    else {
      con0 = con;
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setConstraint(con0);
    setMultiplier_ = false;
  }

  /** \brief Set multiplier vector.

      @param[in]    lvec      is the vector of multiplier vectors
      @param[in]    csampler  is the vector of samplers for the constraints
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setEqualityMultiplier(const std::vector<Teuchos::RCP<Vector<Real> > >          &lvec,
                           const std::vector<Teuchos::RCP<SampleGenerator<Real> > > &csampler,
                           Teuchos::ParameterList                                   &parlist) {
    int nc = static_cast<int>(lvec.size());
    if ( nc != static_cast<int>(csampler.size()) ) {
      throw Exception::NotImplemented(">>> ROL::StochasticProblem::setEqualityMultiplier: Multiplier vector and SampleGenerator vector are not the same size!");
    }
    INPUT_mul_ = lvec;
    std::vector<Teuchos::RCP<Vector<Real> > > lvec0(nc);
    for (int i = 0; i < nc; ++i) {
      if ( csampler[i] != Teuchos::null ) {
        std::vector<Teuchos::RCP<Vector<Real> > > lveci;
        for (int j = 0; j < csampler[i]->numMySamples(); ++j) {
          lveci.push_back(lvec[i]->clone());
          lveci[j]->set(*lvec[i]);
        }
        lvec0[i] = Teuchos::rcp(new DualSimulatedVector<Real>(lveci, csampler[i]->getBatchManager(), csampler[i]));
      }
      else {
        lvec0[i] = lvec[i];
      }
    }
    OptimizationProblem<Real>::setEqualityMultiplier(lvec0);
    setMultiplier_ = true;
  }

  /** \brief Set multiplier vector.

      @param[in]    lvec      is the vector of multiplier vectors
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setEqualityMultiplier(const std::vector<Teuchos::RCP<Vector<Real> > > &lvec,
                           Teuchos::ParameterList                          &parlist) {
    INPUT_mul_ = lvec;
    OptimizationProblem<Real>::setEqualityMultiplier(lvec);
    setMultiplier_ = true;
  }

  /** \brief Set multiplier vector.

      @param[in]    mul       is the multiplier vector
      @param[in]    csampler  is the sampler for the constraint
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setEqualityMultiplier(const Teuchos::RCP<Vector<Real> >          &mul,
                           const Teuchos::RCP<SampleGenerator<Real> > &csampler,
                           Teuchos::ParameterList                     &parlist) {
    INPUT_mul_.clear(); INPUT_mul_.push_back(mul);
    Teuchos::RCP<Vector<Real> > mul0;
    if ( csampler != Teuchos::null ) {
      std::vector<Teuchos::RCP<Vector<Real> > > lveci;
      for (int j = 0; j < csampler[i]->numMySamples(); ++j) {
        lveci.push_back(mul->clone());
        lveci[j]->set(*mul);
      }
      mul0 = Teuchos::rcp(new DualSimulatedVector<Real>(lveci, csampler[i]->getBatchManager(), csampler[i]));
    }
    else {
      mul0 = mul;
    }
    OptimizationProblem<Real>::setEqualityMultiplier(mul0);
    setMultiplier_ = true;
  }

  /** \brief Set multiplier vector.

      @param[in]    mul       is the multiplier vector
      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  void setEqualityMultiplier(const Teuchos::RCP<Vector<Real> > &mul,
                           Teuchos::ParameterList            &parlist) {
    INPUT_mul_.clear(); INPUT_mul_.push_back(mul);
    OptimizationProblem<Real>::setEqualityMultiplier(mul);
    setMultiplier_ = true;
  }

  Teuchos::RCP<Objective<Real> > getObjective(void) {
    if ( !setObjective_ ) {
      throw Exception::NotImplemented(">>> ROL::StochasticProblem::getObjective: No objective inputed!");
    }
    return OptimizationProblem<Real>::getObjective();
  }

  Teuchos::RCP<Vector<Real> > getSolutionVector(void) {
    if ( !setVector_ ) {
      throw Exception::NotImplemented(">>> ROL::StochasticProblem::getSolutionVector: No solution vector inputed!");
    }
    return OptimizationProblem<Real>::getSolutionVector();
  }

  Teuchos::RCP<BoundConstraint<Real> > getBoundConstraint(void) {
    if ( INPUT_bnd_ != Teuchos::null && !setBound_ ) {
      throw Exception::NotImplemented(">>> ROL::StochasticProblem::getBoundConstraint: No bound inputed!");
    }
    return OptimizationProblem<Real>::getBoundConstraint();
  }

  Teuchos::RCP<Vector<Real> > getMultiplier(void) {
    if ( !setMultiplier_ && isCon_ ) {
      throw Exception::NotImplemented(">>> ROL::StochasticProblem::getMultiplier: No multiplier vector inputed!");
    }
    return OptimizationProblem<Real>::getMultiplier();
  }

  /** \brief Returns the statistic from the soluton vector.

      @param[in]    parlist   is the Teuchos::ParameterList with problem specifications
  */
  Real getSolutionStatistic(Teuchos::ParameterList &parlist) {
    try {
      const RiskVector<Real> x = Teuchos::dyn_cast<const RiskVector<Real> >(
        Teuchos::dyn_cast<const Vector<Real> >(*(StochasticProblem<Real>::getSolutionVector())));
      std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      Real val(0);
      if ( type == "Risk Averse" ) {
        Teuchos::ParameterList &list
          = parlist.sublist("SOL").sublist("Risk Measure");
        std::string risk = list.get("Name","CVaR");
        if ( risk == "Mixed-Quantile Quadrangle" ) {
          Teuchos::ParameterList &MQQlist = list.sublist("Mixed-Quantile Quadrangle");
          Teuchos::Array<Real> coeff
            = Teuchos::getArrayFromStringParameter<Real>(MQQlist,"Coefficient Array");
          for (int i = 0; i < coeff.size(); i++) {
            val += coeff[i]*x.getStatistic(i);
          }
        }
        else if ( risk == "Super Quantile Quadrangle" ) {
          SuperQuantileQuadrangle<Real> sqq(parlist);
          val = sqq.computeStatistic(*(StochasticProblem<Real>::getSolutionVector()));
        }
        else if ( risk == "Chebyshev-Kusuoka" ) {
          ChebyshevKusuoka<Real> sqq(parlist);
          val = static_cast<SpectralRisk<Real> >(sqq).computeStatistic(*(StochasticProblem<Real>::getSolutionVector()));
        }
        else if ( risk == "Spectral Risk" ) {
          SpectralRisk<Real> sqq(parlist);
          val = sqq.computeStatistic(*(StochasticProblem<Real>::getSolutionVector()));
        }
        else if ( risk == "Quantile-Radius Quadrangle" ) {
          Real half(0.5);
          val = half*(x.getStatistic(0) + x.getStatistic(1));
        }
        else {
          val = x.getStatistic(0);
        }
      }
      else {
        val = x.getStatistic(0);
      }
      return val;
    }
    catch (std::exception &e) {
      return 0;
    }
  }
};
}
#endif
