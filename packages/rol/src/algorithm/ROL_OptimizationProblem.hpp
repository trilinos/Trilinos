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

#ifndef ROL_OPTIMIZATIONPROBLEM_HPP
#define ROL_OPTIMIZATIONPROBLEM_HPP

#include "Teuchos_ParameterList.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_InteriorPoint.hpp"
#include "ROL_LogBarrierObjective.hpp"
#include "ROL_InequalityConstraint.hpp"
#include "ROL_BoundInequalityConstraint.hpp"
#include "ROL_ObjectiveFromBoundConstraint.hpp"
#include "ROL_RandomVector.hpp"

namespace ROL {

/*
 * Note: We may wish to consider making the get functions private and make Algorithm
 *       a friend of OptimizationProblem as Algorithm is the only class which should 
 *       need these functions and they may return something other than what the user 
 *       expects (penalized instead of raw objective, solution and slack instead of 
 *       solution, etc).
 */

template<class Real>
class OptimizationProblem {

  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type size_type;

private:
  Teuchos::RCP<Objective<Real> >            obj_;
  Teuchos::RCP<Vector<Real> >               sol_;
  Teuchos::RCP<BoundConstraint<Real> >      bnd_;
  Teuchos::RCP<EqualityConstraint<Real> >   con_;
  Teuchos::RCP<InequalityConstraint<Real> > incon_;
  Teuchos::RCP<Vector<Real> >               mul_;
  Teuchos::RCP<Teuchos::ParameterList>      parlist_;

  bool hasSlack_;

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

public:
  virtual ~OptimizationProblem(void) {}

  OptimizationProblem(void)
    : obj_(Teuchos::null), sol_(Teuchos::null), bnd_(Teuchos::null),
      con_(Teuchos::null), mul_(Teuchos::null),
      parlist_(Teuchos::null), hasSlack_(false) {}

  OptimizationProblem(const Teuchos::RCP<Objective<Real> >       &obj,
                      const Teuchos::RCP<Vector<Real> >          &sol,
                      const Teuchos::RCP<BoundConstraint<Real> > &bnd = Teuchos::null,
                      const Teuchos::RCP<Teuchos::ParameterList> &parlist = Teuchos::null)
    : obj_(obj), sol_(sol), bnd_(Teuchos::null), con_(Teuchos::null), mul_(Teuchos::null),
      parlist_(parlist), hasSlack_(false) {
    if ( parlist != Teuchos::null ) {
      if ( bnd != Teuchos::null ) {
        Teuchos::ParameterList &stepList = parlist->sublist("Step");
        std::string step = stepList.get("Type","Trust Region");
        if( step == "Interior Point" ) {
          if ( bnd->isActivated() ) {
            Teuchos::ParameterList &iplist = stepList.sublist("Interior Point");
            Real mu         = iplist.get("Initial Barrier Penalty",1.0);
            Real slack_ival = iplist.get("Initial Slack Variable Value",1.0);
            // Build composite constraint and multipliers
            incon_ = Teuchos::rcp(new BoundInequalityConstraint<Real>(*bnd,*sol));
            con_   = Teuchos::rcp(new InteriorPoint::CompositeConstraint<Real>(incon_));
            Teuchos::RCP<Vector<Real> > lmult1 = sol->dual().clone();
            Teuchos::RCP<Vector<Real> > lmult2 = sol->dual().clone();
            Teuchos::RCP<Vector<Real> > inmul = CreatePartitionedVector(lmult1,lmult2);
            // Create slack variables - fill with parlist value
            Elementwise::Fill<Real> fill(slack_ival);
            Teuchos::RCP<Vector<Real> > slack1 = sol->clone();
            slack1->applyUnary(fill);
            Teuchos::RCP<Vector<Real> > slack2 = sol->clone();
            slack2->applyUnary(fill);
            Teuchos::RCP<Vector<Real> > slack = CreatePartitionedVector(slack1,slack2);
            // Form vector of optimization and slack variables
            sol_ = CreatePartitionedVector(sol,slack);
            // Form partitioned Lagrange multiplier
            mul_ = CreatePartitionedVector(inmul);
            // Create penalty
            Teuchos::RCP<Objective<Real> > barrier
              = Teuchos::rcp( new LogBarrierObjective<Real> );
            obj_ = Teuchos::rcp( new InteriorPoint::PenalizedObjective<Real>(obj,barrier,*sol_,mu) );
          }
          else {
            // Exception
          }
        }
        else { // Not an Interior Point, but have parameters
          bnd_ = bnd;
        }
      }
    }
    else {
      bnd_ = bnd;
    }
  }

  OptimizationProblem(const Teuchos::RCP<Objective<Real> >          &obj,
                      const Teuchos::RCP<Vector<Real> >             &sol,
                      const Teuchos::RCP<EqualityConstraint<Real> > &con,
                      const Teuchos::RCP<Vector<Real> >             &mul,
                      const Teuchos::RCP<Teuchos::ParameterList>    &parlist = Teuchos::null)
    : obj_(obj), sol_(sol), bnd_(Teuchos::null), con_(con), mul_(mul),
      parlist_(parlist), hasSlack_(false) {}

  OptimizationProblem(const Teuchos::RCP<Objective<Real> >          &obj,
                      const Teuchos::RCP<Vector<Real> >             &sol,
                      const Teuchos::RCP<BoundConstraint<Real> >    &bnd,
                      const Teuchos::RCP<EqualityConstraint<Real> > &con,
                      const Teuchos::RCP<Vector<Real> >             &mul,
                      const Teuchos::RCP<Teuchos::ParameterList>    &parlist = Teuchos::null)
    : obj_(obj), sol_(sol), bnd_(Teuchos::null), con_(con), mul_(mul),
      parlist_(parlist), hasSlack_(true) {
    if ( parlist != Teuchos::null ) {
      Teuchos::ParameterList &stepList = parlist->sublist("Step");
      std::string step = stepList.get("Type","Trust Region");
      if ( bnd->isActivated() && step == "Interior Point" ) {
        Teuchos::ParameterList &iplist = stepList.sublist("Interior Point");
        Real mu         = iplist.get("Initial Barrier Penalty",1.0);
        Real slack_ival = iplist.get("Initial Slack Variable Value",1.0);
        // Build composite constraint and multipliers
        incon_ = Teuchos::rcp(new BoundInequalityConstraint<Real>(*bnd,*sol));
        con_   = Teuchos::rcp(new InteriorPoint::CompositeConstraint<Real>(incon_,con));
        Teuchos::RCP<Vector<Real> > lmult1 = sol->clone();
        Teuchos::RCP<Vector<Real> > lmult2 = sol->clone();
        Teuchos::RCP<Vector<Real> > inmul = CreatePartitionedVector(lmult1,lmult2);
        // Create slack variables - fill with parlist value
        Elementwise::Fill<Real> fill(slack_ival);
        Teuchos::RCP<Vector<Real> > slack1 = sol->clone();
        slack1->applyUnary(fill);
        Teuchos::RCP<Vector<Real> > slack2 = sol->clone();
        slack2->applyUnary(fill);
        Teuchos::RCP<Vector<Real> > slack = CreatePartitionedVector(slack1,slack2);
        // Form vector of optimization and slack variables
        sol_ = CreatePartitionedVector(sol,slack);
        // Form partitioned Lagrange multiplier
        mul_ = CreatePartitionedVector(inmul,mul);
        // Create penalty
        Teuchos::RCP<Objective<Real> > barrier
          = Teuchos::rcp( new LogBarrierObjective<Real> );
        obj_ = Teuchos::rcp( new InteriorPoint::PenalizedObjective<Real>(obj,barrier,*sol_,mu) );
      }
      else {
        bnd_ = bnd;
      }
    }
    else {
      bnd_ = bnd;
    }
  }
 

   // For interior points without equality constraint
   OptimizationProblem(const Teuchos::RCP<Objective<Real> > &obj,
                       const Teuchos::RCP<Vector<Real> > &sol,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> > &inmul,
                       const Teuchos::RCP<Teuchos::ParameterList> &parlist )
     : obj_(Teuchos::null), sol_(Teuchos::null),
       con_(Teuchos::null), mul_(Teuchos::null),
       parlist_(Teuchos::null), hasSlack_(true) {

      using InteriorPoint::PenalizedObjective;
      using InteriorPoint::CompositeConstraint;
      using Elementwise::Fill;

      using Teuchos::RCP;  using Teuchos::rcp;

      Teuchos::ParameterList &iplist = parlist->sublist("Interior Point");

      Real mu         = iplist.get("Initial Barrier Penalty",1.0);
      Real slack_ival = iplist.get("Initial Slack Variable Value",1.0);

      con_ = rcp( new CompositeConstraint<Real>(incon) );

      // Create slack variables - fill with parlist value
      RCP<Vector<Real> > slack = inmul->dual().clone();

      Fill<Real> fill(slack_ival);

      slack->applyUnary(fill);

      // Form vector of optimization and slack variables
      sol_ = CreatePartitionedVector(sol,slack);

      // Form partitioned Lagrange multiplier
      mul_ = CreatePartitionedVector(inmul);

      // Create penalty
      RCP<Objective<Real> > barrier = rcp( new LogBarrierObjective<Real> );

      obj_ = rcp( new PenalizedObjective<Real>(obj,barrier,*sol_,mu) );

   }

   // Bound but no equality
   OptimizationProblem(const Teuchos::RCP<Objective<Real> > &obj,
                       const Teuchos::RCP<Vector<Real> > &sol,
                       const Teuchos::RCP<BoundConstraint<Real> > &bnd,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> > &inmul,
                       const Teuchos::RCP<Teuchos::ParameterList> &parlist )
     : obj_(Teuchos::null), sol_(Teuchos::null), bnd_(bnd),
       con_(Teuchos::null), mul_(Teuchos::null),
       parlist_(Teuchos::null), hasSlack_(true) {

      using InteriorPoint::PenalizedObjective;
      using InteriorPoint::CompositeConstraint;
      using Elementwise::Fill;

      using Teuchos::RCP;  using Teuchos::rcp;

      Teuchos::ParameterList &iplist = parlist->sublist("Interior Point");

      Real mu         = iplist.get("Initial Barrier Penalty",1.0);
      Real slack_ival = iplist.get("Initial Slack Variable Value",1.0);

      con_ = rcp( new CompositeConstraint<Real>(incon) );

      // Create slack variables - fill with parlist value
      RCP<Vector<Real> > slack = inmul->dual().clone();

      Fill<Real> fill(slack_ival);

      slack->applyUnary(fill);

      // Form vector of optimization and slack variables
      sol_ = CreatePartitionedVector(sol,slack);

      // Form partitioned Lagrange multiplier
      mul_ = CreatePartitionedVector(inmul);

      // Create penalties
      RCP<Objective<Real> > slack_barrier = rcp( new LogBarrierObjective<Real> );

      RCP<Objective<Real> > bc_barrier = rcp( new ObjectiveFromBoundConstraint<Real>(*bnd,*parlist) );

      obj_ = rcp( new PenalizedObjective<Real>(obj,slack_barrier,bc_barrier,*sol_,mu) );

   }


   // For interior points with equality constraint
   OptimizationProblem(const Teuchos::RCP<Objective<Real> > &obj,
                       const Teuchos::RCP<Vector<Real> > &sol,
                       const Teuchos::RCP<EqualityConstraint<Real> > &eqcon,
                       const Teuchos::RCP<Vector<Real> > &eqmul,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> > &inmul,
                       const Teuchos::RCP<Teuchos::ParameterList> &parlist )
     : obj_(Teuchos::null), sol_(Teuchos::null),
       con_(Teuchos::null), mul_(Teuchos::null),
       parlist_(parlist), hasSlack_(true) {

      using InteriorPoint::PenalizedObjective;
      using InteriorPoint::CompositeConstraint;
      using Elementwise::Fill;

      using Teuchos::RCP;  using Teuchos::rcp;

      Teuchos::ParameterList &iplist = parlist->sublist("Interior Point");

      Real mu         = iplist.get("Initial Barrier Penalty",1.0);
      Real slack_ival = iplist.get("Initial Slack Variable Value",1.0);

      con_ = rcp( new CompositeConstraint<Real>(incon,eqcon) );

      // Create slack variables - fill with parlist value
      RCP<Vector<Real> > slack = inmul->dual().clone();

      Fill<Real> fill(slack_ival);

      slack->applyUnary(fill);

      // Form vector of optimization and slack variables
      sol_ = CreatePartitionedVector(sol,slack);

      // Form partitioned Lagrange multiplier
      mul_ = CreatePartitionedVector(inmul,eqmul);

      // Create penalty
      RCP<Objective<Real> > slack_barrier = rcp( new LogBarrierObjective<Real> );

      obj_ = rcp( new PenalizedObjective<Real>(obj,slack_barrier,*sol_,mu) );

  }

   // Both bound and equality constraint
   OptimizationProblem(const Teuchos::RCP<Objective<Real> > &obj,
                       const Teuchos::RCP<Vector<Real> > &sol,
                       const Teuchos::RCP<BoundConstraint<Real> > &bnd,
                       const Teuchos::RCP<EqualityConstraint<Real> > &eqcon,
                       const Teuchos::RCP<Vector<Real> > &eqmul,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> > &inmul,
                       const Teuchos::RCP<Teuchos::ParameterList> &parlist )
     : obj_(Teuchos::null), sol_(Teuchos::null), bnd_(bnd),
       con_(Teuchos::null), mul_(Teuchos::null),
       parlist_(parlist), hasSlack_(true) {

      using InteriorPoint::PenalizedObjective;
      using InteriorPoint::CompositeConstraint;
      using Elementwise::Fill;

      using Teuchos::RCP;  using Teuchos::rcp;

      Teuchos::ParameterList &iplist = parlist->sublist("Interior Point");

      Real mu         = iplist.get("Initial Barrier Penalty",1.0);
      Real slack_ival = iplist.get("Initial Slack Variable Value",1.0);

      con_ = rcp( new CompositeConstraint<Real>(incon,eqcon) );

      // Create slack variables - fill with parlist value
      RCP<Vector<Real> > slack = inmul->dual().clone();

      Fill<Real> fill(slack_ival);

      slack->applyUnary(fill);

      // Form vector of optimization and slack variables
      sol_ = CreatePartitionedVector(sol,slack);

      // Form partitioned Lagrange multiplier
      mul_ = CreatePartitionedVector(inmul,eqmul);

      // Create penalties
      RCP<Objective<Real> > slack_barrier = rcp( new LogBarrierObjective<Real> );
      RCP<Objective<Real> > bc_barrier = rcp( new ObjectiveFromBoundConstraint<Real>(*bnd) );

      obj_ = rcp( new PenalizedObjective<Real>(obj,slack_barrier,bc_barrier,*sol_,mu) );



  }


  Teuchos::RCP<Objective<Real> > getObjective(void) {
    return obj_;
  }

  void setObjective(const Teuchos::RCP<Objective<Real> > &obj) {
    obj_ = obj;
  }

  Teuchos::RCP<Vector<Real> > getSolutionVector(void) {
    return sol_;
  }

  void setSolutionVector(const Teuchos::RCP<Vector<Real> > &sol) {
    sol_ = sol;
  }

  Teuchos::RCP<BoundConstraint<Real> > getBoundConstraint(void) {
    return bnd_;
  }

  void setBoundConstraint(const Teuchos::RCP<BoundConstraint<Real> > &bnd) {
    bnd_ = bnd;
  }

  Teuchos::RCP<EqualityConstraint<Real> > getEqualityConstraint(void) {
    return con_;
  }

  void setEqualityConstraint(const Teuchos::RCP<EqualityConstraint<Real> > &con) {
    con_ = con;
  }

  Teuchos::RCP<Vector<Real> > getMultiplierVector(void) {
    return mul_;
  }

  void setMultiplierVector(const Teuchos::RCP<Vector<Real> > &mul) {
    mul_ = mul;
  }

  Teuchos::RCP<Teuchos::ParameterList> getParameterList(void) {
    return parlist_;
  }

  void setParameterList( const Teuchos::RCP<Teuchos::ParameterList> &parlist ) {
    parlist_ = parlist;
  }

  virtual std::vector<std::vector<Real> > checkObjectiveGradient( const Vector<Real> &d,
                                                                  const bool printToStream = true,
                                                                  std::ostream & outStream = std::cout,
                                                                  const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                                  const int order = 1 ) {
    if(hasSlack_) {
      Teuchos::RCP<PV> ds = Teuchos::rcp_static_cast<PV>(sol_->clone());
      ds->set(OPT,d);
      RandomizeVector(*(ds->get(SLACK)));
      return obj_->checkGradient(*sol_,*ds,printToStream,outStream,numSteps,order);
    }
    else {
      return obj_->checkGradient(*sol_,d,printToStream,outStream,numSteps,order);
    }
  }

  virtual std::vector<std::vector<Real> > checkObjectiveHessVec( const Vector<Real> &v,
                                                                 const bool printToStream = true,
                                                                 std::ostream & outStream = std::cout,
                                                                 const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                                 const int order = 1 ) {
    if(hasSlack_) {
      Teuchos::RCP<PV> vs = Teuchos::rcp_static_cast<PV>(sol_->clone());
      vs->set(OPT,v);
      RandomizeVector(*(vs->get(SLACK)));
      return obj_->checkHessVec(*sol_,*vs,printToStream,outStream,numSteps,order);     
    } 
    else {
      return obj_->checkHessVec(*sol_,v,printToStream,outStream,numSteps,order);
    }
  }

};
}
#endif
