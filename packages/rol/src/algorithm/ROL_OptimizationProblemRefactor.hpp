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
#include "ROL_BoundConstraint_Partitioned.hpp"
#include "ROL_CompositeConstraint.hpp"
#include "ROL_Objective.hpp"

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

  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type;

private:
  Teuchos::RCP<Objective<Real> >            obj_;
  Teuchos::RCP<Vector<Real> >               sol_;
  Teuchos::RCP<BoundConstraint<Real> >      bnd_;
  Teuchos::RCP<EqualityConstraint<Real> >   con_;
  Teuchos::RCP<Vector<Real> >               mul_;
  Teuchos::RCP<Teuchos::ParameterList>      parlist_;

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

public:
  virtual ~OptimizationProblem(void) {}


  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &sol,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd     = Teuchos::null,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon   = Teuchos::null,
                       const Teuchos::RCP<Vector<Real> >               &eqmul   = Teuchos::null,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon   = Teuchos::null,
                       const Teuchos::RCP<Vector<Real> >               &inmul   = Teuchos::null,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) : 
      obj_(obj), sol_(sol), bnd_(bnd), eqcon_(eqcon), eqmul_(eqcon) {
   
    using Teuchos::RCP; using Teuchos::rcp;

    // If we have an inequality constraint, create slack variables and composite constraint
    if( incon != Teuchos::null && inmul != Teuchos::null ) {

       RCP<Vector<Real> > slack = inmul->dual().clone();

       Real zero(0.0);        

       // Initialize slack variables s_i = max([c(x)]_i,eps)
       incon->value( *slack, *sol, zero );
      
       Real eps = ROL_EPSILON<Real>();
       
       if( parlist != Teuchos::null ) {
         eps = parlist->sublist("General").get("Initial Minimum Slack Value",ROL_EPSILON<Real>());
       }

       slack->applyUnary( Elementwise::ThresholdUpper<Real>(eps) );


       // Upper and lower bounds on slack variables
       RCP<V> lslack = slack->clone();
       lslack->zero();

       RCP<V> uslack = slack->clone();
       uslack->applyUnary(Elementwise::Fill<Real>(ROL_INF<Real>());      

       RCP<BoundConstraint<Real> > bndslack = rcp( new BoundConstraint<Real> (lslack, uslack) );

       // Augment solution vector with slack variables
       sol_ = CreatePartitionedVector(sol,slack);

       // If we have a bound constraint, augment it 
       if( bndcon != Teuchos::null ) {
         bnd_ = rcp( new BoundConstraint_Partitioned<Real>( bnd, bndslack ) );
       }
       else {
         bnd_ = bndslack;
       }


       // If we have an equality constraint, include it
       if( eqcon != Teuchos::null && eqmul != Teuchos::null ) {
         
         con_ = rcp( new CompositeConstraint<Real>( incon, eqcon ) );         

         mul_ = CreatePartitionedVector(inmul,eqmul);   

       }        
       else { // no equality constraint
         
         con_ = rcp( new CompositeConstraint<Real>(incon) );
 
         mul_ = CreatePartitionedVector( inmul );

       }
        
    }  
    
    else { // No inequality constraints


    }

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
