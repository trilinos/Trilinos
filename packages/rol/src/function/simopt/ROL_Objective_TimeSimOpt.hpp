// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_OBJECTIVE_TIMESIMOPT_HPP
#define ROL_OBJECTIVE_TIMESIMOPT_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_VectorWorkspace.hpp"

/** @ingroup func_group
    \class ROL::Objective_TimeSimOpt
    \brief Defines the time-dependent objective function interface for simulation-based optimization.
           Computes time-local contributions of value, gradient, Hessian-vector product etc to a 
           larger composite objective defined over the simulation time. In contrast to other 
           objective classes Objective_TimeSimOpt has a default implementation of value which 
           returns zero, as time-dependent simulation based optimization problems may have an
           objective value which depends only on the final state of the system. 

           This objective interface inherits from ROL_Objective_SimOpt. Though the interface
           takes two simulation space vectors from spaces
           \f$\mathcal{U_o}\times\mathcal{U_n}\f$. The space \f$\mathcal{U_o}\f$ is ``old'' information
           that accounts for the initial condition on the time interval.

 Comments: It may be worthwhile to provide an implementation of OptimizationProblem for TimeSimOpt
           which recognizes this common use case and avoids unnecessary objective calls. 

 NOTE:     As written, this interface is step agnostic, which needs to be changed
           if a final-time cost is desire or if the time-distributed cost is to
           be weighted non-uniformly
           
    ---
*/



namespace ROL {

template<typename Real>
class Objective_TimeSimOpt : public Objective_SimOpt<Real> {
private:

  // Get the end point of the time intervals vector
  template<int I>
  Vector<Real> & getVector(Vector<Real> & x) const  { 
    return *(static_cast<PartitionedVector<Real>&>(x).get(I));
  }
 
  template<int I>
  const Vector<Real> & getVector(const Vector<Real> & x) const { 
    return *(static_cast<const PartitionedVector<Real>&>(x).get(I));
  }

  mutable VectorWorkspace<Real> workspace_;

protected:

  VectorWorkspace<Real>& getVectorWorkspace() const { return workspace_; }

public:

  using Objective_SimOpt<Real>::Objective_SimOpt;

  // Interface functions (to be overloaded)

  /** \brief Update constraint functions.  
                u_old Is the state from the end of the previous time step.
                u_new Is the state from the end of this time step.
                z Is the control variable
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real>& u_old,
                       const Vector<Real>& u_new,
                       const Vector<Real>& z, 
                       bool flag = true, int iter = -1 ) {
    update_1_old( u_old, flag, iter );
    update_1_new( u_new, flag, iter );
    update_2( z, flag, iter );  
  }
};

  /** \brief Update constraint functions with respect to Sim variable.  
                u_old is the state variable
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1_old( const Vector<Real>& u_old, bool flag = true, int iter = -1 ) {}

  /** \brief Update constraint functions with respect to Sim variable.  
                u_new is the state variable
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1_new( const Vector<Real>& u_new, bool flag = true, int iter = -1 ) {}

  /** \brief Update constraint functions with respect to Opt variable.
                z is the control variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) override {}


  /** \brief Compute contribution to objective value from this time step */
  virtual Real value( const Vector<Real>& u_old, const Vector<Real>& u_new, 
                      const Vector<Real>& z, Real& tol ) { return 0; }


  /** \brief Compute contribution to simulation term gradient from this time step */
   virtual void gradient_1_old( Vector<Real>& g, const Vector<Real>& u_old, 
                                Vector<Real>& u_new, const Vector<Real>& z, Real& tol ) {}

  /** \brief Compute contribution to simulation term gradient from this time step */
   virtual void gradient_1_new( Vector<Real>& g, const Vector<Real>& u_old, 
                                Vector<Real>& u_new, const Vector<Real>& z, Real& tol ) {}
  
   
  /** \brief Compute contribution to optimization term gradient from this time step */
  virtual void gradient_2( Vector<Real>& g, const Vector<Real>& u_old, 
                           Vector<Real>& u_new, const Vector<Real>& z, Real& tol ) override {}

  virtual void hessVec_11_old( Vector<Real> &hv, const Vector<Real> &v_old, 
                               const Vector<Real> &u_old,  const Vector<Real>& u_new,
                               const Vector<Real> &z, Real &tol ) {}

  virtual void hessVec_11_new( Vector<Real> &hv, const Vector<Real> &v_new, 
                               const Vector<Real> &u_old,  const Vector<Real>& u_new,
                               const Vector<Real> &z, Real &tol ) {}
     
  // Functions from SimOpt that are overriden
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void update( const Vector<Real>& u, const Vector<Real>& z, 
                       bool flag = true, int iter = -1 ) override {
    update(getVector<0>(u), getVector<1>(u), z, flag,iter);  
  }

  virtual Real value( const Vector<Real>& u, const Vector<Real>& z, 
                      Real& tol ) override {
    return value( getVector<0>(u), getVector<1>(u), z, tol );
  }
  
  virtual void solve( Vector<Real>& c, Vector<Real>& u, const Vector<Real>& z ) override {
    solve( c, getVector<0>(u), getVector<1>(u), z, tol ); 
  }

  virtual void gradient_1( Vector<Real>& g, const Vector<Real>& u, 
                           const Vector<Real>& z, Real& tol ) override {

    auto& u_old = getVector<0>(u); 
    auto& u_new = getVector<1>(u);

    gradient_1_old( g, u_old, u_new, z, tol );

    auto g_new = workspace_.clone(g);
    
    gradient_1_new( *g_new, u_old, u_new, z, tol );

    g.plus(*g_new);
  }

  virtual void gradient_2( Vector<Real>& g, const Vector<Real>& u, 
                           const Vector<Real>& z, Real& tol ) override {
    auto& u_old = getVector<0>(u); 
    auto& u_new = getVector<1>(u);
  
    gradient_2( g, u_old, u_new, z, tol ); 
  }
 
  virtual void hessVec_11( Vector<Real>& hv, const Vector<Real>& v, 
                           const Vector<Real>& u, const Vector<Real>& z, 
                           Real& tol ) override {

    auto& hv_old = getVector<0>(hv);
    auto& hv_new = getVector<1>(hv);
    auto& v_old  = getVector<0>(v);
    auto& v_new  = getVector<1>(v):
    auto& u_old  = getVector<0>(u);
    auto& u_new  = getVector<1>(u);
  
    hessVec_11( hv_old, v_old, u_old, u_new, z, tol );
    hessVec_11( hv_new, v_new, u_old, u_new, z, tol );  
    
  }
                          
   virtual void hessVec_12( Vector<Real>& hv, const Vector<Real>& v, 
                           const Vector<Real>& u, const Vector<Real>& z, 
                           Real& tol ) override { hv.zero(); }
 
   virtual void hessVec_21( Vector<Real>& hv, const Vector<Real>& v, 
                           const Vector<Real>& u, const Vector<Real>& z, 
                           Real& tol ) override { hv.zero(); }

   virtual void hessVec_22( Vector<Real>& hv, const Vector<Real>& v, 
                           const Vector<Real>& u, const Vector<Real>& z, 
                           Real& tol ) override { hv.zero(); }

} // namespace ROL


#endif // ROL_OBJECTIVE_TIMESIMOPT_HPP

