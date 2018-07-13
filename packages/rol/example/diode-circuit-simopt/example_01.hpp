#ifndef ROL_DIODECIRCUIT_HPP
#define ROL_DIODECIRCUIT_HPP

#include "ROL_Algorithm.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

/** \file
    \brief Contains definitions for the diode circuit problem.
    \author Created by T. Takhtaganov, D. Ridzal, D. Kouri
*/


//namespace ROL {
//namespace ZOO {

  /*!
    \brief The diode circuit problem.
    
    The diode circuit problem:
    \f{eqnarray*}{
    \min_{I_S,R_S} \,\, \frac{1}{2}\sum\limits_{n=1}^N (I_n-I_n^{meas})^2 \\
    \text{s.t.}\;\;\begin{cases}c(I_S,R_S,I_1,V^{src}_1)=0\\ \dots \\c(I_S,R_S,I_N,V^{src}_N)=0\end{cases}
    \f}
    where
    \f[c(I_S,R_S,I_n,V^{src}_n)=I_n - I_S\left(\exp\left(\frac{-I_n R_S+V^{src}_n}{V_{th}}\right)-1\right)\f].
  */

template<class Real>
class Constraint_DiodeCircuit : public ROL::Constraint_SimOpt<Real> {
private:
  /// Thermal voltage (constant)
  Real Vth_; 
  /// Vector of source voltages in DC analysis (input) 
  ROL::Ptr<std::vector<Real> > Vsrc_; 
  /// Number of source voltages
  int ns_;
  
  /***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  
  /*!
    \brief Diode equation
    
    Diode equation formula:
      \f$
      I-I_S\left(\exp\left(\frac{V_{src}-IR_S}{V_{th}}\right)-1\right)
      \f$.

      ---
    */
    Real diode(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return I-Is*(exp((Vsrc-I*Rs)/Vth_)-1);
    }
    
    //! Derivative of diode equation wrt I
    Real diodeI(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return 1+Is*exp((Vsrc-I*Rs)/Vth_)*(Rs/Vth_);
    }

    //! Derivative of diode equation wrt Is
    Real diodeIs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return 1-exp((Vsrc-I*Rs)/Vth_);
    }

    //! Derivative of diode equation wrt Rs  
    Real diodeRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return Is*exp((Vsrc-I*Rs)/Vth_)*(I/Vth_);
    }
    
    //! Second derivative of diode equation wrt I^2
    Real diodeII(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return -Is*exp((Vsrc-I*Rs)/Vth_)*(Rs/Vth_)*(Rs/Vth_);
    }

    //! Second derivative of diode equation wrt I and Is
    Real diodeIIs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return exp((Vsrc-I*Rs)/Vth_)*(Rs/Vth_);
    }

    //! Second derivative of diode equation wrt I and Rs
    Real diodeIRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return (Is/Vth_)*exp((Vsrc-I*Rs)/Vth_)*(1-(I*Rs)/Vth_);
    }

    //! Second derivative of diode equation wrt Is^2
    Real diodeIsIs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return 0.0;
    }
    
    //! Second derivative of diode equation wrt Is and Rs
    Real diodeIsRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return exp((Vsrc-I*Rs)/Vth_)*(I/Vth_);
    }

    //! Second derivative of diode equation wrt Rs^2
    Real diodeRsRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      return -Is*exp((Vsrc-I*Rs)/Vth_)*(I/Vth_)*(I/Vth_);
    }

    /*!
      \brief Newton's method with line search

      Solves the diode equation for the current using Newton's method.

      ---
     */
    Real Newton(const Real I, const Real Vsrc, const Real Is, const Real Rs){
      Real EPS = 1.e-16;
      Real TOL = 1.e-13;
      int MAXIT = 200;
      Real IN = I;
      Real fval  = diode(IN,Vsrc,Is,Rs);
      Real dfval = 0.0;
      Real IN_tmp   = 0.0;
      Real fval_tmp = 0.0;
      Real alpha = 1.0;
      for ( int i = 0; i < MAXIT; i++ ) {
	if ( std::abs(fval) < TOL ) {
          // std::cout << "converged with |fval| = " << std::abs(fval) << " and TOL = " << TOL << "\n";
          break;
        }
        dfval = diodeI(IN,Vsrc,Is,Rs);
        if( std::abs(dfval) < EPS ){
          std::cout << "denominator is too small" << std::endl;
          break;
        }
        
        alpha    = 1.0;
        IN_tmp   = IN - alpha*fval/dfval;
        fval_tmp = diode(IN_tmp,Vsrc,Is,Rs);
        while ( std::abs(fval_tmp) >= (1.0-1.e-4*alpha)*std::abs(fval) ) {
          alpha   /= 2.0;
          IN_tmp   = IN - alpha*fval/dfval;
          fval_tmp = diode(IN_tmp,Vsrc,Is,Rs);
          if ( alpha < std::sqrt(EPS) ) { 
            // std::cout << "Step tolerance met\n";
            break;
          }
        }
        IN   = IN_tmp;
        fval = fval_tmp;
    	// if ( i == MAXIT-1){
    	//   std::cout << "did not converge  " << std::abs(fval) << "\n";
    	// }
      }
      return IN;
    }    
    
    //! Solve circuit given optimization parameters Is and Rs
    void solve_circuit(ROL::Vector<Real> &I, const ROL::Vector<Real> &S){
      ROL::Ptr<std::vector<Real> > Ip =
        dynamic_cast<ROL::StdVector<Real>&>(I).getVector();
      ROL::Ptr<const std::vector<Real> > Sp =
        dynamic_cast<const ROL::StdVector<Real>&>(S).getVector();

      int n = Ip->size();
      
      // Using Newton's method      
      Real I0 = 1.e-12; // Initial guess for Newton
      for(int i=0;i<n;i++){
	(*Ip)[i] = Newton(I0,(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
      }
    }    
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:
  
  Constraint_DiodeCircuit(Real Vth, Real Vsrc_min, Real Vsrc_max, Real Vsrc_step){
    Vth_ = Vth;
    ns_ = (Vsrc_max-Vsrc_min)/Vsrc_step + 1;
    Vsrc_ = ROL::makePtr<std::vector<Real>>(ns_,0.0);
    for(int i=0;i<ns_;i++){
      (*Vsrc_)[i] = Vsrc_min+i*Vsrc_step;
    }
  }
  
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol){
    ROL::Ptr<std::vector<Real> > cp = 
      dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0;i<ns_;i++){
      (*cp)[i] = diode((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1]);
    }
  }
  
  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    solve_circuit(u,z);
    value(c,u,z,tol);
  }
  
  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
		       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = 
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0;i<ns_;i++){
      (*jvp)[i] = diodeI((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1])*(*vp)[i];
    }
  }
  
  void applyInverseJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
			      const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = 
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0;i<ns_;i++){
      (*jvp)[i] = (*vp)[i]/diodeI((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1]);
    }
  }
  
  void applyAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = 
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0;i<ns_;i++){
      (*jvp)[i] = diodeI((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1])*(*vp)[i];
    }
  } 
  
  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = 
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0;i<ns_;i++){
      (*jvp)[i] = (*vp)[i]/diodeI((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1]);
    }
  }
  
  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = 
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    jv.zero();
    for (int i=0;i<ns_;i++){
      (*jvp)[i] = diodeIs((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1])*(*vp)[0] + diodeRs((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1])*(*vp)[1];
      }
  }
    
  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = 
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up = 
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = 
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    jv.zero();
    for (int i=0;i<ns_;i++){
      (*jvp)[0] += diodeIs((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1])*(*vp)[i];
      (*jvp)[1] += diodeRs((*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1])*(*vp)[i];
    }
  }
  
  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp =
      dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp =
      dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    ahwv.zero();
    for (int i=0;i<ns_;i++){
      (*ahwvp)[i] = (*vp)[i] * ( diodeII( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] );
    } 
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp =
      dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp =
      dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    ahwv.zero();
    for (int i=0;i<ns_;i++){
      (*ahwvp)[0] += (*vp)[i] * ( diodeIIs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] );
      (*ahwvp)[1] += (*vp)[i] * ( diodeIRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] );
    } 
  }
  
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp =
      dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp =
      dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    ahwv.zero();
    for (int i=0;i<ns_;i++){
      (*ahwvp)[i] = ((*vp)[0] * ( diodeIIs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] ) + (*vp)[1] * ( diodeIRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] ));
    } 
    
  }
  
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp =
      dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp =
      dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    ahwv.zero();
    Real temp1=0.0;
    Real temp2=0.0;
    Real temp3=0.0;
    for (int i=0;i<ns_;i++){
      temp1 += diodeIsIs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i];
      temp2 += diodeIsRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i];
      temp3 += diodeRsRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i];
    }
      // (*ahwvp)[0] += (*vp)[0] * diodeIsIs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] + (*vp)[1] * diodeIsRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i];
      // (*ahwvp)[1] += (*vp)[0] * diodeIsRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i] + (*vp)[1] * diodeRsRs( (*up)[i],(*Vsrc_)[i],(*zp)[0],(*zp)[1] ) * (*wp)[i];
    (*ahwvp)[0] = ((*vp)[0] * temp1 + (*vp)[1] * temp2);
    (*ahwvp)[1] = ((*vp)[0] * temp2 + (*vp)[1] * temp3);
  }
  
};

template<class Real>
class Objective_DiodeCircuit : public ROL::Objective_SimOpt<Real> {
private:
  Real alpha_; // Penalty Parameter
  int ns_,nz_;
  ROL::Ptr<std::vector<Real> > Imeas_;


public: 
  
  Objective_DiodeCircuit(Real alpha, int ns, int nz){
    ns_ = ns;
    nz_ = nz;
    alpha_ = alpha;
    Imeas_ = ROL::makePtr<std::vector<Real>>(ns_,0.0);
    Real temp;
    std::ifstream measurements("measurements.dat");
    //std::cout << "Reading measurements" << std::endl;
    if (measurements.is_open()){
      for (int i=0;i<ns_;i++){
	measurements >> temp;
	measurements >> temp;
	(*Imeas_)[i] = temp;
      }
    }
    measurements.close();
    
  }
    
  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    Real val = 0.0;
    for (int i=0;i<ns_;i++){
      val += ((*up)[i] - (*Imeas_)[i])*((*up)[i] - (*Imeas_)[i]);
    }
    return 0.5*val;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<std::vector<Real> > gup = 
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    // Unwrap x
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE GRADIENT WRT U
    for (int i=0; i<ns_; i++) {
      (*gup)[i] = ((*up)[i]-(*Imeas_)[i]);
    }
    
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<std::vector<Real> > gzp = 
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    // Unwrap x
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE GRADIENT WRT Z
    for (int i=0; i<nz_; i++) {
      (*gzp)[i] = 0.0;
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvup = 
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    // Unwrap v
    ROL::Ptr<const std::vector<Real> > vup =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE GRADIENT WRT U
    for(int i=0; i<ns_; i++){
      (*hvup)[i] = (*vup)[i];
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }
  
  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }
  
  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

};

template<class Real>
class BoundConstraint_DiodeCircuit : public ROL::BoundConstraint<Real> {
private:
  /// Vector of lower bounds
  std::vector<Real> x_lo_;
  /// Vector of upper bounds
  std::vector<Real> x_up_;
  /// Half of the minimum distance between upper and lower bounds
  Real min_diff_;
  /// Scaling for the epsilon margin
  Real scale_;
public:
  BoundConstraint_DiodeCircuit( Real scale, Real lo_Is, Real up_Is, Real lo_Rs, Real up_Rs ){
    x_lo_.push_back(lo_Is);
    x_lo_.push_back(lo_Rs);

    x_up_.push_back(up_Is);
    x_up_.push_back(up_Rs);

    scale_ = scale;
    min_diff_ = 0.5*std::min(x_up_[0]-x_lo_[0],x_up_[1]-x_lo_[1]);
  }
  void project( ROL::Vector<Real> &x ) {
    ROL::Ptr<std::vector<Real> > ex =
      ROL::constPtrCast<std::vector<Real> >((dynamic_cast<ROL::StdVector<Real>&>(x)).getVector());
    (*ex)[0] = std::max(x_lo_[0],std::min(x_up_[0],(*ex)[0]));
    (*ex)[1] = std::max(x_lo_[1],std::min(x_up_[1],(*ex)[1]));
  }
  bool isFeasible( const ROL::Vector<Real> &x ) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    return ((*ex)[0] >= this->x_lo_[0] && (*ex)[1] >= this->x_lo_[1] &&
	    (*ex)[0] <= this->x_up_[0] && (*ex)[1] <= this->x_up_[1]);
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    //epsn *= this->scale_;
    for ( int i = 0; i < 2; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ) {
	(*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    //epsn *= this->scale_;
    for ( int i = 0; i < 2; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn) ) {
	(*ev)[i] = 0.0;
      }
    }
  }
  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    //epsn *= this->scale_;
    for ( int i = 0; i < 2; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ||
	   ((*ex)[i] >= this->x_up_[i]-epsn) ) {
	(*ev)[i] = 0.0;
      }
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<const std::vector<Real> > eg =
      dynamic_cast<const ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    //epsn *= this->scale_;
    for ( int i = 0; i < 2; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ) {
	(*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<const std::vector<Real> > eg =
      dynamic_cast<const ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    //epsn *= this->scale_;
    for ( int i = 0; i < 2; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
	(*ev)[i] = 0.0;
      }
    }
  }
  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<const std::vector<Real> > eg =
      dynamic_cast<const ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    //epsn *= this->scale_;
    for ( int i = 0; i < 2; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ||
	   ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
	(*ev)[i] = 0.0;
      }
    }
  }

  void setVectorToUpperBound( ROL::Vector<Real> &u ) {
    ROL::Ptr<std::vector<Real> > us = ROL::makePtr<std::vector<Real>>(2,0.0);
    us->assign(this->x_up_.begin(),this->x_up_.end());
    ROL::Ptr<ROL::Vector<Real> > up = ROL::makePtr<ROL::StdVector<Real>>(us);
    u.set(*up);
  }

  void setVectorToLowerBound( ROL::Vector<Real> &l ) {
    ROL::Ptr<std::vector<Real> > ls = ROL::makePtr<std::vector<Real>>(2,0.0);
    ls->assign(this->x_lo_.begin(),this->x_lo_.end());
    ROL::Ptr<ROL::Vector<Real> > lp = ROL::makePtr<ROL::StdVector<Real>>(ls);
    l.set(*lp);
  }
};



// /***************************OLD CODE********************************/

//   template<class Real>
//   class Objective_DiodeCircuit : public Objective<Real> {
//   private:
//     /// Thermal voltage (constant)
//     Real Vth_; 
//     /// Vector of measured currents in DC analysis (data)
//     ROL::Ptr<std::vector<Real> > Imeas_;
//     /// Vector of source voltages in DC analysis (input) 
//     ROL::Ptr<std::vector<Real> > Vsrc_; 
//     /// If true, use Lambert-W function to solve circuit, else use Newton's method.
//     bool lambertw_; 
//     /// Percentage of noise to add to measurements; if 0.0 - no noise.
//     Real noise_; 
//     /// If true, use adjoint gradient computation, else compute gradient using sensitivities
//     bool use_adjoint_;
//     /// 0 - use FD(with scaling), 1 - use exact implementation (with second order derivatives), 2 - use Gauss-Newton approximation (first order derivatives only)
//     int use_hessvec_;
  
//   public:

//     /*!
//       \brief A constructor generating data
      
//       Given thermal voltage, minimum and maximum values of source voltages and a step size, values of Is and Rs generates vector of source voltages and solves nonlinear diode equation to  populate the vector of measured currents, which is later used as data. If noise is nonzero, adds random perturbation to data on the order of the magnitude of the components. Sets the flag to use Lambert-W function or Newton's method to solve circuit. Sets the flags to use adjoint gradient computation and one of three Hessian-vector implementations.

//       ---
//      */
//     Objective_DiodeCircuit(Real Vth, Real Vsrc_min, Real Vsrc_max, Real Vsrc_step, Real true_Is, Real true_Rs, bool lambertw, Real noise, bool use_adjoint, int use_hessvec){
//       lambertw_ = lambertw; 
//       Vth_ = Vth;      
//       use_adjoint_ = use_adjoint;
//       use_hessvec_ = use_hessvec;
//       int n = (Vsrc_max-Vsrc_min)/Vsrc_step + 1;
//       Vsrc_ = ROL::makePtr<std::vector<Real>>(n,0.0);
//       Imeas_ = ROL::makePtr<std::vector<Real>>(n,0.0);
//       std::ofstream output ("measurements.dat");
//       Real left = 0.0; Real right = 1.0;
//       if(lambertw_){
// 	// Using Lambert-W function
// 	std::cout << "Generating data using Lambert-W function." << std::endl;
// 	for(int i=0;i<n;i++){
// 	  (*Vsrc_)[i] = Vsrc_min+i*Vsrc_step;
// 	  (*Imeas_)[i] = lambertWCurrent(true_Is,true_Rs,(*Vsrc_)[i]);
// 	  if(noise>0.0){
// 	    (*Imeas_)[i] += noise*pow(10,int(log10((*Imeas_)[i])))*(( (Real)rand() / (Real)RAND_MAX ) * (right - left) + left);
// 	  }
// 	  // Write generated data into file
// 	  if(output.is_open()){
// 	    output << std::setprecision(8) << std::scientific << (*Vsrc_)[i] << "  " << (*Imeas_)[i] << "\n";
// 	  }
// 	}	
//       }
//       else{
// 	// Using Newton's method
// 	std::cout << "Generating data using Newton's method." << std::endl;
// 	for(int i=0;i<n;i++){
// 	  (*Vsrc_)[i] = Vsrc_min+i*Vsrc_step;
// 	  Real I0 = 1.e-12; // initial guess for Newton
// 	  (*Imeas_)[i] = Newton(I0,Vsrc_min+i*Vsrc_step,true_Is,true_Rs);	
// 	  if(noise>0.0){
// 	    (*Imeas_)[i] += noise*pow(10,int(log10((*Imeas_)[i])))*(( (Real)rand() / (Real)RAND_MAX ) * (right - left) + left);
// 	  }
// 	  // Write generated data into file
// 	  if(output.is_open()){
// 	    output << std::setprecision(8) << std::scientific << (*Vsrc_)[i] << "  " << (*Imeas_)[i] << "\n";
// 	  }
// 	}
//       }
	
// 	output.close();
//     }

//     /*!
//       \brief A constructor using data from given file

//       Given thermal voltage and a file with two columns - one for source voltages, another for corresponding currents - populates vectors of source voltages and measured currents. If noise is nonzero, adds random perturbation to data on the order of the magnitude of the components. Sets the flag to use Lambert-W function or Newton's method to solve circuit. Sets the flags to use adjoint gradient computation and one of three Hessian-vector implementations.

//       ---
//     */
//     Objective_DiodeCircuit(Real Vth, std::ifstream& input_file, bool lambertw, Real noise, bool use_adjoint, int use_hessvec){
//       lambertw_ = lambertw; 
//       Vth_ = Vth;     
//       use_adjoint_ = use_adjoint;
//       use_hessvec_ = use_hessvec;
//       std::string line;
//       int dim = 0;
//       for(int k=0;std::getline(input_file,line);++k){dim=k;} // count number of lines
//       input_file.clear(); // reset to beginning of file
//       input_file.seekg(0,std::ios::beg); 
//       Vsrc_ = ROL::makePtr<std::vector<Real>>(dim,0.0);
//       Imeas_ = ROL::makePtr<std::vector<Real>>(dim,0.0);
//       Real Vsrc, Imeas;
//       std::cout << "Using input file to generate data." << "\n";
//       for(int i=0;i<dim;i++){
//         input_file >> Vsrc;
//         input_file >> Imeas;
//         (*Vsrc_)[i] = Vsrc;
//         (*Imeas_)[i] = Imeas;
//       }
//       input_file.close();
//     }

//     //! Change the method for solving the circuit if needed
//     void set_method(bool lambertw){
//       lambertw_ = lambertw;
//     }

//     /*!
//       \brief Diode equation
      
//       Diode equation formula:
//       \f$
//       I-I_S\left(\exp\left(\frac{V_{src}-IR_S}{V_{th}}\right)-1\right)
//       \f$.

//       ---
//     */
//     Real diode(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return I-Is*(exp((Vsrc-I*Rs)/Vth_)-1);
//     }
    
//     //! Derivative of diode equation wrt I
//     Real diodeI(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return 1+Is*exp((Vsrc-I*Rs)/Vth_)*(Rs/Vth_);
//     }

//     //! Derivative of diode equation wrt Is
//     Real diodeIs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return 1-exp((Vsrc-I*Rs)/Vth_);
//     }

//     //! Derivative of diode equation wrt Rs  
//     Real diodeRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return Is*exp((Vsrc-I*Rs)/Vth_)*(I/Vth_);
//     }
    
//     //! Second derivative of diode equation wrt I^2
//     Real diodeII(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return -Is*exp((Vsrc-I*Rs)/Vth_)*(Rs/Vth_)*(Rs/Vth_);
//     }

//     //! Second derivative of diode equation wrt I and Is
//     Real diodeIIs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return exp((Vsrc-I*Rs)/Vth_)*(Rs/Vth_);
//     }

//     //! Second derivative of diode equation wrt I and Rs
//     Real diodeIRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return (Is/Vth_)*exp((Vsrc-I*Rs)/Vth_)*(1-(I*Rs)/Vth_);
//     }

//     //! Second derivative of diode equation wrt Is^2
//     Real diodeIsIs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return 0;
//     }
    
//     //! Second derivative of diode equation wrt Is and Rs
//     Real diodeIsRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return exp((Vsrc-I*Rs)/Vth_)*(I/Vth_);
//     }

//     //! Second derivative of diode equation wrt Rs^2
//     Real diodeRsRs(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       return -Is*exp((Vsrc-I*Rs)/Vth_)*(I/Vth_)*(I/Vth_);
//     }

//     /*!
//       \brief Newton's method with line search

//       Solves the diode equation for the current using Newton's method.

//       ---
//      */
//     Real Newton(const Real I, const Real Vsrc, const Real Is, const Real Rs){
//       Real EPS = 1.e-16;
//       Real TOL = 1.e-13;
//       int MAXIT = 200;
//       Real IN = I;
//       Real fval  = diode(IN,Vsrc,Is,Rs);
//       Real dfval = 0.0;
//       Real IN_tmp   = 0.0;
//       Real fval_tmp = 0.0;
//       Real alpha = 1.0;
//       for ( int i = 0; i < MAXIT; i++ ) {
// 	if ( std::abs(fval) < TOL ) {
//           // std::cout << "converged with |fval| = " << std::abs(fval) << " and TOL = " << TOL << "\n";
//           break;
//         }
//         dfval = diodeI(IN,Vsrc,Is,Rs);
//         if( std::abs(dfval) < EPS ){
//           std::cout << "denominator is too small" << std::endl;
//           break;
//         }
        
//         alpha    = 1.0;
//         IN_tmp   = IN - alpha*fval/dfval;
//         fval_tmp = diode(IN_tmp,Vsrc,Is,Rs);
//         while ( std::abs(fval_tmp) >= (1.0-1.e-4*alpha)*std::abs(fval) ) {
//           alpha   /= 2.0;
//           IN_tmp   = IN - alpha*fval/dfval;
//           fval_tmp = diode(IN_tmp,Vsrc,Is,Rs);
//           if ( alpha < std::sqrt(EPS) ) { 
//             // std::cout << "Step tolerance met\n";
//             break;
//           }
//         }
//         IN   = IN_tmp;
//         fval = fval_tmp;
//     	// if ( i == MAXIT-1){
//     	//   std::cout << "did not converge  " << std::abs(fval) << "\n";
//     	// }
//       }
//       return IN;
//     }    
    
//     /*!
//       \brief Lambert-W function for diodes
      
//       Function      : DeviceSupport::lambertw
//       Purpose       : provides a lambert-w function for diodes and BJT's.
//       Special Notes :
      
//       Purpose.  Evaluate principal branch of Lambert W function at x.
      
//       w = w(x) is the value of Lambert's function.
//       ierr = 0 indicates a safe return.
//       ierr = 1 if x is not in the domain.
//       ierr = 2 if the computer arithmetic contains a bug.
//       xi may be disregarded (it is the error).
      
//       Prototype: void lambertw( Real, Real, int, Real);
      
//       Reference:
//       T.C. Banwell
//       Bipolar transistor circuit analysis using the Lambert W-function,
//       IEEE Transactions on Circuits and Systems I: Fundamental Theory
//       and Applications
      
//       vol. 47, pp. 1621-1633, Nov. 2000.
      
//       Scope         : public
//       Creator       : David Day,  SNL
//       Creation Date : 04/16/02

//       ---
//     */
//     void lambertw(Real x, Real &w, int &ierr, Real &xi){
//       int i=0, maxit = 10;
//       const Real turnpt = -exp(-1.), c1 = 1.5, c2 = .75;
//       Real r, r2, r3, s, mach_eps, relerr = 1., diff;
//       mach_eps = 2.e-15;   // float:2e-7
//       ierr = 0;
      
//       if( x > c1){
// 	w = c2*log(x);
// 	xi = log( x/ w) - w;
//       }
//       else{
// 	if( x >= 0.0){
// 	  w = x;
// 	  if( x == 0. ) return;
// 	  if( x < (1-c2) ) w = x*(1.-x + c1*x*x);
// 	  xi = - w;
// 	}
// 	else{
// 	  if( x >= turnpt){
// 		if( x > -0.2 ){
// 		  w = x*(1.0-x + c1*x*x);
// 		  xi = log(1.0-x + c1*x*x) - w;
// 		}
// 		else{
// 		  diff = x-turnpt;
// 		  if( diff < 0.0 ) diff = -diff;
// 		  w = -1 + sqrt(2.0*exp(1.))*sqrt(x-turnpt);
// 		  if( diff == 0.0 ) return;
// 		  xi = log( x/ w) - w;
// 		}
// 	  }
// 	  else{
// 	    ierr = 1; // x is not in the domain.
// 	    w = -1.0;
// 	    return;
// 	  }
// 	}
//       }
      
//       while( relerr > mach_eps  && i<maxit){
// 	r = xi/(w+1.0);   //singularity at w=-1
// 	r2 = r*r;
// 	r3 = r2*r;
// 	s  = 6.*(w+1.0)*(w+1.0);
// 	w = w * (  1.0 + r + r2/(2.0*( w+1.0)) - (2. * w -1.0)*r3/s  );
// 	if( w * x < 0.0 ) w = -w;
// 	xi = log( x/ w) - w;
	
// 	if( x>1.0 ){
// 	  relerr =  xi / w;
// 	}
// 	else{
// 	  relerr =  xi;
// 	}
// 	if(relerr < 0.0 ) relerr = -relerr;
// 	++i;
//       }
//       if( i == maxit ) ierr = 2;
//     }
    
//     /*!
//       \brief Find currents using Lambert-W function.

//       Reference:
//       T.C. Banwell
//       Bipolar transistor circuit analysis using the Lambert W-function,
//       IEEE Transactions on Circuits and Systems I: Fundamental Theory
//       and Applications
//       vol. 47, pp. 1621-1633, Nov. 2000.

//       ---
//     */
//     Real lambertWCurrent(Real Is, Real Rs, Real Vsrc){
//       Real arg1 = (Vsrc + Is*Rs)/Vth_;
//       Real evd = exp(arg1);
//       Real lambWArg = Is*Rs*evd/Vth_;
//       Real lambWReturn;
//       int ierr;
//       Real lambWError;
//       lambertw(lambWArg, lambWReturn, ierr, lambWError);
//       if(ierr == 1){std::cout << "LambertW error: argument is not in the domain" <<  std::endl; return -1.0;}
//       if(ierr == 2){std::cout << "LambertW error: BUG!" <<  std::endl;}
//       Real Id = -Is+Vth_*(lambWReturn)/Rs;
//       //Real Gd = lambWReturn / ((1 + lambWReturn)*RS);     
//       return Id;
//     }


//     //! Solve circuit given optimization parameters Is and Rs
//     void solve_circuit(Vector<Real> &I, const Vector<Real> &S){
//       ROL::Ptr<std::vector<Real> > Ip =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(I)).getVector());
//       ROL::Ptr<const std::vector<Real> > Sp =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();

//       int n = Ip->size();
      
//       if(lambertw_){
// 	// Using Lambert-W function
// 	Real lambval;
// 	for(int i=0;i<n;i++){
// 	  lambval = lambertWCurrent((*Sp)[0],(*Sp)[1],(*Vsrc_)[i]);
// 	  (*Ip)[i] = lambval;
// 	}
//       }
//       else{	
// 	// Using Newton's method      
// 	Real I0 = 1.e-12; // Initial guess for Newton
// 	for(int i=0;i<n;i++){
// 	  (*Ip)[i] = Newton(I0,(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
// 	}
//       }
//     }

    
//     /*!
//       \brief Evaluate objective function

//       \f$\frac{1}{2}\sum\limits_{i=1}^{N}(I_i-I^{meas}_i)^2\f$
      
//       ---
//      */
//     Real value(const Vector<Real> &S, Real &tol){
//       ROL::Ptr<const std::vector<Real> > Sp =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
//       int n = Imeas_->size();
//       StdVector<Real> I( ROL::makePtr<std::vector<Real>>(n,0.0) );
//       ROL::Ptr<std::vector<Real> > Ip =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(I)).getVector());

//       // Solve state equation
//       this->solve_circuit(I,S);
//       Real val = 0;
      
//       for(int i=0;i<n;i++){
// 	val += ((*Ip)[i]-(*Imeas_)[i])*((*Ip)[i]-(*Imeas_)[i]);
//       }
//       return val/2.0;
//     }
    
//     /*!
//       \brief Solve the adjoint equation
     
//       \f$\lambda_i = \frac{(I^{meas}_i-I_i)}{\frac{\partial c}{\partial I}(I_i,V^{src}_i,I_S,R_S)}\f$

//      ---
//      */
//     void solve_adjoint(Vector<Real> &lambda, const Vector<Real> &I, const Vector<Real> &S){
      
//       ROL::Ptr<std::vector<Real> > lambdap =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(lambda)).getVector());
      
//       ROL::Ptr<const std::vector<Real> > Ip =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(I))).getVector();
//       ROL::Ptr<const std::vector<Real> > Sp =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
      
//       int n = Ip->size();
//       for(int i=0;i<n;i++){
//         (*lambdap)[i] = ((*Imeas_)[i]-(*Ip)[i])/diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
//       }
//     }

//     /*!
//       \brief Solve the sensitivity equation wrt Is
      
//       Computes sensitivity \f[\frac{\partial I}{\partial Is}\f]
      
//       ---
//     */
//     void solve_sensitivity_Is(Vector<Real> &sens, const Vector<Real> &I, const Vector<Real> &S){
      
//       ROL::Ptr<std::vector<Real> > sensp =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(sens)).getVector());
//       ROL::Ptr<const std::vector<Real> > Ip =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(I))).getVector();
//       ROL::Ptr<const std::vector<Real> > Sp =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
      
//       int n = Ip->size();
//       for(int i=0;i<n;i++){
//         (*sensp)[i] = -diodeIs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])/diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
//       }
//     }
    
//     /*!
//       \brief Solve the sensitivity equation wrt Rs
      
//       Computes sensitivity \f[\frac{\partial I}{\partial Rs}\f]
      
//       ---
//     */
//     void solve_sensitivity_Rs(Vector<Real> &sens, const Vector<Real> &I, const Vector<Real> &S){
      
//       ROL::Ptr<std::vector<Real> > sensp =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(sens)).getVector());
//       ROL::Ptr<const std::vector<Real> > Ip =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(I))).getVector();
//       ROL::Ptr<const std::vector<Real> > Sp =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
      
//       int n = Ip->size();
//       for(int i=0;i<n;i++){
//         (*sensp)[i] = -diodeRs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])/diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
//       }
//     }
    
//     //! Compute the gradient of the reduced objective function either using adjoint or using sensitivities
//     void gradient(Vector<Real> &g, const Vector<Real> &S, Real &tol){
      
//       ROL::Ptr<std::vector<Real> > gp =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(g)).getVector());
//       ROL::Ptr<const std::vector<Real> > Sp =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();

//       int n = Imeas_->size();
      
//       StdVector<Real> I( ROL::makePtr<std::vector<Real>>(n,0.0) );
//       ROL::Ptr<std::vector<Real> > Ip =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(I)).getVector());
      
//       // Solve state equation      
//       this->solve_circuit(I,S);
      
//       if(use_adjoint_){      
// 	// Compute the gradient of the reduced objective function using adjoint computation
// 	StdVector<Real> lambda( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	ROL::Ptr<std::vector<Real> > lambdap =
// 	  ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(lambda)).getVector());
	
// 	// Solve adjoint equation
// 	this->solve_adjoint(lambda,I,S);
      
// 	// Compute gradient
// 	(*gp)[0] = 0.0; (*gp)[1] = 0.0;
// 	for(int i=0;i<n;i++){
// 	  (*gp)[0] += diodeIs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])*(*lambdap)[i];
// 	  (*gp)[1] += diodeRs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])*(*lambdap)[i];		
// 	}      
//       }
//       else{
// 	// Compute the gradient of the reduced objective function using sensitivities
// 	StdVector<Real> sensIs( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	StdVector<Real> sensRs( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	// Solve sensitivity equations
// 	this->solve_sensitivity_Is(sensIs,I,S);
// 	this->solve_sensitivity_Rs(sensRs,I,S);
	
// 	ROL::Ptr<std::vector<Real> > sensIsp = ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(sensIs)).getVector());
// 	ROL::Ptr<std::vector<Real> > sensRsp = ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(sensRs)).getVector());
	
// 	// Write sensitivities into file
// 	std::ofstream output ("Sensitivities.dat");
// 	for(int k=0;k<n;k++){
// 	  if(output.is_open()){
// 	    output << std::scientific << (*sensIsp)[k] << " " << (*sensRsp)[k] << "\n";
// 	  }	
// 	}
// 	output.close();
// 	// Compute gradient
// 	(*gp)[0] = 0.0; (*gp)[1] = 0.0;
// 	for(int i=0;i<n;i++){
// 	  (*gp)[0] += ((*Ip)[i]-(*Imeas_)[i])*(*sensIsp)[i];
// 	  (*gp)[1] += ((*Ip)[i]-(*Imeas_)[i])*(*sensRsp)[i];	
// 	}      
//       }
//     }
    

 
//     /*!
//       \brief Compute the Hessian-vector product of the reduced objective function
      
//       Hessian-times-vector computation.

//       ---
//      */
//     void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &S, Real &tol ){
//       if(use_hessvec_==0){
// 	// Use finite-difference approximation
//       	// Modification of parent class function that takes into accout different scale of components
//       	ROL::Ptr<const std::vector<Real> > vp =
//       	  (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(v))).getVector();
//       	ROL::Ptr<const std::vector<Real> > Sp =
//       	  (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
//       	Real gtol = std::sqrt(ROL_EPSILON<Real>()));
	
//       	// Get Step Length                                                                                     
//       	Real h = std::max(1.0,S.norm()/v.norm())*tol;
//       	//Real h = 2.0/(v.norm()*v.norm())*tol;

//       	// Find the scale of componenets of S
//       	Real Is_scale = pow( 10,int( log10( (*Sp)[0] ) ) );                                           
//       	Real Rs_scale = pow( 10,int( log10( (*Sp)[1] ) ) ); 
// 	// Apply scaling
//       	Real h1 = Is_scale*h;
//       	Real h2 = Rs_scale*h;
	
// 	// Compute Gradient at S                                                                               
//       	ROL::Ptr<Vector<Real> > g = S.clone();
//       	this->gradient(*g,S,gtol);

//       	// Compute New Step S + h*v                                                                            
// 	ROL::Ptr<std::vector<Real> > Snewp = ROL::makePtr<std::vector<Real>>(2, 0.0);
//       	ROL::StdVector<Real> Snew(Snewp);
//       	(*Snewp)[0] = (*Sp)[0] + h1*(*vp)[0];
//       	(*Snewp)[1] = (*Sp)[1] + h2*(*vp)[1];
      	
//       	// Compute Gradient at x + h*v                                                                    
//       	hv.zero();
//       	this->gradient(hv,Snew,gtol);

//       	// Compute Newton Quotient                                                                            
//       	hv.axpy(-1.0,*g);
//       	hv.scale(1.0/std::sqrt(h1*h1+h2*h2));
//       }
//       else if(use_hessvec_==1){
// 	ROL::Ptr<std::vector<Real> > hvp =
// 	  ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(hv)).getVector());
// 	ROL::Ptr<const std::vector<Real> > vp =
// 	  (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(v))).getVector();
// 	ROL::Ptr<const std::vector<Real> > Sp =
// 	  (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
	
// 	int n = Imeas_->size();
	
// 	StdVector<Real> I( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	ROL::Ptr<std::vector<Real> > Ip =
// 	  ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(I)).getVector());
	
// 	// Solve state equation      
// 	this->solve_circuit(I,S);
	
	
// 	StdVector<Real> lambda( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	ROL::Ptr<std::vector<Real> > lambdap =
// 	  ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(lambda)).getVector());
	
// 	// Solve adjoint equation
// 	this->solve_adjoint(lambda,I,S);
	
// 	StdVector<Real> w( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	ROL::Ptr<std::vector<Real> > wp =
// 	  ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(w)).getVector());
	
// 	// Solve state sensitivity equation
// 	for(int i=0;i<n;i++){
// 	  (*wp)[i] = ( (*vp)[0] * diodeIs( (*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1] ) + (*vp)[1] * diodeRs( (*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1] ) ) / diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
// 	}
	
// 	StdVector<Real> p( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	ROL::Ptr<std::vector<Real> > pp =
// 	  ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(p)).getVector());
	
// 	// Solve for p
// 	for(int j=0;j<n;j++){
// 	  (*pp)[j] = ( (*wp)[j] + (*lambdap)[j] * diodeII( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] ) * (*wp)[j] - (*lambdap)[j] * diodeIIs( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] ) * (*vp)[0] - (*lambdap)[j] * diodeIRs( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] ) * (*vp)[1] ) / diodeI( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] );
// 	}
	
// 	// Assemble Hessian-vector product
// 	(*hvp)[0] = 0.0;(*hvp)[1] = 0.0;
// 	for(int k=0;k<n;k++){
// 	  (*hvp)[0] += diodeIs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) * (*pp)[k] - (*lambdap)[k] * (*wp)[k] * diodeIIs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) + (*lambdap)[k] * (*vp)[0] * diodeIsIs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) + (*lambdap)[k] * (*vp)[1] * diodeIsRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] );
// 	  (*hvp)[1] += diodeRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) * (*pp)[k] - (*lambdap)[k] * (*wp)[k] * diodeIRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) + (*lambdap)[k] * (*vp)[0] * diodeIsRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) + (*lambdap)[k] * (*vp)[1] * diodeRsRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] );
// 	}
//       }
//       else if(use_hessvec_==2){
// 	//Gauss-Newton approximation
// 	ROL::Ptr<std::vector<Real> > hvp =
//           ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(hv)).getVector());
// 	ROL::Ptr<const std::vector<Real> > vp =
//           (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(v))).getVector();
// 	ROL::Ptr<const std::vector<Real> > Sp =
//           (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(S))).getVector();
	
//         int n = Imeas_->size();

//         StdVector<Real> I( ROL::makePtr<std::vector<Real>>(n,0.0) );
// 	ROL::Ptr<std::vector<Real> > Ip =
//           ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(I)).getVector());

//         // Solve state equation                                                                                
//         this->solve_circuit(I,S);

// 	// Compute sensitivities
// 	StdVector<Real> sensIs( ROL::makePtr<std::vector<Real>>(n,0.0) );
//         StdVector<Real> sensRs( ROL::makePtr<std::vector<Real>>(n,0.0) );
//         // Solve sensitivity equations                                                                          
//         this->solve_sensitivity_Is(sensIs,I,S);
//         this->solve_sensitivity_Rs(sensRs,I,S);
// 	ROL::Ptr<std::vector<Real> > sensIsp = ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(sensIs)).getVector());
// 	ROL::Ptr<std::vector<Real> > sensRsp = ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(sensRs)).getVector());
	
// 	// Compute approximate Hessian
// 	Real H11 = 0.0; Real H12 = 0.0; Real H22 = 0.0;
// 	for(int k=0;k<n;k++){
// 	  H11 += (*sensIsp)[k]*(*sensIsp)[k];
// 	  H12 += (*sensIsp)[k]*(*sensRsp)[k];
// 	  H22 += (*sensRsp)[k]*(*sensRsp)[k];
// 	}
	
// 	// Compute approximate Hessian-times-vector
// 	(*hvp)[0] = H11*(*vp)[0] + H12*(*vp)[1];
// 	(*hvp)[1] = H12*(*vp)[0] + H22*(*vp)[1];
//       }
//       else{
// 	this->ROL::Objective<Real>::hessVec( hv, v, S, tol ); // Use parent class function	
//       }
//     }

//     /*!
//       \brief Generate data to plot objective function

//       Generates a file with three columns - Is value, Rs value, objective value. To plot with gnuplot type:
//       gnuplot;
//       set dgrid3d 100,100;
//       set hidden3d;
//       splot "Objective.dat" u 1:2:3 with lines;

//       ---
//      */
//     void generate_plot(Real Is_lo, Real Is_up, Real Is_step, Real Rs_lo, Real Rs_up, Real Rs_step){
//       ROL::Ptr<std::vector<Real> > S_ptr = ROL::makePtr<std::vector<Real>>(2,0.0);
//       ROL::StdVector<Real> S(S_ptr);
//       std::ofstream output ("Objective.dat");

//       Real Is = 0.0;
//       Real Rs = 0.0;
//       Real val = 0.0;
//       Real tol = 1.e-16;
//       int n = (Is_up-Is_lo)/Is_step + 1;
//       int m = (Rs_up-Rs_lo)/Rs_step + 1;
//       for(int i=0;i<n;i++){
// 	Is = Is_lo + i*Is_step;
// 	for(int j=0;j<m;j++){
// 	  Rs = Rs_lo + j*Rs_step;
// 	  (*S_ptr)[0] = Is;
// 	  (*S_ptr)[1] = Rs;
// 	  val = this->value(S,tol);
// 	  if(output.is_open()){
// 	    output << std::scientific << Is << " " << Rs << " " << val << "\n";
// 	  }
// 	}
//       }
//       output.close();
//     }


//   };
  
//   /*! 
//     \brief Bound constraints on optimization parameters


//     ---
//   */
//   template<class Real>
//   class BoundConstraint_DiodeCircuit : public BoundConstraint<Real> {
//   private:
//     /// Vector of lower bounds
//     std::vector<Real> x_lo_;
//     /// Vector of upper bounds
//     std::vector<Real> x_up_;
//     /// Half of the minimum distance between upper and lower bounds
//     Real min_diff_;
//     /// Scaling for the epsilon margin
//     Real scale_;
//   public:
//     BoundConstraint_DiodeCircuit( Real scale, Real lo_Is, Real up_Is, Real lo_Rs, Real up_Rs ){
//       x_lo_.push_back(lo_Is);
//       x_lo_.push_back(lo_Rs);

//       x_up_.push_back(up_Is);
//       x_up_.push_back(up_Rs);

//       scale_ = scale;
//       min_diff_ = 0.5*std::min(x_up_[0]-x_lo_[0],x_up_[1]-x_lo_[1]);
//     }
//     void project( Vector<Real> &x ) {
//       ROL::Ptr<std::vector<Real> > ex =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(x)).getVector());
//       (*ex)[0] = std::max(x_lo_[0],std::min(x_up_[0],(*ex)[0]));
//       (*ex)[1] = std::max(x_lo_[1],std::min(x_up_[1],(*ex)[1]));
//     }
//     bool isFeasible( const Vector<Real> &x ) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       return ((*ex)[0] >= this->x_lo_[0] && (*ex)[1] >= this->x_lo_[1] &&
//               (*ex)[0] <= this->x_up_[0] && (*ex)[1] <= this->x_up_[1]);
//     }
//     void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       ROL::Ptr<std::vector<Real> > ev =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(v)).getVector());
//       Real epsn = std::min(this->scale_*eps,this->min_diff_);
//       //epsn *= this->scale_;
//       for ( int i = 0; i < 2; i++ ) {
//         if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ) {
//           (*ev)[i] = 0.0;
//         }
//       }
//     }
//     void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       ROL::Ptr<std::vector<Real> > ev =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(v)).getVector());
//       Real epsn = std::min(this->scale_*eps,this->min_diff_);
//       //epsn *= this->scale_;
//       for ( int i = 0; i < 2; i++ ) {
//         if ( ((*ex)[i] >= this->x_up_[i]-epsn) ) {
//           (*ev)[i] = 0.0;
//         }
//       }
//     }
//     void pruneActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       ROL::Ptr<std::vector<Real> > ev =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(v)).getVector());
//       Real epsn = std::min(this->scale_*eps,this->min_diff_);
//       //epsn *= this->scale_;
//       for ( int i = 0; i < 2; i++ ) {
//         if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ||
//              ((*ex)[i] >= this->x_up_[i]-epsn) ) {
//           (*ev)[i] = 0.0;
//         }
//       }
//     }
//     void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       ROL::Ptr<const std::vector<Real> > eg =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(g))).getVector();
//       ROL::Ptr<std::vector<Real> > ev =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(v)).getVector());
//       Real epsn = std::min(this->scale_*eps,this->min_diff_);
//       //epsn *= this->scale_;
//       for ( int i = 0; i < 2; i++ ) {
//         if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ) {
//           (*ev)[i] = 0.0;
//         }
//       }
//     }
//     void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       ROL::Ptr<const std::vector<Real> > eg =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(g))).getVector();
//       ROL::Ptr<std::vector<Real> > ev =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(v)).getVector());
//       Real epsn = std::min(this->scale_*eps,this->min_diff_);
//       //epsn *= this->scale_;
//       for ( int i = 0; i < 2; i++ ) {
//         if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
//           (*ev)[i] = 0.0;
//         }
//       }
//     }
//     void pruneActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
//       ROL::Ptr<const std::vector<Real> > ex =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(x))).getVector();
//       ROL::Ptr<const std::vector<Real> > eg =
//         (dynamic_cast<StdVector<Real> >(const_cast<Vector<Real> &&>(g))).getVector();
//       ROL::Ptr<std::vector<Real> > ev =
//         ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(v)).getVector());
//       Real epsn = std::min(this->scale_*eps,this->min_diff_);
//       //epsn *= this->scale_;
//       for ( int i = 0; i < 2; i++ ) {
//         if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ||
//              ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
//           (*ev)[i] = 0.0;
//         }
//       }
//     }

//     void setVectorToUpperBound( ROL::Vector<Real> &u ) {
//       ROL::Ptr<std::vector<Real> > us = ROL::makePtr<std::vector<Real>>(2,0.0);
//       us->assign(this->x_up_.begin(),this->x_up_.end());
//       ROL::Ptr<ROL::Vector<Real> > up = ROL::makePtr<ROL::StdVector<Real>>(us);
//       u.set(*up);
//     }

//     void setVectorToLowerBound( ROL::Vector<Real> &l ) {
//       ROL::Ptr<std::vector<Real> > ls = ROL::makePtr<std::vector<Real>>(2,0.0);
//       ls->assign(this->x_lo_.begin(),this->x_lo_.end());
//       ROL::Ptr<ROL::Vector<Real> > lp = ROL::makePtr<ROL::StdVector<Real>>(ls);
//       l.set(*lp);
//     }
//   };



//   // template<class Real>
//   // void getDiodeCircuit( ROL::Ptr<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
//   //   // Cast Initial Guess and Solution Vectors                                     
//   //   ROL::Ptr<std::vector<Real> > x0p =
//   //     ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(x0)).getVector());
//   //   ROL::Ptr<std::vector<Real> > xp =
//   //     ROL::constPtrCast<std::vector<Real> >((dynamic_cast<StdVector<Real>&>(x)).getVector());

//   //   int n = xp->size();

//   //   // Resize Vectors                                                                                              
//   //   n = 2;
//   //   x0p->resize(n);
//   //   xp->resize(n);

//   //   // Instantiate Objective Function                                                                              
//   //   obj = ROL::makePtr<Objective_DiodeCircuit<Real>>(0.02585,0.0,1.0,1.e-2);
//   //   //ROL::Objective_DiodeCircuit<Real> obj(0.02585,0.0,1.0,1.e-2);

//   //   // Get Initial Guess
//   //   (*x0p)[0] = 1.e-13;
//   //   (*x0p)[1] = 0.2;
    
//   //   // Get Solution
//   //   (*xp)[0] = 1.e-12;
//   //   (*xp)[1] = 0.25;
    
//   // }


// //} //end namespace ZOO
// //} //end namespace ROL



#endif
