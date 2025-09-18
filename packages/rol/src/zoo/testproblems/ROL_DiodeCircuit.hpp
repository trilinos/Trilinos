#ifndef ROL_DIODECIRCUIT_HPP
#define ROL_DIODECIRCUIT_HPP

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_ScaledStdVector.hpp"

#include <iostream>
#include <fstream>
#include <string>

/** \file
    \brief Contains definitions for the diode circuit problem.
    \author Created by T. Takhtaganov, D. Ridzal, D. Kouri
*/


namespace ROL {
namespace ZOO {

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
class Objective_DiodeCircuit : public Objective<Real> {

  typedef std::vector<Real>            vector;
  typedef Vector<Real>                 V;
  typedef StdVector<Real>              STDV;
  typedef PrimalScaledStdVector<Real>  PSV;
  typedef DualScaledStdVector<Real>    DSV; 
  typedef typename vector::size_type   uint;

private:
  /// Thermal voltage (constant)
  Real Vth_; 
  /// Vector of measured currents in DC analysis (data)
  ROL::Ptr<std::vector<Real> > Imeas_;
  /// Vector of source voltages in DC analysis (input) 
  ROL::Ptr<std::vector<Real> > Vsrc_; 
  /// If true, use Lambert-W function to solve circuit, else use Newton's method.
  bool lambertw_; 
  /// Percentage of noise to add to measurements; if 0.0 - no noise.
  Real noise_; 
  /// If true, use adjoint gradient computation, else compute gradient using sensitivities
  bool use_adjoint_;
  /// 0 - use FD(with scaling),
  /// 1 - use exact implementation (with second order derivatives),
  /// 2 - use Gauss-Newton approximation (first order derivatives only)
  int use_hessvec_;

public:

  /*!
    \brief A constructor generating data
    
    Given thermal voltage, minimum and maximum values of source voltages and
    a step size, values of Is and Rs generates vector of source voltages and
    solves nonlinear diode equation to  populate the vector of measured
    currents, which is later used as data. If noise is nonzero, adds random
    perturbation to data on the order of the magnitude of the components.
    Sets the flag to use Lambert-W function or Newton's method to solve
    circuit. Sets the flags to use adjoint gradient computation and one of
    three Hessian-vector implementations.

    ---
   */
  Objective_DiodeCircuit(Real Vth, Real Vsrc_min, Real Vsrc_max, Real Vsrc_step,
                         Real true_Is, Real true_Rs,
                         bool lambertw, Real noise,
                         bool use_adjoint, int use_hessvec)
    : Vth_(Vth), lambertw_(lambertw), use_adjoint_(use_adjoint), use_hessvec_(use_hessvec) {
    int n  = (Vsrc_max-Vsrc_min)/Vsrc_step + 1;
    Vsrc_  = ROL::makePtr<std::vector<Real>>(n,0.0);
    Imeas_ = ROL::makePtr<std::vector<Real>>(n,0.0);
    std::ofstream output ("Measurements.dat");
    Real left = 0.0, right = 1.0;
    // Generate problem data
    if ( lambertw_ ) {
      std::cout << "Generating data using Lambert-W function." << std::endl;
    }
    else {
      std::cout << "Generating data using Newton's method." << std::endl;
    }
    for ( int i = 0; i < n; i++ ) {
      (*Vsrc_)[i] = Vsrc_min+i*Vsrc_step;
      if (lambertw_) {
        (*Imeas_)[i] = lambertWCurrent(true_Is,true_Rs,(*Vsrc_)[i]);
      }
      else {
        Real I0 = 1.e-12; // initial guess for Newton
        (*Imeas_)[i] = Newton(I0,Vsrc_min+i*Vsrc_step,true_Is,true_Rs);
      }
      if ( noise > 0.0 ) {
        (*Imeas_)[i] += noise*pow(10,(int)log10((*Imeas_)[i]))*random(left, right);
      }
      // Write generated data into file
      if( output.is_open() ) {
        output << std::setprecision(8) << std::scientific << (*Vsrc_)[i] << "  " << (*Imeas_)[i] << "\n";
      }
    }
    output.close();
  }

  /*!
    \brief A constructor using data from given file

    Given thermal voltage and a file with two columns - one for source
    voltages, another for corresponding currents - populates vectors of source
    voltages and measured currents. If noise is nonzero, adds random
    perturbation to data on the order of the magnitude of the components. Sets
    the flag to use Lambert-W function or Newton's method to solve circuit.
    Sets the flags to use adjoint gradient computation and one of three
    Hessian-vector implementations.

    ---
  */
  Objective_DiodeCircuit(Real Vth, std::ifstream& input_file,
                         bool lambertw, Real noise,
                         bool use_adjoint, int use_hessvec)
    : Vth_(Vth), lambertw_(lambertw), use_adjoint_(use_adjoint), use_hessvec_(use_hessvec) {
    std::string line;
    int dim = 0;
    for( int k = 0; std::getline(input_file,line); ++k ) {
      dim = k;
    } // count number of lines
    input_file.clear(); // reset to beginning of file
    input_file.seekg(0,std::ios::beg); 
    Vsrc_  = ROL::makePtr<std::vector<Real>>(dim,0.0);
    Imeas_ = ROL::makePtr<std::vector<Real>>(dim,0.0);
    Real Vsrc, Imeas;
    std::cout << "Using input file to generate data." << "\n";
    for( int i = 0; i < dim; i++ ){
      input_file >> Vsrc;
      input_file >> Imeas;
      (*Vsrc_)[i] = Vsrc;
      (*Imeas_)[i] = Imeas;
    }
    input_file.close();
  }

  //! Change the method for solving the circuit if needed
  void set_method(bool lambertw){
    lambertw_ = lambertw;
  }

  //! Solve circuit given optimization parameters Is and Rs
  void solve_circuit(Vector<Real> &I, const Vector<Real> &S){
    
    ROL::Ptr<vector> Ip = getVector(I);
    ROL::Ptr<const vector> Sp = getVector(S);

    uint n = Ip->size();
    
    if ( lambertw_ ) {
      // Using Lambert-W function
      Real lambval;
      for ( uint i = 0; i < n; i++ ) {
        lambval = lambertWCurrent((*Sp)[0],(*Sp)[1],(*Vsrc_)[i]);
        (*Ip)[i] = lambval;
      }
    }
    else{	
      // Using Newton's method      
      Real I0 = 1.e-12; // Initial guess for Newton
      for ( uint i = 0; i < n; i++ ) {
        (*Ip)[i] = Newton(I0,(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
      }
    }
  }

  /*!
    \brief Evaluate objective function

    \f$\frac{1}{2}\sum\limits_{i=1}^{N}(I_i-I^{meas}_i)^2\f$
    
    ---
   */
  Real value(const Vector<Real> &S, Real &tol){
      
    ROL::Ptr<const vector> Sp = getVector(S);
    uint n = Imeas_->size();
    STDV I( ROL::makePtr<vector>(n,0.0) );
    ROL::Ptr<vector> Ip = getVector(I);

    // Solve state equation
    solve_circuit(I,S);
    Real val = 0;
    
    for ( uint i = 0; i < n; i++ ) {
      val += ((*Ip)[i]-(*Imeas_)[i])*((*Ip)[i]-(*Imeas_)[i]);
    }
    return val/2.0;
  }
    
  //! Compute the gradient of the reduced objective function either using adjoint or using sensitivities
  void gradient(Vector<Real> &g, const Vector<Real> &S, Real &tol){

      
    ROL::Ptr<vector> gp = getVector(g);
    ROL::Ptr<const vector> Sp = getVector(S);
    
    uint n = Imeas_->size();
    
    STDV I( ROL::makePtr<vector>(n,0.0) );
    ROL::Ptr<vector> Ip = getVector(I);
    
    // Solve state equation      
    solve_circuit(I,S);
    
    if ( use_adjoint_ ) {      
      // Compute the gradient of the reduced objective function using adjoint computation
      STDV lambda( ROL::makePtr<vector>(n,0.0) );
      ROL::Ptr<vector> lambdap = getVector(lambda);
      
      // Solve adjoint equation
      solve_adjoint(lambda,I,S);
           
      // Compute gradient
      (*gp)[0] = 0.0; (*gp)[1] = 0.0;
      for ( uint i = 0; i < n; i++ ) {
        (*gp)[0] += diodeIs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])*(*lambdap)[i];
        (*gp)[1] += diodeRs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])*(*lambdap)[i];		
      }      
    }
    else {
      // Compute the gradient of the reduced objective function using sensitivities
      STDV sensIs( ROL::makePtr<vector>(n,0.0) );
      STDV sensRs( ROL::makePtr<vector>(n,0.0) );
      // Solve sensitivity equations
      solve_sensitivity_Is(sensIs,I,S);
      solve_sensitivity_Rs(sensRs,I,S);
      
      ROL::Ptr<vector> sensIsp = getVector(sensIs);
      ROL::Ptr<vector> sensRsp = getVector(sensRs);
      
      // Write sensitivities into file
      std::ofstream output ("Sensitivities.dat");
      for ( uint k = 0; k < n; k++ ) {
        if ( output.is_open() ) {
          output << std::scientific << (*sensIsp)[k] << " " << (*sensRsp)[k] << "\n";
        }
      }
      output.close();
      // Compute gradient
      (*gp)[0] = 0.0; (*gp)[1] = 0.0;
      for( uint i = 0; i < n; i++ ) {
        (*gp)[0] += ((*Ip)[i]-(*Imeas_)[i])*(*sensIsp)[i];
        (*gp)[1] += ((*Ip)[i]-(*Imeas_)[i])*(*sensRsp)[i];	
      }      
    }
  }
  
  /*!
    \brief Compute the Hessian-vector product of the reduced objective function
    
    Hessian-times-vector computation.

    ---
   */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &S, Real &tol ){

          

    if ( use_hessvec_ == 0 ) {
      Objective<Real>::hessVec(hv, v, S, tol);
    }
    else if ( use_hessvec_ == 1 ) {
      ROL::Ptr<vector> hvp = getVector(hv);
      ROL::Ptr<const vector> vp = getVector(v);
      ROL::Ptr<const vector> Sp = getVector(S);
      
      uint n = Imeas_->size();
      
      STDV I( ROL::makePtr<vector>(n,0.0) );
      ROL::Ptr<vector> Ip = getVector(I);
      
      // Solve state equation      
      solve_circuit(I,S);
      
      STDV lambda( ROL::makePtr<vector>(n,0.0) );
      ROL::Ptr<vector> lambdap = getVector(lambda);
      
      // Solve adjoint equation
      solve_adjoint(lambda,I,S);
      
      STDV w( ROL::makePtr<vector>(n,0.0) );
      ROL::Ptr<vector> wp = getVector(w);
      
      // Solve state sensitivity equation
      for ( uint i = 0; i < n; i++ ){
        (*wp)[i] = ( (*vp)[0] * diodeIs( (*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1] )
                   + (*vp)[1] * diodeRs( (*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1] ) )
                   / diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
      }
      
      STDV p( ROL::makePtr<vector>(n,0.0) );
      ROL::Ptr<vector> pp = getVector(p);
      
      // Solve for p
      for ( uint j = 0; j < n; j++ ) {
        (*pp)[j] = ( (*wp)[j] + (*lambdap)[j] * diodeII( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] )
                   * (*wp)[j] - (*lambdap)[j] * diodeIIs( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] )
                   * (*vp)[0] - (*lambdap)[j] * diodeIRs( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] )
                   * (*vp)[1] ) / diodeI( (*Ip)[j],(*Vsrc_)[j],(*Sp)[0],(*Sp)[1] );
      }
      
      // Assemble Hessian-vector product
      (*hvp)[0] = 0.0;(*hvp)[1] = 0.0;
      for ( uint k = 0; k < n; k++ ) {
        (*hvp)[0] += diodeIs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] )* (*pp)[k]
                     - (*lambdap)[k] * (*wp)[k] * diodeIIs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] )
                     + (*lambdap)[k] * (*vp)[0] * diodeIsIs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] )
                     + (*lambdap)[k] * (*vp)[1] * diodeIsRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] );
        (*hvp)[1] += diodeRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] ) * (*pp)[k]
                     - (*lambdap)[k] * (*wp)[k] * diodeIRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] )
                     + (*lambdap)[k] * (*vp)[0] * diodeIsRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] )
                     + (*lambdap)[k] * (*vp)[1] * diodeRsRs( (*Ip)[k],(*Vsrc_)[k],(*Sp)[0],(*Sp)[1] );
      }
    }
    else if ( use_hessvec_ == 2 ) {
      //Gauss-Newton approximation
      ROL::Ptr<vector> hvp = getVector(hv);
      ROL::Ptr<const vector> vp = getVector(v);
      ROL::Ptr<const vector> Sp = getVector(S);
      
      uint n = Imeas_->size();

      STDV I( ROL::makePtr<vector>(n,0.0) );
      ROL::Ptr<vector> Ip = getVector(I);

      // Solve state equation                                                                                
      solve_circuit(I,S);

      // Compute sensitivities
      STDV sensIs( ROL::makePtr<vector>(n,0.0) );
      STDV sensRs( ROL::makePtr<vector>(n,0.0) );

      // Solve sensitivity equations                                                                          
      solve_sensitivity_Is(sensIs,I,S);
      solve_sensitivity_Rs(sensRs,I,S);
      ROL::Ptr<vector> sensIsp = getVector(sensIs);
      ROL::Ptr<vector> sensRsp = getVector(sensRs);
      
      // Compute approximate Hessian
      Real H11 = 0.0; Real H12 = 0.0; Real H22 = 0.0;
      for ( uint k = 0; k < n; k++ ) {
        H11 += (*sensIsp)[k]*(*sensIsp)[k];
        H12 += (*sensIsp)[k]*(*sensRsp)[k];
        H22 += (*sensRsp)[k]*(*sensRsp)[k];
      }
      
      // Compute approximate Hessian-times-vector
      (*hvp)[0] = H11*(*vp)[0] + H12*(*vp)[1];
      (*hvp)[1] = H12*(*vp)[0] + H22*(*vp)[1];
    }
    else {
      ROL::Objective<Real>::hessVec( hv, v, S, tol ); // Use parent class function	
    }
  }

  /*!
    \brief Generate data to plot objective function

    Generates a file with three columns - Is value, Rs value, objective value. To plot with gnuplot type:
    gnuplot;
    set dgrid3d 100,100;
    set hidden3d;
    splot "Objective.dat" u 1:2:3 with lines;

    ---
   */
  void generate_plot(Real Is_lo, Real Is_up, Real Is_step, Real Rs_lo, Real Rs_up, Real Rs_step){
    ROL::Ptr<std::vector<Real> > S_ptr = ROL::makePtr<std::vector<Real>>(2,0.0);
    StdVector<Real> S(S_ptr);
    std::ofstream output ("Objective.dat");

    Real Is = 0.0;
    Real Rs = 0.0;
    Real val = 0.0;
    Real tol = 1.e-16;
    int n = (Is_up-Is_lo)/Is_step + 1;
    int m = (Rs_up-Rs_lo)/Rs_step + 1;
    for ( int i = 0; i < n; i++ ) {
      Is = Is_lo + i*Is_step;
      for ( int j = 0; j < m; j++ ) {
        Rs = Rs_lo + j*Rs_step;
        (*S_ptr)[0] = Is;
        (*S_ptr)[1] = Rs;
        val = value(S,tol);
        if ( output.is_open() ) {
          output << std::scientific << Is << " " << Rs << " " << val << std::endl;
        }
      }
    }
    output.close();
  }

private:

  ROL::Ptr<const vector> getVector( const V& x ) {
    try { 
      return dynamic_cast<const STDV&>(x).getVector();
    }
    catch (std::exception &e) {
      try { 
        return dynamic_cast<const PSV&>(x).getVector();
      }
      catch (std::exception &le) {
        return dynamic_cast<const DSV&>(x).getVector();
      }
    }
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    try {
      return dynamic_cast<STDV&>(x).getVector(); 
    }
    catch (std::exception &e) {
      try {
        return dynamic_cast<PSV&>(x).getVector(); 
      }
      catch (std::exception &le) {
        return dynamic_cast<DSV&>(x).getVector(); 
      }
    }
  }

  Real random(const Real left, const Real right) const {
    return (Real)rand()/(Real)RAND_MAX * (right - left) + left;
  }

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
    return 0;
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
  
  /*!
    \brief Lambert-W function for diodes
    
    Function      : DeviceSupport::lambertw
    Purpose       : provides a lambert-w function for diodes and BJT's.
    Special Notes :
    
    Purpose.  Evaluate principal branch of Lambert W function at x.
    
    w = w(x) is the value of Lambert's function.
    ierr = 0 indicates a safe return.
    ierr = 1 if x is not in the domain.
    ierr = 2 if the computer arithmetic contains a bug.
    xi may be disregarded (it is the error).
    
    Prototype: void lambertw( Real, Real, int, Real);
    
    Reference:
    T.C. Banwell
    Bipolar transistor circuit analysis using the Lambert W-function,
    IEEE Transactions on Circuits and Systems I: Fundamental Theory
    and Applications
    
    vol. 47, pp. 1621-1633, Nov. 2000.
    
    Scope         : public
    Creator       : David Day,  SNL
    Creation Date : 04/16/02

    ---
  */
  void lambertw(Real x, Real &w, int &ierr, Real &xi){
    int i = 0, maxit = 10;
    const Real turnpt = -exp(-1.), c1 = 1.5, c2 = .75;
    Real r, r2, r3, s, mach_eps, relerr = 1., diff;
    mach_eps = 2.e-15;   // float:2e-7
    ierr = 0;
    
    if ( x > c1 ) {
      w  = c2*log(x);
      xi = log( x/ w) - w;
    }
    else {
      if ( x >= 0.0 ) {
        w = x;
        if ( x == 0. ) {
          return;
        }
        if ( x < (1-c2) ) {
          w = x*(1.-x + c1*x*x);
        }
        xi = - w;
      }
      else {
        if ( x >= turnpt ){
      	  if ( x > -0.2 ){
      	    w  = x*(1.0-x + c1*x*x);
      	    xi = log(1.0-x + c1*x*x) - w;
      	  }
      	  else {
      	    diff = x-turnpt;
      	    if ( diff < 0.0 ) {
              diff = -diff;
            }
      	    w = -1 + sqrt(2.0*exp(1.))*sqrt(x-turnpt);
      	    if ( diff == 0.0 ) {
              return;
            }
      	    xi = log( x/ w) - w;
      	  }
        }
        else {
          ierr = 1; // x is not in the domain.
          w = -1.0;
          return;
        }
      }
    }
    
    while ( relerr > mach_eps  && i < maxit ) {
      r  = xi/(w+1.0);   //singularity at w=-1
      r2 = r*r;
      r3 = r2*r;
      s  = 6.*(w+1.0)*(w+1.0);
      w  = w * ( 1.0 + r + r2/(2.0*( w+1.0)) - (2. * w -1.0)*r3/s );
      w  = ((w*x < 0.0) ? -w : w);
      xi = log( x/ w) - w;
      
      relerr = ((x > 1.0) ? xi/w : xi);
      relerr = ((relerr < 0.0) ? -relerr : relerr);
      ++i;
    }
    ierr = ((i == maxit) ? 2 : ierr);
  }
    
  /*!
    \brief Find currents using Lambert-W function.

    Reference:
    T.C. Banwell
    Bipolar transistor circuit analysis using the Lambert W-function,
    IEEE Transactions on Circuits and Systems I: Fundamental Theory
    and Applications
    vol. 47, pp. 1621-1633, Nov. 2000.

    ---
  */
  Real lambertWCurrent(Real Is, Real Rs, Real Vsrc){
    Real arg1        = (Vsrc + Is*Rs)/Vth_;
    Real evd         = exp(arg1);
    Real lambWArg    = Is*Rs*evd/Vth_;
    Real lambWReturn = 0.0;
    Real lambWError  = 0.0;
    int ierr = 0;
    lambertw(lambWArg, lambWReturn, ierr, lambWError);
    if ( ierr == 1 ) {
     std::cout << "LambertW error: argument is not in the domain" <<  std::endl;
     return -1.0;
    }
    if ( ierr == 2 ) {
      std::cout << "LambertW error: BUG!" <<  std::endl;
    }
    Real Id = -Is+Vth_*(lambWReturn)/Rs;
    //Real Gd = lambWReturn / ((1 + lambWReturn)*RS);     
    return Id;
  }
    
  /*!
    \brief Solve the adjoint equation
   
    \f$\lambda_i = \frac{(I^{meas}_i-I_i)}{\frac{\partial c}{\partial I}(I_i,V^{src}_i,I_S,R_S)}\f$

   ---
   */
  void solve_adjoint(Vector<Real> &lambda, const Vector<Real> &I, const Vector<Real> &S){
    
    
    ROL::Ptr<vector> lambdap = getVector(lambda);
    ROL::Ptr<const vector> Ip = getVector(I);
    ROL::Ptr<const vector> Sp = getVector(S);
    
    uint n = Ip->size();
    for ( uint i = 0; i < n; i++ ){
      (*lambdap)[i] = ((*Imeas_)[i]-(*Ip)[i])
                      /diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
    }
  }

  /*!
    \brief Solve the sensitivity equation wrt Is
    
    Computes sensitivity \f[\frac{\partial I}{\partial Is}\f]
    
    ---
  */
  void solve_sensitivity_Is(Vector<Real> &sens, const Vector<Real> &I, const Vector<Real> &S){

    
    ROL::Ptr<vector> sensp = getVector(sens);
    ROL::Ptr<const vector> Ip = getVector(I);
    ROL::Ptr<const vector> Sp = getVector(S);     
    
    uint n = Ip->size();
    for ( uint i = 0; i < n; i++ ) {
      (*sensp)[i] = -diodeIs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])
                    /diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
    }
  }
    
  /*!
    \brief Solve the sensitivity equation wrt Rs
    
    Computes sensitivity \f[\frac{\partial I}{\partial Rs}\f]
    
    ---
  */
  void solve_sensitivity_Rs(Vector<Real> &sens, const Vector<Real> &I, const Vector<Real> &S){
         
    
    ROL::Ptr<vector> sensp = getVector(sens);
    ROL::Ptr<const vector> Ip = getVector(I);
    ROL::Ptr<const vector> Sp = getVector(S);
    
    uint n = Ip->size();
    for ( uint i = 0; i < n; i++ ) {
      (*sensp)[i] = -diodeRs((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1])
                    /diodeI((*Ip)[i],(*Vsrc_)[i],(*Sp)[0],(*Sp)[1]);
    }
  }
};  // class Objective_DiodeCircuit 


  // template<class Real>
  // void getDiodeCircuit( ROL::Ptr<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
  //   // Cast Initial Guess and Solution Vectors                                     
  //   ROL::Ptr<std::vector<Real> > x0p =
  //     ROL::constPtrCast<std::vector<Real> >((dynamic_cast<PrimalScaledStdVector<Real>&>(x0)).getVector());
  //   ROL::Ptr<std::vector<Real> > xp =
  //     ROL::constPtrCast<std::vector<Real> >((dynamic_cast<PrimalScaledStdVector<Real>&>(x)).getVector());

  //   int n = xp->size();

  //   // Resize Vectors                                                                                              
  //   n = 2;
  //   x0p->resize(n);
  //   xp->resize(n);

  //   // Instantiate Objective Function                                                                              
  //   obj = ROL::makePtr<Objective_DiodeCircuit<Real>>(0.02585,0.0,1.0,1.e-2);
  //   //ROL::Objective_DiodeCircuit<Real> obj(0.02585,0.0,1.0,1.e-2);

  //   // Get Initial Guess
  //   (*x0p)[0] = 1.e-13;
  //   (*x0p)[1] = 0.2;
    
  //   // Get Solution
  //   (*xp)[0] = 1.e-12;
  //   (*xp)[1] = 0.25;
    
  // }


} //end namespace ZOO
} //end namespace ROL

#endif
