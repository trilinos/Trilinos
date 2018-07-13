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

#ifndef ROL_BUNDLE_TT_H
#define ROL_BUNDLE_TT_H

#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bundle.hpp"

#include "ROL_Ptr.hpp"
#include "ROL_LAPACK.hpp" 

#include <vector>
#include <limits.h> 
#include <stdint.h> 
#include <float.h> 
#include <math.h> 
#include <algorithm> // TT: std::find

#include "ROL_LinearAlgebra.hpp"

#define EXACT 1
#define TABOO_LIST 1
#define FIRST_VIOLATED 0

/** \class ROL::Bundle_TT
    \brief Provides the interface for and implements a bundle. The semidefinite
           quadratic subproblem is solved using TT algorithm by Antonio
           Frangioni (1996).
*/

namespace ROL {

template<class Real>
class Bundle_TT : public Bundle<Real> {
private: 
  ROL::LAPACK<int, Real> lapack_; // TT

  int QPStatus_;           // QP solver status
  int maxind_;             // maximum integer value
  int entering_;           // index of entering item
  int LiMax_;              // index of max element of diag(L)
  int LiMin_;              // index of min element of diag(L)

  unsigned maxSize_;       // maximum bundle size
  unsigned dependent_;     // number of lin. dependent items in base
  unsigned currSize_;      // current size of base

  bool isInitialized_;
  bool optimal_;           // flag for optimality of restricted solution

  Real rho_;
  Real lhNorm;
  Real ljNorm;
  Real lhz1_;
  Real lhz2_;
  Real ljz1_;
  Real kappa_;             // condition number of matrix L ( >= max|L_ii|/min|L_ii| )
  Real objval_;            // value of objective
  Real minobjval_;         // min value of objective (ever reached)
  Real deltaLh_;           // needed in case dependent row becomes independent
  Real deltaLj_;           // needed in case dependent row becomes independent

  std::vector<int> taboo_; // list of "taboo" items
  std::vector<int> base_;  // base

  LA::Matrix<Real> L_;
  LA::Matrix<Real> Id_;
  LA::Vector<Real> tempv_;
  LA::Vector<Real> tempw1_;
  LA::Vector<Real> tempw2_;
  LA::Vector<Real> lh_;
  LA::Vector<Real> lj_;
  LA::Vector<Real> z1_;
  LA::Vector<Real> z2_;

  Real GiGj(const int i, const int j) const {
    return (Bundle<Real>::subgradient(i)).dot(Bundle<Real>::subgradient(j));
  }

  Real sgn(const Real x) const {
    const Real zero(0), one(1);
    return ((x < zero) ? -one :
           ((x > zero) ?  one : zero));
  }

public:
  Bundle_TT(const unsigned maxSize = 10,
            const Real coeff = 0.0,
            const Real omega = 2.0,
            const unsigned remSize = 2) 
    : Bundle<Real>(maxSize,coeff,omega,remSize),
      maxSize_(maxSize), isInitialized_(false) {
    maxind_ = std::numeric_limits<int>::max();
    Id_.reshape(maxSize_,maxSize_);
    for(unsigned i=0; i<maxSize_; ++i) {
      Id_(i,i) = static_cast<Real>(1);
    }
  }
  
  unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    unsigned iter = 0;
    if (Bundle<Real>::size() == 1) {
      iter = Bundle<Real>::solveDual_dim1(t,maxit,tol);
    }
    else if (Bundle<Real>::size() == 2) {
      iter = Bundle<Real>::solveDual_dim2(t,maxit,tol);
    }
    else {
      iter = solveDual_arbitrary(t,maxit,tol);
    }
    return iter;
  }

/***********************************************************************************************/
/****************** DUAL CUTTING PLANE SUBPROBLEM ROUTINES *************************************/
/***********************************************************************************************/
private:
  void swapRowsL(unsigned ind1, unsigned ind2, bool trans=false) {
    const Real zero(0), one(1);
    if( ind1 > ind2){
      unsigned tmp = ind1;
      ind2 = ind1;
      ind1 = tmp;
    }
    unsigned dd = ind1;
    for (unsigned n=ind1+1; n<=ind2; ++n){
      LA::Matrix<Real> Id_n(LA::Copy,Id_,currSize_,currSize_);
      Id_n(dd,dd) = zero; Id_n(dd,n) = one;
      Id_n(n,dd)  = one;  Id_n(n,n)  = zero;
      LA::Matrix<Real> prod(currSize_,currSize_);
      if( !trans ) {
        prod.multiply(LA::ETransp::NO_TRANS,LA::ETransp::NO_TRANS,one,Id_n,L_,zero);
      }
      else {
        prod.multiply(LA::ETransp::NO_TRANS,LA::ETransp::NO_TRANS,one,L_,Id_n,zero);
      }
      L_ = prod;
      dd++;
    }
  }

  void updateK(void) {
    if (currSize_ <= dependent_) { // L is empty
      kappa_ = static_cast<Real>(1);
    }
    else{
      Real tmpdiagMax = ROL_NINF<Real>();
      Real tmpdiagMin = ROL_INF<Real>();
      for (unsigned j=0;j<currSize_-dependent_;j++){
        if( L_(j,j) > tmpdiagMax ){
          tmpdiagMax = L_(j,j);
          LiMax_ = j;
        }
        if( L_(j,j) < tmpdiagMin ){
          tmpdiagMin = L_(j,j);
          LiMin_ = j;
        }
      }
      kappa_ = tmpdiagMax/tmpdiagMin;
    }
  }

  void addSubgradToBase(unsigned ind, Real delta) {
    // update z1, z2, kappa
    // swap rows if: dependent == 1 and we insert independent row (dependent row is always last)
    //               dependent == 2 and Lj has become independent (Lh still dependent)
    if(dependent_ && (ind == currSize_-1)){
        swapRowsL(currSize_-2,currSize_-1);
        int tmp;
        tmp = base_[currSize_-2];
        base_[currSize_-2] = base_[currSize_-1];
        base_[currSize_-1] = tmp;
        ind--;
    } // end if dependent
    
    L_(ind,ind) = delta;
    
    // update z1 and z2
    unsigned zsize = ind+1;
    z1_.resize(zsize); z2_.resize(zsize);
    z1_[ind] = ( static_cast<Real>(1) - lhz1_ ) / delta;
    z2_[ind] = ( Bundle<Real>::alpha(base_[ind]) - lhz2_ ) / delta;  
    //z2[zsize-1] = ( Bundle<Real>::alpha(entering_) - lhz2_ ) / delta;  
    
    // update kappa
    if(delta > L_(LiMax_,LiMax_)){
      LiMax_ = ind;
      kappa_ = delta/L_(LiMin_,LiMin_);
    }
    if(delta < L_(LiMin_,LiMin_)){
      LiMin_ = ind;
      kappa_ = L_(LiMax_,LiMax_)/delta;
    }
  }

  void deleteSubgradFromBase(unsigned ind, Real tol){
    const Real zero(0), one(1);
    // update L, currSize, base_, z1, z2, dependent, dualVariables_, kappa
    if (ind >= currSize_-dependent_){
      // if dependent > 0, the last one or two rows of L are lin. dependent              
      if (ind < currSize_-1){ // eliminate currSize_-2 but keep currSize_-1
        // swap the last row with the second to last
        swapRowsL(ind,currSize_-1);
        base_[ind] = base_[currSize_-1];                
#if( ! EXACT )
        lhNorm = ljNorm; // new last row is lh
#endif
      }
 
      dependent_--;
      currSize_--;
      L_.reshape(currSize_,currSize_); // the row to be eliminated is the last row
      base_.resize(currSize_);

      // note: z1, z2, kappa need not be updated
      return;
    } // end if dependent item

    /* currently L_B is lower trapezoidal
            
             | L_1  0  0   |
       L_B = | l    d  0   |
             | Z    v  L_2 | 
            
       Apply Givens rotations to transform it to
            
       | L_1  0  0    |
       | l    d  0    |
       | Z    0  L_2' |
            
       then delete row and column to obtain factorization of L_B' with B' = B/{i}
             
       L_B' = | L_1  0    |
              | Z    L_2' |
            
    */
    for (unsigned j=ind+1; j<currSize_-dependent_; ++j){
      Real ai = L_(j,ind);
      if (std::abs(ai) <= tol*currSize_) { // nothing to do
        continue;
      }
      Real aj = L_(j,j);
      Real d, Gc, Gs;
      // Find Givens
      // Anderson's implementation
      if (std::abs(aj) <= tol*currSize_){ // aj is zero
        Gc = zero;
        d  = std::abs(ai);
        Gs = -sgn(ai); // Gs = -sgn(ai)
      }
      else if ( std::abs(ai) > std::abs(aj) ){
        Real t = aj/ai;
        Real sgnb = sgn(ai);
        Real u = sgnb * std::sqrt(one + t*t );
        Gs = -one/u;
        Gc = -Gs*t;
        d  = ai*u;
      }
      else{
        Real t = ai/aj;
        Real sgna = sgn(aj);
        Real u = sgna * std::sqrt(one + t*t );
        Gc = one/u;
        Gs = -Gc*t;
        d  = aj*u;
      }

      // // "Naive" implementation
      // d  = hypot(ai,aj);
      // Gc = aj/d;
      // Gs = -ai/d;
              
      L_(j,j) = d; L_(j,ind) = zero; 
      // apply Givens to columns i,j of L
      for (unsigned h=j+1; h<currSize_; ++h){
        Real tmp1 = L_(h,ind);
        Real tmp2 = L_(h,j);
        L_(h,ind) = Gc*tmp1 + Gs*tmp2;
        L_(h,j) = Gc*tmp2 - Gs*tmp1;
      }
      // apply Givens to z1, z2
      Real tmp1 = z1_[ind];
      Real tmp2 = z1_[j];
      Real tmp3 = z2_[ind];
      Real tmp4 = z2_[j];
      z1_[ind] = Gc*tmp1 + Gs*tmp2;
      z1_[j] = Gc*tmp2 - Gs*tmp1;
      z2_[ind] = Gc*tmp3 + Gs*tmp4;
      z2_[j] = Gc*tmp4 - Gs*tmp3;
    }// end loop over j

    if(dependent_){
      deltaLh_ = L_(currSize_-dependent_,ind);  // h = currSize_ - dependent
      if( dependent_ > 1 )                 // j = currSize_ - 1, h = currSize_ - 2
        deltaLj_ = L_(currSize_-1,ind);
    }
            
    // shift rows and columns of L by exchanging i-th row with next row and i-th column with next column until the row to be deleted is the last, then deleting last row and column
    swapRowsL(ind,currSize_-1);
    swapRowsL(ind,currSize_-1,true);
    L_.reshape(currSize_-1,currSize_-1);

    // delete i-th item from z1 and z2
    // note: z1 and z2 are of size currSize_-dependent
    unsigned zsize = currSize_-dependent_;
    for( unsigned m=ind; m<zsize; m++ ){
      z1_[m] = z1_[m+1];
      z2_[m] = z2_[m+1];
    }
    z1_.resize(zsize-1);
    z2_.resize(zsize-1);

    // update base
    base_.erase(base_.begin()+ind);
    currSize_--; // update size

    // update kappa
    updateK();
            
    if(dependent_){
      // if some previously dependent item have become independent
      // recompute deltaLh
      Real ghNorm = GiGj(base_[currSize_-dependent_],base_[currSize_-dependent_]);
      Real lhnrm(0); // exact lhNorm
#if( EXACT)
      for (unsigned ii=0; ii<currSize_-dependent_; ++ii){
        lhnrm += L_(currSize_-dependent_,ii)*L_(currSize_-dependent_,ii);
      }
      deltaLh_ = std::abs(ghNorm - lhnrm);
#else
      Real sgn1 = sgn(deltaLh_);
      Real sgn2 = sgn(dletaLj);
      Real signum = sgn1 * sgn2; // sgn( deltaLh ) * sgn ( deltaLj );
      deltaLh_ = std::abs( ghNorm - lhNorm + deltaLh_ * deltaLh_);
#endif
              
      if( std::sqrt(deltaLh_) > tol*kappa_*std::max(static_cast<Real>(1),ghNorm) ){ // originally had deltaLh without sqrt 
        unsigned newind = currSize_-dependent_;
        dependent_--;
        // get the last row of L
        lh_.size(newind); // initialize to zeros;
        lhz1_ = zero;
        lhz2_ = zero;
        for (unsigned ii=0; ii<newind; ++ii){
          lh_[ii] = L_(newind,ii);
          lhz1_ += lh_[ii]*z1_[ii];
          lhz2_ += lh_[ii]*z2_[ii];
        }
        deltaLh_ = std::sqrt(deltaLh_);
        addSubgradToBase(newind,deltaLh_);                

        if(dependent_){ // dependent was 2
#if( ! EXACT )
          Real gjNorm = GiGj(base_[currSize_-1],base_[currSize_-1]);
          ljNorm     -= deltaLj_ * deltaLj_;
          lhNorm      = gjNorm;
          deltaLj_    = std::abs(gjNorm - ljNorm);
          if ( signum < 0 )
            deltaLj_ = - std::sqrt( deltaLj_ );
          else
            deltaLj_ = std::sqrt( deltaLj_ );
#else
          // recompute deltaLj
          Real gjTgh = GiGj(base_[currSize_-1],base_[currSize_-2]);
          Real ljTlh = 0.0;
          for (unsigned ii=0;ii<currSize_;ii++){
            ljTlh += L_(currSize_-1,ii)*L_(currSize_-2,ii);
          }
          deltaLj_ = (gjTgh - ljTlh) / deltaLh_;
#endif
          L_(currSize_-1,currSize_-2) = deltaLj_;
        }
      } // end if deltaLh > 0

      if (dependent_ > 1){ // deltaLh is 0 but deltaLj is not
        // recompute deltaLj
        Real gjNorm = GiGj(base_[currSize_-1],base_[currSize_-1]);
        Real ljnrm = zero; // exact ljNorm
#if( EXACT )
        for (unsigned ii=0; ii<currSize_; ++ii) {
          ljnrm += L_(currSize_-1,ii)*L_(currSize_-1,ii);
        }
        deltaLj_ = std::abs(gjNorm - ljnrm);
#else
        deltaLj_ = std::abs( gjNorm - ljNorm + deltaLj_ * deltaLj_);
#endif
                
        if( std::sqrt(deltaLj_) > tol*kappa_*std::max(static_cast<Real>(1),gjNorm) ){ // originally had deltaLj without sqrt
          unsigned newind = currSize_-1;
          dependent_--;
          // get the last row of L
          lj_.size(newind-1); // initialize to zeros;
          Real ljz1_ = zero;
          Real ljTz2 = zero;
          for (unsigned ii=0;ii<newind-1;ii++){
            lj_[ii] = L_(newind,ii);
            ljz1_ += lj_[ii]*z1_[ii];
            ljTz2 += lj_[ii]*z2_[ii];
          }
          deltaLj_ = std::sqrt(deltaLj_);
          addSubgradToBase(newind,deltaLj_);        
#if( EXACT )
          deltaLh_ = GiGj(base_[currSize_-2],base_[currSize_-1]);
          for (unsigned ii=0;ii<currSize_-1;ii++){
            deltaLh_ -= L_(currSize_-2,ii)*L_(currSize_-1,ii);
          }
          deltaLh_ /= deltaLj_;
#else
          if ( signum < 0) {
            deltaLh_ = - std::sqrt( deltaLh_ );
          }
          else {
            deltaLh_ = std::sqrt ( deltaLh_ );
          }
#endif
          L_(currSize_-1,currSize_-2) = deltaLh_;
        } // end if deltaLj > 0
      } // end if ( dependent > 1 )
    } // end if(dependent)
  }// end deleteSubgradFromBase()
  
  // TT: solving triangular system for TT algorithm
  void solveSystem(int size, char tran, LA::Matrix<Real> &L, LA::Vector<Real> &v){
    int info;
    if( L.numRows()!=size )
      std::cout << "Error: Wrong size matrix!" << std::endl;
    else if( v.numRows()!=size )
      std::cout << "Error: Wrong size vector!" << std::endl;
    else if( size==0 )
      return;
    else{
      //std::cout << L_.stride() << ' ' << size << std::endl;
      lapack_.TRTRS( 'L', tran, 'N', size, 1, L.values(), L.stride(), v.values(), v.stride(), &info );
    }
  }

  // TT: check that inequality constraints are satisfied for dual variables
  bool isFeasible(LA::Vector<Real> &v, const Real &tol){
    bool feas = true;
    for(int i=0;i<v.numRows();i++){
      if(v[i]<-tol){
        feas = false;
      }
    }
    return feas;
  }

  unsigned solveDual_TT(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    const Real zero(0), half(0.5), one(1);
    Real z1z2(0), z1z1(0);
    QPStatus_ = 1; // normal status
    entering_ = maxind_;

    // cold start
    optimal_   = true; 
    dependent_ = 0; 
    rho_       = ROL_INF<Real>(); // value of rho = -v
    currSize_  = 1;               // current base size
    base_.clear();
    base_.push_back(0);           // initial base
    L_.reshape(1,1);
    L_(0,0) = std::sqrt(GiGj(0,0));
    Bundle<Real>::resetDualVariables();
    Bundle<Real>::setDualVariable(0,one);
    tempv_.resize(1);
    tempw1_.resize(1);
    tempw2_.resize(1);
    lh_.resize(1);
    lj_.resize(1);
    z1_.resize(1); z2_.resize(1);
    z1_[0]     = one/L_(0,0);
    z2_[0]     = Bundle<Real>::alpha(0)/L_(0,0);
    LiMax_     = 0;                    // index of max element of diag(L)
    LiMin_     = 0;                    // index of min element of diag(L)
    kappa_     = one;                  // condition number of matrix L ( >= max|L_ii|/min|L_ii| )
    objval_    = ROL_INF<Real>();      // value of objective
    minobjval_ = ROL_INF<Real>();      // min value of objective (ever reached)
    
    unsigned iter;
    //-------------------------- MAIN LOOP --------------------------------//
    for (iter=0;iter<maxit;iter++){
      //---------------------- INNER LOOP -----------------------//
      while( !optimal_ ){
        switch( dependent_ ){
        case(0): // KT system admits unique solution
          {
            /*
              L = L_B'
            */
            z1z2    = z1_.dot(z2_);
            z1z1    = z1_.dot(z1_); 
            rho_    = ( one + z1z2/t )/z1z1;
            tempv_  = z1_; tempv_.scale(rho_);
            tempw1_ = z2_; tempw1_.scale(one/t);
            tempv_ -= tempw1_;
            solveSystem(currSize_,'T',L_,tempv_); // tempv contains solution
            optimal_ = true;
            break;
          }
        case(1): 
          {
            /*
              L = | L_B'   0 | \ currSize
                  | l_h^T  0 | /
            */
            LA::Matrix<Real> LBprime( LA::Copy,L_,currSize_-1,currSize_-1);     
            lh_.size(currSize_-1); // initialize to zeros;
            lhz1_ = zero;
            lhz2_ = zero;
            for(unsigned i=0; i<currSize_-1; ++i){
              Real tmp = L_(currSize_-1,i);
              lhz1_ += tmp*z1_(i);
              lhz2_ += tmp*z2_(i); 
              lh_[i] = tmp;
            }
            // Test for singularity
            if (std::abs(lhz1_-one) <= tol*kappa_){ 
              // tempv is an infinite direction
              tempv_ = lh_;  
              solveSystem(currSize_-1,'T',LBprime,tempv_);
              tempv_.resize(currSize_);   // add last entry
              tempv_[currSize_-1] = -one;
              optimal_ = false;
            }
            else{
              // system has (unique) solution
              rho_ = ( (Bundle<Real>::alpha(base_[currSize_-1]) - lhz2_)/t ) / ( one - lhz1_ );
              z1z2 = z1_.dot(z2_);
              z1z1 = z1_.dot(z1_); 
              Real tmp = ( one + z1z2 / t - rho_ * z1z1 )/( one - lhz1_ );
              tempw1_ = z1_; tempw1_.scale(rho_);
              tempw2_ = z2_; tempw2_.scale(one/t);
              tempw1_ -= tempw2_;
              tempw2_ = lh_; tempw2_.scale(tmp);
              tempw1_ -= tempw2_;
              solveSystem(currSize_-1,'T',LBprime,tempw1_);
              tempv_ = tempw1_;
              tempv_.resize(currSize_);
              tempv_[currSize_-1] = tmp;
              optimal_ = true;
            }
            break;
          } // case(1)
        case(2):
          {
            /*     | L_B'  0 0 | \
               L = | l_h^T 0 0 | | currSize
                   | l_j^T 0 0 | /
            */
            LA::Matrix<Real> LBprime( LA::Copy,L_,currSize_-2,currSize_-2 );
               lj_.size(currSize_-2); // initialize to zeros;
            lh_.size(currSize_-2); // initialize to zeros;
            ljz1_ = zero;
            lhz1_ = zero;
            for(unsigned i=0; i<currSize_-2; ++i){
              Real tmp1 = L_(currSize_-1,i);
              Real tmp2 = L_(currSize_-2,i);
              ljz1_ += tmp1*z1_(i);
              lhz1_ += tmp2*z1_(i); 
              lj_[i] = tmp1;
              lh_[i] = tmp2;
            }
            if(std::abs(ljz1_-one) <= tol*kappa_){ 
              // tempv is an infinite direction
              tempv_ = lj_;    
              solveSystem(currSize_-2,'T',LBprime,tempv_);
              tempv_.resize(currSize_);   // add two last entries
              tempv_[currSize_-2] = zero;
              tempv_[currSize_-1] = -one;
            }
            else{
              // tempv is an infinite direction
              Real mu = ( one - lhz1_ )/( one - ljz1_ );
              tempw1_ = lj_; tempw1_.scale(-mu);
              tempw1_ += lh_;
              solveSystem(currSize_-2,'T',LBprime,tempw1_);
              tempv_ = tempw1_;
              tempv_.resize(currSize_);
              tempv_[currSize_-2] = -one;
              tempv_[currSize_-1] = mu;
            }
            optimal_ = false;
          }// case(2)
        } // end switch(dependent_)

        // optimal is true if tempv is a solution, otherwise, tempv is an infinite direction
        if (optimal_){
          // check feasibility (inequality constraints)
          if (isFeasible(tempv_,tol*currSize_)){
            // set dual variables to values in tempv
            Bundle<Real>::resetDualVariables();
            for (unsigned i=0; i<currSize_; ++i){
              Bundle<Real>::setDualVariable(base_[i],tempv_[i]); 
            }
          }
          else{
            // w_B = /bar{x}_B - x_B
            for (unsigned i=0; i<currSize_; ++i){
              tempv_[i] -= Bundle<Real>::getDualVariable(base_[i]);
            }
            optimal_ = false;
          }
        } // if(optimal)
        else{ // choose the direction of descent
          if ( ( entering_ < maxind_ ) && ( Bundle<Real>::getDualVariable(entering_) == zero ) ){
            if ( tempv_[currSize_-1] < zero ) // w_h < 0
              tempv_.scale(-one);
          }
          else{ // no i such that dualVariables_[i] == 0
            Real sp(0);
            for( unsigned kk=0; kk<currSize_; ++kk)
              sp += tempv_[kk]*Bundle<Real>::alpha(base_[kk]);
            if ( sp > zero )
              tempv_.scale(-one);
          }
        }// end else ( not optimal )

        if(!optimal_){
          // take a step in direction tempv (possibly infinite)
          Real myeps = tol*currSize_;
          Real step  = ROL_INF<Real>();
          for (unsigned h=0; h<currSize_; ++h){
            if ( (tempv_[h] < -myeps) && (-Bundle<Real>::getDualVariable(base_[h])/tempv_[h] < step) )
              if ( (Bundle<Real>::getDualVariable(base_[h]) > myeps)
                || (Bundle<Real>::getDualVariable(base_[h]) < -myeps) ) // otherwise, consider it 0
                step = -Bundle<Real>::getDualVariable(base_[h])/tempv_[h];
          }

          if (step <= zero || step == ROL_INF<Real>()){
            QPStatus_ = -1; // invalid step
            return iter;
          }
          
          for (unsigned i=0; i<currSize_; ++i)
            Bundle<Real>::setDualVariable(base_[i],Bundle<Real>::getDualVariable(base_[i]) + step * tempv_[i]);          
        }// if(!optimal)
        
        //------------------------- ITEMS ELIMINATION ---------------------------//
        
        // Eliminate items with 0 multipliers from base
        bool deleted = optimal_;
        for (unsigned baseitem=0; baseitem<currSize_; ++baseitem){
          if (Bundle<Real>::getDualVariable(base_[baseitem]) <= tol){
            deleted = true;
            
#if( TABOO_LIST )
            // item that just entered shouldn't exit; if it does, mark it as taboo
            if( base_[baseitem] == entering_ ){
              taboo_.push_back(entering_);
              entering_ = maxind_;
            }
#endif
            
            // eliminate item from base; 
            deleteSubgradFromBase(baseitem,tol);

          } // end if(dualVariables_[baseitem] < tol)
        } // end loop over baseitem 
          
        if(!deleted){ // nothing deleted and not optimal
          QPStatus_ = -2; // loop
          return iter;
        }
      } // end inner loop
      
      Real newobjval(0), Lin(0), Quad(0); // new objective value
      for (unsigned i=0; i<currSize_; ++i){
        Lin += Bundle<Real>::alpha(base_[i])*Bundle<Real>::getDualVariable(base_[i]);
      }
      if (rho_ == ROL_NINF<Real>()){
        Quad = -Lin/t;
        newobjval = -half*Quad;
      }
      else{
        Quad = rho_ - Lin/t;
        newobjval = half*(rho_ + Lin/t);
      }

#if( TABOO_LIST )
      // -- test for strict decrease -- // 
      // if item didn't provide decrease, move it to taboo list ...
      if( ( entering_ < maxind_ ) && ( objval_ < ROL_INF<Real>() ) ){
        if( newobjval >= objval_ - std::max( tol*std::abs(objval_), ROL_EPSILON<Real>() ) ){
          taboo_.push_back(entering_);
        } 
      }
#endif

      objval_ = newobjval;

      // if sufficient decrease obtained
      if ( objval_ + std::max( tol*std::abs(objval_), ROL_EPSILON<Real>() ) <= minobjval_ ){
        taboo_.clear(); // reset taboo list
        minobjval_ = objval_;
      }

      //---------------------- OPTIMALITY TEST -------------------------//
      if ( (rho_ >= ROL_NINF<Real>()) && (objval_ <= ROL_NINF<Real>()) ) // if current x (dualVariables_) is feasible
        break;
        
      entering_  = maxind_;
      Real minro = - std::max( tol*currSize_*std::abs(objval_), ROL_EPSILON<Real>() );
#if ( ! FIRST_VIOLATED )
      Real diff  = ROL_NINF<Real>(), olddiff = ROL_NINF<Real>();
#endif

      for (unsigned bundleitem=0; bundleitem<Bundle<Real>::size(); ++bundleitem){ // loop over items in bundle
      //for (int bundleitem=size_-1;bundleitem>-1;bundleitem--){ // loop over items in bundle (in reverse order)
        if ( std::find(taboo_.begin(),taboo_.end(),bundleitem) != taboo_.end() ){
          continue; // if item is taboo, move on
        }

        if ( std::find(base_.begin(),base_.end(),bundleitem) == base_.end() ){
          // base does not contain bundleitem
          Real ro = zero;
          for (unsigned j=0;j<currSize_;j++){
            ro += Bundle<Real>::getDualVariable(base_[j]) * GiGj(bundleitem,base_[j]);
          }
          ro += Bundle<Real>::alpha(bundleitem)/t;


          if (rho_ >= ROL_NINF<Real>()){
            ro = ro - rho_; // note: rho = -v 
          }
          else{
            ro         = ROL_NINF<Real>();
            minobjval_ = ROL_INF<Real>();
            objval_    = ROL_INF<Real>();
          }
          
          if (ro < minro){                       
#if ( FIRST_VIOLATED )
            entering_ = bundleitem;
            break; // skip going through rest of constraints; alternatively, could look for "most violated"
#else
            diff = minro - ro;
            if ( diff > olddiff ){
              entering_ = bundleitem;
              olddiff = diff;
            }
#endif
          }
          
        } // end if item not in base
      }// end of loop over items in bundle

      //----------------- INSERTING ITEM ------------------------// 
      if (entering_ < maxind_){ // dual constraint is violated
        optimal_ = false;
        Bundle<Real>::setDualVariable(entering_,zero);
        base_.push_back(entering_);
        // construct new row of L_B
        unsigned zsize = currSize_ - dependent_; // zsize is the size of L_Bprime (current one)
        lh_.size(zsize); // initialize to zeros;
        lhz1_ = zero;
        lhz2_ = zero;
        for (unsigned i=0; i<zsize; ++i){
          lh_[i] = GiGj(entering_,base_[i]);
        }
        LA::Matrix<Real> LBprime( LA::Copy,L_,zsize,zsize);      
        solveSystem(zsize,'N',LBprime,lh_); // lh = (L_B^{-1})*(G_B^T*g_h)
        for (unsigned i=0; i<zsize; ++i){
          lhz1_ += lh_[i]*z1_[i];
          lhz2_ += lh_[i]*z2_[i];
        }

        Real nrm = lh_.dot(lh_);
        Real delta = GiGj(entering_,entering_) - nrm; // delta squared
#if( ! EXACT )
        if(dependent_)
          ljNorm = nrm; // adding second dependent
        else
          lhNorm = nrm; // adding first dependent
#endif

        currSize_++; // update base size
        
        L_.reshape(currSize_,currSize_);
        zsize = currSize_ - dependent_; // zsize is the size of L_Bprime (new one)
        for (unsigned i=0; i<zsize-1; ++i){
          L_(currSize_-1,i) = lh_[i];
        }

        Real deltaeps = tol*kappa_*std::max(one,lh_.dot(lh_)); 
        if ( delta > deltaeps ){ // new row is independent
          // add subgradient to the base
          unsigned ind = currSize_-1;
          delta = std::sqrt(delta);
          addSubgradToBase(ind,delta);
        }
        else if(delta < -deltaeps){
          dependent_++;
          QPStatus_ = 0; // negative delta
          return iter;
        }
        else{ // delta zero
          dependent_++;
        }
      } // end if(entering_ < maxind_)
      else{ // no dual constraint violated
        if( objval_ - std::max( tol*std::abs( objval_ ), ROL_EPSILON<Real>() ) > minobjval_ ){ // check if f is as good as minf
          QPStatus_ = -3; // loop
          return iter;
        }
      }

      if(optimal_)
        break; 
    } // end main loop

    taboo_.clear();
    return iter;
  }// end solveDual_TT()

  unsigned solveDual_arbitrary(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    Real mytol = tol;
    unsigned outermaxit = 20;
    bool increase = false, decrease = false;
    unsigned iter = 0;
    for ( unsigned it=0; it < outermaxit; ++it ){
      iter += solveDual_TT(t,maxit,mytol);
      if ( QPStatus_ == 1 ) {
        break;
      }
      else if ( QPStatus_ == -2  || QPStatus_ == -3 ) {
        mytol /= static_cast<Real>(10);
        decrease = true;
      }
      else {
        mytol *= static_cast<Real>(10);
        increase = true;
      }
      if ( (mytol > static_cast<Real>(1e-4))
        || (mytol < static_cast<Real>(1e-16)) ){
        break;
      }
      if ( increase && decrease ) {
        break;
      }
    }// end outer loop
    return iter;
  }

}; // class Bundle_TT

} // namespace ROL

#endif

