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

#include "Teuchos_RCP.hpp"

#include <vector>
#include <limits.h> 
#include <stdint.h> 
#include <float.h> 
#include <math.h> 
#include <algorithm> // TT: std::find

#include "Teuchos_SerialDenseMatrix.hpp" 
#include "Teuchos_SerialDenseVector.hpp" 
#include "Teuchos_LAPACK.hpp" 

#define DEBUG_TT 0
#define EXACT 1
#define TABOO_LIST 1
#define FIRST_VIOLATED 0

/** \class ROL::Bundle_TT
    \brief Provides the interface for and implements a bundle. The semidefinite quadratic subproblem is solved using TT algorithm by Antonio Frangioni (1996).
*/

namespace ROL {

template<class Real>
class Bundle_TT : public Bundle<Real> {

private: 
  unsigned maxSize_;
  Teuchos::LAPACK<int, Real> lapack_; // TT

public:
  Bundle_TT (const unsigned maxSize = 10, const Real coeff = 0.0, const unsigned remSize = 2) 
    : Bundle<Real>(maxSize,coeff,remSize), maxSize_(maxSize) {
    INF = std::numeric_limits<double>::max();
    maxind = std::numeric_limits<int>::max();
    Id.reshape(maxSize_,maxSize_);
    for(unsigned i=0;i<maxSize_;i++)
      Id(i,i) = 1.0;
  }

  Real GiTGj(const int i, const int j){
    return (this->subgradient(i)).dot(this->subgradient(j));
  }

/***********************************************************************************************/
/****************** DUAL CUTTING PLANE SUBPROBLEM ROUTINES *************************************/
/***********************************************************************************************/
private:
  int QPStatus_;
  Real INF;
  int maxind;
  int entering_; // index of entering item
  std::vector<int> taboo_; // list of "taboo" items
  bool optimal_; // flag for optimality of restricted solution
  unsigned dependent_; // number of lin. dependent items in base
  Real rho_;
  unsigned currSize_; // current size of base
  std::vector<int> base; // base
  Teuchos::SerialDenseMatrix<int, Real> L;
  Teuchos::SerialDenseMatrix<int, Real> Id;
  Teuchos::SerialDenseVector<int, Real> tempv;
  Teuchos::SerialDenseVector<int, Real> tempw1;
  Teuchos::SerialDenseVector<int, Real> tempw2;
  Teuchos::SerialDenseVector<int, Real> lh;
  Teuchos::SerialDenseVector<int, Real> lj;
  Teuchos::SerialDenseVector<int, Real> z1;
  Teuchos::SerialDenseVector<int, Real> z2;
  Real lhNorm, ljNorm, z1Tz2, z1Tz1, lhTz1, lhTz2, ljTz1;
  int LiMax; // index of max element of diag(L)
  int LiMin; // index of min element of diag(L)
  Real kappa; // condition number of matrix L ( >= max|L_ii|/min|L_ii| )
  Real Quad, Lin; // quadratic and linear terms of objective
  Real objval; // value of objective
  Real minobjval; // min value of objective (ever reached)
  Real deltaLh, deltaLj; // needed in case dependent row becomes independent

  void swapRowsL(unsigned ind1, unsigned ind2, bool trans=false){
    if( ind1 > ind2){
      unsigned tmp = ind1;
      ind2 = ind1;
      ind1 = tmp;
    }
    unsigned dd = ind1;
    for (unsigned n=ind1+1; n<=ind2; n++){
      Teuchos::SerialDenseMatrix<int, Real> Id_n(Teuchos::Copy,Id,currSize_,currSize_);
      Id_n(dd,dd) = 0; Id_n(dd,n) = 1.0; 
      Id_n(n,dd) = 1.0; Id_n(n,n) = 0;
      Teuchos::SerialDenseMatrix<int, Real> prod(currSize_,currSize_);
      if( !trans )
	prod.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,Id_n,L,0.0);
      else
	prod.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,L,Id_n,0.0);
      L = prod;
      dd++;
    }
  }

  void updateK(){
    if (currSize_ <= dependent_) { // L is empty
      kappa = 1.0;
    }
    else{
      Real tmpdiagMax = -INF;
      Real tmpdiagMin = INF;
      for (unsigned j=0;j<currSize_-dependent_;j++){
	if( L(j,j) > tmpdiagMax ){
	  tmpdiagMax = L(j,j);
	  LiMax = j;
	}
	if( L(j,j) < tmpdiagMin ){
	  tmpdiagMin = L(j,j);
	  LiMin = j;
	}
      }
      kappa = tmpdiagMax/tmpdiagMin;
    }
  }

  void addSubgradToBase(unsigned ind, Real delta){
    // update z1, z2, kappa
    // swap rows if: dependent == 1 and we insert independent row (dependent row is always last)
    //               dependent == 2 and Lj has become independent (Lh still dependent)
    if(dependent_ && (ind == currSize_-1)){
	swapRowsL(currSize_-2,currSize_-1);
	int tmp;
	tmp = base[currSize_-2];
	base[currSize_-2] = base[currSize_-1];
	base[currSize_-1] = tmp;
	ind--;
#if( DEBUG_TT )	  
	std::cout << "Swapped last two rows of L " << std::endl;
#endif
    } // end if dependent
    
    L(ind,ind) = delta;
    
    // update z1 and z2
    unsigned zsize = ind+1;
    z1.resize(zsize); z2.resize(zsize);
    z1[ind] = ( 1.0 - lhTz1 ) / delta;
    z2[ind] = ( this->alpha(base[ind]) - lhTz2 ) / delta;  
    //z2[zsize-1] = ( this->alpha(entering_) - lhTz2 ) / delta;  
    
    // update kappa
    if(delta > L(LiMax,LiMax)){
      LiMax = ind;
      kappa = delta/L(LiMin,LiMin);
    }
    if(delta < L(LiMin,LiMin)){
      LiMin = ind;
      kappa = L(LiMax,LiMax)/delta;
    }
  }

  void deleteSubgradFromBase(unsigned ind, Real tol){
    // update L, currSize, base, z1, z2, dependent, dualVariables_, kappa
    if (ind >= currSize_-dependent_){
#if( DEBUG_TT )
      std::cout << "Eliminating dependent item " << base[ind] << std::endl;
#endif
      // if dependent > 0, the last one or two rows of L are lin. dependent	      
      if (ind < currSize_-1){ // eliminate currSize_-2 but keep currSize_-1
#if( DEBUG_TT )
	std::cout << "Eliminating Lh but keeping Lj" << std::endl;
#endif
	// swap the last row with the second to last
	swapRowsL(ind,currSize_-1);
	base[ind] = base[currSize_-1];		
#if( ! EXACT )
	lhNorm = ljNorm; // new last row is lh
#endif
      }
 
      dependent_--;
      currSize_--;
      L.reshape(currSize_,currSize_); // the row to be eliminated is the last row
      base.resize(currSize_);

#if( DEBUG_TT )
      std::cout << "New base = " << std::endl;
      for (unsigned kk=0;kk<base.size();kk++){
	std::cout << base[kk] << std::endl;
      }
      std::cout << "\n";
      std::cout << "New L " << std::endl;
      std::cout << L << std::endl;
      std::cout << "\n";
#endif
      // note: z1, z2, kappa need not be updated
      return;
    } // end if dependent item

#if( DEBUG_TT )
    std::cout << "Eliminating INdependent item " << base[baseitem] << std::endl;
#endif
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
    for (unsigned j=ind+1;j<currSize_-dependent_;j++){
      Real ai = L(j,ind);
      if (std::abs(ai) <= tol*currSize_) // nothing to do
	continue;
      Real aj = L(j,j);
      Real d, Gc, Gs;
      // Find Givens
      // Anderson's implementation
      if (std::abs(aj) <= tol*currSize_){ // aj is zero
	Gc = 0.0; d = std::abs(ai); Gs = (ai < 0) - (ai > 0); // Gs = -sgn(ai)
      }
      else if ( std::abs(ai) > std::abs(aj) ){
	Real t = aj/ai;
	Real sgnb = (ai > 0) - (ai < 0);
	Real u = sgnb * std::sqrt(1.0 + t*t );
	Gs = -1.0/u;
	Gc = -Gs*t;
	d  = ai*u;
      }
      else{
	Real t = ai/aj;
	Real sgna = (aj > 0) - (aj < 0);
	Real u = sgna * std::sqrt(1.0 + t*t );
	Gc = 1.0/u;
	Gs = -Gc*t;
	d  = aj*u;
      }

      // // "Naive" implementation
      // d  = hypot(ai,aj);
      // Gc = aj/d;
      // Gs = -ai/d;
	      
      L(j,j) = d; L(j,ind) = 0.0; 
      // apply Givens to columns i,j of L
      for (unsigned h=j+1;h<currSize_;h++){
	Real tmp1 = L(h,ind);
	Real tmp2 = L(h,j);
	L(h,ind) = Gc*tmp1 + Gs*tmp2;
	L(h,j) = Gc*tmp2 - Gs*tmp1;
      }
      // apply Givens to z1, z2
      Real tmp1 = z1[ind];
      Real tmp2 = z1[j];
      Real tmp3 = z2[ind];
      Real tmp4 = z2[j];
      z1[ind] = Gc*tmp1 + Gs*tmp2;
      z1[j] = Gc*tmp2 - Gs*tmp1;
      z2[ind] = Gc*tmp3 + Gs*tmp4;
      z2[j] = Gc*tmp4 - Gs*tmp3;
    }// end loop over j

#if( DEBUG_TT )
    std::cout << "After Givens: L,z1,z2 " << std::endl;
    std::cout << L << std::endl;
    std::cout << z1 << std::endl;
    std::cout << z2 << std::endl;
#endif

    if(dependent_){
      deltaLh = L(currSize_-dependent_,ind);  // h = currSize_ - dependent
      if( dependent_ > 1 )                 // j = currSize_ - 1, h = currSize_ - 2
	deltaLj = L(currSize_-1,ind);
    }
	    
    // shift rows and columns of L by exchanging i-th row with next row and i-th column with next column until the row to be deleted is the last, then deleting last row and column
    swapRowsL(ind,currSize_-1);
    swapRowsL(ind,currSize_-1,true);
    L.reshape(currSize_-1,currSize_-1);

    // delete i-th item from z1 and z2
    // note: z1 and z2 are of size currSize_-dependent
    unsigned zsize = currSize_-dependent_;
    for( unsigned m=ind; m<zsize; m++ ){
      z1[m] = z1[m+1];
      z2[m] = z2[m+1];
    }
    z1.resize(zsize-1);
    z2.resize(zsize-1);

#if( DEBUG_TT )
    std::cout << "After elimination: L,z1,z2 " << std::endl;
    std::cout << L << std::endl;
    std::cout << z1 << std::endl;
    std::cout << z2 << std::endl;
#endif

    // update base
    base.erase(base.begin()+ind);

#if( DEBUG_TT )
    std::cout << "New base = " << std::endl;
    for (unsigned kk=0;kk<base.size();kk++){
      std::cout << base[kk] << std::endl;
    }
#endif
	    
    currSize_--; // update size

    // update kappa
    updateK();
	    
    if(dependent_){
      // if some previously dependent item have become independent
#if( DEBUG_TT )
      std::cout << "deltaLh = " << deltaLh << std::endl;
#endif
      // recompute deltaLh
      Real ghNorm = GiTGj(base[currSize_-dependent_],base[currSize_-dependent_]);
      Real lhnrm = 0.0; // exact lhNorm
#if( EXACT)
      for (unsigned ii=0;ii<currSize_-dependent_;ii++){
	lhnrm += L(currSize_-dependent_,ii)*L(currSize_-dependent_,ii);
      }
      deltaLh = std::abs(ghNorm - lhnrm);
#else
      Real sgn1 = (deltaLh > 0) ? 1 : ((deltaLh < 0) ? -1 : 0);
      Real sgn2 = (deltaLj > 0) ? 1 : ((deltaLj < 0) ? -1 : 0);
      Real signum = sgn1 * sgn2; // sgn( deltaLh ) * sgn ( deltaLj );
      deltaLh = std::abs( ghNorm - lhNorm + deltaLh * deltaLh);
#endif
	      
#if( DEBUG_TT )
      std::cout << "ghNorm = " << ghNorm << std::endl;
      std::cout << "lhNorm (exact) = " << lhnrm << std::endl;
      std::cout << "lhNorm = " << lhNorm << std::endl;
      std::cout << "deltaLh = " << std::sqrt(deltaLh) << std::endl;
      std::cout << "kappa = " << kappa << std::endl;
#endif

      if( std::sqrt(deltaLh) > tol*kappa*std::max(1.0,ghNorm) ){ // originally had just deltaLh (without sqrt) 
#if( DEBUG_TT )
	std::cout << "Lh has become lin. INdependent" << std::endl;
#endif
	unsigned newind = currSize_-dependent_;
	dependent_--;
	// get the last row of L
	lh.size(newind); // initialize to zeros;
	lhTz1 = 0.0;
	lhTz2 = 0.0;
	for (unsigned ii=0;ii<newind;ii++){
	  lh[ii] = L(newind,ii);
	  lhTz1 += lh[ii]*z1[ii];
	  lhTz2 += lh[ii]*z2[ii];
	}
	deltaLh = std::sqrt(deltaLh);
	addSubgradToBase(newind,deltaLh);		

	if(dependent_){ // dependent was 2
#if( ! EXACT )
	  Real gjNorm = GiTGj(base[currSize_-1],base[currSize_-1]);
	  ljNorm -= deltaLj * deltaLj;
	  lhNorm = gjNorm;
		  
	  deltaLj = std::abs(gjNorm - ljNorm);
	  if ( signum < 0 )
	    deltaLj = - std::sqrt( deltaLj );
	  else
	    deltaLj = std::sqrt( deltaLj );
#else
	  // recompute deltaLj
	  Real gjTgh = GiTGj(base[currSize_-1],base[currSize_-2]);
	  Real ljTlh = 0.0;
	  for (unsigned ii=0;ii<currSize_;ii++){
	    ljTlh += L(currSize_-1,ii)*L(currSize_-2,ii);
	  }
	  deltaLj = (gjTgh - ljTlh) / deltaLh;
#endif

	  L(currSize_-1,currSize_-2) = deltaLj;
	}

#if( DEBUG_TT )
	std::cout << "Updated L, z1, z2: " << std::endl;
	std::cout << L  << std::endl;
	std::cout << z1 << std::endl;
	std::cout << z2 << std::endl;
	std::cout << "kappa = " << kappa << std::endl;
#endif
		
      } // end if deltaLh > 0

      if (dependent_ > 1){ // deltaLh is 0 but deltaLj is not
	// recompute deltaLj
	Real gjNorm = GiTGj(base[currSize_-1],base[currSize_-1]);
	Real ljnrm = 0.0; // exact ljNorm
#if( EXACT )
	for (unsigned ii=0;ii<currSize_;ii++){
	  ljnrm += L(currSize_-1,ii)*L(currSize_-1,ii);
	}
	deltaLj = std::abs(gjNorm - ljnrm);
#else
	deltaLj = std::abs( gjNorm - ljNorm + deltaLj * deltaLj);
#endif
		
#if( DEBUG_TT )
	std::cout << "gjNorm = " << gjNorm << std::endl;
	std::cout << "ljNorm (exact) = " << ljnrm << std::endl;
	std::cout << "ljNorm = " << ljNorm << std::endl; 
	std::cout << "deltaLj = " << std::sqrt(deltaLj) << std::endl;
	std::cout << "kappa = " << kappa << std::endl;
#endif
		
	if( std::sqrt(deltaLj) > tol*kappa*std::max(1.0,gjNorm) ){ // originally just had deltaLj (without sqrt)
#if( DEBUG_TT )
	  std::cout << "Lj has become lin. INdependent" << std::endl;
#endif
	  unsigned newind = currSize_-1;
	  dependent_--;
	  // get the last row of L
	  lj.size(newind-1); // initialize to zeros;
	  Real ljTz1 = 0.0;
	  Real ljTz2 = 0.0;
	  for (unsigned ii=0;ii<newind-1;ii++){
	    lj[ii] = L(newind,ii);
	    ljTz1 += lj[ii]*z1[ii];
	    ljTz2 += lj[ii]*z2[ii];
	  }
	  deltaLj = std::sqrt(deltaLj);
	  addSubgradToBase(newind,deltaLj);	
#if( EXACT )
	  deltaLh = GiTGj(base[currSize_-2],base[currSize_-1]);
	  for (unsigned ii=0;ii<currSize_-1;ii++){
	    deltaLh -= L(currSize_-2,ii)*L(currSize_-1,ii);
	  }
	  deltaLh /= deltaLj;
#else
	  if ( signum < 0)
	    deltaLh = - std::sqrt( deltaLh );
	  else
	    deltaLh = std::sqrt ( deltaLh );
#endif
		  
	  L(currSize_-1,currSize_-2) = deltaLh;
		  
	} // end if deltaLj > 0
      } // end if ( dependent > 1 )

    } // end if(dependent)

  }// end deleteSubgradFromBase()
  
  Real evaluateObjective(std::vector<Real> &g, const std::vector<Real> &x, const Real t) const {
    this->gx_->zero(); this->eG_->zero();
    for (unsigned i = 0; i < this->size(); i++) {
      // Compute Gx using Kahan's compensated sum
      this->tG_->set(*this->gx_);
      this->yG_->set(*this->eG_); this->yG_->axpy(x[i],this->subgradient(i));
      this->gx_->set(*this->tG_); this->gx_->plus(*this->yG_);
      this->eG_->set(*this->tG_); this->eG_->axpy(-1.0,*this->gx_); this->eG_->plus(*this->yG_);
    }
    Real Hx = 0.0, val = 0.0, err = 0.0, tmp = 0.0, y = 0.0;
    for (unsigned i = 0; i < this->size(); i++) {
      // Compute < g_i, Gx >
      Hx   = this->gx_->dot(this->subgradient(i));
      // Add to the objective function value using Kahan's compensated sum
      tmp  = val;
      y    = x[i]*(0.5*Hx + this->alpha(i)/t) + err;
      val  = tmp + y;
      err  = (tmp - val) + y;
      // Add gradient component
      g[i] = Hx + this->alpha(i)/t;
    }
    return val;
  }

  unsigned solveDual_dim1(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    this->setDualVariables(0,1.0);
    return 0;
  }

  unsigned solveDual_dim2(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    this->gx_->set(this->subgradient(0)); 
    //gx_ = this->subgradient(0).clone();
    this->gx_->axpy(-1.0,this->subgradient(1));
    Real diffg  = this->gx_->dot(*this->gx_);
    if ( std::abs(diffg) > ROL_EPSILON<Real>() ) {
      Real diffa  = (this->alpha(0)-this->alpha(1))/t;
      Real gdiffg = this->subgradient(1).dot(*this->gx_);
      this->setDualVariables(0,std::min(1.0,std::max(0.0,-(gdiffg+diffa)/diffg)));
      this->setDualVariables(1,1.0 - this->getDualVariables(0));
    }
    else {
      if ( std::abs(this->alpha(0)-this->alpha(1)) > ROL_EPSILON<Real>() ) {
        if ( this->alpha(0) < this->alpha(1) ) {
          this->setDualVariables(0,1.0); this->setDualVariables(1,0.0);
        }
        else if ( this->alpha(0) > this->alpha(1) ) {
          this->setDualVariables(0,0.0); this->setDualVariables(1,1.0);
        }
      }
      else {
        this->setDualVariables(0,0.5); this->setDualVariables(1,0.5);
      }
    }
    return 0;
  }
  
  // TT: solving triangular system for TT algorithm
  void solveSystem(int size, char tran, Teuchos::SerialDenseMatrix<int,Real> &L, Teuchos::SerialDenseVector<int,Real> &v){
    int info;
    if( L.numRows()!=size )
      std::cout << "Error: Wrong size matrix!" << std::endl;
    else if( v.numRows()!=size )
      std::cout << "Error: Wrong size vector!" << std::endl;
    else if( size==0 )
      return;
    else{
      //std::cout << L.stride() << ' ' << size << std::endl;
      lapack_.TRTRS( 'L', tran, 'N', size, 1, L.values(), L.stride(), v.values(), v.stride(), &info );
    }
  }

  // TT: check that inequality constraints are satisfied for dual variables
  bool isFeasible(Teuchos::SerialDenseVector<int,Real> &v, const Real &tol){
    bool feas = true;
    for(int i=0;i<v.numRows();i++){
      if(v[i]<-tol){
	feas = false;
      }
    }
    return feas;
  }

  unsigned solveDual_TT(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
#if( DEBUG_TT )
    std::cout << "Calling solveDual_TT" << std::endl;
    std::cout << "t = " << t << std::endl;
    std::cout << "maxit = " << maxit << std::endl;
    std::cout << "tol = " << tol << std::endl;
#endif
    QPStatus_ = 1; // normal status
    entering_ = maxind;

    // cold start
    optimal_ = true; 
    dependent_ = 0; 
    rho_ = INF; // value of rho = -v
    currSize_ = 1; // current base size
    base.clear();
    base.push_back(0); // initial base
    L.reshape(1,1);
    L(0,0) = std::sqrt(GiTGj(0,0));
    this->resetDualVariables();
    this->setDualVariables(0,1.0);
    tempv.resize(1);
    tempw1.resize(1);
    tempw2.resize(1);
    lh.resize(1);
    lj.resize(1);
    z1.resize(1); z2.resize(1);
    z1[0] = 1.0/L(0,0);
    z2[0] = this->alpha(0)/L(0,0);
    LiMax = 0; // index of max element of diag(L)
    LiMin = 0; // index of min element of diag(L)
    kappa = 1.0; // condition number of matrix L ( >= max|L_ii|/min|L_ii| )
    objval = INF; // value of objective
    minobjval = INF; // min value of objective (ever reached)
    
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
	    z1Tz2 = z1.dot(z2);
	    z1Tz1 = z1.dot(z1); 
	    rho_ = ( 1 + z1Tz2/t )/z1Tz1;
	    tempv = z1; tempv.scale(rho_);
	    tempw1 = z2; tempw1.scale(1/t);
	    tempv -= tempw1;
	    solveSystem(currSize_,'T',L,tempv); // tempv contains solution
	    optimal_ = true;
#if( DEBUG_TT )
	    std::cout << "In case 0" << std::endl;
	    std::cout << "rho_ = " << rho_ << std::endl;
	    std::cout << "Solution tempv = \n" << tempv << std::endl;
#endif
	    break;
	  }
	case(1): 
	  {
	    /*
	      L = | L_B'   0 | \ currSize
	          | l_h^T  0 | /
	    */
	    Teuchos::SerialDenseMatrix<int,Real> LBprime( Teuchos::Copy,L,currSize_-1,currSize_-1);     
	    lh.size(currSize_-1); // initialize to zeros;
	    lhTz1 = 0.0;
	    lhTz2 = 0.0;
	    for(unsigned i=0;i<currSize_-1;i++){
	      Real tmp = L(currSize_-1,i);
	      lhTz1 += tmp*z1(i);
	      lhTz2 += tmp*z2(i); 
	      lh[i] = tmp;
	    }
#if( DEBUG_TT )
	    bool singular = false;
#endif
	    // Test for singularity
	    if (std::abs(lhTz1-1.0) <= tol*kappa){ 
	      // tempv is an infinite direction
#if( DEBUG_TT )
	      singular = true;
#endif
	      tempv = lh;  
	      solveSystem(currSize_-1,'T',LBprime,tempv);
	      tempv.resize(currSize_);   // add last entry
	      tempv[currSize_-1] = -1.0;
	      optimal_ = false;
	    }
	    else{
	      // system has (unique) solution
	      rho_ = ( (this->alpha(base[currSize_-1]) - lhTz2)/t ) / ( 1.0 - lhTz1 );
	      z1Tz2 = z1.dot(z2);
	      z1Tz1 = z1.dot(z1); 
	      Real tmp = ( 1.0 + z1Tz2 / t - rho_ * z1Tz1 )/( 1.0 - lhTz1 );
	      tempw1 = z1; tempw1.scale(rho_);
	      tempw2 = z2; tempw2.scale(1.0/t);
	      tempw1 -= tempw2;
	      tempw2 = lh; tempw2.scale(tmp);
	      tempw1 -= tempw2;
	      solveSystem(currSize_-1,'T',LBprime,tempw1);
	      tempv = tempw1;
	      tempv.resize(currSize_);
	      tempv[currSize_-1] = tmp;
	      optimal_ = true;
	    }
#if( DEBUG_TT )
	    std::cout << "In case 1" << std::endl;
	    if (!singular){
	      std::cout << "rho_ = " << rho_ << std::endl;
	      std::cout << "Solution tempv = \n" << tempv << std::endl;
	    }
	    else
	      std::cout << "Direction tempv = \n" << tempv << std::endl;
#endif
	    break;
	  } // case(1)
	case(2):
	  {
	    /*     | L_B'  0 0 | \
	       L = | l_h^T 0 0 | | currSize
	           | l_j^T 0 0 | /
	    */
	    Teuchos::SerialDenseMatrix<int,Real> LBprime( Teuchos::Copy,L,currSize_-2,currSize_-2 );
   	    lj.size(currSize_-2); // initialize to zeros;
	    lh.size(currSize_-2); // initialize to zeros;
	    ljTz1 = 0.0;
	    lhTz1 = 0.0;
	    for(unsigned i=0;i<currSize_-2;i++){
	      Real tmp1 = L(currSize_-1,i);
	      Real tmp2 = L(currSize_-2,i);
	      ljTz1 += tmp1*z1(i);
	      lhTz1 += tmp2*z1(i); 
	      lj[i] = tmp1;
	      lh[i] = tmp2;
	    }
	    if(std::abs(ljTz1-1.0) <= tol*kappa){ 
	      // tempv is an infinite direction
	      tempv = lj;    
	      solveSystem(currSize_-2,'T',LBprime,tempv);
	      tempv.resize(currSize_);   // add two last entries
	      tempv[currSize_-2] = 0.0;
	      tempv[currSize_-1] = -1.0;
	    }
	    else{
	      // tempv is an infinite direction
	      Real mu = ( 1.0 - lhTz1 )/( 1.0 - ljTz1 );
	      tempw1 = lj; tempw1.scale(-mu);
	      tempw1 += lh;
	      solveSystem(currSize_-2,'T',LBprime,tempw1);
	      tempv = tempw1;
	      tempv.resize(currSize_);
	      tempv[currSize_-2] = -1.0;
	      tempv[currSize_-1] = mu;
	    }
	    optimal_ = false;
#if( DEBUG_TT )
	    std::cout << "In case 2" << std::endl;
	    std::cout << "Direction tempv = \n" << tempv << std::endl;
#endif
	  }// case(2)
	} // end switch(dependent_)

	// optimal is true if tempv is a solution, otherwise, tempv is an infinite direction
	if (optimal_){
	  // check feasibility (inequality constraints)
	  if (isFeasible(tempv,tol*currSize_)){
#if( DEBUG_TT )
	    std::cout << "Solution tempv is feasible" << std::endl;
#endif
	    // set dual variables to values in tempv
	    this->resetDualVariables();
	    for (unsigned i=0;i<currSize_;i++){
	      this->setDualVariables(base[i],tempv[i]); 
	    }
	  }
	  else{
#if( DEBUG_TT )
	    std::cout << "Solution tempv is NOT feasible" << std::endl;
#endif
	    // w_B = /bar{x}_B - x_B
	    for (unsigned i=0;i<currSize_;i++){
	      tempv[i] -= this->getDualVariables(base[i]);
	    }
	    optimal_ = false;
	  }
	} // if(optimal)
	else{ // choose the direction of descent
	  if ( ( entering_ < maxind ) && ( this->getDualVariables(entering_) == 0.0 ) ){
	    if ( tempv[currSize_-1] < 0 ) // w_h < 0
	      tempv.scale(-1.0);
	  }
	  else{ // no i such that dualVariables_[i] == 0
	    Real sp = 0.0;
	    for( unsigned kk=0;kk<currSize_;kk++)
	      sp += tempv[kk]*this->alpha(base[kk]);
	    if ( sp > 0 )
	      tempv.scale(-1.0);
	  }
	}// end else ( not optimal )

	if(!optimal_){
	  // take a step in direction tempv (possibly infinite)
	  Real myeps = tol*currSize_;
	  Real step = INF;
	  for (unsigned h=0;h<currSize_;h++){
	    if ( (tempv[h] < -myeps) && (-this->getDualVariables(base[h])/tempv[h] < step) )
	      if ( (this->getDualVariables(base[h]) > myeps) || (this->getDualVariables(base[h]) < -myeps) ) // otherwise, consider it 0
		step = -this->getDualVariables(base[h])/tempv[h];
#if( DEBUG_TT )	  	      
	    std::cout << "h = " << h << " tempv[h] = " << tempv[h] << " dualV[base[h]] = " << this->getDualVariables(base[h]) << std::endl;
#endif
	  }

#if( DEBUG_TT )	  
	  std::cout << "Taking step of size " << step << std::endl;
#endif

	  if (step <= 0 || step == INF){
#if( DEBUG_TT )
	    std::cout << "Invalid step!" << std::endl;
#endif
	    QPStatus_ = -1; // invalid step
	    return iter;
	  }
	  
	  for (unsigned i=0;i<currSize_;i++)
	    this->setDualVariables(base[i],this->getDualVariables(base[i]) + step * tempv[i]);	  
	}// if(!optimal)
	
	//------------------------- ITEMS ELIMINATION ---------------------------//
	
	// Eliminate items with 0 multipliers from base
	bool deleted = optimal_;
	for (unsigned baseitem=0;baseitem<currSize_;baseitem++){
	  if (this->getDualVariables(base[baseitem]) <= tol){
	    deleted = true;
	    
#if( TABOO_LIST )
	    // item that just entered shouldn't exit; if it does, mark it as taboo
	    if( base[baseitem] == entering_ ){
#if( DEBUG_TT )
	      std::cout << "Blocking " << entering_ << " because it just entered" << std::endl;
#endif
	      taboo_.push_back(entering_);
	      entering_ = maxind;
	    }
#endif
	    
	    // eliminate item from base; 
	    deleteSubgradFromBase(baseitem,tol);

	  } // end if(dualVariables_[baseitem] < tol)
	} // end loop over baseitem 
	  
	if(!deleted){ // nothing deleted and not optimal
#if( DEBUG_TT )	  
	  std::cout << "Returning because nothing deleted and not optimal" << std::endl;
#endif
	  QPStatus_ = -2; // loop
	  return iter;
	}
      } // end inner loop
      
      Real newobjval; // new objective value
      Lin = 0.0;
      for (unsigned i=0;i<currSize_;i++){
	Lin += this->alpha(base[i])*this->getDualVariables(base[i]);
      }
      
      if (rho_ == -INF){
	Quad = -Lin/t;
	newobjval = - Quad/2;
      }
      else{
	Quad = rho_ - Lin/t;
	newobjval = (rho_ + Lin/t)/2;
      }

#if( DEBUG_TT )	  
      std::cout << "New Obj value = " << newobjval << std::endl;
#endif

#if( TABOO_LIST )
      // -- test for strict decrease -- // 
      // if item didn't provide decrease, move it to taboo list ...
      if( ( entering_ < maxind ) && ( objval < INF ) ){
	if( newobjval >= objval - std::max( tol*std::abs(objval), 1.e-16 ) ){
#if( DEBUG_TT )
	  std::cout << "Blocking " << entering_ << " because it didn't provide decrease" << std::endl;
#endif
	  taboo_.push_back(entering_);
	} 
      }
#endif

      objval = newobjval;

      // if sufficient decrease obtained
      if ( objval + std::max( tol*std::abs(objval), 1.e-16 ) <= minobjval ){
	taboo_.clear(); // reset taboo list
#if( DEBUG_TT )
	std::cout << "Taboo list has been reset." << std::endl;
#endif
	minobjval = objval;
      }

      //---------------------- OPTIMALITY TEST -------------------------//
      if ( (rho_ >= -INF) && (objval <= -INF) ) // if current x (dualVariables_) is feasible
	break;
	
      entering_ = maxind;
      Real minro = - std::max( tol*currSize_*std::abs(objval),1.e-16 );
#if ( ! FIRST_VIOLATED )
      Real diff = -INF, olddiff = -INF;
#endif

      for (unsigned bundleitem=0;bundleitem<this->size();bundleitem++){ // loop over items in bundle
      //for (int bundleitem=size_-1;bundleitem>-1;bundleitem--){ // loop over items in bundle (in reverse order)
	if ( std::find(taboo_.begin(),taboo_.end(),bundleitem) != taboo_.end() ){
#if( DEBUG_TT )	  
	  std::cout << "Item " << bundleitem << " is blocked." << std::endl;
#endif	  
	  continue; // if item is taboo, move on
	}

	if ( std::find(base.begin(),base.end(),bundleitem) == base.end() ){
	  // base does not contain bundleitem
#if( DEBUG_TT )	  
	  std::cout << "Base does not contain index " << bundleitem << std::endl;
#endif
	  Real ro = 0.0;
	  for (unsigned j=0;j<currSize_;j++){
	    ro += this->getDualVariables(base[j]) * GiTGj(bundleitem,base[j]);
	  }
	  ro += this->alpha(bundleitem)/t;

#if( DEBUG_TT )	  
	  std::cout << "RO = " << ro << std::endl;
#endif	  

	  if (rho_ >= -INF){
	    ro = ro - rho_; // note: rho = -v 
	  }
	  else{
	    ro = -INF;
	    minobjval = INF;
	    objval = INF;
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

#if( DEBUG_TT )	  
	    std::cout << "entering_ = " << entering_ << std::endl;
#endif

      //----------------- INSERTING ITEM ------------------------// 
      if (entering_ < maxind){ // dual constraint is violated
#if( DEBUG_TT )	  
	std::cout << "Inserting " << entering_ << std::endl;
#endif
	optimal_ = false;
	this->setDualVariables(entering_,0.0);
	base.push_back(entering_);
#if( DEBUG_TT )	  
	std::cout << "Current base = " << std::endl;
	for (unsigned k=0;k<base.size();k++){
	  std::cout << base[k] << std::endl;
	}
	std::cout << "dependent_ = " << dependent_ << std::endl;
#endif
	
	// construct new row of L_B
	unsigned zsize = currSize_ - dependent_; // zsize is the size of L_Bprime (current one)
	lh.size(zsize); // initialize to zeros;
	lhTz1 = 0.0;
	lhTz2 = 0.0;
	for (unsigned i=0;i<zsize;i++){
	  lh[i] = GiTGj(entering_,base[i]);
	}
	Teuchos::SerialDenseMatrix<int,Real> LBprime( Teuchos::Copy,L,zsize,zsize);      
	solveSystem(zsize,'N',LBprime,lh); // lh = (L_B^{-1})*(G_B^T*g_h)
	for (unsigned i=0;i<zsize;i++){
	  lhTz1 += lh[i]*z1[i];
	  lhTz2 += lh[i]*z2[i];
	}

	Real nrm = lh.dot(lh);
	Real delta = GiTGj(entering_,entering_) - nrm; // delta squared
#if( DEBUG_TT )	  
	std::cout << "GiTGj = " << GiTGj(entering_,entering_) << std::endl;
	std::cout << "lh_dot_lh = " << nrm << std::endl;
	std::cout << "delta = " << delta << std::endl;
#endif	
	
#if( ! EXACT )
	if(dependent_)
	  ljNorm = nrm; // adding second dependent
	else
	  lhNorm = nrm; // adding first dependent
#endif

	currSize_++; // update base size
	
	L.reshape(currSize_,currSize_);
	zsize = currSize_ - dependent_; // zsize is the size of L_Bprime (new one)
	for (unsigned i=0;i<zsize-1;i++){
	  L(currSize_-1,i) = lh[i];
	}

	Real deltaeps = tol*kappa*std::max(1.0,lh.dot(lh)); 
#if( DEBUG_TT )	  
	std::cout << "kappa = " << kappa << std::endl;
	std::cout << "deltaeps = " << deltaeps << std::endl;
#endif
	if ( delta > deltaeps ){ // new row is independent
	  // add subgradient to the base
	  unsigned ind = currSize_-1;
	  delta = std::sqrt(delta);
	  addSubgradToBase(ind,delta);
	}
	else if(delta < -deltaeps){
#if( DEBUG_TT )	  
	  std::cout << "NEGATIVE delta!" << std::endl;
#endif
	  dependent_++;
	  QPStatus_ = 0; // negative delta
	  return iter;
	}
	else{ // delta zero
	  dependent_++;
	}

#if( DEBUG_TT )	  
	std::cout << "Current L = " << std::endl;
	std::cout << L << std::endl;
	std::cout << "Current z1 = " << std::endl;
	std::cout << z1 << std::endl;
	std::cout << "Current z2 = " << std::endl;
	std::cout << z2 << std::endl;	  
#endif
      } // end if(entering_ < maxind)
      
      else{ // no dual constraint violated
	if( objval - std::max( tol*std::abs( objval ), 1.e-16 ) > minobjval ){ // check if f is as good as minf
#if( DEBUG_TT )	  
	  std::cout << "Returning because no dual constraint violated and f cannot reach min value " << minobjval << std::endl;
#endif
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

public:
  
  unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    unsigned iter = 0;
    if (this->size() == 1) {
      iter = solveDual_dim1(t,maxit,tol);
    }
    else if (this->size() == 2) {
      iter = solveDual_dim2(t,maxit,tol);
    }
    else {
      Real mytol = tol;
      unsigned outermaxit = 20;
      bool increase = false, decrease = false;
      iter = 0;
      for ( unsigned it=0; it < outermaxit; it++ ){
	iter += solveDual_TT(t,maxit,mytol);
	if ( QPStatus_ == 1 )
	  break;
	else if ( QPStatus_ == -2  || QPStatus_ == -3 ){
	  mytol /= 10.0;
	  decrease = true;
	}
	else {
	  mytol *= 10.0;
	  increase = true;
	}
	if ( (mytol > 1.e-4) || (mytol < 1.e-16) ){
	  break;
	}
	if ( increase && decrease ){
	  break;
	}
      }// end outer loop
#if ( DEBUG_TT )
      std::cout << "SolveDual returned after " << iter << " iterations with status " << QPStatus_ << " and solution" << std::endl;
      std::vector<Real> sol;
      for(unsigned i=0;i<this->size();i++){
	sol.push_back(this->getDualVariables(i));
      	std::cout << "x[" << i << "] = " << sol[i] << std::endl;
      }
      std::cout << std::endl;
      std::vector<Real> g(this->size(),0.0);
      Real val = evaluateObjective(g,sol,t);
      std::cout << "and objective value = " << val << std::endl; 
      std::cout << "Checking constraints" << std::endl;
      bool success = checkPrimary(g,sol,t);
      std::cout << success << std::endl;
#endif
    }
    return iter;
  }

  bool checkPrimary(std::vector<Real> &g, const std::vector<Real> &x, const Real t) const {
    bool success = true;
    this->gx_->zero(); this->eG_->zero();
    for (unsigned i = 0; i < this->size(); i++) {
      // Compute Gx using Kahan's compensated sum
      this->tG_->set(*this->gx_);
      this->yG_->set(*this->eG_); this->yG_->axpy(x[i],this->subgradient(i));
      this->gx_->set(*this->tG_); this->gx_->plus(*this->yG_);
      this->eG_->set(*this->tG_); this->eG_->axpy(-1.0,*this->gx_); this->eG_->plus(*this->yG_);
    }
    Real Hx = 0.0, v = 0.0, err = 0.0, tmp = 0.0, y = 0.0;
    for (unsigned i = 0; i < this->size(); i++) {
      // Compute < g_i, Gx > = - < g_i, d >
      Hx   = this->gx_->dot(this->subgradient(i));
      // Add to the objective function value using Kahan's compensated sum
      tmp  = v;
      y    = x[i]*(t*Hx + this->alpha(i)) + err;
      v    = tmp + y;
      err  = (tmp - v) + y;
      // Add gradient component
      g[i] = - t*Hx - this->alpha(i);
    }
    v *= -1.0;
    Real myeps = 1.e-8; 
    for (unsigned i = 0; i < this->size(); i++) {
      if ( g[i] > v + myeps ){
  	std::cout << "Constraint " << i << " is violated!: g[" << i << "] = " << g[i] << ", v = " << v  << std::endl;
  	success = false;
      }
      else if ( g[i] < v - myeps ){
  	std::cout << "Constraint " << i << " is inactive" << std::endl;
      }
      else{
  	std::cout << "Constraint " << i << " is active" << std::endl;
      }
    }
    return success;
  }

}; // class Bundle_TT

} // namespace ROL

#endif

