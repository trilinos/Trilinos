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

#ifndef ROL_OPTIMIZATIONPROBLEMREFACTOR_HPP
#define ROL_OPTIMIZATIONPROBLEMREFACTOR_HPP

#include "ROL_BoundConstraint_Partitioned.hpp"
#include "ROL_CompositeConstraint.hpp"
#include "ROL_SlacklessObjective.hpp"
#include "ROL_RandomVector.hpp"

namespace ROL {

/* Represents optimization problems in Type-EB form 
 */

template<class Real>
class OptimizationProblem {

  typedef Vector<Real>                       V;
  typedef BoundConstraint<Real>              BND;
  typedef CompositeConstraint<Real>          CCON;
  typedef EqualityConstraint<Real>           EQCON;
  typedef InequalityConstraint<Real>         INCON;
  typedef Objective<Real>                    OBJ;
  typedef PartitionedVector<Real>            PV;
  typedef SlacklessObjective<Real>           SLOBJ;

  typedef Elementwise::AbsoluteValue<Real>   ABS;
  typedef Elementwise::Fill<Real>            FILL;

  typedef typename PV::size_type  size_type;

private:

  Teuchos::RCP<OBJ>      ORIGINAL_obj_;
  Teuchos::RCP<V>        ORIGINAL_sol_;
  Teuchos::RCP<BND>      ORIGINAL_bnd_;
  Teuchos::RCP<EQCON>    ORIGINAL_econ_;
  Teuchos::RCP<V>        ORIGINAL_emul_;
  Teuchos::RCP<INCON>    ORIGINAL_icon_;
  Teuchos::RCP<V>        ORIGINAL_imul_;

  Teuchos::RCP<OBJ>      obj_;
  Teuchos::RCP<V>        sol_;
  Teuchos::RCP<BND>      bnd_;
  Teuchos::RCP<EQCON>    con_;
  Teuchos::RCP<V>        mul_;

  EProblem problemType_;

  bool isInitialized_;

protected:
  void initialize(const Teuchos::RCP<Objective<Real> >            &obj,
                  const Teuchos::RCP<Vector<Real> >               &x,
                  const Teuchos::RCP<BoundConstraint<Real> >      &bnd,
                  const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon,
                  const Teuchos::RCP<Vector<Real> >               &le,
                  const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                  const Teuchos::RCP<Vector<Real> >               &li ) {
    using Teuchos::RCP; using Teuchos::rcp;
    if (!isInitialized_) {
      // If we have an inequality constraint
      if( incon != Teuchos::null ) {

        Real tol = std::sqrt(ROL_EPSILON<Real>());      

        // Create slack variables s = |c_i(x)|
        RCP<V> s = li->dual().clone();
        incon->value(*s,*x,tol);
        s->applyUnary(ABS()); 

        sol_ = CreatePartitionedVector(x,s);

        RCP<BND> xbnd, sbnd; 

        RCP<V> sl = s->clone();
        RCP<V> su = s->clone();

        sl->applyUnary( FILL(0.0) );
        su->applyUnary( FILL(ROL_INF<Real>()) );

        sbnd = rcp( new BND(sl,su) );    
  
        // Create a do-nothing bound constraint for x if we don't have one
        if( bnd == Teuchos::null ) {
          xbnd = rcp( new BND(*x) );
        }
        else { // Otherwise use the given bound constraint on x
          xbnd = bnd;  
        }
       
        // Create a partitioned bound constraint on the optimization and slack variables 
        bnd_ = CreateBoundConstraint_Partitioned(xbnd,sbnd);

        // Create partitioned lagrange multiplier and composite constraint
        if( eqcon == Teuchos::null ) {
          mul_ = CreatePartitionedVector(li);
          con_ = rcp( new CCON(incon,x) );
        }
        else {
          mul_ = CreatePartitionedVector(li,le);
          con_ = rcp( new CCON(incon,eqcon,x) );
        }

        obj_ = rcp( new SLOBJ(obj) );
      }
      else {  // There is no inequality constraint
   
        obj_ = obj;
        sol_ = x;
        mul_ = le;
        bnd_ = bnd;
        con_ = eqcon;
      }

      if( con_ == Teuchos::null ) {    // Type-U or Type-B
        if( bnd_ == Teuchos::null || !bnd_->isActivated() ) {  // Type-U
          problemType_ = TYPE_U;        
        }
        else { // Type-B
          problemType_ = TYPE_B; 
        }
      }
      else { // Type-E or Type-EB
        if( bnd_ == Teuchos::null || !bnd_->isActivated() ) { // Type-E
          problemType_ = TYPE_E;     
        }
        else { // Type-EB
          problemType_ = TYPE_EB; 
        }
      }
      isInitialized_ = true;
    }
  }

public:
  virtual ~OptimizationProblem(void) {}

  // Complete option constructor [1]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon,
                       const Teuchos::RCP<Vector<Real> >               &le,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li )
    : ORIGINAL_obj_(obj), ORIGINAL_sol_(x), ORIGINAL_bnd_(bnd),
      ORIGINAL_econ_(eqcon), ORIGINAL_emul_(le),
      ORIGINAL_icon_(incon), ORIGINAL_imul_(li),
      isInitialized_(false) {
    initialize(obj,x,bnd,eqcon,le,incon,li);
  }

  // No inequality constructor [2]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,   
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon, 
                       const Teuchos::RCP<Vector<Real> >               &le ) :    
     OptimizationProblem( obj, x, bnd, eqcon, le, Teuchos::null, Teuchos::null ) { } 

  // No equality constructor [3]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li ): 
     OptimizationProblem( obj, x, bnd, Teuchos::null, Teuchos::null, incon, li ) { } 

  // No bound constuctor [4]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon, 
                       const Teuchos::RCP<Vector<Real> >               &le,   
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li ) :
     OptimizationProblem( obj, x, Teuchos::null, eqcon, le, incon, li ) {}

  // No inequality or equality [5]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >     &bnd ) :
     OptimizationProblem( obj, x, bnd, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null ) { } 
  // No inequality or bound [6]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon,
                       const Teuchos::RCP<Vector<Real> >               &le ) :
     OptimizationProblem( obj, x, Teuchos::null, eqcon, le, Teuchos::null, Teuchos::null ) { } 

  // No equality or bound [7]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li ) :
     OptimizationProblem( obj, x, Teuchos::null, Teuchos::null, Teuchos::null, incon, li ) { } 

  // Unconstrained problem [8]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x ) :
     OptimizationProblem( obj, x, Teuchos::null, Teuchos::null, Teuchos::null, 
                          Teuchos::null, Teuchos::null ) { } 

  /* Set methods */

  virtual void setObjective(const Teuchos::RCP<Objective<Real> > &obj) {
    isInitialized_ = false;
    ORIGINAL_obj_ = obj;
  }

  virtual void setSolutionVector(const Teuchos::RCP<Vector<Real> > &sol) {
    isInitialized_ = false;
    ORIGINAL_sol_ = sol;
  }

  virtual void setBoundConstraint(const Teuchos::RCP<BoundConstraint<Real> > &bnd) {
    isInitialized_ = false;
    ORIGINAL_bnd_ = bnd;
  }

  virtual void setEqualityConstraint(const Teuchos::RCP<EqualityConstraint<Real> > &con) {
    isInitialized_ = false;
    ORIGINAL_econ_ = con;
  }

  virtual void setEMultiplierVector(const Teuchos::RCP<Vector<Real> > &mul) {
    isInitialized_ = false;
    ORIGINAL_emul_ = mul;
  }

  virtual void setInequalityConstraint(const Teuchos::RCP<InequalityConstraint<Real> > &con) {
    isInitialized_ = false;
    ORIGINAL_icon_ = con;
  }

  virtual void setIMultiplierVector(const Teuchos::RCP<Vector<Real> > &mul) {
    isInitialized_ = false;
    ORIGINAL_imul_ = mul;
  }

  /* Get methods */

  Teuchos::RCP<Objective<Real> > getObjective(void) {
    initialize(ORIGINAL_obj_,ORIGINAL_sol_,ORIGINAL_bnd_,
               ORIGINAL_econ_,ORIGINAL_emul_,
               ORIGINAL_icon_,ORIGINAL_imul_);
    return obj_;
  }

  Teuchos::RCP<Vector<Real> > getSolutionVector(void) {
    initialize(ORIGINAL_obj_,ORIGINAL_sol_,ORIGINAL_bnd_,
               ORIGINAL_econ_,ORIGINAL_emul_,
               ORIGINAL_icon_,ORIGINAL_imul_);
    return sol_;
  }

  Teuchos::RCP<BoundConstraint<Real> > getBoundConstraint(void) {
    initialize(ORIGINAL_obj_,ORIGINAL_sol_,ORIGINAL_bnd_,
               ORIGINAL_econ_,ORIGINAL_emul_,
               ORIGINAL_icon_,ORIGINAL_imul_);
    return bnd_;
  }

  Teuchos::RCP<EqualityConstraint<Real> > getEqualityConstraint(void) {
    initialize(ORIGINAL_obj_,ORIGINAL_sol_,ORIGINAL_bnd_,
               ORIGINAL_econ_,ORIGINAL_emul_,
               ORIGINAL_icon_,ORIGINAL_imul_);
    return con_;
  }

  Teuchos::RCP<Vector<Real> > getMultiplierVector(void) {
    initialize(ORIGINAL_obj_,ORIGINAL_sol_,ORIGINAL_bnd_,
               ORIGINAL_econ_,ORIGINAL_emul_,
               ORIGINAL_icon_,ORIGINAL_imul_);
    return mul_;
  }

  EProblem getProblemType(void) {
    return problemType_;
  }

  // Check derivatives, and consistency 
  void check( std::ostream &outStream = std::cout, const int numSteps = ROL_NUM_CHECKDERIV_STEPS, const int order = 1 ) {

    Teuchos::RCP<V> x = sol_->clone();
    Teuchos::RCP<V> y = sol_->clone();
    Teuchos::RCP<V> u = sol_->clone();
    Teuchos::RCP<V> v = sol_->clone();

    RandomizeVector(*x);
    RandomizeVector(*y);
    RandomizeVector(*u);
    RandomizeVector(*v);

    outStream << "\nPerforming OptimizationProblem diagnostics.\n\n";

    outStream << "Checking vector operations in optimization vector space X." << std::endl;
    x->checkVector(*y,*u,true,outStream);
 
    outStream << "Checking objective function." << std::endl;

    obj_->checkGradient(*x,*v,true,outStream,numSteps,order);                     outStream << std::endl;
    obj_->checkHessVec(*x,*u,true,outStream,numSteps,order);                      outStream << std::endl;
    obj_->checkHessSym(*x,*u,*v,true,outStream);                                  outStream << std::endl;
    
    if(con_ != Teuchos::null) {
      Teuchos::RCP<V> c = mul_->dual().clone();
      Teuchos::RCP<V> l = mul_->clone();
      Teuchos::RCP<V> w = mul_->clone();    
      Teuchos::RCP<V> q = mul_->clone();    

      RandomizeVector(*c);
      RandomizeVector(*l);
      RandomizeVector(*w);
      RandomizeVector(*q);   

      outStream << "Checking vector operations in constraint multiplier space C*." << std::endl;
      l->checkVector(*q,*w,true,outStream);

      outStream << "Checking equality constraint." << std::endl;
      con_->checkApplyJacobian(*x,*v,*c,true,outStream,numSteps,order);             outStream << std::endl;
      con_->checkAdjointConsistencyJacobian(*l,*u,*x,true,outStream);               outStream << std::endl;
      con_->checkApplyAdjointHessian(*x,*l,*v,*u,true,outStream,numSteps,order);    outStream << std::endl;  
    }

  }
  
  

}; // class OptimizationProblem

}  // namespace ROL

#endif // ROL_OPTIMIZATIONPROBLEMREFACTOR_HPP
