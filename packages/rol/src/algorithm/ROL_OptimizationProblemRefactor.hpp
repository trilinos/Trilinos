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

namespace ROL {
namespace Refactor {

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

  Teuchos::RCP<OBJ>      obj_;
  Teuchos::RCP<V>        sol_;
  Teuchos::RCP<BND>      bnd_;
  Teuchos::RCP<EQCON>    con_;
  Teuchos::RCP<V>        mul_;

  Teuchos::RCP<Teuchos::ParameterList>  parlist_;

public:
  virtual ~OptimizationProblem(void) {}

  // Complete option constructor [1]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon,
                       const Teuchos::RCP<Vector<Real> >               &le,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) : 
      parlist_(parlist) { 
   
    using Teuchos::RCP; using Teuchos::rcp; 

    // If we have an inequality constraint
    if( incon != Teuchos::null ) {

      Real tol = 0;      

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
        con_ = rcp( new CCON(incon) );
      }
      else {
        mul_ = CreatePartitionedVector(li,le);
        con_ = rcp( new CCON(incon,eqcon) );
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
  }

  // No inequality constructor [2]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,   
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon, 
                       const Teuchos::RCP<Vector<Real> >               &le,    
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, bnd, eqcon, le, Teuchos::null, Teuchos::null, parlist ) { } 

  // No equality constructor [3]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, bnd, Teuchos::null, Teuchos::null, incon, li, parlist ) { } 

  // No bound constuctor [4]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon, 
                       const Teuchos::RCP<Vector<Real> >               &le,   
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li,   
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, Teuchos::null, eqcon, le, incon, li, parlist ) {}

  // No inequality or equality [5]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<BoundConstraint<Real> >      &bnd,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, bnd, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, parlist ) { } 
  // No inequality or bound [6]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<EqualityConstraint<Real> >   &eqcon,
                       const Teuchos::RCP<Vector<Real> >               &le,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, Teuchos::null, eqcon, le, Teuchos::null, Teuchos::null, parlist ) { } 

  // No equality or bound [7]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<Vector<Real> >               &li,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, Teuchos::null, Teuchos::null, Teuchos::null, incon, li, parlist ) { } 

  // Unconstrained problem [8]
  OptimizationProblem( const Teuchos::RCP<Objective<Real> >            &obj,
                       const Teuchos::RCP<Vector<Real> >               &x,
                       const Teuchos::RCP<Teuchos::ParameterList>      &parlist = Teuchos::null ) :
     OptimizationProblem( obj, x, Teuchos::null, Teuchos::null, Teuchos::null, 
                          Teuchos::null, Teuchos::null, parlist ) { } 

  /* Get/Set methods */

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

}; // class OptimizationProblem

}  // namespace Refactor

}  // namespace ROL

#endif // ROL_OPTIMIZATIONPROBLEMREFACTOR_HPP
