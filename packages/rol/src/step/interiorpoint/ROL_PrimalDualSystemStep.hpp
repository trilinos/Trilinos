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

#ifndef ROL_PRIMALDUALSYSTEMSTEP_H
#define ROL_PRIMALDUALSYSTEMSTEP_H

#include "ROL_NewtonKrylovStep.hpp"
#include "ROL_PrimalDualInteriorPointOperator.hpp"
#include "ROL_SchurComplememt.hpp"

/** @ingroup step_group
    \class ROL::PrimalDualSystemStep
    \brief Provides the interface to compute approximate
           solutions to 2x2 block systems arising from primal-dual 
           interior point methods

           Note that as we do not need an additional Lagrange multiplier
           for the primal dual system, the vector expected to be passed
           in its place is the primal-dual residual
 */

namespace ROL { 

template<class Real> 
class PrimalDualSystemStep : public Step<Real> {

  typedef Vector<Real>             V;
  typedef PartitionedVector<Real>  PV;
  typedef Objective<Real>          OBJ;
  typedef BoundConstraint<Real>    BND;
  typedef Constraint<Real> CON;
  typedef AlgorithmState<Real>     AS;
  typedef SchurComplement<Real>    SCHUR;

  typedef PrimalDualInteriorPointBlock11  OP11;
  typedef PrimalDualInteriorPointBlock12  OP12;
  typedef PrimalDualInteriorPointBlock21  OP21;
  typedef PrimalDualInteriorPointBlock22  OP22;


private:
  
  // Block indices
  static const size_type OPT   = 0;
  static const size_type EQUAL = 1;
  static const size_type LOWER = 2;
  static const size_type UPPER = 3;  

  // Super block indices
  static const size_type OPTMULT = 0;  // Optimization and equality multiplier components
  static const size_type BNDMULT = 1;  // Bound multiplier components

  ROL::Ptr<Secant<Real> > secant_;
  ROL::Ptr<Krylov<Real> > krylov_;
  ROL::Ptr<V> scratch1_;           // scratch vector 
  ROL::Ptr<V> scratch_; 

  ROL::Ptr<OP11> A_;
  ROL::Ptr<OP12> B_;
  ROL::Ptr<OP21> C_;
  ROL::Ptr<OP22> D_;

  ROL::Ptr<SCHUR> schur_; // Allows partial decoupling of (x,lambda) and (zl,zu)
  ROL::Ptr<OP>    op_;    // Solve fully coupled system

  int iterKrylov_; ///< Number of Krylov iterations (used for inexact Newton)
  int flagKrylov_; ///< Termination flag for Krylov method (used for inexact Newton)
  int verbosity_;  ///< Verbosity level

  bool useSecantPrecond_;
  bool useSchurComplement_;

  

  // Repartition (x,lambda,zl,zu) as (xlambda,z) = ((x,lambda),(zl,zu))
  ROL::Ptr<PV> repartition( V &x ) {
     
    PV &x_pv = dynamic_cast<PV&>(x);
    ROL::Ptr<V> xlambda = CreatePartitionedVector(x_pv.get(OPT),x_pv.get(EQUAL));  
    ROL::Ptr<V> z = CreatePartitionedVector(x_pv.get(LOWER),x_pv.get(UPPER));  
 
    ROL::Ptr<V> temp[] = {xlambda,z};

    return ROL::makePtr<PV( std::vector<ROL::Ptr<V> >>(temp,temp+2) );

  }

  // Repartition (x,lambda,zl,zu) as (xlambda,z) = ((x,lambda),(zl,zu))
  ROL::Ptr<const PV> repartition( const V &x ) {
    const PV &x_pv = dynamic_cast<const PV&>(x);
    ROL::Ptr<const V> xlambda = CreatePartitionedVector(x_pv.get(OPT),x_pv.get(EQUAL));  
    ROL::Ptr<const V> z = CreatePartitionedVector(x_pv.get(LOWER),x_pv.get(UPPER));  

    ROL::Ptr<const V> temp[] = {xlambda,z};

    return ROL::makePtr<PV( std::vector<ROL::Ptr<const V> >>(temp,temp+2) );
         
  }

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;


  PrimalDualSystemStep( Teuchos::ParameterList &parlist, 
                        const ROL::Ptr<Krylov<Real> > &krylov,
                        const ROL::Ptr<Secant<Real> > &secant,
                        ROL::Ptr<V> &scratch1 ) : Step<Real>(),
    krylov_(krylov), secant_(secant), scratch1_(scratch1), schur_(ROL::nullPtr),
    op_(ROL::nullPtr), useSchurComplement_(false) {

    PL &iplist = parlist.sublist("Step").sublist("Primal Dual Interior Point");
    PL &syslist = iplist.sublist("System Solver");

    useSchurComplement_ = syslist.get("Use Schur Complement",false);    
     
  }
 
  PrimalDualSystemStep( Teuchos::ParameterList &parlist,
                        ROL::Ptr<V> &scratch1_ ) : Step<Real>() {
    PrimalDualSystemStep(parlist,ROL::nullPtr,ROL::nullPtr,scratch1); 
  }

  void initialize( V &x, const V &g, V &res, const V &c,
                   OBJ &obj, CON &con, BND &bnd, AS &algo_state ) {

    Step<Real>::initialize(x,g,res,c,obj,con,bnd,algo_state);
 
     
     
    ;

    ROL::Ptr<OBJ> pObj = ROL::makePtrFromRef(obj);
    ROL::Ptr<CON> pCon = ROL::makePtrFromRef(con);
    ROL::Ptr<BND> pBnd = ROL::makePtrFromRef(bnd);
 
    ROL::Ptr<PV> x_pv = repartition(x);

    ROL::Ptr<V> xlambda = x_pv->get(OPTMULT);
    ROL::Ptr<V> z = x_pv->get(BNDMULT);
 
    A_ = ROL::makePtr<OP11>( pObj, pCon, *xlambda, scratch1_ );
    B_ = ROL::makePtr<OP12>( );
    C_ = ROL::makePtr<OP21>( *z );
    D_ = ROL::makePtr<OP22>( pBnd, *xlambda );

    if( useSchurComplement_ ) {
      schur_ = ROL::makePtr<SCHUR>(A_,B_,C_,D_,scratch1_);
    } 
    else {
      op_ = BlockOperator2<Real>(A_,B_,C_,D_);
    }
  }

  void compute( V &s, const V &x, const V &res, OBJ &obj, CON &con, 
                BND &bnd, AS &algo_state ) {

    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();


    if( useSchurComplement_ ) {
      
      ROL::Ptr<const PV> x_pv = repartition(x);
      ROL::Ptr<const PV> res_pv = repartition(res);
      ROL::Ptr<PV> s_pv = repartition(s);


      // Decouple (x,lambda) from (zl,zu) so that s <- L

      ROL::Ptr<V> sxl   = s_pv->get(OPTMULT);
      ROL::Ptr<V> sz    = s_pv->get(BNDMULT);
 
      

    }
    else {

    }

  }

  void update( V &x, V &res, const V &s, OBJ &obj, CON &con, 
               BND &bnd, AS &algo_state ) {

    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();

    
  }
  

};

} // namespace ROL

#endif  // ROL_PRIMALDUALSYSTEMSTEP_H
