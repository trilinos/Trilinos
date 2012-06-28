/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#ifndef MULTIPRECCG_HPP_
#define MULTIPRECCG_HPP_

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_MatrixIO.hpp>

#include <iostream>
#include <functional>

#ifdef HAVE_TPETRA_QD
# include <qd/qd_real.h>
#endif

#define XFORM1     TPETRA_UNARY_TRANSFORM
#define XFORM2     TPETRA_BINARY_TRANSFORM
#define XFORM3     TPETRA_TERTIARY_TRANSFORM

#define REDUCE1(in, gexp) \
  Tpetra::RTI::reduce( *in,                                                           \
    Tpetra::RTI::reductionGlob<Tpetra::RTI::ZeroOp<decltype((in)->meanValue())>>(     \
                               [=]( decltype((in)->meanValue()) in ){ return gexp; }, \
                               std::plus<decltype((in)->meanValue())>() ) )

#define REDUCE2(in1, in2, gexp) \
  Tpetra::RTI::reduce( *in1, *in2,                                                      \
    Tpetra::RTI::reductionGlob<Tpetra::RTI::ZeroOp<decltype((in1)->meanValue())>>(      \
                               [=]( decltype((in1)->meanValue()) in1,                   \
                                    decltype((in2)->meanValue()) in2 ) { return gexp; },\
                               std::plus<decltype((in1)->meanValue())>() ) )

#define XFORM2RED(out,in,texp,gexp) \
  Tpetra::RTI::binary_pre_transform_reduce( *out, *in,                                  \
    Tpetra::RTI::reductionGlob<Tpetra::RTI::ZeroOp<decltype((out)->meanValue())>>(      \
                               [=]( decltype((out)->meanValue()) out,                   \
                                    decltype((in)->meanValue()) in )  { return texp; }, \
                               [=]( decltype((out)->meanValue()) out,                   \
                                    decltype((in)->meanValue()) in ) { return gexp; },  \
                               std::plus<decltype((out)->meanValue())>() ) )

#define XFORM3RED(out,in2,in3, texp, gexp)                                          \
  Tpetra::RTI::tertiary_pre_transform_reduce( *out, *in2, *in3,                     \
    Tpetra::RTI::reductionGlob<Tpetra::RTI::ZeroOp<decltype((out)->meanValue())>>(  \
                               [=]( decltype((out)->meanValue()) out,               \
                                    decltype((in2)->meanValue()) in2,               \
                                    decltype((in3)->meanValue()) in3 )              \
                                    { return texp; },                               \
                               [=]( decltype((out)->meanValue()) out,               \
                                    decltype((in2)->meanValue()) in2,               \
                                    decltype((in3)->meanValue()) in3 )              \
                                    { return gexp; },                               \
                               std::plus<decltype((out)->meanValue())>() ) )

namespace TpetraExamples {

  using Teuchos::RCP;
  using Tpetra::createVector;
  using Teuchos::ParameterList;

  namespace IRTRdetails {

    struct trivial_fpu_fix {
      void fix() {}
      void unfix() {}
    };
#ifdef HAVE_TPETRA_QD
    struct nontrivial_fpu_fix {
      unsigned int old_cw;
      void fix()   {fpu_fix_start(&old_cw);}
      void unfix() {fpu_fix_end(&old_cw);}
    };
#endif
    // implementations
    template <class T> struct fpu_fix : trivial_fpu_fix {};
#ifdef HAVE_TPETRA_QD
    template <> struct fpu_fix<qd_real> : nontrivial_fpu_fix {};
    template <> struct fpu_fix<dd_real> : nontrivial_fpu_fix {};
#endif

  }

  /***************************************************
  *   Implicit Riemannian Trust-Region Eigensolver
  *   Baker, Absil, Gallivan 2008
  ***************************************************/

  template <class S, class LO, class GO, class Node>
  S diagPrecond(const RCP<const Tpetra::Vector<S,LO,GO,Node> > &D,
                const RCP<const Tpetra::Vector<S,LO,GO,Node> > &x,
                const RCP<const Tpetra::Vector<S,LO,GO,Node> > &iDx,
                const S &xiDx,
                const RCP<const Tpetra::Vector<S,LO,GO,Node> > &r,
                const RCP<      Tpetra::Vector<S,LO,GO,Node> > &z)
  {
    /*
       Solve 
          P D P z = r
       for solution
          r = P_{u,x} inv(D) r
            = (I - u inv(x'*u) x') inv(D) r
            = z - u inv(x'*u) x'*z
       for u = inv(D) x = iDx
           z = inv(D) r
       
       for our purposes,
    1)        z = r
    2) fuse   z  \= D        
        and   xz  = x'*z     
    3) fuse   z = z - iDx ( xt / xiDx )
        and   return z'*r
    */
    XFORM2(z,r,  r);
    const S alpha = XFORM3RED(z,D,x,  z/D, x*z )  /  xiDx;
    return XFORM3RED(z,iDx,r,  z - iDx*alpha, z*r );
  }

  template <class S, class LO, class GO, class Node>
  struct tCG_workspace {
    RCP<Tpetra::Vector<S,LO,GO,Node>> r, z, iDx, d, Hd;
  };

  /***************************************************
  *   Trust-region Model Subproblem CG Solver
  *   Given x, compute eta that minimizes the
  *   quadratic Taylor expansion of f around x, 
  *   subject to  rho_x(eta) >= rho_prime, 
  *   or, equiv., eta'*eta <= ONE/rho_prime - ONE
  ***************************************************/
  //! \brief Model trust-region subproblem CG solver
  template <class S, class LO, class GO, class Node>
  void truncateCG(const RCP<Teuchos::FancyOStream> &out, int verbose, 
                  const RCP<const Tpetra::Operator<S,LO,GO,Node> > &A,
                  const RCP<const Tpetra::Vector<S,LO,GO,Node> > &diag,
                  const RCP<const Tpetra::Vector<S,LO,GO,Node> > &x,
                  const RCP<const Tpetra::Vector<S,LO,GO,Node> > &Ax,
                  const RCP<const Tpetra::Vector<S,LO,GO,Node> > &grad,
                  const S &lx, 
                  const RCP<Tpetra::Vector<S,LO,GO,Node> > &eta,
                  const tCG_workspace<S,LO,GO,Node> &workspace,
                  double rho_prime, double kappa, double theta, int maxItersCG) 
  {
    typedef Teuchos::ScalarTraits<S> ST;
    Teuchos::OSTab tab(out,2);
    const S ZERO = ST::zero();
    const S  ONE = ST::one();
    const S   D2 = ONE/rho_prime - ONE;
    decltype(eta) r = workspace.r,
                  z = workspace.z,
                iDx = workspace.iDx,
                  d = workspace.d,
                 Hd = workspace.Hd;
    //////////////////////////////////////////////
    // init precond. compute:
    //     iDx = inv(D)*x 
    //    xiDx = x'*iDx
    const S xiDx = XFORM3RED(iDx,diag,x,  x/diag,  x*iDx );
    //////////////////////////////////////////////
    // init solver
    XFORM1( eta,      ZERO );                     // eta = 0
    S r_r = XFORM2RED( r,grad,   grad,  r*r );    // r = grad
    S e_e = ZERO;
    S z_r = diagPrecond<S,LO,GO,Node>(diag,x,iDx,xiDx,r,z);
    XFORM2( d,z,  -z );                           // d = -z
    const S rnrm0 = ST::squareroot(r_r);
    const S kconv = rnrm0*kappa;
    const S tconv = ST::pow(rnrm0, (theta+1));
    const S convtol = TEUCHOS_MAX(TEUCHOS_MIN(kconv,tconv), lx*ST::eps());
    std::string convstr;
    if (kconv < tconv) convstr = "non-local convergence";
    else               convstr = "local convergence";
    //////////////////////////////////////////////
    // iterate
    if (verbose > 2) *out << "starting tCG with |res| =  " << rnrm0 << std::endl;
    if (verbose > 2) *out << "targeting tol = " << convtol << std::endl; // DEBUG
    int k = 0;
    std::string reason;
    while (true)
    {
      // apply Hessian Hd = A*d - lx*d
      A->apply(*d,*Hd);
      const S s1 = lx;
      decltype(Hd) v1 = Hd;
      decltype(x)  v2 = x;
      decltype(d)  v3 = d;
      // something is wrong with this form: it crashes GCC 4.5.1 and GCC 4.6.1; i don't know whether it's me or them
      //const S xHd = XFORM3RED(Hd,x,d,  Hd -  lx*d,    x*Hd );
      const S xHd = XFORM3RED(v1,v2,v3,  v1 -  s1*v3,    v2*v1 );
      const S dHd = XFORM3RED(Hd,d,x,  Hd - xHd*x,    d*Hd );
      if (verbose > 9) *out << "curvature: " << dHd << std::endl; // DEBUG
      // compute step size
      const S alpha = z_r / dHd;
      if (verbose > 9) *out << "alpha: " << alpha << std::endl;
      // anticipate e'*e for this step size
      const S e_d = REDUCE2(eta,d,   eta*d),
              d_d = REDUCE1(d,       d*d);
      const S new_e_e = e_e + 2.0*alpha*e_d + alpha*alpha*d_d;

      // check for truncation: negative curvature or implicit trust-region boundary
      if (dHd <= ZERO || new_e_e > D2) 
      {
        // find truncated step to boundary
        if (verbose > 9) *out << "e_d: " << e_d << std::endl; // DEBUG
        if (verbose > 9) *out << "d_d: " << d_d << std::endl; // DEBUG
        if (verbose > 9) *out << "e_e: " << e_e << std::endl; // DEBUG
        const S alpha_prime = (-e_d + ST::squareroot(e_d*e_d + d_d*(D2-e_e))) / d_d;
        if (verbose > 9) *out << "alpha_prime: " << alpha_prime << std::endl;
        XFORM2(eta,d,  eta + d*alpha_prime );
        if (dHd <= ZERO) reason="negative curvature";
        else              reason="exceeded trust-region";
        //
        break;
      }

      // compute new optimal point
      XFORM2(eta,d,  eta + d*alpha );
      e_e = new_e_e;    

      // update the residual
      const S x_r = XFORM3RED(r,Hd,x,    r + Hd*alpha,  x*r );
      // re-tangentialize r
      r_r = XFORM2RED(r,x,  r - x*x_r,  r*r );
      // compute norm
      S rnrm = ST::squareroot(r_r);

      // check stopping criteria
      if (rnrm <= convtol) {reason = convstr; break;}

      // save the old z'*r
      const S zold_rold = z_r;
      // z = diagPrec(r)
      z_r = diagPrecond<S,LO,GO,Node>(diag,x,iDx,xiDx,r,z);
      // compute new search direction
      const S beta = z_r / zold_rold;
      XFORM2(d,z,   -z + d*beta );

      k += 1;
      if (k >= maxItersCG) {reason="maximum number of inner iterations";break;}
    }
    if (verbose > 1) {
      *out << "leaving tCG after " << k << " iterations: " << reason << std::endl;
    }
  }

  //! \brief Templated IRTR solver
  template <class S, class LO, class GO, class Node>      
  S IRTR_singleprec(const RCP<Teuchos::FancyOStream> &out, ParameterList &params,
                    const RCP<const Tpetra::CrsMatrix<S,LO,GO,Node> > &A,
                    const RCP<Tpetra::Vector<S,LO,GO,Node> > &x)
  {
    IRTRdetails::fpu_fix<S> ff; ff.fix();
    Teuchos::OSTab tab(out,2);
    std::string methName("IRTR_singleprec<"+Teuchos::TypeNameTraits<S>::name()+">");
    typedef Tpetra::Vector<S,LO,GO,Node>   Vec;
    typedef Tpetra::Operator<S,LO,GO,Node>  Op;
    typedef Teuchos::ScalarTraits<S>        ST;
    // get objects from level database
    const int    maxIters   = params.get<int   >("maxIters",  100);
    const int    maxItersCG = params.get<int   >("maxItersCG",A->getDomainMap()->getGlobalNumElements());
    const int    verbose    = params.get<int   >("verbose",0);
    const double tolerance  = params.get<double>("tolerance");
    const double kappa      = params.get<double>("kappa",0.1);
    const double theta      = params.get<double>("theta",2.0);
    const double rho_prime  = params.get<double>("rho_prime",.1);
    if (verbose > 10) *out << params << std::endl;
    if (verbose > 10) *out << "maxIters: " << maxIters << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(rho_prime >= 1.0, std::runtime_error, "IRTR: rho_prime must be less than 1");
    RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer(methName);
    // x, ax, gradfx, eta
    auto map = A->getDomainMap();
    auto Ax = createVector<S>(map),
          r = createVector<S>(map),
        eta = createVector<S>(map),
       diag = createVector<S>(map);
    tCG_workspace<S,LO,GO,Node> workspace;
    workspace.r   = createVector<S>(map);
    workspace.z   = createVector<S>(map);
    workspace.iDx = createVector<S>(map);
    workspace.d   = createVector<S>(map);
    workspace.Hd  = createVector<S>(map);
    A->getLocalDiagCopy(*diag);
    timer->start(); 
    //////////////////////////////
    // init the solver
    S x_x = REDUCE1(x,         x*x );                           // x'*x
    S xnrm  = ST::squareroot(x_x);
    XFORM1( x,  x/xnrm);                                        // normalize x
    A->apply(*x,*Ax);                                           // Ax = A*x
    S lambda = REDUCE2(x,Ax,  x*Ax );                           // x'*A*x
    S rr = XFORM3RED(r,x,Ax,  Ax - x*lambda,  r*r );            // r = (I - x*x')Ax = Ax - x*lambda and compute its norm
    S rnrm = ST::squareroot(rr);
    if (verbose) {
      *out << "Beginning " << methName << " with λ = " << lambda << std::endl;
    }
    //////////////////////////////
    // iterate
    int k=0;
    while (k < maxIters && rnrm > tolerance) {
      //////////////////////////////
      // solve the trust-region subproblem for update eta
      truncateCG<S,LO,GO,Node>(out,verbose,A,diag,
                               x,Ax,r,lambda,
                               eta, workspace,
                               rho_prime,kappa,theta,maxItersCG);
      //////////////////////////////
      // update iterate, compute new eigenvalue
      x_x = XFORM2RED(x,eta,   x + eta,  x*x );                 // x = x+eta and x'*x
      xnrm = ST::squareroot(x_x);
      if (verbose > 2) *out << "|x+eta| = " << xnrm << std::endl;
      XFORM1(x,  x/xnrm);                                       // normalize x
      A->apply(*x,*Ax);                                         // Ax = A*x
      lambda = REDUCE2(x,Ax,    x*Ax );                         // lambda = x'*A*x
      rr = XFORM3RED(r,x,Ax,  Ax - x*lambda, r*r );             // r = A*x - x*lambda and r'*r
      rnrm = ST::squareroot(rr) / lambda;
      k += 1;
      if (verbose) {
        *out << "After iteration " << k << ", λ = " << lambda << ", |Ax-xλ|/|λ| = " << rnrm << std::endl;
      }
    }
    timer->stop(); 
    ff.unfix();
    if (verbose) {
      *out << "Leaving " << methName << " after " << k << " iterations." << std::endl;
    }
    return lambda;
  }

} // end of namespace TpetraExamples

#endif // MULTIPRECCG_HPP_
