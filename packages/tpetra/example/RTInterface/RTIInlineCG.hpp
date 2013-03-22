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

#include <iostream>
#include <functional>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>


#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>

// no gpu support for RTI
#undef HAVE_KOKKOSCLASSIC_THRUST
#include <Tpetra_HybridPlatform.hpp>

#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_RTIOp.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

namespace Tpetra {
  namespace RTI {
    // specialization for pair
    template <class T1, class T2>
    class ZeroOp<std::pair<T1,T2>> {
      public:
      static inline std::pair<T1,T2> identity() {
        return std::make_pair( Teuchos::ScalarTraits<T1>::zero(), 
                               Teuchos::ScalarTraits<T2>::zero() );
      }
    };
  }
}

namespace TpetraExamples {

  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::Time;
  using Teuchos::null;
  using std::binary_function;
  using std::pair;
  using std::make_pair;
  using std::plus;
  using std::multiplies;

  template <class T1, class T2, class Op>
  class pair_op : public binary_function<pair<T1,T2>,pair<T1,T2>,pair<T1,T2>> {
  private:
    Op op_;
  public:
    pair_op(Op op) : op_(op) {}
    inline pair<T1,T2> operator()(const pair<T1,T2>& a, const pair<T1,T2>& b) const {
      return make_pair(op_(a.first,b.first),op_(a.second,b.second));
    }
  };

  template <class T1, class T2, class Op>
  pair_op<T1,T2,Op> make_pair_op(Op op) { return pair_op<T1,T2,Op>(op); }

  //! Class to demonstrate a simple diagonal matrix.
  template <class S>
  class DiagKernel : public Tpetra::RTI::detail::StdOpKernel<S>
  {
    protected:
      int N_;
      S kappa_;
    public:
      DiagKernel(int N, const S & kappa) : N_(N), kappa_(kappa) {}
      inline void execute(int i) { 
        this->_vec_inout[i] = this->_alpha * ((kappa_-1.0) * (S)(i)/(S)(N_-1) + 1.0) * this->_vec_in2[i] + this->_beta * this->_vec_inout[i];
      }
  };



  //! \brief Recursive, self-preconditioning flexible CG.
  template <class S, class LO, class GO, class Node>      
  void RTICG(const RCP<const Tpetra::Operator<S,LO,GO,Node>> &A, RCP<Tpetra::Vector<S,LO,GO,Node>> b, 
             const RCP<Teuchos::FancyOStream> &out, ParameterList &plist)
  {
    using Teuchos::as;
    using Tpetra::RTI::ZeroOp;
    typedef Tpetra::Vector<S ,LO,GO,Node> Vector;
    typedef Tpetra::Operator<S,LO,GO,Node>    Op;
    typedef Teuchos::ScalarTraits<S>          ST;
    // get objects from level database
    const int numIters = plist.get<int>("numIters");
    const S tolerance  = plist.get<double>("tolerance");
    const int verbose  = plist.get<int>("verbose");
    RCP<Time> timer = Teuchos::TimeMonitor::getNewTimer(
                        "CG<"+Teuchos::TypeNameTraits<S>::name()+">"
                      );
    auto r = Tpetra::createVector<S>(A->getRangeMap()),
         p = Tpetra::createVector<S>(A->getRangeMap()),
        Ap = Tpetra::createVector<S>(A->getRangeMap()); 
    auto x = b;

    Teuchos::OSTab tab(out);

    if (verbose) {
      *out << "Beginning CG<" << Teuchos::TypeNameTraits<S>::name() << ">" << std::endl;
    }
    int k;
    { // begin timer scop
      Teuchos::TimeMonitor lcltimer(*timer);
      // initial guess is x=0, so initial residual is r = b - A*x = b
      const S r2 = TPETRA_BINARY_PRETRANSFORM_REDUCE(
                      r, b,                                         // fused: 
                      b,                                            //      : r = x  
                      r*r, ZeroOp<S>, plus<S>() );                  //      : sum r'*r
      // forget b; refer to it as x now
      b = null;
      // b comes in, x goes out. now we're done with b, so zero the solution.
      TPETRA_UNARY_TRANSFORM( x,  ST::zero() );                     // x = 0
      S rr = r2;
      TPETRA_BINARY_TRANSFORM( p, r,    r );                        // p = r
      ///////////////////////////////////
      for (k=0; k<numIters; ++k) 
      {
        A->apply(*p,*Ap);                                           // Ap = A*p
        S pAp = TPETRA_REDUCE2( p, Ap,     
                                p*Ap, ZeroOp<S>, plus<S>() );       // p'*Ap
        const S alpha = rr / pAp;
        TPETRA_BINARY_TRANSFORM( x,    p,  x + alpha*p  );          // x = x + alpha*p
        S rrold = rr;
        rr = TPETRA_BINARY_PRETRANSFORM_REDUCE(
                               r, Ap,                               // fused:
                               r - alpha*Ap,                        //      : r - alpha*Ap
                               r*r, ZeroOp<S>, plus<S>() );         //      : sum r'*r
        if (verbose > 1) *out << "|res|/|res_0|: " << ST::squareroot(rr/r2) 
                              << std::endl;
        if (rr/r2 < tolerance*tolerance) {
          if (verbose) {
            *out << "Convergence detected!" << std::endl;
          }
          break;
        }
        const S beta = rr / rrold;
        TPETRA_BINARY_TRANSFORM( p, r,   r + beta*p );               // p = z + beta*p
      }
    } // end timer scop
    if (verbose) {
      *out << "Leaving recursiveFPCG<" << Teuchos::TypeNameTraits<S>::name() 
           << "> after " << k << " iterations." << std::endl;
    }
  }

} // end of namespace TpetraExamples

template <class S>
class CGDriver {
  public:
  // input
  Teuchos::RCP<Teuchos::FancyOStream>     out; 
  Teuchos::RCP<Teuchos::ParameterList> params;
  // output
  bool                             testPassed;

  template <class Node> 
  void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) 
  {
    using std::pair;
    using std::make_pair;
    using std::plus;
    using std::endl;
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::null;
    using TpetraExamples::make_pair_op;
    using Tpetra::RTI::reductionGlob;
    using Tpetra::RTI::ZeroOp;
    using Tpetra::RTI::binary_pre_transform_reduce;
    using Tpetra::RTI::binary_transform;

    // Static types
    typedef int                                    LO;
    typedef int                                    GO;
    typedef Tpetra::Operator<S,LO,GO,Node>   Operator;
    typedef Tpetra::Vector<S,LO,GO,Node>       Vector;
    const int        N = params->get<int>("size");
    const double kappa = params->get<double>("kappa");

    *out << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;

    // create the operator
    *out << "Building problem of size " << N << " with condition number " << kappa << std::endl;
    RCP<const Operator> A;
    {
      auto map = Tpetra::createUniformContigMapWithNode<LO,GO,Node>(N,comm,node);
      A = Tpetra::RTI::kernelOp<S>( TpetraExamples::DiagKernel<S>(N,kappa), map );
    }

    testPassed = true;

    // choose a solution, compute a right-hand-side
    auto x = Tpetra::createVector<S>(A->getDomainMap()),
         b = Tpetra::createVector<S>(A->getDomainMap()),
      xhat = Tpetra::createVector<S>(A->getDomainMap());
    x->randomize();
    A->apply(*x,*b);
    TPETRA_BINARY_TRANSFORM(xhat, b,    b); // xhat = b

    // call the solve
    TpetraExamples::RTICG(A,xhat,out,*params);

    // check that residual is as requested
    {
      auto bhat = Tpetra::createVector<S>(A->getDomainMap());
      A->apply(*xhat,*bhat);
      // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
      typedef ZeroOp<pair<S,S>> ZERO;
      const pair<S,S> nrms = TPETRA_BINARY_PRETRANSFORM_REDUCE(                                               
                               bhat, b,                                                                       
                               b - bhat,                                                                        
                               make_pair(bhat*bhat, b*b), ZERO, (make_pair_op<S,S>(plus<S>()))
                             );                                                                               
      //auto nrms = binary_pre_transform_reduce(*bhat, *b, 
      //                                        reductionGlob<ZeroOp<pair<S,S>>>( 
      //                                          [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
      //                                          [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
      //                                          make_pair_op<S,S>(plus<S>())) );
      const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
              bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
      // check that residual is as requested
      *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
      const double tolerance = params->get<double>("tolerance");
      // give a little slack
      if (enrm / bnrm > 5*tolerance) testPassed = false;
    }
         
    // print timings
    Teuchos::TimeMonitor::summarize( *out );
  }
};
