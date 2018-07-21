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
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>


#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>

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

  /// \brief Operator representing a diagonal matrix.
  /// \tparam S The "Scalar" type; the type of the entries of the
  ///   input and output vectors.
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
  void
  RTICG (const RCP<const Tpetra::Operator<S,LO,GO,Node> >& A,
         RCP<Tpetra::Vector<S,LO,GO,Node> > b,
         const RCP<Teuchos::FancyOStream>& out,
         ParameterList& plist)
  {
    using Tpetra::RTI::ZeroOp;
    using Teuchos::as;
    using std::endl;
    typedef Tpetra::Vector<S ,LO,GO,Node> Vector;
    typedef Tpetra::Operator<S,LO,GO,Node>    Op;
    typedef Teuchos::ScalarTraits<S>          ST;

    Teuchos::FancyOStream& os = *out;
    Teuchos::OSTab tab0 (os);

    // get objects from level database
    const int numIters = plist.get<int>("numIters");
    const S tolerance  = plist.get<double>("tolerance");
    const int verbose  = plist.get<int>("verbose");
    RCP<Time> timer = Teuchos::TimeMonitor::getNewTimer(
                        "CG<"+Teuchos::TypeNameTraits<S>::name()+">"
                      );
    auto r = Tpetra::createVector<S> (A->getRangeMap ());
    auto p = Tpetra::createVector<S> (A->getRangeMap ());
    auto Ap = Tpetra::createVector<S> (A->getRangeMap ());
    auto x = b;

    if (verbose) {
      os << "Beginning CG<" << Teuchos::TypeNameTraits<S>::name () << ">" << std::endl;
    }
    int k;
    { // begin timer scope
      Teuchos::TimeMonitor lcltimer (*timer);
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

      if (verbose) {
        os << "Finished initial setup; now iterating" << endl;
      }

      ///////////////////////////////////
      for (k = 0; k < numIters; ++k) {
        Teuchos::OSTab tab1 (os);
        os << "Iteration " << k << ":" << endl;
        Teuchos::OSTab tab2 (os);

        A->apply (*p, *Ap);                                         // Ap = A*p
        S pAp = TPETRA_REDUCE2( p, Ap,
                                p*Ap, ZeroOp<S>, plus<S>() );       // p'*Ap
        const S alpha = rr / pAp;
        TPETRA_BINARY_TRANSFORM( x,    p,  x + alpha*p  );          // x = x + alpha*p
        S rrold = rr;
        rr = TPETRA_BINARY_PRETRANSFORM_REDUCE(
                               r, Ap,                               // fused:
                               r - alpha*Ap,                        //      : r - alpha*Ap
                               r*r, ZeroOp<S>, plus<S>() );         //      : sum r'*r
        if (verbose) { // if (verbose > 1)
          os << "|res|/|res_0|: " << ST::squareroot (rr / r2) << endl;
        }
        if (rr/r2 < tolerance*tolerance) {
          if (verbose) {
            os << "Convergence detected! rr/r2 = " << rr/r2
               << " < tolerance^2 = " << tolerance*tolerance << endl;
          }
          break;
        }
        const S beta = rr / rrold;
        TPETRA_BINARY_TRANSFORM( p, r,   r + beta*p );               // p = z + beta*p
      }
    } // end timer scope
    if (verbose) {
      os << "Leaving recursiveFPCG<" << Teuchos::TypeNameTraits<S>::name()
         << "> after " << k << " iterations." << endl;
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

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Vector<>::node_type Node;

  void
  run (const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
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
    typedef Teuchos::ScalarTraits<S> STS;
    typedef typename STS::magnitudeType magnitude_type;

    // Static types
    typedef Tpetra::Operator<S>   Operator;
    typedef Tpetra::Vector<S>       Vector;
    typedef Tpetra::Map<> map_type;

    Teuchos::FancyOStream& os = *out;
    Teuchos::OSTab tab0 (os);

    os << "CG test:" << endl;
    Teuchos::OSTab tab1 (os);

    const int        N = params->get<int> ("size");
    const double kappa = params->get<double> ("kappa");
    const double tolerance = params->get<double> ("tolerance");
    {
      os << "Problem characteristics:" << endl;
      Teuchos::OSTab tab2 (os);
      os << "Node type: " << Node::name () << endl
         << "Number of processes: " << comm->getSize () << endl
         << "Problem size: " << N << endl
         << "Condition number: " << kappa << endl
         << "Required tolerance: " << tolerance << endl;
    }

    // create the operator
    RCP<const Operator> A;
    {
      typedef Tpetra::global_size_t GST;
      const GST numGlobal = static_cast<GST> (N);
      const GO indexBase = 0;
      RCP<const map_type> map = rcp (new map_type (numGlobal, indexBase, comm));
      A = Tpetra::RTI::kernelOp<S> (TpetraExamples::DiagKernel<S> (N, kappa), map);
    }
    os << "Created operator" << endl;

    testPassed = true;

    // Choose a random exact solution and compute a right-hand-side from it.
    auto x = Tpetra::createVector<S> (A->getDomainMap ());
    auto b = Tpetra::createVector<S> (A->getRangeMap ());
    auto xhat = Tpetra::createVector<S> (A->getRangeMap());
    x->randomize ();
    A->apply (*x, *b);
    TPETRA_BINARY_TRANSFORM(xhat, b, b); // xhat := b
    os << "Chose exact solution and computed right-hand side" << endl;
    {
      Teuchos::OSTab tab2 (os);
      os << "Solution characteristics:" << endl;
      Teuchos::OSTab tab3 (os);
      os << "Exact solution 2-norm: " << x->norm2 () << endl
         << "Right-hand side 2-norm: " << b->norm2 () << endl
         << "Right-hand side copy 2-norm (should be same as right-hand side 2-norm): "
         << xhat->norm2 () << endl;
    }

    // call the solve
    TpetraExamples::RTICG (A, xhat, out, *params);
    os << "Finished CG solve" << endl;
    {
      Teuchos::OSTab tab2 (os);
      os << "Check solution using Tpetra::Vector methods:" << endl;
      Teuchos::OSTab tab3 (os);
      os << "Computed solution 2-norm: " << xhat->norm2 () << endl;

      // FIXME (mfh 23 Feb 2014) I'm not sure I can trust the
      // operator's implementation of apply() to be correct for alpha
      // != 1 or beta != 0.  If it is, we could just use that and omit
      // the update() line below it.
      auto r = Tpetra::createVector<S> (b->getMap ());
      A->apply (*xhat, *r); // r := A * xhat
      r->update (STS::one (), *b, -STS::one ()); // r := 1*b - A*xhat
      const magnitude_type r_norm = r->norm2 ();
      const magnitude_type b_norm = b->norm2 ();

      os << "Absolute residual 2-norm: " << r_norm << endl
         << "Right-hand side 2-norm: " << b_norm << endl
         << "Relative residual 2-norm: " << (r_norm / b_norm) << endl;
    }

    // check that residual is as requested
    {
      auto bhat = Tpetra::createVector<S>(A->getDomainMap());
      A->apply (*xhat, *bhat);
      // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
      typedef ZeroOp<pair<S,S> > ZERO;
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
      const S enrm = STS::squareroot (nrms.first);
      const S bnrm = STS::squareroot (nrms.second);

      Teuchos::OSTab tab2 (os);
      os << "Check solution using RTI:" << endl;
      {
        Teuchos::OSTab tab2 (os);
        os << "Absolute residual 2-norm: " << enrm << endl
           << "Right-hand side 2-norm: " << bnrm << endl
           << "Relative residual 2-norm: " << enrm / bnrm << endl;
      }
      // give a little slack
      if (enrm / bnrm > 5*tolerance) {
        testPassed = false;
      }
    }

    // print timings
    Teuchos::TimeMonitor::summarize (*out);
  }
};
