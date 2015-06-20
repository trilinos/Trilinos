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

#ifndef IRTR_DRIVER_HPP
#define IRTR_DRIVER_HPP

#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_Version.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include "IRTR_solvers.hpp"


#define XFORM_REDUCE TPETRA_BINARY_PRETRANSFORM_REDUCE

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


template <class S, class Sinner = S>
class IRTR_Driver {
  public:
  // input
  Teuchos::RCP<Teuchos::FancyOStream>     out;
  Teuchos::RCP<Teuchos::ParameterList> params;
  std::string                      matrixFile;
  // output
  bool                             testPassed;

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type Node;

  void
  run (const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  {
    using std::plus;
    using std::endl;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Tpetra::RTI::ZeroOp;

    // Static types
    typedef Tpetra::Map<>            Map;
    typedef Tpetra::CrsMatrix<S>     CrsMatrix;
    typedef Tpetra::Vector<S>        Vector;
    typedef Teuchos::ScalarTraits<S> ST;

    IRTRdetails::fpu_fix<S> ff; ff.fix();

    *out << "Running test with Node==" << Teuchos::TypeNameTraits<Node>::name ()
         << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;

    // read the matrix
    RCP<CrsMatrix> A;
    RCP<const Map> rowMap = null;
    RCP<ParameterList> fillParams = parameterList();
    // must preserve the local graph in order to do convert() calls later
    if (Teuchos::TypeTraits::is_same<S,Sinner>::value) {
      fillParams->set("Preserve Local Graph",false);
    }
    else {
      fillParams->set("Preserve Local Graph",true);
    }
    RCP<Node> node = rcp (new Node ());
    Tpetra::Utils::readHBMatrix (matrixFile, comm, node, A, rowMap, fillParams);
    rowMap = A->getRowMap ();

    testPassed = true;

    // compute an inital vector
    auto x = Tpetra::createVector<S> (rowMap);
    x->randomize ();

    // call the solve
    S lambda = TpetraExamples::IRTR<Sinner,S,LO,GO,Node> (out, *params, A, x);

    // check that residual is as requested
    {
      auto r = Tpetra::createVector<S> (rowMap);
      A->apply (*x, *r);
      // compute A*x - x*lambda, while simultaneously computing |A*x - x*lambda|
      const S r_r =
        XFORM_REDUCE(r, x,                          // fused:
                     r - x*lambda,                  //      : r = r - x*lambda = A*x - x*lambda
                     r*r, ZeroOp<S>, plus<S>() );   //      : sum r'*r
      const S rnrm = Teuchos::ScalarTraits<S>::squareroot (r_r);
      // check that residual is as requested
      *out << "|A*x - x*lambda|/|lambda|: " << rnrm / ST::magnitude (lambda) << endl;
      const double tolerance = params->get<double> ("tolerance");
      if (rnrm / lambda > tolerance) {
        testPassed = false;
      }
    }

    ff.unfix();

    // print timings
    Teuchos::TimeMonitor::summarize( *out );
  }
};

#endif // IRTR_DRIVER_HPP
