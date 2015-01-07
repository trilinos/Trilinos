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

#ifndef MULTIPREC_DRIVER_HPP
#define MULTIPREC_DRIVER_HPP

#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_Version.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include "MultiPrecCG.hpp"

template <class MPStack>
class MultiPrecDriver {
  public:
  // input
  Teuchos::RCP<Teuchos::FancyOStream>     out; 
  Teuchos::RCP<Teuchos::ParameterList> params;
  std::string                      matrixFile;
  bool                             unfusedTest;
  // output
  bool                             testPassed;

  template <class Node> 
  void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) 
  {
    using std::pair;
    using std::make_pair;
    using std::plus;
    using std::endl;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using TpetraExamples::make_pair_op;
    using Tpetra::RTI::reductionGlob;
    using Tpetra::RTI::ZeroOp;
    using Tpetra::RTI::binary_pre_transform_reduce;
    using Tpetra::RTI::binary_transform;

    // Static types
    typedef typename MPStack::type   S;
    typedef int                     LO;
    typedef int                     GO;
    typedef Tpetra::Map<LO,GO,Node>               Map;
    typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;
    typedef Tpetra::Vector<S,LO,GO,Node>       Vector;

    *out << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;

    // read the matrix
    RCP<CrsMatrix> A;
    RCP<const Map> rowMap = null;
    RCP<ParameterList> fillParams = parameterList();
    fillParams->set("Preserve Local Graph",true);
    // must preserve the local graph in order to do convert() calls later
    Tpetra::Utils::readHBMatrix(matrixFile,comm,node,A,rowMap,fillParams);
    rowMap = A->getRowMap();

    // init the solver stack
    TpetraExamples::RFPCGInit<S,LO,GO,Node> init(A);
    RCP<ParameterList> db = Tpetra::Ext::initStackDB<MPStack>(*params,init);

    testPassed = true;

    // choose a solution, compute a right-hand-side
    auto x = Tpetra::createVector<S>(rowMap),
         b = Tpetra::createVector<S>(rowMap);
    x->randomize();
    A->apply(*x,*b);
    {
      // init the rhs
      auto bx = db->get<RCP<Vector>>("bx");
      binary_transform( *bx, *b, [](S, S bi) {return bi;}); // bx = b
    }

    // call the solve
    TpetraExamples::recursiveFPCG<MPStack,LO,GO,Node>(out,*db);

    // check that residual is as requested
    {
      auto xhat = db->get<RCP<Vector>>("bx"),
           bhat = Tpetra::createVector<S>(rowMap);
      A->apply(*xhat,*bhat);
      // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
      auto nrms = binary_pre_transform_reduce(*bhat, *b, 
                                              reductionGlob<ZeroOp<pair<S,S>>>( 
                                                [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
                                                [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
                                                make_pair_op<S,S>(plus<S>())) );
      const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
              bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
      // check that residual is as requested
      *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
      const double tolerance = db->get<double>("tolerance");
      if (MPStack::bottom) {
        // give a little slack
        if (enrm / bnrm > 5*tolerance) testPassed = false;
      }
      else {
        if (enrm / bnrm > tolerance) testPassed = false;
      }
    }

    // 
    // solve again, with the unfused version, just for timings purposes
    if (unfusedTest) 
    {
      // init the rhs
      auto bx = db->get<RCP<Vector>>("bx");
      binary_transform( *bx, *b, [](S, S bi) {return bi;}); // bx = b
      // call the solve
      TpetraExamples::recursiveFPCGUnfused<MPStack,LO,GO,Node>(out,*db);
      //
      // test the result
      auto xhat = db->get<RCP<Vector>>("bx"),
           bhat = Tpetra::createVector<S>(rowMap);
      A->apply(*xhat,*bhat);
      // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
      auto nrms = binary_pre_transform_reduce(*bhat, *b, 
                                              reductionGlob<ZeroOp<pair<S,S>>>( 
                                                [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
                                                [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
                                                make_pair_op<S,S>(plus<S>())) );
      const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
              bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
      // check that residual is as requested
      *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
      const double tolerance = db->get<double>("tolerance");
      if (MPStack::bottom) {
        // give a little slack
        if (enrm / bnrm > 5*tolerance) testPassed = false;
      }
      else {
        if (enrm / bnrm > tolerance) testPassed = false;
      }
    }    
         
         
    // print timings
    Teuchos::TimeMonitor::summarize( *out );
  }
};

#endif // MULTIPREC_DRIVER_HPP
