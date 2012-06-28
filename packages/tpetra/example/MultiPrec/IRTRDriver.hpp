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


template <class S>
class IRTR_Driver {
  public:
  // input
  Teuchos::RCP<Teuchos::FancyOStream>     out; 
  Teuchos::RCP<Teuchos::ParameterList> params;
  std::string                      matrixFile;
  // output
  bool                             testPassed;

  template <class Node> 
  void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) 
  {
    using std::plus;
    using std::endl;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Tpetra::RTI::ZeroOp;

    // Static types
    typedef int                     LO;
    typedef int                     GO;
    typedef Tpetra::Map<LO,GO,Node>               Map;
    typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;
    typedef Tpetra::Vector<S,LO,GO,Node>       Vector;

    *out << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;

    // read the matrix
    RCP<CrsMatrix> A;
    Tpetra::Utils::readHBMatrix(matrixFile,comm,node,A);
    RCP<const Map> rowMap = A->getRowMap();

    testPassed = true;

    // compute an inital vector
    auto x = Tpetra::createVector<S>(rowMap);
    x->randomize();

    // call the solve
    S lambda = TpetraExamples::IRTR_singleprec<S,LO,GO,Node>(out,*params,A,x);

    // check that residual is as requested
    {
      auto r = Tpetra::createVector<S>(rowMap);
      A->apply(*x,*r);
      // compute A*x - x*lambda, while simultaneously computing |A*x - x*lambda|
      const S r_r = XFORM_REDUCE(r, x,                          // fused: 
                                 r - x*lambda,                  //      : r = r - x*lambda = A*x - x*lambda
                                 r*r, ZeroOp<S>, plus<S>() );   //      : sum r'*r
      const S rnrm = Teuchos::ScalarTraits<S>::squareroot(r_r);
      // check that residual is as requested
      *out << "|A*x - x*lambda|/|lambda|: " << rnrm / lambda << endl;
      const double tolerance = params->get<double>("tolerance");
      if (rnrm / lambda > tolerance) testPassed = false;
    }

    // print timings
    Teuchos::TimeMonitor::summarize( *out );
  }
};

#endif // IRTR_DRIVER_HPP
