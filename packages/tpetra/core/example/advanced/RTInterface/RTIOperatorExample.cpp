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

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Tpetra_Version.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_RTIOp.hpp"

#include <iostream>
#include <iomanip>
#include <functional>

/** \file RTIOperatorExample.cpp
    \brief An example of the Tpetra::Operator interface for Tpetra::RTI.
 */


namespace TpetraExamples {

  //! Class to demonstrate a simple, finite-different type tridiagonal stencil.
  template <class S>
  class FDStencil : public Tpetra::RTI::detail::StdOpKernel<S>
  {
    protected:
      int _myImageID, _numImages;
      int _n;
      S _sub, _diag, _super;
    public:
      FDStencil(int myImageID, int numImages, int n, const S &sub, const S &diag, const S &super) : _myImageID(myImageID), _numImages(numImages), _n(n), _sub(sub), _diag(diag), _super(super) {}
      inline void execute(int i) const
      {
        S res = _diag * this->_vec_in2[i];
        if (i >  0 || _myImageID != 0           ) res +=   _sub * this->_vec_in2[i-1];
        if (i < _n || _myImageID != _numImages-1) res += _super * this->_vec_in2[i+1];
        if (this->_alpha == Teuchos::ScalarTraits<S>::one() && this->_beta == Teuchos::ScalarTraits<S>::zero()) {
          this->_vec_inout[i] = res;
        }
        else {
          this->_vec_inout[i] = this->_alpha * res + this->_beta * this->_vec_inout[i];
        }
      }
  };

  //! Class to demonstrate a simple diagonal matrix.
  template <class S>
  class ScaleKernel : public Tpetra::RTI::detail::StdOpKernel<S>
  {
  protected:
    S _gamma;
    const S zero_; // Teuchos::ScalarTraits::{zero,one} are not device...
    const S one_;  // ...functions, so we construct zero and one on the host.

  public:
    ScaleKernel (const S & gamma) :
      _gamma (gamma),
      zero_ (Teuchos::ScalarTraits<S>::zero ()),
      one_ (Teuchos::ScalarTraits<S>::one ())
    {}

    inline void execute (const int i) const {
      if (this->_alpha == one_ && this->_beta == zero_) {
        this->_vec_inout[i] = _gamma * this->_vec_in2[i];
      }
      else {
        this->_vec_inout[i] = this->_alpha * _gamma * this->_vec_in2[i] +
          this->_beta * this->_vec_inout[i];
      }
    }
  };

}

int main(int argc, char *argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  using Teuchos::TimeMonitor;
  using TpetraExamples::FDStencil;
  using TpetraExamples::ScaleKernel;

  //
  // Get the default communicator and node
  //
  auto comm = Tpetra::getDefaultComm();
  const int myImageID = comm->getRank();
  const int numImages = comm->getSize();

  //
  // Get example parameters from command-line processor
  //
  bool verbose = (myImageID==0);
  int numGlobal_user = 100*comm->getSize();
  int numTimeTrials = 3;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("global-size",&numGlobal_user,"Global test size.");
  cmdp.setOption("num-time-trials",&numTimeTrials,"Number of trials in timing loops.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  //
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl;
    std::cout << "Comm info: " << *comm;
    std::cout << std::endl;
  }

  //
  // Create a simple map for domain and range
  //
  Tpetra::global_size_t numGlobalRows = numGlobal_user;
  auto map = Tpetra::createUniformContigMapWithNode<int,int>(numGlobalRows, comm);
  // const size_t numLocalRows = map->getNodeNumElements();
  auto x = Tpetra::createVector<double>(map),
       y = Tpetra::createVector<double>(map);

  // Create a simple diagonal operator using lambda function
  auto fTwoOp = Tpetra::RTI::binaryOp<double>( [](double /*y*/, double xx) { return 2.0 * xx; } , map );
  // y = 3*fTwoOp*x + 2*y = 3*2*1 + 2*1 = 8
  x->putScalar(1.0);
  y->putScalar(1.0);
  fTwoOp->apply( *x, *y, Teuchos::NO_TRANS, 3.0, 2.0 );
  // check that y == eights
  double norm = y->norm1();
  if (verbose) {
    std::cout << "Tpetra::RTI::binaryOp" << std::endl
              << "norm(y) result == " << std::setprecision(2) << std::scientific << norm
              << ", expected value is " << numGlobalRows * 8.0 << std::endl;
  }

  // Create the same diagonal operator using a Kokkos kernel
  auto kTwoOp = Tpetra::RTI::kernelOp<double>( ScaleKernel<double>(2.0), map );
  // y = 0.5*kTwop*x + 0.75*y = 0.5*2*1 + 0.75*8 = 7
  kTwoOp->apply( *x, *y, Teuchos::NO_TRANS, 0.5, 0.75 );
  // check that y == sevens
  norm = y->norm1();
  if (verbose) {
    std::cout << "Tpetra::RTI::kernelOp" << std::endl
              << "norm(y) result == " << std::setprecision(2) << std::scientific << norm
              << ", expected value is " << numGlobalRows * 7.0 << std::endl;
  }

  //
  // Create a finite-difference stencil using a Kokkos kernel and non-trivial maps
  //
  decltype(map) colmap;
  if (numImages > 1) {
    Teuchos::Array<int>           colElements;
    Teuchos::ArrayView<const int> rowElements = map->getNodeElementList();
    // This isn't the most efficient Map/Import layout, but it makes for a very straight-forward kernel
    if (myImageID != 0) colElements.push_back( map->getMinGlobalIndex() - 1 );
    colElements.insert(colElements.end(), rowElements.begin(), rowElements.end());
    if (myImageID != numImages-1) colElements.push_back( map->getMaxGlobalIndex() + 1 );
    colmap = Tpetra::createNonContigMapWithNode<int,int>(colElements(), comm);
  }
  else {
    colmap = map;
  }
  auto importer = createImport(map,colmap);
  // Finite-difference kernel = tridiag(-1, 2, -1)
  FDStencil<double> kern(myImageID, numImages, map->getNodeNumElements(), -1.0, 2.0, -1.0);
  auto FDStencilOp = Tpetra::RTI::kernelOp<double>( kern, map, map, importer );
  // x = ones(), FD(x) = [1 zeros() 1]
  auto timeFDStencil = TimeMonitor::getNewTimer("FD RTI Stencil");
  {
    TimeMonitor lcltimer(*timeFDStencil);
    for (int l=0; l != numTimeTrials; ++l) {
      FDStencilOp->apply( *x, *y );
    }
  }
  norm = y->norm1();
  if (verbose) {
    std::cout << std::endl
              << "TpetraExamples::FDStencil" << std::endl
              << "norm(y) result == " << std::setprecision(2) << std::scientific << norm
              << ", expected value is " << 2.0 << std::endl;
  }

  //
  // Create a finite-difference stencil using a CrsMatrix
  //
  auto FDMatrix = Tpetra::createCrsMatrix<double>(map);
  for (int r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
    if (r == map->getMinAllGlobalIndex()) {
      FDMatrix->insertGlobalValues(r, Teuchos::tuple<int>(r,r+1), Teuchos::tuple<double>(2.0,-1.0));
    }
    else if (r == map->getMaxAllGlobalIndex()) {
      FDMatrix->insertGlobalValues(r, Teuchos::tuple<int>(r-1,r), Teuchos::tuple<double>(-1.0,2.0));
    }
    else {
      FDMatrix->insertGlobalValues(r, Teuchos::tuple<int>(r-1,r,r+1), Teuchos::tuple<double>(-1.0,2.0,-1.0));
    }
  }
  FDMatrix->fillComplete();
  auto timeFDMatrix = TimeMonitor::getNewTimer("FD CrsMatrix");
  {
    TimeMonitor lcltimer(*timeFDMatrix);
    for (int l=0; l != numTimeTrials; ++l) {
      FDMatrix->apply(*x, *y);
    }
  }

  //
  // Print timings
  //
  if (verbose) {
    std::cout << std::endl;
    TimeMonitor::summarize( std::cout );
  }

  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}

/** \example RTIOperatorExample.cpp
    Demonstrate using Tpetra::RTI Operator wrappers.
  */

