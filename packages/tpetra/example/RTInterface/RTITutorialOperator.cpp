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

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_RTIOp.hpp>

#include <iostream>
#include <iomanip>

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
  public:
    ScaleKernel(const S & gamma) : _gamma(gamma) {}
    inline void execute(int i) { 
      if (this->_alpha == Teuchos::ScalarTraits<S>::one() && this->_beta == Teuchos::ScalarTraits<S>::zero()) {
        this->_vec_inout[i] = _gamma * this->_vec_in2[i];
      }
      else {
        this->_vec_inout[i] = this->_alpha * _gamma * this->_vec_in2[i] + this->_beta * this->_vec_inout[i];
      }
    }
};

int main(int argc, char *argv[])
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  // 
  // Get the default communicator and node
  //
  auto &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  auto comm = platform.getComm();
  auto node = platform.getNode();
  const int myImageID = comm->getRank();
  const int numImages = comm->getSize();
  const bool verbose = (myImageID==0);

  // 
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl;
    std::cout << "Comm info: " << *comm;
    std::cout << "Node type: " << Teuchos::typeName(*node) << std::endl;
    std::cout << std::endl;
  }

  // 
  // Create a simple map for domain and range
  // 
  Tpetra::global_size_t numGlobalRows = 1000*numImages;
  auto map = Tpetra::createUniformContigMapWithNode<int,int>(numGlobalRows, comm, node);
  auto x = Tpetra::createVector<double>(map),
       y = Tpetra::createVector<double>(map);

  // Create a simple diagonal operator using lambda function
  auto fTwoOp = Tpetra::RTI::binaryOp<double>( [](double /*y*/, double x) { return 2.0 * x; } , map );
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
    colmap = Tpetra::createNonContigMapWithNode<int,int>(colElements(), comm, node);
  }
  else {
    colmap = map;
  }
  auto importer = createImport(map,colmap);
  // Finite-difference kernel = tridiag(-1, 2, -1)
  FDStencil<double> kern(myImageID, numImages, map->getNodeNumElements(), -1.0, 2.0, -1.0);
  auto FDStencilOp = Tpetra::RTI::kernelOp<double>( kern, map, map, importer );
  // x = ones(), FD(x) = [1 zeros() 1]
  FDStencilOp->apply( *x, *y );
  norm = y->norm1();
  if (verbose) {
    std::cout << std::endl
              << "TpetraExamples::FDStencil" << std::endl
              << "norm(y) result == " << std::setprecision(2) << std::scientific << norm 
              << ", expected value is " << 2.0 << std::endl;
  }

  std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  return 0;
}
