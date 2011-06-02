#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_Version.hpp"
#include "Tpetra_DefaultPlatform.hpp"
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

}

int main(int argc, char *argv[])
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  using Teuchos::TimeMonitor;
  using TpetraExamples::FDStencil;
  using TpetraExamples::ScaleKernel;

  // 
  // Get the default communicator and node
  //
  auto &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  auto comm = platform.getComm();
  auto node = platform.getNode();
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
    std::cout << "Node type: " << Teuchos::typeName(*node) << std::endl;
    std::cout << std::endl;
  }

  // 
  // Create a simple map for domain and range
  // 
  Tpetra::global_size_t numGlobalRows = numGlobal_user;
  auto map = Tpetra::createUniformContigMapWithNode<int,int>(numGlobalRows, comm, node);
  // const size_t numLocalRows = map->getNodeNumElements();
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

