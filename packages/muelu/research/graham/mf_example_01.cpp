// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// Matrix-free example 01: Xpetra operator that generates 
// a simple tridiagonal finite difference poisson matrix
//
// This solves the problem f''(x) = h^2*(-pi^2)*sin(pi*x)
// on the domain [0,1] with a finite difference method.
// This supplies n global DOFs, and breaks it into
// contiguous segments depending on the number of MPI ranks.


// STL includes
#include <iostream>
#include <vector>
#include <set>
#include <stdio.h>
#include <random>

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Kokkos include
#include "Kokkos_Core.hpp"
//#include "kokkosTools.hpp"

//Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "MatrixMarket_Tpetra.hpp"

//Xpetra includes
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Operator.hpp"
#include "Xpetra_Import.hpp"
#include "Xpetra_Export.hpp"
#include "Xpetra_Utils.hpp"

// Belos includes
#include "BelosTpetraAdapter.hpp"
#include "BelosTpetraOperator.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBiCGStabSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosPCPGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockStochasticCGSolMgr.hpp"
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosTFQMRSolMgr.hpp"

// MueLu headers
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_Factory.hpp"

#include "mpi.h"

// Kokkos typedefs
typedef Kokkos::Serial HostExec;
typedef Kokkos::HostSpace HostMem;
typedef Kokkos::Device<HostExec,HostMem> HostDevice;
typedef Kokkos::Compat::KokkosSerialWrapperNode HostNode;
#ifdef DREAM_USE_CUDA
  typedef Kokkos::Cuda DeviceExec;
  typedef Kokkos::CudaSpace DeviceMem;
  typedef Kokkos::Compat::KokkosCudaWrapperNode DeviceNode;
#else
  typedef Kokkos::Serial DeviceExec;
  typedef Kokkos::HostSpace DeviceMem;
  typedef Kokkos::Compat::KokkosSerialWrapperNode DeviceNode;
#endif
typedef Kokkos::Device<DeviceExec,DeviceMem> DeviceDevice;

#define PI 3.141592653589793238463

// some things like multivectors are "2D views" but only appear as 1D in practice, so we macro a print statement
#define PRINT_VIEW2_LINEAR(view)                                                               \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << ") = [" << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i)                                      \
  for(unsigned int j=0; j<view.extent(1); ++j)                                    \
    std::cout << view(i,j) << " ";                                                \
std::cout << "]" << std::endl;

/**
 * This class defines an operator corresponding to the 
 * [-1 2 -1] tridiagonal matrix usually obtained from 
 * a 1D finite difference discretization.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TridiagonalOperator : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:
  // Xpetra::Operator subclasses should always define typedefs according to Xpetra
  typedef typename Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::node_type            node_type;
  typedef typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>                    MV;
  typedef typename Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>                                    map_type;
  typedef typename Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>                                 import_type;
  typedef typename Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node>                                 export_type;
public:
  /** Constructor
   * \param[in] n The number of global DOFs
   * \param[in] comm The Teuchos::Comm for the object
   */
  TridiagonalOperator(const GlobalOrdinal n,
                      const Teuchos::RCP<const Teuchos::Comm<int> > comm)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::invalid_argument, "TridiagonalOperator constructor: The input Comm object must be nonnull.");
    
    const int my_rank = comm->getRank();
    const int num_procs = comm->getSize();

    // Construct a default map and let it choose how DOFs are divided
    // Note: This assumes the map constructor is generating contiguous local DOFs
    const GlobalOrdinal index_base = 0;
    opMap_ = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Xpetra::UseTpetra, n, index_base, comm);
    LocalOrdinal nlocal = opMap_->getLocalNumElements();

    // Ghosting: procs 0,1,...,n_p-1 are ordered left to right on [0,1]
    if(my_rank > 0)
      ++nlocal;
    if(my_rank < num_procs - 1)
      ++nlocal;

    // Construct a list of columns where this process has nonzero elements  
    // For this tridiagonal matrix, this is firstRowItOwns-1:lastRowItOwns+1
    std::vector<GlobalOrdinal> indices;
    indices.reserve(nlocal);
    if(my_rank > 0)
      indices.push_back(opMap_->getMinGlobalIndex() - 1);
    for(GlobalOrdinal i = opMap_->getMinGlobalIndex(); i <= opMap_->getMaxGlobalIndex(); ++i)
      indices.push_back(i);
    if(my_rank < num_procs - 1)
      indices.push_back(opMap_->getMaxGlobalIndex() + 1);
    Teuchos::ArrayView<const GlobalOrdinal> element_list(indices);

    // column Map for handling the redistribution
    const GlobalOrdinal num_global_elements = n + 2*(num_procs - 1);
    redistMap_ = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Xpetra::UseTpetra, num_global_elements, element_list, index_base, comm);

    // import object that describes how data will be redistributed
    importer_ = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(opMap_, redistMap_);
  };

  // Required since we inherit from Xpetra::Operator
  // Destructor
  virtual ~TridiagonalOperator() {}
  
  /**
   * \brief Compute Y := alpha Op X + beta Y.
   * \param[in] X Vector to apply the operator to
   * \param[in] Y Vector to update
   * \param[in] mode Transpose mode (unused since this is symmetric)
   * \param[in] alpha Magnitude of Op*X term
   * \param[in] beta Magnitude of Y term
   */
  void
  apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
        Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
  {
    // Setup: get comms, ranks, procs
    Teuchos::RCP<const Teuchos::Comm<int> > comm = opMap_->getComm();
    const int my_rank = comm->getRank();
    const int num_procs = comm->getSize();
    const size_t num_vecs = X.getNumVectors();
    const LocalOrdinal numlocrows = static_cast<LocalOrdinal>(X.getLocalLength());
    
    // Make a temporary multivector for holding the redistributed data and then redistribute
    Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> redistDataX = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(redistMap_, num_vecs);
    redistDataX->doImport(X, *importer_, Xpetra::INSERT);

    // Get a view of the multivector
    auto KokkosViewX = redistDataX->getDeviceLocalView(Xpetra::Access::ReadOnly);
    auto KokkosViewY = Y.getDeviceLocalView(Xpetra::Access::ReadWrite);

    // Perform the matvec with the locally owned data
    // For each column
    for(size_t c = 0; c < num_vecs; ++c) {
      LocalOrdinal offset;
      // On ranks greater than 0, we need the ghosted X values for the computation
      // Y[0,c] = beta*Y[0,c] + alpha*(-colViewX[0] + 2*colViewX[1] - colViewX[2])
      if(my_rank > 0) {
        KokkosViewY(0,c) = beta*KokkosViewY(0, c) + alpha*(-KokkosViewX(0, c) + 2*KokkosViewX(1, c) - KokkosViewX(2, c));
        offset = 0;
      }
      // On rank 0, we only have two entries in the first row
      // Y[0,c] = beta*Y[0,c] + alpha*(2*colViewX[1] - colViewX[2])
      else {
        KokkosViewY(0,c) = beta*KokkosViewY(0, c) + alpha*(2*KokkosViewX(0, c) - KokkosViewX(1, c));
        offset = 1;
      }
      // For all other rows, we need the full stencil
      // Y[r,c] = beta*Y[r,c] + alpha*(-colViewX[r-offset] + 2*colViewX[r+1-offset] - colViewX[r+2-offset])
      for(LocalOrdinal r = 1; r < numlocrows - 1; ++r) {
        const Scalar newVal = beta*KokkosViewY(r, c) + 
                                   alpha*(-KokkosViewX(r-offset, c) + 2*KokkosViewX(r+1-offset, c) - KokkosViewX(r+2-offset, c));
        KokkosViewY(r,c) = newVal;
      }
      // On ranks other than the last rank, we need the ghosted X values for the computation
      // Y[numlocrows-1,c] = beta*Y[numlocrows-1,c] + alpha*(-colViewX[numlocrows-1-offset] + 2*colViewX[numlocrows-offset]
      //                   - colViewX[numlocrows+1-offset])
      if(my_rank < num_procs - 1) {
        const Scalar newVal = beta*KokkosViewY(numlocrows-1, c) + 
                                   alpha*(-KokkosViewX(numlocrows-1-offset, c) + 2*KokkosViewX(numlocrows-offset, c)
                                          - KokkosViewX(numlocrows+1-offset, c));
        KokkosViewY(numlocrows-1,c) = newVal;
      }
      // On the last rank, we only have two entries in the last row
      // Y[numlocrows-1,c] = beta*Y[numlocrows-1,c] + alpha*(-colViewX[numlocrows-1-offset] + 2*colViewX[numlocrows-offset])
      else {
        const Scalar newVal = beta*KokkosViewY(numlocrows-1, c) + 
          alpha*(-KokkosViewX(numlocrows-1-offset, c) + 2*KokkosViewX(numlocrows-offset, c));
        KokkosViewY(numlocrows-1,c) = newVal;
      }
    }
  }

  //! Returns the Xpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const { return opMap_; }

  //! Returns the Xpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const { return opMap_; }

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const { return true; }

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & R) const {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    R.update(STS::one(), B, STS::zero()); // R = 1*B + 0*R
    this->apply(X, R, Teuchos::NO_TRANS, -STS::one(), STS::one()); // R = R - Op*X
  }

private:
  Teuchos::RCP<const map_type> opMap_, redistMap_;
  Teuchos::RCP<const import_type> importer_;
};


/**
 * This class defines an operator corresponding to the 
 * [1 2 1] interpolation stencil. For interpolation to
 * behave correctly, each MPI rank should have a
 * number of DOFs that is divisible by 3, as the 
 * coarsening ratio is 3->1. Communication is ignored.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MFProlongatorOperator : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:
  // Xpetra::Operator subclasses should always define typedefs according to Xpetra
  typedef typename Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::node_type            node_type;
  typedef typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>                    MV;
  typedef typename Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>                                    map_type;
  typedef typename Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>                                 import_type;
  typedef typename Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node>                                 export_type;
public:
  /** Constructor
   * \param[in] n The number of global DOFs
   * \param[in] comm The Teuchos::Comm for the object
   */
  MFProlongatorOperator(const Teuchos::RCP<const map_type> fine_map)
  {
    const GlobalOrdinal n = fine_map->getGlobalNumElements();
    const LocalOrdinal n_local = fine_map->getLocalNumElements();
    const Teuchos::RCP<const Teuchos::Comm<int> > comm = fine_map->getComm();

    TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::invalid_argument, "MFProlongatorOperator constructor: The input Comm object must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION(n_local % 3 != 0, std::invalid_argument, "MFProlongatorOperator constructor: The number of local DOFs is not divisible by 3.");

    // Construct a default map and let it choose how DOFs are divided
    // Note: This assumes the map constructor is generating contiguous local DOFs
    const GlobalOrdinal index_base = 0;
    rangeMap_ = fine_map;
    domainMap_ = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Xpetra::UseTpetra, n/3, index_base, comm);

    TEUCHOS_TEST_FOR_EXCEPTION(domainMap_->getLocalNumElements() != ((size_t) n_local)/3, std::invalid_argument, "MFProlongatorOperator constructor: The operator is not coarsening by 3.");
  };

  // Required since we inherit from Xpetra::Operator
  // Destructor
  virtual ~MFProlongatorOperator() {}
  
  /**
   * \brief Compute Y := alpha Op X + beta Y.
   * \param[in] X Vector to apply the operator to
   * \param[in] Y Vector to update
   * \param[in] mode Transpose mode (unused since this is symmetric)
   * \param[in] alpha Magnitude of Op*X term
   * \param[in] beta Magnitude of Y term
   */
  void
  apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
        Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
  {
    // // Setup: get comms, ranks, procs
    // Teuchos::RCP<const Teuchos::Comm<int> > comm = domainMap_->getComm();
    // const int my_rank = comm->getRank();
    // const int num_procs = comm->getSize();
    // const size_t num_vecs = X.getNumVectors();
    // const LocalOrdinal numlocrows = static_cast<LocalOrdinal>(X.getLocalLength());
    
    // // Make a temporary multivector for holding the redistributed data and then redistribute
    // Teuchos::RCP<MV> redistDataX = Teuchos::rcp(new MV(redistMap_, num_vecs));
    // redistDataX->doImport(X, *importer_, Xpetra::INSERT);

    // // Get a view of the multivector
    // auto KokkosViewX = redistDataX->getLocalView<HostDevice>(Xpetra::Access::ReadOnly);
    // auto KokkosViewY = Y.getLocalView<HostDevice>(Xpetra::Access::ReadWrite);

    // // Perform the matvec with the locally owned data
    // // For each column
    // for(size_t c = 0; c < num_vecs; ++c) {
    //   LocalOrdinal offset;
    //   // On ranks greater than 0, we need the ghosted X values for the computation
    //   // Y[0,c] = beta*Y[0,c] + alpha*(-colViewX[0] + 2*colViewX[1] - colViewX[2])
    //   if(my_rank > 0) {
    //     KokkosViewY(0,c) = beta*KokkosViewY(0, c) + alpha*(-KokkosViewX(0, c) + 2*KokkosViewX(1, c) - KokkosViewX(2, c));
    //     offset = 0;
    //   }
    //   // On rank 0, we only have two entries in the first row
    //   // Y[0,c] = beta*Y[0,c] + alpha*(2*colViewX[1] - colViewX[2])
    //   else {
    //     KokkosViewY(0,c) = beta*KokkosViewY(0, c) + alpha*(2*KokkosViewX(0, c) - KokkosViewX(1, c));
    //     offset = 1;
    //   }
    //   // For all other rows, we need the full stencil
    //   // Y[r,c] = beta*Y[r,c] + alpha*(-colViewX[r-offset] + 2*colViewX[r+1-offset] - colViewX[r+2-offset])
    //   for(LocalOrdinal r = 1; r < numlocrows - 1; ++r) {
    //     const Scalar newVal = beta*KokkosViewY(r, c) + 
    //                                alpha*(-KokkosViewX(r-offset, c) + 2*KokkosViewX(r+1-offset, c) - KokkosViewX(r+2-offset, c));
    //     KokkosViewY(r,c) = newVal;
    //   }
    //   // On ranks other than the last rank, we need the ghosted X values for the computation
    //   // Y[numlocrows-1,c] = beta*Y[numlocrows-1,c] + alpha*(-colViewX[numlocrows-1-offset] + 2*colViewX[numlocrows-offset]
    //   //                   - colViewX[numlocrows+1-offset])
    //   if(my_rank < num_procs - 1) {
    //     const Scalar newVal = beta*KokkosViewY(numlocrows-1, c) + 
    //                                alpha*(-KokkosViewX(numlocrows-1-offset, c) + 2*KokkosViewX(numlocrows-offset, c)
    //                                       - KokkosViewX(numlocrows+1-offset, c));
    //     KokkosViewY(numlocrows-1,c) = newVal;
    //   }
    //   // On the last rank, we only have two entries in the last row
    //   // Y[numlocrows-1,c] = beta*Y[numlocrows-1,c] + alpha*(-colViewX[numlocrows-1-offset] + 2*colViewX[numlocrows-offset])
    //   else {
    //     const Scalar newVal = beta*KokkosViewY(numlocrows-1, c) + 
    //       alpha*(-KokkosViewX(numlocrows-1-offset, c) + 2*KokkosViewX(numlocrows-offset, c));
    //     KokkosViewY(numlocrows-1,c) = newVal;
    //   }
    // }
  }

  //! Returns the Xpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const { return domainMap_; }

  //! Returns the Xpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const { return rangeMap_; }

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const { return true; }

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & R) const {
    //throw Exceptions::RuntimeError("Interface not supported");
  }

private:
  Teuchos::RCP<const map_type> domainMap_, rangeMap_, redistMap_;
  Teuchos::RCP<const import_type> importer_;
};


int main(int argc, char* argv[])
{
  // MPI initialization using Teuchos
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  // Xpetra nodes call Kokkos::execution_space::initialize if the execution
  // space is not initialized, but they don't call Kokkos::initialize.
  // Teuchos::GlobalMPISession captures its command-line arguments for later
  // use that Xpetra takes advantage of.
  //
  // We call Kokkos::initialize() after MPI so that MPI has the chance to bind
  // processes correctly before Kokkos touches things.
  Kokkos::initialize(argc, argv);
  // Boilerplate  MPI/Kokkos initialization
  Teuchos::RCP<Teuchos::MpiComm<int>> comm = Teuchos::rcp(new typename Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  const int my_rank = comm->getRank();
  const int num_procs = comm->getSize();
  {
    // Necessary typedefs
    using SC = double;
    using LO = int;
    using GO = long long;
    using NO = Tpetra::MultiVector<>::node_type;
    //using map_type = Xpetra::Map<>; // unused
    using MV = Xpetra::MultiVector<SC,LO,GO,NO>;

    // Set a command line processor and parse it
    Teuchos::CommandLineProcessor clp(false);
    int n = 300;
    int max_iterations = 1000;
    double tol = 1e-10;
    bool do_multigrid = true;
    bool show_timer_summary = false;
    bool belos_verbose = false;
    bool show_kokkos = false; 
    bool print_RHS_and_solution = false;
    bool dump_matrix_market = false;
    // problem config
    clp.setOption("n", &n, "Size of the n-by-n operator (default: 300)");
    clp.setOption("verbose", "no-verbose", &belos_verbose, "Use Belos verbose output (default: false)");
    // linear solver config
    clp.setOption("max-its", &max_iterations, "Maximum number of Belos iterations (default: 1000)");
    clp.setOption("tol", &tol, "Belos convergence tolerance (default: 1e-10)");
    // muelu-specific config
    clp.setOption("multigrid", "no-multigrid", &do_multigrid, "Use MueLu to solve this via multigrid (default: true)");
    // other print statements and dumps
    clp.setOption("show-timer-summary", "no-timer-summary", &show_timer_summary, "Show a summary of the timer information (default: false)");
    clp.setOption("show-kokkos", "no-kokkos", &show_kokkos, "Show the Kokkos configuration (default: false) [currently does nothing]");
    clp.setOption("print-RHS-and-solution", "no-print-RHS-and-solution", &print_RHS_and_solution, "Print the RHS vector and the solution vector to verify the answer to the problem (default: false)");
    clp.setOption("dump-matrix-market", "no-dump-matrix-market", &dump_matrix_market, "Dump the solution and right-hand side to MatrixMarket format (default: false)");

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: 
        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: 
        return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
        break;
    }
    if(my_rank == 0)
      std::cout << "Running example-01 with n=" << n << " verbose=" << belos_verbose << " config=" << show_kokkos << "..." << std::endl; 

    // print configuration details if needed
    Kokkos::Serial::print_configuration(std::cout, true/*details*/);
    //Kokkos::OpenMP::print_configuration(std::cout, true/*details*/);
    //std::cout << "OpenMP Max Threads = " << omp_get_max_threads() << std::endl;
    //Kokkos::Cuda::print_configuration(std::cout, true/*details*/);
    //Kokkos::Experimental::HIP::print_configuration(std::cout, true/*details*/);

    // Create the operator
    Teuchos::RCP<TridiagonalOperator<SC,LO,GO,NO>> matrix = Teuchos::rcp(new TridiagonalOperator<SC,LO,GO,NO>(n, comm));

    // A useful class name means excruciating template arguments... but I want to avoid headaches
    std::cout << "Creating the vectors..." << std::endl;
    
    // Construct the right-hand side
    Teuchos::RCP<MV> rhs = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(matrix->getRangeMap(),1);
    rhs->putScalar(0.0);

    // Construct the initial guess (seedrandom is not always reproducible across machines, see Github)
    Teuchos::RCP<MV> solution = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(matrix->getDomainMap(),1);
    Teuchos::ScalarTraits<SC>::seedrandom(314159);
    solution->randomize();
    // solution->putScalar(1.0) is an alternative if the random seed isn't reproducible

    // some computations to make the RHS match the example
    // f''(x) = h^2*(-pi^2)*sin(pi*x)
    {
      // divide domain up according to MPI ranks
      const SC dx = 1.0/num_procs;
      const SC x_left = dx*my_rank;
      const SC x_right = dx*(my_rank+1);
      const size_t n_local = matrix->getRangeMap()->getLocalNumElements();
      const SC h = (x_right - x_left)/(n_local-1);

      // fill the RHS appropriately
      auto rhs_2d = rhs->getHostLocalView(Xpetra::Access::OverwriteAll);
      auto rhs_1d = Kokkos::subview (rhs_2d, Kokkos::ALL(), 0);
      SC x = x_left;
      for(size_t i=0; i<n_local; ++i) {
        rhs_1d(i) = 4*h*h*PI*PI*sin(2*PI*x);
        x += h;
      }
    }

    if(!do_multigrid) {
      // TODO: Artifact of the Tpetra->Xpetra switch
      // std::cout << "Creating Belos solver..." << std::endl;

      // // Create the solver
      // std::string solver_name = "Tpetra CG";
      // Teuchos::RCP<Belos::SolverManager<SC, Tpetra::MultiVector<SC,LO,GO,NO>, Tpetra::Operator<SC,LO,GO,NO>> solver;
      // Belos::SolverFactory<SC, Tpetra::MultiVector<SC,LO,GO,NO>, Tpetra::Operator<SC,LO,GO,NO>> factory;
      // solver = factory.create(solver_name, Teuchos::null);

      // // Add settings to the solver
      // Teuchos::RCP<Teuchos::ParameterList> belos_settings = Teuchos::parameterList("Belos");
      // belos_settings->set("Verbosity", belos_verbose ? 1 : 0);
      // belos_settings->set("Output Style", 1);
      // belos_settings->set("Maximum Iterations",    max_iterations);
      // //belos_settings->set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
      // //belos_settings->set("Output Frequency",      1);
      // //belos_settings->set("Output Style",          Belos::Brief);
      // solver->setParameters(belos_settings);

      // // Define the linear problem
      // Teuchos::RCP<Belos::LinearProblem<SC, MV, OP>> lp = 
      //   Teuchos::rcp(new Belos::LinearProblem<SC, MV, OP>(matrix, solution, rhs));
      // lp->setProblem();

      // // Solve the problem
      // std::cout << "Solving the problem..." << std::endl;
      // solver->setProblem(lp);
      // const Belos::ReturnType belos_result = solver->solve();

      // if(my_rank == 0) {
      //   std::cout << "Belos solver wrapper results: "
      //             << (belos_result == Belos::Converged ? "Converged" : "Unconverged")
      //             << std::endl
      //             << "Number of iterations: " << solver->getNumIters()
      //             << std::endl;
      // }
    } else {
      std::cout << "Generating multigrid objects..." << std::endl;
      
      // generate coarse matrix-free operator
      Teuchos::RCP<TridiagonalOperator<SC,LO,GO,NO>> coarse_matrix = Teuchos::rcp(new TridiagonalOperator<SC,LO,GO,NO>(n/3, comm));
      Teuchos::RCP<MFProlongatorOperator<SC,LO,GO,NO>> P = Teuchos::rcp(new MFProlongatorOperator<SC,LO,GO,NO>(matrix->getDomainMap()));

      // create MueLu hierarchy and levels
      Teuchos::RCP<MueLu::Hierarchy<SC,LO,GO,NO>> H = Teuchos::rcp(new typename MueLu::Hierarchy<SC,LO,GO,NO>());
      Teuchos::RCP<MueLu::Level> fineLevel = H->GetLevel(0);
      fineLevel->Set("A", matrix);
      
      H->AddNewLevel();
      Teuchos::RCP<MueLu::Level> coarseLevel = H->GetLevel(1);
      coarseLevel->Set("P", P);
      //coarseLevel->Set("A", coarse_matrix);

      Teuchos::ParameterList params;
      params.set("coarse: max size", 1);
      params.set("max levels", 2);
      params.set("transpose: use implicit", true);

      Teuchos::RCP<MueLu::HierarchyManager<SC,LO,GO,NO>> mueLuFactory = Teuchos::rcp(new MueLu::ParameterListInterpreter<SC,LO,GO,NO>(params,matrix->getDomainMap()->getComm()));
      //H->setlib(matrix->getDomainMap()->lib());
      H->SetProcRankVerbose(matrix->getDomainMap()->getComm()->getRank());
      mueLuFactory->SetupHierarchy(*H);

      std::cout << "Finished hierarchy!" << std::endl;

      Teuchos::ParameterList status;
      //status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
      if (comm->getRank() == 0) {
        std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
        status.print(std::cout, Teuchos::ParameterList::PrintOptions().indent(2));
      }

      //H->Iterate(*rhs, *solution, max_iterations);
    }

    // output the RHS and solution for validation
    // (sleep my_rank is very hacky here, but keeps prints contiguous)
    if(print_RHS_and_solution) {
      sleep(my_rank);
      auto rhs_2d = rhs->getHostLocalView(Xpetra::Access::ReadOnly);
      PRINT_VIEW2_LINEAR(rhs_2d)
      auto solution_2d = solution->getHostLocalView(Xpetra::Access::ReadOnly);
      PRINT_VIEW2_LINEAR(solution_2d) 
    }
    
    if(dump_matrix_market) {
      //Xpetra::MatrixMarket::Writer<vector_type>::writeDenseFile("example_01_solution.mm", *solution);
      //Xpetra::MatrixMarket::Writer<vector_type>::writeDenseFile("example_01_rhs.mm", *rhs);
    }

  }
  Kokkos::finalize();
  return 0;
}
