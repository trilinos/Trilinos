// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NodeT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DIScalarLAIMED. IN Node EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NodeT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GlobalOrdinalODS OR SERVICES; LocalOrdinalSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

// Class implementing our problem
#include "twoD_diffusion_problem_tpetra.hpp"

// solver
#include "Ifpack2_Factory.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Utilities
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"

//MueLu headers
#include <MueLu.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>  // => This header defines Belos::MueLuOp
#include <BelosLinearProblem.hpp>

//Tpetra node conversion
#include <Kokkos_SerialNode.hpp>
#include <Kokkos_TPINode.hpp>
#include <Kokkos_ThrustGPUNode.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "KokkosCore_config.h"
#include "Kokkos_Cuda.hpp"

#include <iostream>

using Teuchos::ArrayRCP;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::parameterList;

using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;
using Tpetra::Operator;
using Ifpack2::Preconditioner;

// Krylov preconditioning approaches
enum Prec { ILU, CHEBY, MG };
const int num_prec = 3;
const Prec prec_values[] = {ILU, CHEBY, MG};
const char *prec_names[] = { "ILU", "Chebyshev", "Multigrid"};
template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node,
          typename LocalMatOps,
          typename CloneNode,
          typename CloneLocalMatOps>
Scalar MGclone_and_solve(
  RCP< Xpetra::Matrix <Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mueluJ,
  RCP< const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > f,
  RCP< const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > dx,
  RCP< MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > H,
  ParameterList& nodeParams,
  RCP<ParameterList> belosParams,
  bool symmetric)
{
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode, CloneLocalMatOps> Clone_Mat;
  typedef Map<LocalOrdinal,GlobalOrdinal,CloneNode> Clone_Map;
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode> Clone_MV;
  typedef Belos::OperatorT<Clone_MV> Clone_OP;
  typedef Belos::LinearProblem<Scalar,Clone_MV,Clone_OP> BLinProb;
  typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode, CloneLocalMatOps> Clone_Hierarchy;
  typedef Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, CloneNode, CloneLocalMatOps> CloneBelos_MueLuOperator;

  //Convert J to clone node type
  RCP<ParameterList> plClone = parameterList();
  RCP<CloneNode> node_clone = rcp(new CloneNode(nodeParams));
  RCP<Clone_Mat> J_clone = Xpetra::clone(*mueluJ, node_clone);

  //Clone the hierarchy
  RCP<Clone_Hierarchy> H_clone = H->template clone< CloneNode, CloneLocalMatOps >(node_clone);

  //Convert f, dx to clone node type
  RCP<Clone_MV> f_clone = f->clone(node_clone);
  RCP<Clone_MV> dx_clone = dx->clone(node_clone);
  dx_clone->putScalar(Scalar(0.0));
  // Define Operator and Preconditioner
  RCP<Clone_OP> OP_clone = rcp(new Belos::XpetraOp<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode,CloneLocalMatOps>(J_clone));  // Turns a Xpetra::Matrix object into a Belos operator
  const RCP<const CloneBelos_MueLuOperator> M_clone = rcp(new CloneBelos_MueLuOperator(H_clone)); // Turns a MueLu::Hierarchy object into a Belos operator

  //Create problem for clone node
  RCP< BLinProb > problem = rcp(new BLinProb(OP_clone, dx_clone,f_clone));
  problem->setRightPrec(M_clone);
  problem->setProblem();

  // Create solver for clone node type
  RCP<Belos::SolverManager<Scalar,Clone_MV,Clone_OP> > solver;
  if (symmetric)
    solver =
      rcp(new Belos::PseudoBlockCGSolMgr<Scalar,Clone_MV,Clone_OP>(
            problem, belosParams));
  else
    solver =
      rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,Clone_MV,Clone_OP>(
            problem, belosParams));

  // Solve linear system for clone node
  solver->solve();

  // Compute norm of difference
  ParameterList serial_params;
  RCP<Node> node = rcp(new Node(serial_params));
  RCP<MV> dx_cpu = dx_clone->clone(node);
  dx_cpu->update(1.0, *dx, -1.0);
  Scalar norm;
  dx_cpu->norm2(Teuchos::arrayView(&norm,1));

  return norm;
}


template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node,
          typename CloneNode>
Scalar clone_and_solve(
  RCP< const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > J,
  RCP< const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > f,
  RCP< const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > dx,
  RCP< Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node> > M,
  ParameterList& nodeParams,
  RCP<ParameterList> belosParams,
  ParameterList& precParams,

  bool symmetric)
{
  typedef CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Mat;
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode> Clone_Mat;
  typedef Map<LocalOrdinal,GlobalOrdinal,CloneNode> Clone_Map;
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode> Clone_MV;
  typedef Operator<Scalar,LocalOrdinal,GlobalOrdinal,CloneNode> Clone_OP;
  typedef Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, CloneNode> Clone_Prec;
  typedef Belos::LinearProblem<Scalar,Clone_MV,Clone_OP> BLinProb;

  //Convert J to clone node type
  RCP<ParameterList> plClone = parameterList();
  RCP<CloneNode> node_clone = rcp(new CloneNode(nodeParams));
  RCP<const Clone_Mat> J_clone = J->clone(node_clone, plClone);
  J_clone->print(std::cout);

  //Clone preconditioner
  Ifpack2::Factory factory;
  RCP<Clone_Prec> M_clone = factory.clone<Mat, Clone_Mat>(M, J_clone, precParams);

  //Convert f, dx to clone node type
  RCP<Clone_MV> f_clone = f->clone(node_clone);
  RCP<Clone_MV> dx_clone = dx->clone(node_clone);
  dx_clone->putScalar(Scalar(0.0));
  //Create problem for clone node
  RCP<BLinProb> problem = rcp(new BLinProb(J_clone, dx_clone, f_clone));
  problem->setRightPrec(M_clone);
  problem->setProblem();

  // Create solver for clone node type
  RCP<Belos::SolverManager<Scalar,Clone_MV,Clone_OP> > solver;
  if (symmetric)
    solver =
      rcp(new Belos::PseudoBlockCGSolMgr<Scalar,Clone_MV,Clone_OP>(
            problem, belosParams));
  else
    solver =
      rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,Clone_MV,Clone_OP>(
            problem, belosParams));

  // Solve linear system for clone node
  solver->solve();

  // Compute norm of difference
  ParameterList serial_params;
  RCP<Node> node = rcp(new Node(serial_params));
  RCP<MV> dx_cpu = dx_clone->clone(node);
  dx_cpu->update(1.0, *dx, -1.0);
  Scalar norm;
  dx_cpu->norm2(Teuchos::arrayView(&norm,1));

  return norm;
}

int main(int argc, char *argv[]) {
  typedef double Scalar;
  typedef double MeshScalar;
  typedef double BasisScalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef KokkosClassic::TPINode TPINode;
  typedef KokkosClassic::ThrustGPUNode GPUNode;
  typedef KokkosClassic::DefaultKernels<void, LocalOrdinal, Node>::SparseOps LocalMatOps;
  typedef KokkosClassic::DefaultKernels<void, LocalOrdinal, TPINode>::SparseOps TPILocalMatOps;
  typedef KokkosClassic::DefaultKernels<void, LocalOrdinal, GPUNode>::SparseOps GPULocalMatOps;
  typedef Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Belos::OperatorT<MV> BelosOP;
  typedef Belos::LinearProblem<Scalar,MV,OP> BLinProb;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node> Map;
  typedef Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Belos_MueLuOperator;
  // Start up MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    // Create comm
    RCP<const Teuchos::Comm<LocalOrdinal> > comm =
      Teuchos::DefaultComm<int>::getComm();
    int my_rank = comm->getRank();

    // Parse command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example solves a diffusion problem with different node types.\n");
    bool gpu = true;
    CLP.setOption("gpu", "no-gpu", &gpu, "Enable GPU solve");
    bool tpi = true;
    CLP.setOption("tpi", "no-tpi", &tpi, "Enable TPI solve");
    LocalOrdinal n = 16;
    CLP.setOption("num_mesh", &n, "Number of mesh points in each direction");
    bool symmetric = true;
    CLP.setOption("symmetric", "unsymmetric", &symmetric,
                  "Make matrix symmetric by eliminating Dirichlet BCs");
    int num_threads = 0;
    CLP.setOption("num_threads", &num_threads,
                  "Number of threads for TPI node");
    int ranks_per_node = 1;
    CLP.setOption("ranks_per_node", &ranks_per_node,
                  "Number of MPI ranks per node");
    int gpu_ranks_per_node = 1;
    CLP.setOption("gpu_ranks_per_node", &gpu_ranks_per_node,
                  "Number of MPI ranks per node for GPUs");
    int device_offset = 0;
    CLP.setOption("device_offset", &device_offset,
                  "Offset for attaching MPI ranks to CUDA devices");
    int chebyDeg = 5;
    CLP.setOption("chebyDeg", &chebyDeg, "Chebyshev degree");

    Prec precMethod = MG;
    CLP.setOption("prec_method", &precMethod,
                  num_prec, prec_values, prec_names,
                  "Preconditioner method");


    CLP.parse( argc, argv );
    if (my_rank == 0)
      std::cout << "Summary of command line options:" << std::endl
                << "\tgpu                = " << gpu << std::endl
                << "\ttpi                = " << tpi << std::endl
                << "\tnum_mesh           = " << n << std::endl
                << "\tsymmetric          = " << symmetric << std::endl
                << "\tnum_threads        = " << num_threads << std::endl
                << "\tranks_per_node     = " << ranks_per_node << std::endl
                << "\tgpu_ranks_per_node = " << gpu_ranks_per_node << std::endl
                << "\tdevice_offset      = " << device_offset << std::endl
                << "\tpreconditioner     = " << prec_names[precMethod] << std::endl
                << "\tChebyshev poly deg = " << chebyDeg << std::endl;

    // Create application
    typedef twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,Node> problem_type;
    LocalOrdinal num_KL = 2;
    BasisScalar s = 0.1, mu = 0.2;
    bool log_normal = false;
    // There seems to be a problem with the symmetric formulation
    // but CG appears to work with the unsymmetric form anyway
    RCP<problem_type> model =
      rcp(new problem_type(comm, n, num_KL, s, mu, log_normal, false));

    // Create vectors and operators
    typedef problem_type::Tpetra_Vector Tpetra_Vector;
    typedef problem_type::Tpetra_CrsMatrix Tpetra_CrsMatrix;
    typedef Tpetra::MatrixMarket::Writer<Tpetra_CrsMatrix> Writer;
    RCP<const Tpetra_Vector> p = model->get_p_init(0);
    RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(model->get_x_map());
    x->putScalar(0.0);
    RCP<Tpetra_Vector> f = Tpetra::createVector<Scalar>(model->get_f_map());
    RCP<Tpetra_Vector> dx = Tpetra::createVector<Scalar>(model->get_x_map());
    RCP<Tpetra_CrsMatrix> J = model->create_W();

    // Evaluate model
    model->computeResidual(*x, *p, *f);
    model->computeJacobian(*x, *p, *J);

    //Create preconditioners
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    RCP<Tprec> M;
    Ifpack2::Factory factory;
    ParameterList precParams;

    //Create RILUK preconditioner using Prec factory
    if(precMethod==ILU){
      std::string prec_name = "RILUK";
      precParams.set("fact: iluk level-of-fill", 1);
      precParams.set("fact: iluk level-of-overlap", 0);
      M = factory.create<Tpetra_CrsMatrix>(prec_name, J);
      M->setParameters(precParams);
      M->initialize();
      M->compute();
    }

    //Create Chebyshev preconditioner using Prec factory
    else if(precMethod==CHEBY){
      M = factory.create<Tpetra_CrsMatrix>("CHEBYSHEV", J);
      precParams.set("chebyshev: degree", (LocalOrdinal) chebyDeg);
      precParams.set("chebyshev: ratio eigenvalue", (Scalar) 7);
      M->setParameters(precParams);
      M->initialize();
      M->compute();
    }

    RCP<ParameterList> belosParams = rcp(new ParameterList);
    belosParams->set("Num Blocks", 20);
    belosParams->set("Convergence Tolerance", 1e-12);
    belosParams->set("Maximum Iterations", 1000);
    belosParams->set("Verbosity", 33);
    belosParams->set("Output Style", 1);
    belosParams->set("Output Frequency", 10);

    if(precMethod == MG){
      //Create MueLu precondtioner with coasre grid solver and smoother = Chebyshev
      // Turns a Tpetra::CrsMatrix into a MueLu::Matrix
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mueluA_ = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(J));
      RCP<Xpetra::Matrix <Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mueluA  = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(mueluA_));
      RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > H = rcp(new MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(mueluA));
      H->setDefaultVerbLevel(Teuchos::VERB_MEDIUM);
      MueLu::FactoryManager<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> FM;
      //Setup smoother (chebyshev)
      std::string ifpackType;
      Teuchos::ParameterList ifpackList;
      ifpackType = "CHEBYSHEV";
      ifpackList.set("chebyshev: degree", (LocalOrdinal) chebyDeg);
      ifpackList.set("chebyshev: ratio eigenvalue", (Scalar) 7);
      RCP<MueLu::SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > smootherPrototype = rcp(new MueLu::Ifpack2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(ifpackType, ifpackList));
      FM.SetFactory("Smoother", rcp(new MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(smootherPrototype)));

      //Setup coarse grid smoother (chebyshev)
      Teuchos::ParameterList coarsestSmooList;
      RCP<MueLu::SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > coarsestSmooProto = rcp( new MueLu::Ifpack2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("RILUK",coarsestSmooList) );
      FM.SetFactory("CoarseSolver", rcp(new MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(smootherPrototype)));

      //Setup hierarchy with startlevel=0 and max levels = 10
      H->Setup(FM, 0, 10);
      RCP<BelosOP> belosOp =
        rcp(new Belos::XpetraOp<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(mueluA));
      const RCP<const Belos_MueLuOperator> belosPrec =
        rcp(new Belos_MueLuOperator(H));
      RCP< Belos::LinearProblem<Scalar, MV, BelosOP > > MGproblem =
        rcp(new Belos::LinearProblem<Scalar, MV, BelosOP >(belosOp, dx,f));
      MGproblem->setRightPrec(belosPrec);
      MGproblem->setProblem();

      //Create solver
      RCP<Belos::SolverManager<Scalar,MV,BelosOP> > MGsolver;
      if (symmetric)
        MGsolver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,MV,BelosOP>(MGproblem,
                                                                belosParams));
      else
        MGsolver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV,BelosOP>(MGproblem,
                                                                 belosParams));
      // Solve linear system
      if (my_rank == 0)
        std::cout << "Solving with default node..." << std::endl;
      MGsolver->solve();
      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();

      Scalar tpi_norm = 0.0, gpu_norm = 0.0;

      // Solve linear system with TPI node type
      if (tpi) {
        ParameterList node_params;
        node_params.set("Verbose", 1);
        node_params.set("Num Threads", num_threads);

        if (my_rank == 0)
            std::cout << "Solving with TPI node..." << std::endl;

        tpi_norm =
           MGclone_and_solve<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps,TPINode,TPILocalMatOps>(
                mueluA, f, dx, H, node_params, belosParams, symmetric);
        if (my_rank == 0)
          std::cout << "\nNorm of serial node soln - tpi node soln = "
                  << tpi_norm << std::endl;

        Teuchos::TimeMonitor::summarize(std::cout);
        Teuchos::TimeMonitor::zeroOutTimers();
      }

     // Solve linear system with GPU node type
    if (gpu) {
      if (my_rank == 0)
        std::cout << "Solving with GPU node..." << std::endl;

      // Compute whether we are a CPU or GPU rank, and GPU device ID
      // The first gpu_ranks_per_node ranks are associated with GPUs
      // (may not be the best choice performance-wise)
      int num_ranks = comm->getSize();
      int num_node = num_ranks / ranks_per_node;
      int node_rank = num_node == 1 ? my_rank : my_rank % num_node;
      bool gpu_rank = node_rank < gpu_ranks_per_node;
      int num_device; cudaGetDeviceCount(&num_device);
      int device_id = node_rank + device_offset;
      TEUCHOS_TEST_FOR_EXCEPTION(
        num_node*ranks_per_node != num_ranks, std::logic_error,
        "ranks_per_node does not evenly divide num_ranks");
      TEUCHOS_TEST_FOR_EXCEPTION(
        gpu_ranks_per_node > num_device, std::logic_error,
        "gpu_ranks_per_node cannot exceed number of GPU devices");
      TEUCHOS_TEST_FOR_EXCEPTION(
        gpu_rank && device_id > num_device, std::logic_error,
        "Invalid device ID " << device_id << ".  You probably are trying" <<
        " to run with too many MPI ranks");

      ParameterList node_params;
      if (gpu_rank) {
         std::cout << "MPI rank " << my_rank
                   << ":  Attached to GPU " << device_id
                   << std::endl;

        node_params.set("Verbose", 1);
        node_params.set("Device Number", device_id);

#if TPETRA_USE_KOKKOS_DISTOBJECT && defined(KOKKOS_HAVE_CUDA)
        // Initialize Cuda
        if (!Kokkos::Cuda::is_initialized())
          Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_id) );
#endif

        gpu_norm =
          MGclone_and_solve<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps,GPUNode,GPULocalMatOps>(
            mueluA, f, dx, H, node_params, belosParams, symmetric);

#if TPETRA_USE_KOKKOS_DISTOBJECT && defined(KOKKOS_HAVE_CUDA)
        Kokkos::Cuda::finalize();
#endif
      }
      else {
        // Note for the non-GPU ranks, we still have to clone since new
        // matrices are created which call fillComplete() (which involves
        // communication)
         MGclone_and_solve<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps,Node,LocalMatOps>(
            mueluA, f, dx, H, node_params, belosParams, symmetric);
      }

      if (my_rank == 0)
        std::cout << "\nNorm of serial node soln - gpu node soln = "
                  << gpu_norm << std::endl;

      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();
      }

    //Determine if example passed
    bool passed = false;
    if (gpu_norm < Scalar(1e-12) && tpi_norm < Scalar(1e-12))
      passed = true;
    if (my_rank == 0) {
      if (passed)
        std::cout << "Example Passed!" << std::endl;
      else
        std::cout << "Example Failed!" << std::endl;
    }

    }
    else {
      RCP< BLinProb > problem = rcp(new BLinProb(J, dx, f));
      problem->setRightPrec(M);
      problem->setProblem();

      //Create solver
      RCP<Belos::SolverManager<Scalar,MV,OP> > solver;
      if (symmetric)
        solver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,MV,OP>(problem,
                                                                belosParams));
      else
        solver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV,OP>(problem,
                                                                 belosParams));

      // Solve linear system
      if (my_rank == 0)
        std::cout << "Solving with default node..." << std::endl;
      solver->solve();
      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();

      Scalar tpi_norm = 0.0, gpu_norm = 0.0;

      // Solve linear system with TPI node type
      if (tpi) {
        ParameterList node_params;
        node_params.set("Verbose", 1);
        node_params.set("Num Threads", num_threads);

        if (my_rank == 0)
          std::cout << "Solving with TPI node..." << std::endl;

        tpi_norm =
          clone_and_solve<Scalar,LocalOrdinal,GlobalOrdinal,Node,TPINode>(
            J, f, dx, M, node_params, belosParams, precParams, symmetric);
        if (my_rank == 0)
          std::cout << "\nNorm of serial node soln - tpi node soln = "
                  << tpi_norm << std::endl;

      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // Solve linear system with GPU node type
    if (gpu) {
      if (my_rank == 0)
        std::cout << "Solving with GPU node..." << std::endl;

      // Compute whether we are a CPU or GPU rank, and GPU device ID
      // The first gpu_ranks_per_node ranks are associated with GPUs
      // (may not be the best choice performance-wise)
      int num_ranks = comm->getSize();
      int num_node = num_ranks / ranks_per_node;
      int node_rank = num_node == 1 ? my_rank : my_rank % num_node;
      bool gpu_rank = node_rank < gpu_ranks_per_node;
      int num_device; cudaGetDeviceCount(&num_device);
      int device_id = node_rank + device_offset;
      TEUCHOS_TEST_FOR_EXCEPTION(
        num_node*ranks_per_node != num_ranks, std::logic_error,
        "ranks_per_node does not evenly divide num_ranks");
      TEUCHOS_TEST_FOR_EXCEPTION(
        gpu_ranks_per_node > num_device, std::logic_error,
        "gpu_ranks_per_node cannot exceed number of GPU devices");
      TEUCHOS_TEST_FOR_EXCEPTION(
        gpu_rank && device_id > num_device, std::logic_error,
        "Invalid device ID " << device_id << ".  You probably are trying" <<
        " to run with too many MPI ranks");

      ParameterList node_params;
      if (gpu_rank) {
         std::cout << "MPI rank " << my_rank
                   << ":  Attached to GPU " << device_id
                   << std::endl;

        node_params.set("Verbose", 1);
        node_params.set("Device Number", device_id);

#if TPETRA_USE_KOKKOS_DISTOBJECT && defined(KOKKOS_HAVE_CUDA)
        // Initialize Cuda
        if (!Kokkos::Cuda::is_initialized())
          Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_id) );
#endif

        gpu_norm =
          clone_and_solve<Scalar,LocalOrdinal,GlobalOrdinal,Node,GPUNode>(
            J, f, dx, M, node_params, belosParams, precParams, symmetric);

#if TPETRA_USE_KOKKOS_DISTOBJECT && defined(KOKKOS_HAVE_CUDA)
        Kokkos::Cuda::finalize();
#endif
      }
      else {
        // Note for the non-GPU ranks, we still have to clone since new
        // matrices are created which call fillComplete() (which involves
        // communication)
        gpu_norm =
          clone_and_solve<Scalar,LocalOrdinal,GlobalOrdinal,Node,Node>(
            J, f, dx, M, node_params, belosParams, precParams, symmetric);
      }

      if (my_rank == 0)
        std::cout << "\nNorm of serial node soln - gpu node soln = "
                  << gpu_norm << std::endl;

      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    //Determine if example passed
    bool passed = false;
    if (gpu_norm < Scalar(1e-12) && tpi_norm < Scalar(1e-12))
        passed = true;
    if (my_rank == 0) {
      if (passed)
        std::cout << "Example Passed!" << std::endl;
      else
        std::cout << "Example Failed!" << std::endl;
    }

   }
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (std::string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
  }

}
