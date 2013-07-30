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

//Tpetra node conversion
#include <Kokkos_SerialNode.hpp>
#include <Kokkos_TPINode.hpp>
#include <Kokkos_ThrustGPUNode.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  typedef double Scalar;
  typedef double MeshScalar;
  typedef double BasisScalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Belos::LinearProblem<Scalar,MV,OP> BLinProb;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode> Mat_TPINode;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::ThrustGPUNode> Mat_GPUNode;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,KokkosClassic::ThrustGPUNode> GPU_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::ThrustGPUNode> GPU_MV;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode> TPI_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode> TPI_MV;
  typedef Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, KokkosClassic::ThrustGPUNode> GPUPrec;
  typedef Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, KokkosClassic::TPINode> TPIPrec;

  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  // Start up MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    // Create comm
    RCP<const Teuchos::Comm<LocalOrdinal> > comm =
      Teuchos::DefaultComm<int>::getComm();
    int MyPID = comm->getRank();

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
    int device_offset = 0;
    CLP.setOption("device_offset", &device_offset,
                  "Offset for attaching MPI ranks to CUDA devices");
    CLP.parse( argc, argv );
    if (MyPID == 0)
      std::cout << "Summary of command line options:" << std::endl
                << "\tgpu           = " << gpu << std::endl
                << "\ttpi           = " << tpi << std::endl
                << "\tnum_mesh      = " << n << std::endl
                << "\tsymmetric     = " << symmetric << std::endl
                << "\tnum_threads   = " << num_threads << std::endl
                << "\tdevice_offset = " << device_offset << std::endl;

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

    //Create RILUK preconditioner using Prec factory
    ParameterList precParams;
    std::string prec_name = "RILUK";
    precParams.set("fact: iluk level-of-fill", 1);
    precParams.set("fact: iluk level-of-overlap", 0);
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    RCP<Tprec> M;
    Ifpack2::Factory factory;
    M = factory.create<Tpetra_CrsMatrix>(prec_name, J);
    M->setParameters(precParams);
    M->initialize();
    M->compute();

    //Create Chebyshev preconditioner using Prec factory
    Teuchos::ParameterList chevprecParams;
    Teuchos::RCP<Tprec> M_chev;
    M_chev = factory.create<Tpetra_CrsMatrix>("CHEBYSHEV", J);
    M_chev->setParameters(chevprecParams);
    M_chev->initialize();
    M_chev->compute();

    // Setup Belos solver
    RCP<ParameterList> belosParams = rcp(new ParameterList);
    belosParams->set("Num Blocks", 20);
    belosParams->set("Convergence Tolerance", 1e-12);
    belosParams->set("Maximum Iterations", 1000);
    belosParams->set("Verbosity", 33);
    belosParams->set("Output Style", 1);
    belosParams->set("Output Frequency", 10);

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
    std::cout << "Solving with default node..." << std::endl;
    solver->solve();
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

    Scalar tpi_norm = 0.0, gpu_norm = 0.0;

    // Solve linear system with TPI node type
    if (tpi) {

      //Convert J to TPI node type
      RCP<ParameterList> plClone = parameterList();
      ParameterList pl_tpi;
      pl_tpi.set("Verbose", 1);
      pl_tpi.set("Num Threads", num_threads);
      RCP<KokkosClassic::TPINode> tpi_node =
        rcp(new KokkosClassic::TPINode(pl_tpi));
      const RCP<const Mat_TPINode> J_TPI = J->clone(tpi_node, plClone);
      J_TPI->print(std::cout);

      //Clone RILUK preconditioner
      RCP<TPIPrec> M_TPI =
        factory.clone<Tpetra_CrsMatrix, Mat_TPINode>(M, J_TPI);

      //Clone Chebyshev preconditioner
      RCP<TPIPrec>  M_chev_tpi =
        factory.clone<Tpetra_CrsMatrix, Mat_TPINode>(M_chev,  J_TPI);

      //Convert f, dx to TPI node type
      RCP<TPI_MV> f_tpi = f->clone(tpi_node);
      RCP<TPI_MV> dx_tpi = dx->clone(tpi_node);

      //Create problem for TPI node
      typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode> TPI_OP;
      typedef Belos::LinearProblem<Scalar,TPI_MV,TPI_OP> TPI_BLinProb;

      RCP <TPI_BLinProb> tpi_problem =
        rcp(new TPI_BLinProb(J_TPI, dx_tpi, f_tpi));
      tpi_problem->setRightPrec(M_TPI);
      tpi_problem->setProblem();

      // Create solve for TPI node type
      RCP<Belos::SolverManager<Scalar,TPI_MV,TPI_OP> > tpi_solver;
      if (symmetric)
        tpi_solver =
          rcp(new Belos::PseudoBlockCGSolMgr<Scalar,TPI_MV,TPI_OP>(
                tpi_problem, belosParams));
      else
        tpi_solver =
          rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,TPI_MV,TPI_OP>(
                tpi_problem, belosParams));

      // Solve linear system for TPI node
      std::cout << "Solving with TPI node..." << std::endl;
      tpi_solver->solve();
      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();

      ParameterList pl_serial;
      RCP<KokkosClassic::SerialNode> serialnode =
        rcp(new KokkosClassic::SerialNode(pl_serial));
      RCP<MV> dx_tpioncpu = dx_tpi->clone(serialnode);
      dx_tpioncpu->update(1.0, *dx, -1.0);
      dx_tpioncpu->norm2(Teuchos::arrayView(&tpi_norm,1));
      if (MyPID == 0)
        std::cout << "\nNorm of serial node soln - tpi node soln = " << tpi_norm << std::endl;
    }

    // Solve linear system with GPU node type
    if (gpu) {

      // Compute CUDA device ID
      int num_device; cudaGetDeviceCount(&num_device);
      int device_id = MyPID % num_device + device_offset;
      TEUCHOS_TEST_FOR_EXCEPTION(
        device_id > num_device, std::logic_error,
        "Invalid device ID " << device_id << ".  You probably are trying" <<
        " to run with too many MPI ranks");

      //Convert J to GPU node type
      RCP<ParameterList> plClone = parameterList();
      ParameterList pl_gpu;
      pl_gpu.set("Verbose", 1);
      pl_gpu.set("Device Number", device_id);
      RCP<KokkosClassic::ThrustGPUNode> thrustnode =
        rcp(new KokkosClassic::ThrustGPUNode(pl_gpu));
      const RCP<const Mat_GPUNode> J_GPU = J->clone(thrustnode, plClone);
      J_GPU->print(std::cout);

      //Clone RILUK preconditioner
      RCP<GPUPrec> M_GPU =
        factory.clone<Tpetra_CrsMatrix, Mat_GPUNode>(M, J_GPU);

      //Clone Chebyshev preconditioner
      RCP<GPUPrec>  M_chev_gpu =
        factory.clone<Tpetra_CrsMatrix, Mat_GPUNode>(M_chev, J_GPU);

      //Convert f, dx to GPU node type
      RCP<GPU_MV> f_gpu = f->clone(thrustnode);
      RCP<GPU_MV> dx_gpu = dx->clone(thrustnode);

      //Create problem for GPU node
      typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::ThrustGPUNode> GPU_OP;
      typedef Belos::LinearProblem<Scalar,GPU_MV,GPU_OP> GPU_BLinProb;
      RCP <GPU_BLinProb> gpu_problem =
        rcp(new GPU_BLinProb(J_GPU, dx_gpu, f_gpu));
      gpu_problem->setRightPrec(M_GPU);
      gpu_problem->setProblem();

      //Create solver for GPU node
      RCP<Belos::SolverManager<Scalar,GPU_MV,GPU_OP> > gpu_solver;
      if (symmetric)
        gpu_solver =
          rcp(new Belos::PseudoBlockCGSolMgr<Scalar,GPU_MV,GPU_OP>(
                gpu_problem, belosParams));
      else
        gpu_solver =
          rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,GPU_MV,GPU_OP>(
                gpu_problem, belosParams));

      // Solve linear system for GPU node
      std::cout << "Solving with GPU node..." << std::endl;
      gpu_solver->solve();
      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();

      ParameterList pl_serial;
      RCP<KokkosClassic::SerialNode> serialnode =
        rcp(new KokkosClassic::SerialNode(pl_serial));
      RCP<MV> dx_gpuoncpu = dx_gpu->clone(serialnode);
      dx_gpuoncpu->update(1.0, *dx, -1.0);
      dx_gpuoncpu->norm2(Teuchos::arrayView(&gpu_norm,1));
      if (MyPID == 0)
        std::cout << "\nNorm of serial node soln - gpu node soln = " << gpu_norm << std::endl;
    }

    //Determine if example passed
    bool passed = false;
    if (gpu_norm < Scalar(1e-12) && tpi_norm < Scalar(1e-12))
        passed = true;
    if (MyPID == 0) {
      if (passed)
        std::cout << "Example Passed!" << std::endl;
      else
        std::cout << "Example Failed!" << std::endl;
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
