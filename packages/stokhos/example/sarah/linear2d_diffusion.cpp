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

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// solver
#include "Ifpack2_Factory.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Timing utilities
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
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  LocalOrdinal n = 16;               // spatial discretization (per dimension)
  LocalOrdinal num_KL = 2;           // number of KL terms
  LocalOrdinal p = 3;                // polynomial order
  BasisScalar mu = 0.2;              // mean of exponential random field
  BasisScalar s = 0.1;               // std. dev. of exponential r.f.
  bool nonlinear_expansion = false;  // nonlinear expansion of diffusion coeff
                                     // (e.g., log-normal)
  bool symmetric = false;            // use symmetric formulation

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  LocalOrdinal MyPID;

  try {

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create a communicator for Epetra objects
    RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< RCP<const Stokhos::OneDOrthogPolyBasis<LocalOrdinal,BasisScalar> > > bases(num_KL); 
    for (LocalOrdinal i=0; i<num_KL; i++)
      bases[i] = rcp(new Stokhos::LegendreBasis<LocalOrdinal,BasisScalar>(p,true));
    RCP<const Stokhos::CompletePolynomialBasis<LocalOrdinal,BasisScalar> > basis = 
      rcp(new Stokhos::CompletePolynomialBasis<LocalOrdinal,BasisScalar>(bases,
		     1e-12));
    LocalOrdinal sz = basis->size();
    RCP<Stokhos::Sparse3Tensor<LocalOrdinal,BasisScalar> > Cijk;
    if (nonlinear_expansion)
      Cijk = basis->computeTripleProductTensor();
    else
      Cijk = basis->computeLinearTripleProductTensor();
    RCP<Stokhos::OrthogPolyExpansion<LocalOrdinal,BasisScalar> > expansion = 
      rcp(new Stokhos::AlgebraicOrthogPolyExpansion<LocalOrdinal,BasisScalar>(basis,
									 Cijk));
    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    LocalOrdinal num_spatial_procs = -1;
    if (argc > 1)
      num_spatial_procs = std::atoi(argv[1]);
    ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    // parallelParams.set("Rebalance Stochastic Graph", true);
    // Teuchos::ParameterList& isorropia_params = 
    //   parallelParams.sublist("Isorropia");
    // isorropia_params.set("Balance objective", "nonzeros");
    RCP<Stokhos::ParallelData> sg_parallel_data =
      rcp(new Stokhos::ParallelData(basis, Cijk, globalComm,
					     parallelParams));
    RCP<const EpetraExt::MultiComm> sg_comm = 
      sg_parallel_data->getMultiComm();
    RCP<const Epetra_Comm> app_comm = 
      sg_parallel_data->getSpatialComm();

    // Create Teuchos::Comm from Epetra_Comm
    RCP< Teuchos::Comm<int> > teuchos_app_comm;
#ifdef HAVE_MPI
    RCP<const Epetra_MpiComm> app_mpi_comm = 
      Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(app_comm);
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > raw_mpi_comm = 
      Teuchos::opaqueWrapper(app_mpi_comm->Comm());
    teuchos_app_comm = rcp(new Teuchos::MpiComm<int>(raw_mpi_comm));
#else
    teuchos_app_comm = rcp(new Teuchos::SerialComm<int>());
#endif

    // Create application
    typedef twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,Node> problem_type;
    RCP<problem_type> model = 
      rcp(new problem_type(teuchos_app_comm, n, num_KL, s, mu, 
			   nonlinear_expansion, symmetric));
    

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
    
    //Convert J to different node types

    //typedefs 
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Belos::LinearProblem<Scalar,MV,OP> BLinProb;

    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;

    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode> Mat_TPINode;    
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode> Mat_GPUNode;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode> GPU_Map;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode> GPU_MV;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::TPINode> TPI_Map;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode> TPI_MV;
    typedef Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::ThrustGPUNode> GPUPrec;
    typedef Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::TPINode> TPIPrec;


    //Convert J to TPI node type
    RCP<ParameterList> plClone = parameterList();
    ParameterList pl;
    pl.set<LocalOrdinal>("Verbose", 1);
    RCP<Kokkos::TPINode> tpi_node = rcp(new Kokkos::TPINode(pl));
    const RCP<const Mat_TPINode> J_TPI = J->clone(tpi_node, plClone);
    J_TPI->print(cout);

    //Convert J to GPU node type
    RCP<Kokkos::ThrustGPUNode> thrustnode = rcp(new Kokkos::ThrustGPUNode(pl));
    const RCP<const Mat_GPUNode> J_GPU = J->clone(thrustnode, plClone);
    J_GPU->print(cout);

    //Clone RILUK preconditioner
    RCP<GPUPrec> M_GPU = factory.clone<Tpetra_CrsMatrix, Mat_GPUNode>(M);
    RCP<TPIPrec> M_TPI = factory.clone<Tpetra_CrsMatrix, Mat_TPINode>(M);

    //Create Chebyshev preconditioner using Prec factory
    Teuchos::ParameterList chevprecParams;
    Teuchos::RCP<Tprec> M_chev;
    M_chev = factory.create<Tpetra_CrsMatrix>("CHEBYSHEV", J);
    M_chev->setParameters(chevprecParams);
    M_chev->initialize();
    M_chev->compute();
    

    //Clone Chebyshev preconditioner
    RCP<GPUPrec>  M_chev_gpu = factory.clone<Tpetra_CrsMatrix, Mat_GPUNode>(M_chev, J_GPU);
    RCP<TPIPrec>  M_chev_tpi = factory.clone<Tpetra_CrsMatrix, Mat_TPINode>(M_chev,  J_TPI);

    // Print nitial residual norm
    Scalar norm_f;
    f->norm2(Teuchos::arrayView(&norm_f,1));
    if (MyPID == 0)
      std::cout << "\nInitial residual norm = " << norm_f << std::endl;

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

    //Convert f, dx to TPI node type
    RCP<TPI_MV> f_tpi = f->clone(tpi_node);
    RCP<TPI_MV> dx_tpi = dx->clone(tpi_node);

    //Convert f, dx to GPU node type
    RCP<GPU_MV> f_gpu = f->clone(thrustnode);
    RCP<GPU_MV> dx_gpu = dx->clone(thrustnode);

    //Create problem for GPU node
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode> GPU_OP;
    typedef Belos::LinearProblem<Scalar,GPU_MV,GPU_OP> GPU_BLinProb;
    RCP <GPU_BLinProb> gpu_problem = rcp(new GPU_BLinProb(J_GPU, dx_gpu, f_gpu));
    gpu_problem->setRightPrec(M_GPU);
    gpu_problem->setProblem(); 
 
    //Create problem for TPI node
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode> TPI_OP;
    typedef Belos::LinearProblem<Scalar,TPI_MV,TPI_OP> TPI_BLinProb;

    RCP <TPI_BLinProb> tpi_problem = rcp(new TPI_BLinProb(J_TPI, dx_tpi, f_tpi));
    tpi_problem->setRightPrec(M_TPI);
    tpi_problem->setProblem(); 

    //Create solver
    RCP<Belos::SolverManager<Scalar,MV,OP> > solver;
    if (symmetric)
      solver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,MV,OP>(problem,
								belosParams));
    else
      solver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV,OP>(problem,
								 belosParams));

    // Create solve for TPI node type
    RCP<Belos::SolverManager<Scalar,TPI_MV,TPI_OP> > tpi_solver;
    if (symmetric)
      tpi_solver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,TPI_MV,TPI_OP>(tpi_problem,
                                                                belosParams));
    else
      tpi_solver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,TPI_MV,TPI_OP>(tpi_problem,
                                                                   belosParams));
								   
    //Create solver for GPU node
    RCP<Belos::SolverManager<Scalar,GPU_MV,GPU_OP> > gpu_solver;
    if (symmetric)
      gpu_solver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,GPU_MV,GPU_OP>(gpu_problem,
                                                                belosParams));
    else
      gpu_solver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,GPU_MV,GPU_OP>(gpu_problem,
                                                                   belosParams));

    // Solve linear system
    std::cout << "Solving with default node..." << std::endl;
    solver->solve();
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

    // Solve linear system for GPU node
    std::cout << "Solving with GPU node..." << std::endl;
    gpu_solver->solve();
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

    // Solve linear system for TPI node
    std::cout << "Solving with TPI node..." << std::endl;
    tpi_solver->solve();
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

    // Update x
    x->update(-1.0, *dx, 1.0);

    // Compute new residual & response function
    RCP<Tpetra_Vector> g = Tpetra::createVector<Scalar>(model->get_g_map(0));
    f->putScalar(0.0);
    model->computeResidual(*x, *p, *f);
    model->computeResponse(*x, *p, *g);

    // Print final residual norm
    f->norm2(Teuchos::arrayView(&norm_f,1));
    if (MyPID == 0)
      std::cout << "\nFinal residual norm = " << norm_f << std::endl;

    

    }


  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
