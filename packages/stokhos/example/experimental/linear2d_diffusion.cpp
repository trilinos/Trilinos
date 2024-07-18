// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  typedef double Scalar;
  typedef double MeshScalar;
  typedef double BasisScalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef Tpetra::Map<>::node_type Node;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  LocalOrdinal n = 32;               // spatial discretization (per dimension)
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

    // Create preconditioner
    Teuchos::ParameterList precParams;
    std::string prec_name = "RILUK";
    precParams.set("fact: iluk level-of-fill", 1);
    precParams.set("fact: iluk level-of-overlap", 0);
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    Teuchos::RCP<Tprec> M;
    Ifpack2::Factory factory;
    M = factory.create<Tpetra_CrsMatrix>(prec_name, J);
    M->setParameters(precParams);

    // Evaluate model
    model->computeResidual(*x, *p, *f);
    model->computeJacobian(*x, *p, *J);
    M->initialize();
    M->compute();

    // Print initial residual norm
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
    belosParams->set("Output Frequency", 1);
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Belos::OperatorTraits<Scalar,MV,OP> BOPT;
    typedef Belos::MultiVecTraits<Scalar,MV> BMVT;
    typedef Belos::LinearProblem<Scalar,MV,OP> BLinProb;
    RCP< BLinProb > problem = rcp(new BLinProb(J, dx, f));
    problem->setRightPrec(M);
    problem->setProblem();
    RCP<Belos::SolverManager<Scalar,MV,OP> > solver;
    if (symmetric)
      solver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,MV,OP>(problem,
								belosParams));
    else
      solver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV,OP>(problem,
								   belosParams));

    // Solve linear system
    Belos::ReturnType ret = solver->solve();

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

    // Print response
    std::cout << "\nResponse =      " << std::endl;
    Writer::writeDense(std::cout, g);

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

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
