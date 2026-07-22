// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "Stokhos_Sacado.hpp"
#include "Stokhos_Ifpack2.hpp"

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
#include "kokkos_pce_specializations.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Utilities
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Scalar types
#include "linear2d_diffusion_scalar_types.hpp"

// MueLu includes
#include "Stokhos_MueLu.hpp"
#include "Stokhos_MueLu_QR_Interface_decl.hpp"
#include "Stokhos_MueLu_QR_Interface_def.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;
//#include <MueLu_UseShortNames.hpp>

#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_Map.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_SaPFactory.hpp"

// Random field types
enum SG_RF { UNIFORM, LOGNORMAL };
const int num_sg_rf = 2;
const SG_RF sg_rf_values[] = { UNIFORM, LOGNORMAL };
const char *sg_rf_names[] = { "Uniform", "Log-Normal" };

// Krylov methods
enum Krylov_Method { GMRES, CG };
const int num_krylov_method = 2;
const Krylov_Method krylov_method_values[] = { GMRES, CG };
const char *krylov_method_names[] = { "GMRES", "CG" };

// Multigrid preconditioning
enum Multigrid_Smoother { CHEBYSHEV, SGS };
const int num_multigrid_smoother = 2;
const Multigrid_Smoother multigrid_smoother_values[] = { CHEBYSHEV, SGS };
const char *multigrid_smoother_names[] = { "Chebyshev", "SGS" };

// Preconditioning approaches
enum SG_Prec { NONE, MEAN, STOCHASTIC };
const int num_sg_prec = 3;
const SG_Prec sg_prec_values[] = { NONE, MEAN, STOCHASTIC };
const char *sg_prec_names[] = { "None",
                                "Mean-Based", 
                                "Stochastic" };

// Stochastic division approaches
enum SG_Div { DIRECT, SPD_DIRECT, MEAN_DIV, QUAD, CGD };
const int num_sg_div = 5;
const SG_Div sg_div_values[] = { DIRECT, SPD_DIRECT, MEAN_DIV, QUAD, CGD };
const char *sg_div_names[] = { "Direct",
                               "SPD-Direct",
                               "Mean-Based", 
                               "Quadrature",
                               "CG"};

// Stochastic division preconditioner approaches
enum SG_DivPrec { NO, DIAG, JACOBI, GS, SCHUR };
const int num_sg_divprec = 5;
const SG_DivPrec sg_divprec_values[] = {NO, DIAG, JACOBI, GS, SCHUR};
const char *sg_divprec_names[] = { "None",
                                   "Diag",
                                   "Jacobi",
                                   "GS", 
                                   "Schur"};


// Option for Schur complement precond: full or diag D
enum Schur_option { full, diag };
const int num_schur_option = 2;
const Schur_option Schur_option_values[] = { full, diag };
const char *schur_option_names[] = { "full", "diag"};

// Full matrix or linear matrix (pb = dim + 1 ) used for preconditioner
enum Prec_option { whole, linear};
const int num_prec_option = 2;
const Prec_option Prec_option_values[] = { whole, linear };
const char *prec_option_names[] = { "full", "linear"};

// #define _GNU_SOURCE 1
// #include <fenv.h>

template<typename ordinal_type, typename value_type, typename Storage>
//void returnScalarAsDenseMatrix(Scalar const &inval,
void returnScalarAsDenseMatrix(Sacado::PCE::OrthogPoly<value_type,Storage> const &inval,
                               Teuchos::RCP<Teuchos::SerialDenseMatrix<ordinal_type,value_type> > & denseEntry,
                               Teuchos::RCP<Stokhos::Sparse3Tensor<ordinal_type,value_type> > const &Cijk)
{
    Stokhos::OrthogPolyApprox<ordinal_type, value_type> val= inval.getOrthogPolyApprox();
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;
    ordinal_type pb = val.size();
    const value_type* cv = val.coeff();

    denseEntry->putScalar(0.0);
    typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
    typename Cijk_type::k_iterator k_end = Cijk->k_end();
    if (pb < Cijk->num_k())
      k_end = Cijk->find_k(pb);
    value_type cijk;
    ordinal_type i,j,k;
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      k = index(k_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); j_it != Cijk->j_end(k_it); ++j_it) {
         j = index(j_it);
         for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it); i_it != Cijk->i_end(j_it); ++i_it) {
           i = index(i_it);
           cijk = value(i_it);
           (*denseEntry)(i,j) += cijk*cv[k];
         }
      }
    }
}

typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Xpetra_Matrix;
typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Xpetra_Map;

template<typename ordinal_type,typename value_type>
void PrintMatrix(Teuchos::FancyOStream &fos, Teuchos::RCP<Xpetra_Matrix> const &A,
                Teuchos::RCP<Stokhos::Sparse3Tensor<ordinal_type, value_type> > const & Cijk,
                Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > const & basis)
{
  ordinal_type sz = basis->size();
  Teuchos::RCP<Teuchos::SerialDenseMatrix<ordinal_type,value_type> > denseEntry = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
             sz, sz));
    size_t maxLength = A->getLocalMaxNumRowEntries();
    size_t NumEntries;
    Scalar val;
    Teuchos::Array<ordinal_type> Indices(maxLength);
    Teuchos::Array<Scalar> Values(maxLength);
    Teuchos::RCP<const Xpetra_Map> colMap = A->getColMap();
    for (ordinal_type i = 0 ; i < Teuchos::as<ordinal_type>(A->getLocalNumRows()); ++i) {
      A->getLocalRowCopy(i, Indices(), Values(), NumEntries);
      fos << "++++++++++++++" << std::endl << "row " << A->getRowMap()->getGlobalElement(i) << ": ";
      fos << "  col ids: ";
      for (size_t ii=0; ii<NumEntries; ++ii) fos << colMap->getGlobalElement(Indices[ii]) << " ";
      fos << std::endl << "++++++++++++++" << std::endl;
      for (size_t k=0; k< NumEntries; ++k) {
        val = Values[k];
        Teuchos::OSTab tab1(fos);
        fos << std::endl << "col " << colMap->getGlobalElement(Indices[k]) << std::endl;
        returnScalarAsDenseMatrix(val,denseEntry,Cijk);
        //TODO tab thing
        Teuchos::OSTab tab2(fos);
        denseEntry->print(fos);
      }
    }
}

int main(int argc, char *argv[]) {
  typedef double MeshScalar;
  typedef double BasisScalar;
  typedef Tpetra::Map<>::node_type Node;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::ParameterList;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  LocalOrdinal MyPID;

  try {

    // Create a communicator for Epetra objects
    RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs an interlaced stochastic Galerkin solvers.\n");

    int n = 32;
    CLP.setOption("num_mesh", &n, "Number of mesh points in each direction");

    // multigrid specific options
    int minAggSize = 1;
    CLP.setOption("min_agg_size", &minAggSize, "multigrid aggregate size");
    Multigrid_Smoother smoother_type = CHEBYSHEV;
    CLP.setOption("smoother_type", &smoother_type, 
                  num_multigrid_smoother, multigrid_smoother_values, multigrid_smoother_names, 
                  "Multigrid smoother types");
    int smootherSweeps = 3;
    CLP.setOption("smoother_sweeps", &smootherSweeps, "# multigrid smoother sweeps");
    bool plainAgg=false;
    CLP.setOption("plain_aggregation", "smoothed_aggregation", &plainAgg, "multigrid prolongator smoothing");
    LocalOrdinal nsSize=1;
    CLP.setOption("nullspace_size", &nsSize, "multigrid nullspace dimension");
    int maxAMGLevels=10;
    CLP.setOption("max_multigrid_levels", &maxAMGLevels, "maximum # levels in multigrid preconditioner");

    bool printTimings=true;
    CLP.setOption("timings", "notimings", &printTimings, "print timing summary");


    bool symmetric = false;
    CLP.setOption("symmetric", "unsymmetric", &symmetric, "Symmetric discretization");

    int num_spatial_procs = -1;
    CLP.setOption("num_spatial_procs", &num_spatial_procs, "Number of spatial processors (set -1 for all available procs)");

    SG_RF randField = UNIFORM;
    CLP.setOption("rand_field", &randField, 
                  num_sg_rf, sg_rf_values, sg_rf_names,
                  "Random field type");

    double mu = 0.2;
    CLP.setOption("mean", &mu, "Mean");

    double s = 0.1;
    CLP.setOption("std_dev", &s, "Standard deviation");

    int num_KL = 2;
    CLP.setOption("num_kl", &num_KL, "Number of KL terms");

    int order = 3;
    CLP.setOption("order", &order, "Polynomial order");

    bool normalize_basis = true;
    CLP.setOption("normalize", "unnormalize", &normalize_basis, 
                  "Normalize PC basis");

    Krylov_Method solver_method = GMRES;
    CLP.setOption("solver_method", &solver_method, 
                  num_krylov_method, krylov_method_values, krylov_method_names, 
                  "Krylov solver method");

    SG_Prec prec_method = STOCHASTIC;
    CLP.setOption("prec_method", &prec_method, 
                  num_sg_prec, sg_prec_values, sg_prec_names,
                  "Preconditioner method");

    SG_Div division_method = DIRECT;
    CLP.setOption("division_method", &division_method, 
                  num_sg_div, sg_div_values, sg_div_names,
                  "Stochastic division method");

    SG_DivPrec divprec_method = NO;
    CLP.setOption("divprec_method", &divprec_method,
                  num_sg_divprec, sg_divprec_values, sg_divprec_names,
                  "Preconditioner for division method");
    Schur_option schur_option = diag;
    CLP.setOption("schur_option", &schur_option,
                  num_schur_option, Schur_option_values, schur_option_names,
                  "Schur option");
    Prec_option prec_option = whole;
    CLP.setOption("prec_option", &prec_option,
                  num_prec_option, Prec_option_values, prec_option_names,
                  "Prec option");


    double solver_tol = 1e-12;
    CLP.setOption("solver_tol", &solver_tol, "Outer solver tolerance");

    double div_tol = 1e-6;
    CLP.setOption("div_tol", &div_tol, "Tolerance in Iterative Solver");
    
    int prec_level = 1;
    CLP.setOption("prec_level", &prec_level, "Level in Schur Complement Prec 0->Solve A0u0=g0 with division; 1->Form 1x1 Schur Complement");

    int max_it_div = 50;
    CLP.setOption("max_it_div", &max_it_div, "Maximum # of Iterations in Iterative Solver for Division");

    bool equilibrate = true; //JJH 8/26/12 changing to true to match ETP example
    CLP.setOption("equilibrate", "noequilibrate", &equilibrate,
                  "Equilibrate the linear system");

    bool printHierarchy = false;
    CLP.setOption("print_hierarchy", "noprint_Hierarchy", &printHierarchy,
                  "Print matrices in multigrid hierarchy");


    CLP.parse( argc, argv );

    if (MyPID == 0) {
      std::cout << "Summary of command line options:" << std::endl
                << "\tnum_mesh           = " << n << std::endl
                << "\tsymmetric          = " << symmetric << std::endl
                << "\tnum_spatial_procs  = " << num_spatial_procs << std::endl
                << "\trand_field         = " << sg_rf_names[randField] 
                << std::endl
                << "\tmean               = " << mu << std::endl
                << "\tstd_dev            = " << s << std::endl
                << "\tnum_kl             = " << num_KL << std::endl
                << "\torder              = " << order << std::endl
                << "\tnormalize_basis    = " << normalize_basis << std::endl
                << "\tsolver_method      = " << krylov_method_names[solver_method] << std::endl
                << "\tprec_method        = " << sg_prec_names[prec_method]    << std::endl
                << "\tdivision_method    = " << sg_div_names[division_method]     << std::endl
                << "\tdiv_tol            = " << div_tol << std::endl
                << "\tdiv_prec           = " << sg_divprec_names[divprec_method]      << std::endl
                << "\tprec_level         = " << prec_level << std::endl
                << "\tmax_it_div         = " << max_it_div << std::endl
                << "\t~~~ multigrid options ~~~" << std::endl
                << "\tsmoother_type      = " << multigrid_smoother_names[smoother_type] << std::endl
                << "\tsmoother_sweeps    = " << smootherSweeps << std::endl
                << "\tplain_aggregation  = " << plainAgg << std::endl
                << "\tmax_multigrid_levels = " << maxAMGLevels << std::endl
                << "\tnullspace_size     = " << nsSize << std::endl
                << "\tmin_agg_size       = " << minAggSize << std::endl;
    }
    bool nonlinear_expansion = false;
    if (randField == UNIFORM)
      nonlinear_expansion = false;
    else if (randField == LOGNORMAL)
      nonlinear_expansion = true;

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< RCP<const Stokhos::OneDOrthogPolyBasis<LocalOrdinal,BasisScalar> > > bases(num_KL); 
    for (LocalOrdinal i=0; i<num_KL; i++)
      if (randField == UNIFORM)
        bases[i] = rcp(new Stokhos::LegendreBasis<LocalOrdinal,BasisScalar>(order, normalize_basis));
      else if (randField == LOGNORMAL)
        bases[i] = rcp(new Stokhos::HermiteBasis<int,double>(order, normalize_basis));
    RCP<const Stokhos::CompletePolynomialBasis<LocalOrdinal,BasisScalar> > basis = 
      rcp(new Stokhos::CompletePolynomialBasis<LocalOrdinal,BasisScalar>(bases, 1e-12));
    LocalOrdinal sz = basis->size();
    RCP<Stokhos::Sparse3Tensor<LocalOrdinal,BasisScalar> > Cijk = 
      basis->computeTripleProductTensor();
    RCP<const Stokhos::Quadrature<int,double> > quad = 
      rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    RCP<ParameterList> expn_params = Teuchos::rcp(new ParameterList);
    if (division_method == MEAN_DIV) {
      expn_params->set("Division Strategy", "Mean-Based");
      expn_params->set("Use Quadrature for Division", false);
    }
    else if (division_method == DIRECT) {
      expn_params->set("Division Strategy", "Dense Direct");
      expn_params->set("Use Quadrature for Division", false);
    }
    else if (division_method == SPD_DIRECT) {
      expn_params->set("Division Strategy", "SPD Dense Direct");
      expn_params->set("Use Quadrature for Division", false);
    }
    else if (division_method == CGD) {
      expn_params->set("Division Strategy", "CG");
      expn_params->set("Use Quadrature for Division", false);
    }

    else if (division_method == QUAD) {
      expn_params->set("Use Quadrature for Division", true);
    }

    if (divprec_method == NO)
         expn_params->set("Prec Strategy", "None");
    else if (divprec_method == DIAG)
         expn_params->set("Prec Strategy", "Diag");
    else if (divprec_method == JACOBI)
         expn_params->set("Prec Strategy", "Jacobi");
    else if (divprec_method == GS)
         expn_params->set("Prec Strategy", "GS");
    else if (divprec_method == SCHUR)
         expn_params->set("Prec Strategy", "Schur");

    if (schur_option == diag)
        expn_params->set("Schur option", "diag");
    else
        expn_params->set("Schur option", "full");
    if (prec_option == linear)
        expn_params->set("Prec option", "linear");


    if (equilibrate)
      expn_params->set("Equilibrate", 1);
    else
      expn_params->set("Equilibrate", 0); 
    expn_params->set("Division Tolerance", div_tol);
    expn_params->set("prec_iter", prec_level);
    expn_params->set("max_it_div", max_it_div);

    RCP<Stokhos::OrthogPolyExpansion<LocalOrdinal,BasisScalar> > expansion = 
      rcp(new Stokhos::QuadOrthogPolyExpansion<LocalOrdinal,BasisScalar>(
            basis, Cijk, quad, expn_params));

    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    // parallelParams.set("Rebalance Stochastic Graph", true);
    // Teuchos::ParameterList& isorropia_params = 
    //   parallelParams.sublist("Isorropia");
    // isorropia_params.set("Balance objective", "nonzeros");
    RCP<Stokhos::ParallelData> sg_parallel_data =
      rcp(new Stokhos::ParallelData(basis, Cijk, globalComm, parallelParams));
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
    //Xpetra matrices
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Xpetra_CrsMatrix;
    typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Xpetra_MultiVector;
    typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Xpetra_MultiVectorFactory;
    typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Xpetra_TpetraCrsMatrix;
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> Xpetra_CrsMatrixWrap;
    typedef Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> Belos_MueLuOperator;
    //MueLu typedefs
    typedef MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> MueLu_Hierarchy;
    typedef MueLu::SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherPrototype;
    typedef MueLu::TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> TrilinosSmoother;
    typedef MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherFactory;
    typedef MueLu::FactoryManager<Scalar,LocalOrdinal,GlobalOrdinal,Node> FactoryManager;
  
    //feenableexcept(FE_ALL_EXCEPT);
    //feenableexcept(FE_INVALID | FE_DIVBYZERO);

    RCP<Tpetra_Vector> p = Tpetra::createVector<Scalar>(model->get_p_map(0));
    RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(model->get_x_map());
    x->putScalar(0.0);
    RCP<Tpetra_Vector> f = Tpetra::createVector<Scalar>(model->get_f_map());
    RCP<Tpetra_Vector> dx = Tpetra::createVector<Scalar>(model->get_x_map());
    RCP<Tpetra_CrsMatrix> J = model->create_W();
    RCP<Tpetra_CrsMatrix> J0;
    if (prec_method == MEAN)
      J0 = model->create_W();

    // Set PCE expansion of p
    p->putScalar(0.0);
    ArrayRCP<Scalar> p_view = p->get1dViewNonConst();
    for (ArrayRCP<Scalar>::size_type i=0; i<p_view.size(); i++) {
      p_view[i].reset(expansion);
      p_view[i].copyForWrite();
    }
    Array<double> point(num_KL, 1.0);
    Array<double> basis_vals(sz);
    basis->evaluateBases(point, basis_vals);
    if (order > 0) {
      for (int i=0; i<num_KL; i++) {
        p_view[i].term(i,1) = 1.0 / basis_vals[i+1];
      }
    }

    // Create preconditioner
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    RCP<Belos_MueLuOperator> M;
    RCP<MueLu_Hierarchy> H;
    RCP<Xpetra_CrsMatrix> xcrsJ = rcp(new Xpetra_TpetraCrsMatrix(J));
    RCP<Xpetra_Matrix> xopJ = rcp(new Xpetra_CrsMatrixWrap(xcrsJ));
    RCP<Xpetra_Matrix> xopJ0;
    if (prec_method != NONE) {
      ParameterList precParams;
      std::string prec_name = "RILUK";
      precParams.set("fact: iluk level-of-fill", 1);
      precParams.set("fact: iluk level-of-overlap", 0);
      //Ifpack2::Factory factory;
      if (prec_method == MEAN) {
        RCP<Xpetra_CrsMatrix> xcrsJ0 = rcp(new Xpetra_TpetraCrsMatrix(J0));
        xopJ0 = rcp(new Xpetra_CrsMatrixWrap(xcrsJ0));
        //M = factory.create<Tpetra_CrsMatrix>(prec_name, J0);
      } else if (prec_method == STOCHASTIC) {
        xopJ0 = xopJ;
        //M = factory.create<Tpetra_CrsMatrix>(prec_name, J);
      }
      H = rcp(new MueLu_Hierarchy(xopJ0));
      M = rcp(new Belos_MueLuOperator(H));
      //M->setParameters(precParams);
      if (nsSize==-1) nsSize=sz;
      RCP<Xpetra_MultiVector> Z = Xpetra_MultiVectorFactory::Build(xcrsJ->getDomainMap(), nsSize);
      size_t n = Z->getLocalLength();
      for (LocalOrdinal j=0; j<nsSize; ++j) {
        ArrayRCP<Scalar> col = Z->getDataNonConst(j);
        for (size_t i=0; i<n; ++i) {
          col[i].reset(expansion);
          col[i].copyForWrite();
          col[i].fastAccessCoeff(j) = 1.0;
        }
      }
      H->GetLevel(0)->Set("Nullspace", Z);
      //RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      //fos->setOutputToRootOnly(-1);
      //Z->describe(*fos);
    }

    // Evaluate model
    model->computeResidual(*x, *p, *f);
    model->computeJacobian(*x, *p, *J);


    // Compute mean for mean-based preconditioner
    if (prec_method == MEAN) {
      size_t nrows = J->getLocalNumRows();
      ArrayView<const LocalOrdinal> indices;
      ArrayView<const Scalar> values;
      J0->resumeFill();
      for (size_t i=0; i<nrows; i++) {
        J->getLocalRowView(i, indices, values);
        Array<Scalar> values0(values.size());
        for (LocalOrdinal j=0; j<values.size(); j++)
          values0[j] = values[j].coeff(0);
        J0->replaceLocalValues(i, indices, values0);
      }
      J0->fillComplete();
    }


    // compute preconditioner
    if (prec_method != NONE) {
      //M->initialize();
      //M->compute();

      //override MueLu defaults via factory manager
      RCP<FactoryManager> fm = rcp( new FactoryManager() );;

      //smoother
      ParameterList smootherParamList;
      RCP<SmootherPrototype> smooPrototype;
      switch(smoother_type) {
        case CHEBYSHEV: {
          smootherParamList.set("chebyshev: degree", smootherSweeps);
          smootherParamList.set("chebyshev: ratio eigenvalue", (double) 20);
          Scalar minusOne=-1.0;
          smootherParamList.set("chebyshev: max eigenvalue", minusOne);
          smootherParamList.set("chebyshev: min eigenvalue", (double) 1.0);
          smootherParamList.set("chebyshev: zero starting solution", true);
          smooPrototype     = rcp( new TrilinosSmoother("CHEBYSHEV", smootherParamList) );
          break;
          }

        case SGS:
        default:
          smootherParamList.set("relaxation: sweeps", smootherSweeps);
          smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
          smooPrototype     = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
          break;
      }

      RCP<SmootherFactory>   smooFact      = rcp( new SmootherFactory(smooPrototype) );
      fm->SetFactory("Smoother", smooFact);

      // coarse level solve
      // TODO until KLU in Amesos2 is fully templated, use incomplete factorization as coarsest level solve
      ParameterList coarseParamList;
      coarseParamList.set("fact: level-of-fill", 0);
      RCP<SmootherPrototype> coarsePrototype     = rcp( new TrilinosSmoother("ILUT", coarseParamList) );
      //RCP<SmootherPrototype> coarsePrototype     = rcp( new TrilinosSmoother("RILUK", coarseParamList) );
      RCP<SmootherFactory>   coarseSolverFact      = rcp( new SmootherFactory(coarsePrototype, Teuchos::null) );
      fm->SetFactory("CoarseSolver", coarseSolverFact);
      //fm->SetFactory("CoarseSolver", smooFact);

      //allow for larger aggregates
      typedef MueLu::CoupledAggregationFactory<LocalOrdinal,GlobalOrdinal,Node>
      MueLu_CoupledAggregationFactory;
      RCP<MueLu_CoupledAggregationFactory> aggFact = rcp(new MueLu_CoupledAggregationFactory());
      aggFact->SetMinNodesPerAggregate(minAggSize);
      fm->SetFactory("Aggregates", aggFact);

      //turn off damping
      typedef MueLu::SaPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> MueLu_SaPFactory;
      if (plainAgg) {
        RCP<MueLu_SaPFactory> sapFactory = rcp(new MueLu_SaPFactory);
        sapFactory->SetDampingFactor( (Scalar) 0.0 );
        fm->SetFactory("P", sapFactory);
      }

      H->Setup(*fm,0,maxAMGLevels); //start at level 0, at most 3 levels

      if (printHierarchy)
      {
        //FIXME #levels
        RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        int numLevels = H->GetNumLevels();
        for (int i=0; i<numLevels; ++i) {
          RCP<Xpetra_Matrix> A = H->GetLevel(i)->Get<RCP<Xpetra_Matrix> >("A");
          *fos << "================\n" << "matrix A, multigrid level " << i << "\n================" << std::endl;
          PrintMatrix<LocalOrdinal,BasisScalar>(*fos,A,Cijk,basis);
        }
      }
    }

    // Setup Belos solver
    RCP<ParameterList> belosParams = rcp(new ParameterList);
    belosParams->set("Flexible Gmres", false);
    belosParams->set("Num Blocks", 500);//20
    belosParams->set("Convergence Tolerance", solver_tol);
    belosParams->set("Maximum Iterations", 1000);
    belosParams->set("Verbosity", 33);
    belosParams->set("Output Style", 1);
    belosParams->set("Output Frequency", 1);
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Belos::OperatorT<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > OP;
    typedef Belos::OperatorTraits<Scalar,MV,OP> BOPT;
    typedef Belos::MultiVecTraits<Scalar,MV> BMVT;
    typedef Belos::MultiVecTraits<double,MV> BTMVT;
    typedef Belos::LinearProblem<double,MV,OP> BLinProb;
    typedef Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> BXpetraOp;
    RCP<OP> belosJ = rcp(new BXpetraOp(xopJ)); // Turns an Xpetra::Matrix object into a Belos operator
    RCP< BLinProb > problem = rcp(new BLinProb(belosJ, dx, f));
    if (prec_method != NONE)
      problem->setRightPrec(M);
    problem->setProblem();
    RCP<Belos::SolverManager<double,MV,OP> > solver;
    if (solver_method == CG)
      solver = rcp(new Belos::PseudoBlockCGSolMgr<double,MV,OP>(problem, belosParams));
    else if (solver_method == GMRES)
      solver = rcp(new Belos::BlockGmresSolMgr<double,MV,OP>(problem, belosParams));
    
    // Print initial residual norm
    std::vector<double> norm_f(1);
    //BMVT::MvNorm(*f, norm_f);
    BTMVT::MvNorm(*f, norm_f);
    if (MyPID == 0)
      std::cout << "\nInitial residual norm = " << norm_f[0] << std::endl;

    // Solve linear system
    Belos::ReturnType ret = solver->solve();

    if (MyPID == 0) {
      if (ret == Belos::Converged)
        std::cout << "Solver converged!" << std::endl;
      else
        std::cout << "Solver failed to converge!" << std::endl;
    }

    // Update x
    x->update(-1.0, *dx, 1.0);
    Writer::writeDenseFile("stochastic_solution.mm", x);

    // Compute new residual & response function
    RCP<Tpetra_Vector> g = Tpetra::createVector<Scalar>(model->get_g_map(0));
    f->putScalar(0.0);
    model->computeResidual(*x, *p, *f);
    model->computeResponse(*x, *p, *g);

    // Print final residual norm
    //BMVT::MvNorm(*f, norm_f);
    BTMVT::MvNorm(*f, norm_f);
    if (MyPID == 0)
      std::cout << "\nFinal residual norm = " << norm_f[0] << std::endl;

    // Expected results for num_mesh=32
    double g_mean_exp = 0.172988;      // expected response mean
    double g_std_dev_exp = 0.0380007;  // expected response std. dev.
    double g_tol = 1e-6;               // tolerance on determining success
    if (n == 8) {
      g_mean_exp = 1.327563e-01;
      g_std_dev_exp = 2.949064e-02;
    }

    double g_mean = g->get1dView()[0].mean();
    double g_std_dev = g->get1dView()[0].standard_deviation();
    std::cout << std::endl;
    std::cout << "Response Mean =      " << g_mean << std::endl;
    std::cout << "Response Std. Dev. = " << g_std_dev << std::endl;
    bool passed = false;
    if (norm_f[0] < 1.0e-10 &&
        std::abs(g_mean-g_mean_exp) < g_tol &&
        std::abs(g_std_dev - g_std_dev_exp) < g_tol)
      passed = true;
    if (MyPID == 0) {
      if (passed)
        std::cout << "Example Passed!" << std::endl;
      else{
        std::cout << "Example Failed!" << std::endl;
        std::cout << "Expected Response Mean      = "<< g_mean_exp << std::endl;
        std::cout << "Expected Response Std. Dev. = "<< g_std_dev_exp << std::endl;
      }
    }

    }

    if (printTimings) {
      Teuchos::TimeMonitor::summarize(std::cout);
      Teuchos::TimeMonitor::zeroOutTimers();
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
