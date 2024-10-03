// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ModelEvaluator implementing our problem
#include "twoD_diffusion_ME.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// AztecOO solver
#include "AztecOO.h"

// Stokhos Stochastic Galerkin
#include "Stokhos.hpp"

// Timing utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Block utilities
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_RowMatrixOut.h"

// Krylov methods
enum Krylov_Method { GMRES, CG };
const int num_krylov_method = 2;
const Krylov_Method krylov_method_values[] = { GMRES, CG };
const char *krylov_method_names[] = { "GMRES", "CG" };

// Random field types
enum SG_RF { UNIFORM, CC_UNIFORM, RYS, LOGNORMAL };
const int num_sg_rf = 4;
const SG_RF sg_rf_values[] = { UNIFORM, CC_UNIFORM, RYS, LOGNORMAL };
const char *sg_rf_names[] = { "Uniform", "CC-Uniform", "Rys", "Log-Normal" };

// Quadrature types
enum SG_Quad { TENSOR, SPARSE_GRID };
const int num_sg_quad = 2;
const SG_Quad sg_quad_values[] = { TENSOR, SPARSE_GRID };
const char *sg_quad_names[] = { "tensor", "sparse-grid" };

// Sparse grid growth rules
enum SG_GROWTH { SLOW_RESTRICTED, MODERATE_RESTRICTED, UNRESTRICTED };
const int num_sg_growth = 3;
const SG_GROWTH sg_growth_values[] = {
  SLOW_RESTRICTED, MODERATE_RESTRICTED, UNRESTRICTED };
const char *sg_growth_names[] = {
  "Slow Restricted", "Moderate Restricted", "Unrestricted" };

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Array;
  using Teuchos::CommandLineProcessor;
  using Teuchos::TimeMonitor;
  using Teuchos::ParameterList;
  using EpetraExt::BlockUtility;
  using EpetraExt::ModelEvaluator;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> Comm;
#ifdef HAVE_MPI
  Comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  Comm = rcp(new Epetra_SerialComm);
#endif

  int MyPID = Comm->MyPID();

  try {

    // Setup command line options
    CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs a stochastic collocation method.\n");

    int n = 32;
    CLP.setOption("num_mesh", &n, "Number of mesh points in each direction");

    SG_RF randField = UNIFORM;
    CLP.setOption("rand_field", &randField,
                   num_sg_rf, sg_rf_values, sg_rf_names,
                  "Random field type");

    double mean = 0.2;
    CLP.setOption("mean", &mean, "Mean");

    double sigma = 0.1;
    CLP.setOption("std_dev", &sigma, "Standard deviation");

    double weightCut = 1.0;
    CLP.setOption("weight_cut", &weightCut, "Weight cut");

    int num_KL = 2;
    CLP.setOption("num_kl", &num_KL, "Number of KL terms");

    int p = 3;
    CLP.setOption("order", &p, "Polynomial order");

    bool normalize_basis = true;
    CLP.setOption("normalize", "unnormalize", &normalize_basis,
                  "Normalize PC basis");

    Krylov_Method krylov_method = CG;
    CLP.setOption("krylov_method", &krylov_method,
                  num_krylov_method, krylov_method_values, krylov_method_names,
                  "Krylov method");

    double tol = 1e-12;
    CLP.setOption("tol", &tol, "Solver tolerance");

#ifdef HAVE_STOKHOS_DAKOTA
    SG_Quad quad_method = SPARSE_GRID;
#else
    SG_Quad quad_method = TENSOR;
#endif
    CLP.setOption("quadrature", &quad_method,
                   num_sg_quad, sg_quad_values, sg_quad_names,
                  "Quadrature type");

    SG_GROWTH sg_growth = MODERATE_RESTRICTED;
    CLP.setOption("sg_growth", &sg_growth,
                   num_sg_growth, sg_growth_values, sg_growth_names,
                  "Sparse grid growth rule");

    CLP.parse( argc, argv );

    if (MyPID == 0) {
      std::cout << "Summary of command line options:" << std::endl
                << "\tnum_mesh        = " << n << std::endl
                << "\trand_field      = " << sg_rf_names[randField] << std::endl
                << "\tmean            = " << mean << std::endl
                << "\tstd_dev         = " << sigma << std::endl
                << "\tnum_kl          = " << num_KL << std::endl
                << "\torder           = " << p << std::endl
                << "\tnormalize_basis = " << normalize_basis << std::endl
                << "\tkrylov_method   = " << krylov_method_names[krylov_method]
                << std::endl
                << "\ttol             = " << tol << std::endl
                << "\tquadrature      = " << sg_quad_names[quad_method]
                << std::endl
                << "\tsg_growth       = " << sg_growth_names[sg_growth]
                << std::endl;
    }

    bool nonlinear_expansion = false;
    if (randField == LOGNORMAL)
      nonlinear_expansion = true;

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total Collocation Calculation Time");

    // Create Stochastic Galerkin basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL);
    for (int i=0; i<num_KL; i++) {
      RCP<Stokhos::OneDOrthogPolyBasis<int,double> > b;
      if (randField == UNIFORM) {
        b = rcp(new Stokhos::LegendreBasis<int,double>(p, normalize_basis));
      }
      else if (randField == CC_UNIFORM) {
        b =
          rcp(new Stokhos::ClenshawCurtisLegendreBasis<int,double>(
                p, normalize_basis, true));
      }
      else if (randField == RYS) {
        b = rcp(new Stokhos::RysBasis<int,double>(
                  p, weightCut, normalize_basis));
      }
      else if (randField == LOGNORMAL) {
        b = rcp(new Stokhos::HermiteBasis<int,double>(
                  p, normalize_basis));
      }
      bases[i] = b;
    }
    RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis =
      rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Create sparse grid
    RCP<Stokhos::Quadrature<int,double> > quad;
    if (quad_method == SPARSE_GRID) {
#ifdef HAVE_STOKHOS_DAKOTA
      int sparse_grid_growth = Pecos::MODERATE_RESTRICTED_GROWTH;
      if (sg_growth == SLOW_RESTRICTED)
        sparse_grid_growth = Pecos::SLOW_RESTRICTED_GROWTH;
      else if (sg_growth == MODERATE_RESTRICTED)
        sparse_grid_growth = Pecos::MODERATE_RESTRICTED_GROWTH;
      else if (sg_growth == UNRESTRICTED)
        sparse_grid_growth = Pecos::UNRESTRICTED_GROWTH;
      quad = rcp(new Stokhos::SparseGridQuadrature<int,double>(
                   basis,p,1e-12,sparse_grid_growth));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Sparse grids require building with Dakota support!");
#endif
    }
    else if (quad_method == TENSOR) {
      quad = rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    }
    const Array< Array<double> >& quad_points = quad->getQuadPoints();
    const Array<double>& quad_weights = quad->getQuadWeights();
    int nqp = quad_weights.size();

    // Create application
    twoD_diffusion_ME model(Comm, n, num_KL, sigma, mean,
                            basis, nonlinear_expansion);
    RCP<Epetra_Vector> p = rcp(new Epetra_Vector(*(model.get_p_map(0))));
    RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*(model.get_x_map())));
    RCP<Epetra_Vector> f = rcp(new Epetra_Vector(*(model.get_f_map())));
    RCP<Epetra_CrsMatrix> A =
      rcp_dynamic_cast<Epetra_CrsMatrix>(model.create_W());

    // Build commuted multipoint map and operator
    RCP<const Epetra_Map> spatial_map = model.get_x_map();
    Epetra_LocalMap stochastic_map(nqp, 0, *Comm);
    RCP<const Epetra_Map> product_map =
      rcp(BlockUtility::GenerateBlockMap(stochastic_map, *spatial_map, *Comm),
          false);
    const Epetra_CrsGraph& spatial_graph = A->Graph();
    Epetra_CrsGraph stochastic_graph(Copy, stochastic_map, 1, true);
    for (int i=0; i<nqp; ++i)
      stochastic_graph.InsertGlobalIndices(i, 1, &i);
    stochastic_graph.FillComplete();
    RCP<const Epetra_CrsGraph> product_graph =
      rcp(BlockUtility::GenerateBlockGraph(stochastic_graph, spatial_graph, *Comm));
    Epetra_CrsMatrix A_mp(Copy, *product_graph);
    Epetra_Vector f_mp(*product_map);
    Epetra_Vector x_mp(*product_map);

    if (MyPID == 0) {
      std::cout << "spatial size = " << spatial_map->NumGlobalElements()
                << " quadrature size = " << nqp
                << " multipoint size = " << product_map->NumGlobalElements()
                << std::endl;
    }

    // Loop over colloction points and assemble global operator, RHS
    int max_num_entries = A->MaxNumEntries();
    Array<int> block_indices(max_num_entries);
    double *values;
    int *indices;
    int num_entries;
    for (int qp=0; qp<nqp; qp++) {
      TEUCHOS_FUNC_TIME_MONITOR("Compute MP Matrix, RHS");

      // Set parameters
      for (int i=0; i<num_KL; i++)
        (*p)[i] = quad_points[qp][i];

      // Set in/out args
      ModelEvaluator::InArgs inArgs = model.createInArgs();
      ModelEvaluator::OutArgs outArgs = model.createOutArgs();
      inArgs.set_p(0, p);
      inArgs.set_x(x);
      outArgs.set_f(f);
      outArgs.set_W(A);

      // Evaluate model at collocation point
      x->PutScalar(0.0);
      model.evalModel(inArgs, outArgs);

      // Put f, A into global objects
      const int num_rows = spatial_map->NumMyElements();
      for (int local_row=0; local_row<num_rows; ++local_row) {
        const int global_row = spatial_map->GID(local_row);
        f_mp[nqp*local_row+qp] = (*f)[local_row];
        A->ExtractMyRowView(local_row, num_entries, values, indices);
        for (int j=0; j<num_entries; ++j) {
          int local_col = indices[j];
          int global_col = A->GCID(local_col);
          block_indices[j] = nqp*global_col+qp;
        }
        A_mp.ReplaceGlobalValues(
          nqp*global_row+qp, num_entries, values, &block_indices[0]);
      }

    }
    f_mp.Scale(-1.0);
    A_mp.FillComplete();

    //EpetraExt::RowMatrixToMatrixMarketFile("A_mp.mm", A_mp);

    // Setup preconditioner
    ParameterList precParams;
    precParams.set("default values", "SA");
    precParams.set("ML output", 10);
    precParams.set("max levels",5);
    precParams.set("increasing or decreasing","increasing");
    precParams.set("aggregation: type", "Uncoupled");
    precParams.set("smoother: type","ML symmetric Gauss-Seidel");
    precParams.set("smoother: sweeps",2);
    precParams.set("smoother: pre or post", "both");
    precParams.set("coarse: max size", 200);
    precParams.set("PDE equations",nqp);
#ifdef HAVE_ML_AMESOS
    precParams.set("coarse: type","Amesos-KLU");
#else
    precParams.set("coarse: type","Jacobi");
#endif
    RCP<ML_Epetra::MultiLevelPreconditioner> M_mp;
    {
      TEUCHOS_FUNC_TIME_MONITOR("Preconditioner Setup");
      M_mp = rcp(new ML_Epetra::MultiLevelPreconditioner(A_mp, precParams));
    }

    // Setup AztecOO solver
    AztecOO aztec;
    if (krylov_method == GMRES)
      aztec.SetAztecOption(AZ_solver, AZ_gmres);
    else if (krylov_method == CG)
      aztec.SetAztecOption(AZ_solver, AZ_cg);
    aztec.SetAztecOption(AZ_precond, AZ_none);
    aztec.SetAztecOption(AZ_kspace, 100);
    aztec.SetAztecOption(AZ_conv, AZ_r0);
    aztec.SetAztecOption(AZ_output, 10);
    aztec.SetUserOperator(&A_mp);
    aztec.SetPrecOperator(M_mp.get());
    aztec.SetLHS(&x_mp);
    aztec.SetRHS(&f_mp);

    // Solve linear system
    {
      TEUCHOS_FUNC_TIME_MONITOR("System Solve");
      aztec.Iterate(1000, tol);
    }

    // Compute responses
    RCP<Epetra_Vector> g = rcp(new Epetra_Vector(*(model.get_g_map(0))));
    Epetra_Vector g2(*(model.get_g_map(0)));
    Epetra_Vector g_mean(*(model.get_g_map(0)));
    Epetra_Vector g_var(*(model.get_g_map(0)));
    for (int qp=0; qp<nqp; qp++) {
      TEUCHOS_FUNC_TIME_MONITOR("Compute Responses");

      // Set parameters
      for (int i=0; i<num_KL; i++)
        (*p)[i] = quad_points[qp][i];

      // Extract x
      const int num_rows = spatial_map->NumMyElements();
      for (int local_row=0; local_row<num_rows; ++local_row)
        (*x)[local_row] = x_mp[nqp*local_row+qp];

      // Set in/out args
      ModelEvaluator::InArgs inArgs = model.createInArgs();
      ModelEvaluator::OutArgs outArgs = model.createOutArgs();
      inArgs.set_p(0, p);
      inArgs.set_x(x);
      outArgs.set_g(0,g);

      // Evaluate model at collocation point
      model.evalModel(inArgs, outArgs);

      // Sum contributions to mean and variance
      g2.Multiply(1.0, *g, *g, 0.0);
      g_mean.Update(quad_weights[qp], *g, 1.0);
      g_var.Update(quad_weights[qp], g2, 1.0);
    }
    g2.Multiply(1.0, g_mean, g_mean, 0.0);
    g_var.Update(-1.0, g2, 1.0);

    // Compute standard deviations
    for (int i=0; i<g_var.MyLength(); i++)
      g_var[i] = std::sqrt(g_var[i]);

    std::cout.precision(16);
    std::cout << "\nResponse Mean =      " << std::endl << g_mean << std::endl;
    std::cout << "Response Std. Dev. = " << std::endl << g_var << std::endl;

    }

    TimeMonitor::summarize(std::cout);
    TimeMonitor::zeroOutTimers();

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
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
