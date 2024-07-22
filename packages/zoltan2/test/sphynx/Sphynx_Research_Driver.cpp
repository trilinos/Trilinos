// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_KokkosCompat_DefaultNode.hpp"
#include "Zoltan2_PartitioningProblem.hpp"

#include "Zoltan2_SphynxProblem.hpp"
#include "Zoltan2_XpetraCrsGraphAdapter.hpp"

#include "Teuchos_TimeMonitor.hpp" 
#include "Teuchos_StackedTimer.hpp"
#include "readMatrixFromBinaryFile.hpp"

#include <Galeri_MultiVectorTraits.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <MatrixMarket_Tpetra.hpp>

/////////////////////////////////////////////////////////////////////////////
// This is a driver with many available options that can be used to test
// a variety of features of Sphynx.
// Note: This is research code. We do not guarantee it is without bugs. 
/////////////////////////////////////////////////////////////////////////////

/* -------------------------------------------------------------------------
 * Function: buildCrsMatrix
 * Purpose:  When the user does not input a matrix file, buildCrsMatrix will
 * 	     construct a Brick3D matrix in either 1, 2, or 3 dimensions. The
 * 	     Brick3D matrix is a 27-point difference stencil for the Laplace
 * 	     operator on a hex mesh.
 * ------------------------------------------------------------------------- */

template <typename lno_t, typename gno_t, typename scalar_t, typename nod_t>
int buildCrsMatrix(int xdim, int ydim, int zdim, std::string problemType,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
    Teuchos::RCP<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, nod_t>> &M_)
{
  /* Print size of the mesh being constructed */
  if (comm->getRank() == 0){
    std::cout << "Create matrix with " << problemType;
    std::cout << " (a " << xdim;
    if (zdim > 0)
      std::cout << " x " << ydim << " x " << zdim << " ";
    else if (ydim > 0)
      std::cout << " x"  << ydim << " x 1 ";
    else
      std::cout << "x 1 x 1 ";
    std::cout << " mesh)" << std::endl;
  }

  Teuchos::CommandLineProcessor tclp;
  Galeri::Xpetra::Parameters<gno_t> params(tclp, xdim, ydim, zdim, problemType);

  Teuchos::RCP<const Tpetra::Map<lno_t, gno_t> > map =
    Teuchos::rcp(new Tpetra::Map<lno_t, gno_t>(params.GetNumGlobalElements(), 0, comm));

  /* Build the Brick3D matrix */
  try{
    Teuchos::RCP<Galeri::Xpetra::Problem<Tpetra::Map<lno_t, gno_t>,
      Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, nod_t>,
      Tpetra::MultiVector<scalar_t, lno_t, gno_t> > > Pr=
        Galeri::Xpetra::BuildProblem<scalar_t, lno_t, gno_t,
      Tpetra::Map<lno_t, gno_t>,
      Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, nod_t>,
      Tpetra::MultiVector<scalar_t, lno_t, gno_t, nod_t> >
        (params.GetMatrixType(), map, params.GetParameterList());

    M_ = Pr->BuildMatrix();
  }
  catch (std::exception &e) {    // Probably not enough memory
    if(comm->getRank() == 0) std::cout << "Error returned from Galeri " << e.what() << std::endl;
    exit(-1);
  }
  if (M_.is_null())
    return 1;
  else
    return 0;
}

template <typename adapter_type>
  void 
compute_edgecut(Teuchos::RCP<adapter_type> &adapter,
    Zoltan2::PartitioningSolution<adapter_type> &solution )
{
  typedef typename adapter_type::user_t graph_type;
  typedef typename graph_type::global_ordinal_type GO;
  typedef typename graph_type::local_ordinal_type LO;
  typedef typename graph_type::node_type NO;
  typedef typename adapter_type::part_t PT;

  using ordinal_view_t = Kokkos::View<GO*, typename NO::device_type>;
  using part_view_t = Kokkos::View<PT*, typename NO::device_type>;

  auto graph = adapter->getUserGraph();
  auto rowMap = graph->getRowMap();
  auto colMap = graph->getColMap();

  size_t numLclRows = rowMap->getLocalNumElements();
  size_t numGblRows = rowMap->getGlobalNumElements();
  size_t numLclCols = colMap->getLocalNumElements();

  ordinal_view_t colLocalToGlobal(Kokkos::view_alloc("colLocalToGlobal", Kokkos::WithoutInitializing), numLclCols);
  auto colMapHost = Kokkos::create_mirror_view (Kokkos::HostSpace (), colLocalToGlobal);
  for(size_t i = 0; i < numLclCols; ++i)
    colMapHost[i] = colMap->getGlobalElement(i);
  Kokkos::deep_copy (colLocalToGlobal, colMapHost);

  ordinal_view_t rowLocalToGlobal(Kokkos::view_alloc("rowLocalToGlobal", Kokkos::WithoutInitializing), numLclRows);
  auto rowMapHost = Kokkos::create_mirror_view (Kokkos::HostSpace (), rowLocalToGlobal);
  for(size_t i = 0; i < numLclRows; ++i)
    rowMapHost[i] = rowMap->getGlobalElement(i);
  Kokkos::deep_copy (rowLocalToGlobal, rowMapHost);

  part_view_t localParts("localParts", numGblRows);
  part_view_t globalParts("globalParts", numGblRows);
  auto localPartsHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), localParts);

  auto parts = solution.getPartListView();
  for(size_t i = 0; i < numLclRows; i++){
    GO gi = rowMap->getGlobalElement(i);
    localPartsHost(gi) = parts[i];
  }
  Kokkos::deep_copy(localParts, localPartsHost);

  auto comm = graph->getComm();
  Teuchos::reduceAll<int, PT> (*comm, Teuchos::REDUCE_SUM, numGblRows, localParts.data(), globalParts.data());

  auto rowPtr = graph->getLocalGraphHost().row_map;
  auto colInd = graph->getLocalGraphHost().entries;

  size_t localtotalcut = 0, totalcut = 0;

  using execution_space = typename NO::device_type::execution_space;
  using range_policy = Kokkos::RangePolicy<execution_space, Kokkos::IndexType<LO>>;
  Kokkos::parallel_reduce("Compute cut", range_policy(0, numLclRows),
      KOKKOS_LAMBDA(const LO i, size_t &cut){

      const GO gRid = rowLocalToGlobal(i);
      const PT gi = globalParts(gRid);

      const size_t start = rowPtr(i);
      const size_t end = rowPtr(i+1);
      for(size_t j = start; j < end; ++j) {

      const GO gCid = colLocalToGlobal(colInd(j)); 
      PT gj = globalParts(gCid);
      if(gi != gj)
      cut += 1;
      }
      }, localtotalcut);

  Teuchos::reduceAll (*comm, Teuchos::REDUCE_SUM, 1, &localtotalcut, &totalcut);

  // compute imbalance
  auto rowPtr_h = Kokkos::create_mirror_view(rowPtr);
  Kokkos::deep_copy(rowPtr_h, rowPtr);
  int nparts = (int)solution.getTargetGlobalNumberOfParts();

  size_t *partw = new size_t[nparts];
  size_t *partc = new size_t[nparts];

  size_t *gpartw = new size_t[nparts];
  size_t *gpartc = new size_t[nparts];

  for(int i = 0; i < nparts; i++){
    partw[i] = 0; partc[i] = 0;
    gpartw[i] = 0; gpartc[i] = 0;
  }

  for(size_t i = 0; i < numLclRows; i++){
    partw[parts[i]] += rowPtr_h(i+1) - rowPtr_h(i) - 1;
    partc[parts[i]] ++;
  }

  Teuchos::reduceAll (*comm, Teuchos::REDUCE_SUM, nparts, partw, gpartw); 
  Teuchos::reduceAll (*comm, Teuchos::REDUCE_SUM, nparts, partc, gpartc); 

  size_t maxc = 0, totc = 0;
  size_t maxw = 0, totw = 0;

  for(int i = 0; i < nparts; i++){
    if(gpartw[i] > maxw)
      maxw = gpartw[i];
    if(gpartc[i] > maxc)
      maxc = gpartc[i];
    totw += gpartw[i];
    totc += gpartc[i];
  }

  double imbw = (double)maxw/((double)totw/nparts);
  double imbc = (double)maxc/((double)totc/nparts);

  if(comm->getRank() == 0) {

    std::cout << "\n\n************************************************" << std::endl;
    std::cout << "                              EDGECUT: " << totalcut <<  std::endl;
    std::cout << "                       MAX/AVG WEIGHT: " << imbw <<  std::endl;
    std::cout << "                        MAX/AVG COUNT: " << imbc <<  std::endl;
    std::cout << "************************************************\n\n" << std::endl;

  }

}

int main(int narg, char *arg[]) 
{

  Tpetra::ScopeGuard tpetraScope (&narg, &arg);
  {

    const Teuchos::RCP<const Teuchos::Comm<int>> pComm= Tpetra::getDefaultComm();

    int me = pComm->getRank();

    // Parameters
    int nparts = 64;     
    int max_iters = 1000;
    int block_size = -1;
    int rand_seed =1;
    int resFreq = 0;
    int orthoFreq = 0;
    std::string matrix_file = "";
    std::string eigensolve = "LOBPCG"; 
    bool parmetis = false;
    bool pulp = false;

    int verbosity = 1;

    std::string ptype = "";
    std::string prec = "jacobi";
    std::string init = "random";
    double tol = 1e-1;

    // Echo the command line
    if (me == 0)  {
      for (int i = 0; i < narg; i++)
        std::cout << arg[i] << " ";
      std::cout << std::endl;
    }

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("matrix_file",&matrix_file,
        "Path and filename of the matrix to be read.");
    cmdp.setOption("nparts",&nparts,
        "Number of global parts desired in the resulting partition.");
    cmdp.setOption("rand_seed",&rand_seed,
        "Seed for the random multivector.");
    cmdp.setOption("max_iters",&max_iters,
        "Maximum iters for the eigensolver.");
    cmdp.setOption("block_size",&block_size,
        "Block size.");
    cmdp.setOption("verbosity", &verbosity,
        "Verbosity level");
    cmdp.setOption("parmetis", "sphynx", &parmetis,
        "Whether to use parmetis.");
    cmdp.setOption("pulp", "sphynx", &pulp,
        "Whether to use pulp.");
    cmdp.setOption("prec", &prec,
        "Prec type to use.");
    cmdp.setOption("eigensolve", &eigensolve,
        "Eigensolver to use: LOBPCG, BlockDavidson, GeneralizedDavidson, BlockKrylovSchur or randomized.");
    cmdp.setOption("prob", &ptype,
        "Problem type to use. Options are combinatorial, normalized or generalized.");
    cmdp.setOption("tol", &tol,
        "Tolerance to use.");
    cmdp.setOption("init",  &init,
        "Sphynx Initial guess. Options: random or constants. Default: random if randomized solver is used.");
    cmdp.setOption("resFreq",  &resFreq,
        "(For randomized) Specify how often to check the residuals. Orthogonalization of the basis is also done.");
    cmdp.setOption("orthoFreq",  &orthoFreq,
        "(For randomized) Specify how often to orthogonalize the basis.");

    if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return -1;

    // Print the most essential options (not in the MyPL parameters later)
    if (me == 0 ){
      std::cout << std::endl << "--------------------------------------------------" << matrix_file << std::endl;
      std::cout << "| Sphynx Parameters                              |" << matrix_file << std::endl;
      std::cout << "--------------------------------------------------" << matrix_file << std::endl;
      std::cout << "  matrix file = " << matrix_file << std::endl;
      std::cout << "  nparts      = " << nparts << std::endl;
      std::cout << "  verbosity   = " << verbosity << std::endl;
      std::cout << "  parmetis    = " << parmetis << std::endl;
      std::cout << "  pulp        = " << pulp << std::endl;
      std::cout << "  prec        = " << prec << std::endl;
      std::cout << "  eigensolver = " << eigensolve << std::endl;
      std::cout << "  prob        = " << ptype << std::endl;
      std::cout << "  tol         = " << tol << std::endl;
      std::cout << "  init        = " << init << std::endl;
      std::cout << "  resFreq     = " << resFreq << std::endl;
      std::cout << "  orthoFreq   = " << orthoFreq << std::endl;
      std::cout << "--------------------------------------------------" << matrix_file << std::endl << std::endl;
    }

    using scalar_type = Tpetra::Details::DefaultTypes::scalar_type;
    using local_ordinal_type = Tpetra::Details::DefaultTypes::local_ordinal_type;
    using global_ordinal_type = Tpetra::Details::DefaultTypes::global_ordinal_type;
    using node_type = Tpetra::Details::DefaultTypes::node_type;

    using crs_matrix_type = Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;  
    using adapter_type = Zoltan2::XpetraCrsGraphAdapter<typename crs_matrix_type::crs_graph_type>;
    using solution_type = Zoltan2::PartitioningSolution<adapter_type>;  

    // Read the input matrix
    Teuchos::RCP<adapter_type> adapter;
    Teuchos::RCP<crs_matrix_type> tmatrix;

    // Set the random seed and hope it goes through to Tpetra.
    std::srand(rand_seed);

    // Read in user-specified matrix or create default Brick3D matrix (100^3)
    std::string mtx = ".mtx", mm = ".mm", lc = ".largestComp", lc2 = ".bin";
    if(std::equal(lc.rbegin(), lc.rend(), matrix_file.rbegin()) || std::equal(lc2.rbegin(), lc2.rend(), matrix_file.rbegin())) {
      tmatrix  = readMatrixFromBinaryFile<crs_matrix_type>(matrix_file, pComm, true, verbosity>0);
      if(me == 0 && verbosity > 1) std::cout << "Used Seher's reader for Largest Comp." << std::endl;
    }
    else if(std::equal(mtx.rbegin(), mtx.rend(), matrix_file.rbegin()) || std::equal(mm.rbegin(), mm.rend(), matrix_file.rbegin())) {
      typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
      reader_type r;
      tmatrix = r.readSparseFile(matrix_file, pComm);
      if(me == 0 && verbosity > 1) std::cout << "Used standard Matrix Market reader." << std::endl;
    }
    else {									// Build Brick3D matrix
      /* If the user entered something that was not an .mtx or .largestComp file,
       * check if it is a mesh size. If invalid input, default to 100^3 Brick3D.
       */
      int meshdim = 100;
      if(matrix_file != "")
      {
        std::string::const_iterator it = matrix_file.begin();
        while(it != matrix_file.end() && std::isdigit(*it)) ++it;
        if(it == matrix_file.end())
        {
          meshdim = std::stoi(matrix_file);
        }
        else {
          if(me == 0) std::cout << "Invalid matrix file entered. Reverting to default matrix." << std::endl;
        }
      } 
      int ierr = buildCrsMatrix<local_ordinal_type, global_ordinal_type, scalar_type>
        (meshdim, meshdim, meshdim, "Brick3D", pComm, tmatrix);
      if (ierr != 0) std::cout << "Error! Brick3D failed to build!" << std::endl;
    }

    if(me == 0 && verbosity > 0) std::cout << "Done with reading/creating the matrix." << std::endl;

    Teuchos::RCP<const Tpetra::Map<> > map = tmatrix->getMap();

    adapter = Teuchos::rcp(new adapter_type(tmatrix->getCrsGraph(), 1));
    adapter->setVertexWeightIsDegree(0);

    // Set the Sphynx parameters
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::RCP<Teuchos::ParameterList> sphynxParams(new Teuchos::ParameterList);
    params->set("num_global_parts", nparts);

    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;  
    stacked_timer = Teuchos::rcp(new Teuchos::StackedTimer("SphynxDriver"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
    if(parmetis || pulp) {

      params->set("partitioning_approach", "partition");
      params->set("imbalance_tolerance", 1.01);
      if(parmetis) {
        params->set("algorithm", "parmetis");
        params->set("imbalance_tolerance", 1.01);
      }
      else {
        params->set("algorithm", "pulp");
        params->set("pulp_vert_imbalance", 1.01);
      }

      using problem_type = Zoltan2::PartitioningProblem<adapter_type>;
      Teuchos::RCP<problem_type> problem;
      pComm->barrier();
      {
        Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::All"));
        {
          Teuchos::TimeMonitor t2(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Problem"));
          problem = Teuchos::rcp(new problem_type(adapter.getRawPtr(), params.getRawPtr(), Tpetra::getDefaultComm()));
        }
        {
          Teuchos::TimeMonitor t3(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Solve"));
          problem->solve();
        }
      }
      pComm->barrier();

      solution_type solution = problem->getSolution();
      compute_edgecut<adapter_type>(adapter, solution);

    } else {

      sphynxParams->set("sphynx_verbosity", verbosity);
      sphynxParams->set("sphynx_max_iterations", max_iters);
      sphynxParams->set("sphynx_ortho_freq", orthoFreq);
      sphynxParams->set("sphynx_res_freq", resFreq);
      sphynxParams->set("sphynx_skip_preprocessing", true);
      sphynxParams->set("sphynx_eigensolver", eigensolve);
      if (block_size > 0) sphynxParams->set("sphynx_block_size", block_size);
      if (ptype != "") sphynxParams->set("sphynx_problem_type", ptype);
      if (init != "") sphynxParams->set("sphynx_initial_guess", init);
      if (prec != "") sphynxParams->set("sphynx_preconditioner_type", prec);
      if (tol != -1) sphynxParams->set("sphynx_tolerance", tol);

      using problem_type = Zoltan2::SphynxProblem<adapter_type>; //We found sphynx
      Teuchos::RCP<problem_type> problem;
      pComm->barrier();
      {
        Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::All"));
        {
          Teuchos::TimeMonitor t2(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Problem"));
          problem = Teuchos::rcp(new problem_type(adapter.get(), params.get(), sphynxParams, Tpetra::getDefaultComm()));
        }
        {
          if(me == 0) std::cout << eigensolve << " will be used to solve the partitioning problem." << std::endl;
          problem->solve();
        }
        pComm->barrier();
      }
      solution_type solution = problem->getSolution();
      compute_edgecut<adapter_type>(adapter, solution);
    }
    stacked_timer->stopBaseTimer();       

    Teuchos::RCP<Teuchos::FancyOStream> fancy2 = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out2 = *fancy2;
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(out2, pComm, options);

    Teuchos::TimeMonitor::summarize();

  } //End Tpetra scope guard
  return 0;
} 
