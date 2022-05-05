#include "Teuchos_CommandLineProcessor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Zoltan2_PartitioningProblem.hpp"

#include "Zoltan2_SphynxProblem.hpp"
#include "Zoltan2_XpetraCrsGraphAdapter.hpp"

#include "Teuchos_TimeMonitor.hpp" 
#include "Teuchos_StackedTimer.hpp"
#include "ReadMatrixFromFile.hpp"

#include <Galeri_MultiVectorTraits.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <MatrixMarket_Tpetra.hpp>

template <typename lno_t, typename gno_t, typename scalar_t, typename nod_t>
int buildCrsMatrix(int xdim, int ydim, int zdim, std::string problemType,
                   const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
		   Teuchos::RCP<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, nod_t>> &M_)
{
  if (comm->getRank() == 0){
    std::cout << "Create matrix with " << problemType;
    std::cout << " (and " << xdim;
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
    std::cout << "Error returned from Galeri " << e.what() << std::endl;
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

  size_t numLclRows = rowMap->getNodeNumElements();
  size_t numGblRows = rowMap->getGlobalNumElements();
  size_t numLclCols = colMap->getNodeNumElements();

  
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
  for(LO i = 0; i < numLclRows; i++){
    
    GO gi = rowMap->getGlobalElement(i);
    localPartsHost(gi) = parts[i];
  }
  Kokkos::deep_copy(localParts, localPartsHost);
  
  auto comm = graph->getComm();
  Teuchos::reduceAll<int, PT> (*comm, Teuchos::REDUCE_SUM, numGblRows, localParts.data(), globalParts.data());

  auto rowPtr = graph->getLocalGraph().row_map;
  auto colInd = graph->getLocalGraph().entries;

  
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

  for(LO i = 0; i < numLclRows; i++){
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

template <typename adapter_type>
void 
compute_edgecut_old(Teuchos::RCP<adapter_type> &adapter,
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

  size_t numLclRows = rowMap->getNodeNumElements();
  size_t numGblRows = rowMap->getGlobalNumElements();
  size_t numLclCols = colMap->getNodeNumElements();

  
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

  
  part_view_t localParts(Kokkos::view_alloc("localParts", Kokkos::WithoutInitializing), numGblRows);
  //part_view_t localParts("localParts", numGblRows);
  part_view_t globalParts("globalParts", numGblRows);
  auto localPartsHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), localParts);
  
  auto parts = solution.getPartListView();
  for(LO i = 0; i < numLclRows; i++){
    
    GO gi = rowMap->getGlobalElement(i);
    localPartsHost(gi) = parts[i];
  }
  Kokkos::deep_copy(localParts, localPartsHost);
  
  auto comm = graph->getComm();
  Teuchos::reduceAll<int, PT> (*comm, Teuchos::REDUCE_SUM, numGblRows, localParts.data(), globalParts.data());

  auto rowPtr = graph->getLocalGraph().row_map;
  auto colInd = graph->getLocalGraph().entries;

  
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

  if(comm->getRank() == 0) {
    std::cout << "\n\n************************************************" << std::endl;
    std::cout << "                              EDGECUT: " << totalcut <<  std::endl;
    std::cout << "************************************************\n\n" << std::endl;

  }
  
}


int main(int narg, char *arg[]) 
{

  Tpetra::ScopeGuard tpetraScope (&narg, &arg);
  {

    const Teuchos::RCP<const Teuchos::Comm<int>> pComm= Tpetra::getDefaultComm();

    int me = pComm->getRank();
    int np = pComm->getSize();

    // Parameters
    int nparts = 64;     
    std::string matrix_file = "";
    std::string vector_file = "";
    bool parmetis = false;
    bool pulp = false;
    
    int verbosity = 0;
    
    std::string ptype = "";
    std::string prec = "";
    std::string init = "";
    double tol = -1;
    
    // Echo the command line
    if (me == 0)  {
      for (int i = 0; i < narg; i++)
	std::cout << arg[i] << " ";
      std::cout << std::endl;
    }

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("matrix_file",&matrix_file,
		   "Path and filename of the matrix to be read.");
    cmdp.setOption("vector_file",&vector_file,
		   "Path and filename of the vector to be read.");
    cmdp.setOption("nparts",&nparts,
		   "Number of global parts desired in the resulting partition.");
    cmdp.setOption("verbosity", &verbosity,
		   "Verbosity level");
    cmdp.setOption("parmetis", "sphynx", &parmetis,
		   "Whether to use parmetis.");
    cmdp.setOption("pulp", "sphynx", &pulp,
		   "Whether to use pulp.");
    cmdp.setOption("prec", &prec,
		   "Prec type to use.");
    cmdp.setOption("prob", &ptype,
		   "Problem type to use.");
    cmdp.setOption("tol", &tol,
		   "Tolerance to use.");
    cmdp.setOption("init",  &init,
		   "Initial guess.");

    if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }


    // Print the most essential options (not in the MyPL parameters later)
    if (me==0){
      std::cout << "matrix file = " << matrix_file << std::endl;
      std::cout << "vector file = " << vector_file << std::endl;
      std::cout << "nparts = " << nparts << std::endl;
      std::cout << "verbosity = " << verbosity << std::endl;
      std::cout << "parmetis = " << parmetis << std::endl;
      std::cout << "pulp = " << pulp << std::endl;
      std::cout << "prec = " << prec << std::endl;
      std::cout << "prob = " << ptype << std::endl;
      std::cout << "tol = " << tol << std::endl;
      std::cout << "init = " << init << std::endl;

    }

    using scalar_type = double;
    using local_ordinal_type = int;
    using global_ordinal_type = long long;
    //using node_type = Kokkos::Compat::KokkosCudaWrapperNode;
    //using node_type = Kokkos::Compat::KokkosHIPWrapperNode;
    //using node_type = Kokkos::Compat::KokkosSerialWrapperNode; 
    using node_type = Kokkos::Compat::KokkosOpenMPWrapperNode; 
    
    using crs_matrix_type = Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;  
    using adapter_type = Zoltan2::XpetraCrsGraphAdapter<typename crs_matrix_type::crs_graph_type>;
    using solution_type = Zoltan2::PartitioningSolution<adapter_type>;  


    // Read the input matrix
    Teuchos::RCP<adapter_type> adapter;
    Teuchos::RCP<crs_matrix_type> tmatrix;



    std::string mtx = ".mtx", lc = ".largestComp";
    if(std::equal(lc.rbegin(), lc.rend(), matrix_file.rbegin())) {
      tmatrix  = readMatrixFromFile<crs_matrix_type>(matrix_file, pComm, true, verbosity>0);
    }
    else if(std::equal(mtx.rbegin(), mtx.rend(), matrix_file.rbegin())) {
      typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
      reader_type r;
      tmatrix = r.readSparseFile(matrix_file, pComm);
    }
    else {
      int meshdim = 100;
      if(matrix_file == "200")
    	meshdim = 200;
      else if(matrix_file == "400")
    	meshdim = 400;
      int ierr = buildCrsMatrix<local_ordinal_type, global_ordinal_type, scalar_type>
    	(meshdim, meshdim, meshdim, "Brick3D", pComm, tmatrix);
      //Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile(matrix_file+".mtx", tmatrix);
    }
    if(me == 0)
      std::cout << "Done with reading/creating the matrix." << std::endl;
  Teuchos::RCP<const Tpetra::Map<> > map = tmatrix->getMap();
  using ST = double;
  using MultiVector  = Tpetra::MultiVector<ST>;
    typedef Belos::MultiVecTraits<ST,MultiVector>    MVT;
  Teuchos::RCP<MultiVector > V = Teuchos::rcp(new MultiVector(map, 10));
  V->randomize();
  std::cout << "Randomized some pretend evects." << std::endl;
   //MVT::MvPrint(*V,std::cout);
   // TODO Insert MultiVector Reader 
  /*if(vector_file ~= ""){
  else{
  }
  }*/
    adapter = Teuchos::rcp(new adapter_type(tmatrix->getCrsGraph(), 1));
    adapter->setVertexWeightIsDegree(0);
    
    // Set the parameters
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
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
	  Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Problem"));
	  problem = Teuchos::rcp(new problem_type(adapter.getRawPtr(), params.getRawPtr(), Tpetra::getDefaultComm()));
	}
	{
	  Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Solve"));
	  problem->solve();
	}
      }
      pComm->barrier();
      
      solution_type solution = problem->getSolution();
      compute_edgecut<adapter_type>(adapter, solution);
	
    }
    else {
   
      params->set("sphynx_verbosity", verbosity);
      params->set("sphynx_skip_preprocessing", true);
      if (ptype != "") params->set("sphynx_problem_type", ptype);
      if (init != "") params->set("sphynx_initial_guess", init);
      if (prec != "") params->set("sphynx_preconditioner_type", prec);
      if (tol != -1) params->set("sphynx_tolerance", tol);
      
      using problem_type = Zoltan2::SphynxProblem<adapter_type>; //We found sphynx
      Teuchos::RCP<problem_type> problem;
      pComm->barrier();
      {
  	Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::All"));
	{
	  Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Problem"));
	  problem = Teuchos::rcp(new problem_type(adapter, params, Tpetra::getDefaultComm()));
	}//TODO: Change constructor problem_type to take in vector_file *if it exists*//
	{
	  Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Partitioning::Solve"));
	  problem->solve(V);
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
