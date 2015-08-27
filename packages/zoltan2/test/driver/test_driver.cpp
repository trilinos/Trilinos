// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// taking headers from existing driver template
// will keep or remove as needed
#include <UserInputForTests.hpp>
#include <AdapterForTests.hpp>
#include <Zoltan2_ComparisonHelper.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolutionQuality.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Zoltan2_Parameters.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <TpetraExt_MatrixMatrix_def.hpp>


#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <queue>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::XMLObject;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::map;
using std::pair;
using std::exception;
using std::ostringstream;
using std::queue;

#define ERRMSG(msg) if (rank == 0){ cerr << "FAIL: " << msg << endl; }
#define EXC_ERRMSG(msg, e) \
if (rank==0){ cerr << "FAIL: " << msg << endl << e.what() << endl;}

// temporary methods for debugging and leanring
void readXML(const XMLObject &xml, const string &title);
void readPList(const ParameterList &plist,
               const string &title,
               bool doc = false,
               bool unused = false);

typedef Zoltan2::MetricValues<zscalar_t> metric_t; // typedef metric_type

void xmlToModelPList(const Teuchos::XMLObject &xml, Teuchos::ParameterList & plist)
{
  // This method composes a plist for the problem
  Teuchos::XMLParameterListReader reader;
  plist = reader.toParameterList(xml);
  
  //  Get list of valid Zoltan2 Parameters
  // Zoltan 2 parameters appear in the input file
  // Right now we have default values stored in
  // the parameter list, we would like to apply
  // the options specified by the user in their
  // input file
  Teuchos::ParameterList zoltan2Parameters;
  Zoltan2::createAllParameters(zoltan2Parameters);
  
  if (plist.isSublist("Zoltan2Parameters")) {
    // Apply user specified zoltan2Parameters
    ParameterList &sub = plist.sublist("Zoltan2Parameters");
    zoltan2Parameters.setParameters(sub);
  }
  
  zoltan2Parameters.set("compute_metrics", "true");
  // update zoltan 2 parameters
  //        plist.remove("Zoltan2Parameters");
  //        plist.set("Zoltan2Parameters", zoltan2Parameters);
}


void getParameterLists(const string &inputFileName,
                       queue<ParameterList> &problems,
                       queue<ParameterList> &comparisons,
                       const RCP<const Teuchos::Comm<int> > & comm)
{
  int rank = comm->getRank();
  // return a parameter list of problem definitions
  // and a parameter list for solution comparisons
  Teuchos::FileInputSource inputSource(inputFileName);
  cout << "input file source: " << inputFileName << endl;
  XMLObject xmlInput;
  
  // Try to get xmlObject from inputfile
  try{
    xmlInput = inputSource.getObject();
  }
  catch(exception &e)
  {
    EXC_ERRMSG("Test Driver error: reading", e); // error reading input
  }
  
  // get the parameter lists for each model
  for(int i = 0; i < xmlInput.numChildren(); i++)
  {
    ParameterList plist;
    xmlToModelPList(xmlInput.getChild(i), plist);
    if(plist.name() == "Comparison") comparisons.emplace(plist);
    else problems.emplace(plist);
  }
  
}

bool MinMaxTest(const metric_t & metric,
                const Teuchos::ParameterList & metricPlist,
                ostringstream &msg)
{
  // run a comparison of min and max agains a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metric.getName() + " test";
  if (metricPlist.isParameter("lower"))
  {
    double min = metricPlist.get<double>("lower");
    
    if(metric.getMinImbalance() < min)
    {
      msg << test_name << " FAILED: Minimum imbalance per part, "
      << metric.getMinImbalance() <<
      ", less than specified allowable minimum, " << min;
      pass = false;
    }
  }
  
  if(metricPlist.isParameter("upper" ) && pass != false) {
    double max = metricPlist.get<double>("upper");
    if (metric.getMaxImbalance() > max)
    {
      msg << test_name << " FAILED: Maximum imbalance per part, "
      << metric.getMaxImbalance() <<
      ", greater than specified allowable maximum, " << max;
      pass = false;
    }
    
  }
  
  if(pass){
    msg << test_name << " PASSED.";
    pass = true;
  }
  
  return pass;
}


template<typename T>
void writePartionSolution(const T * part, size_t N, const RCP<const Teuchos::Comm<int> > & comm)
{
  std::ofstream file;
  char title[256];
  
  static const string path = "/Users/davidson/trilinosall/trilinosrepo/packages/zoltan2/test/driver/pamgen_mesh_data";
  
  sprintf(title, "/partition_%d", comm->getRank());
  file.open(path + string(title));
  for (size_t i = 0; i < N; i++) {
    file << part[i] << "\n";
  }
  
  file.close();
  
}

void run(const UserInputForTests &uinput,
         const ParameterList &problem_parameters,
         RCP<ComparisonHelper> & comparison_helper,
         const RCP<const Teuchos::Comm<int> > & comm)
{
  // Major steps in running a problem in zoltan 2
  // 1. get an input adapter
  // 2. construct the problem
  // 3. solve the problem
  // 4. analyze metrics
  // 5. clean up
  
  typedef AdapterForTests::base_adapter_t base_t;
  typedef AdapterForTests::basic_id_t basic_id_t; // basic_identifier_type
  typedef AdapterForTests::xpetra_mv_adapter xpetra_mv_t; // xpetra_mv_type
  typedef AdapterForTests::xcrsGraph_adapter xcrsGraph_t;
  typedef AdapterForTests::xcrsMatrix_adapter xcrsMatrix_t;
  typedef AdapterForTests::basic_vector_adapter basic_vector_t;
  typedef AdapterForTests::pamgen_adapter_t pamgen_t;

  typedef Zoltan2::Problem<base_t> problem_t;
  typedef Zoltan2::PartitioningProblem<base_t> partioning_problem_t; // base abstract type
  typedef Zoltan2::PartitioningProblem<basic_id_t> basic_problem_t; // basic id problem type
  typedef Zoltan2::PartitioningProblem<xpetra_mv_t> xpetra_mv_problem_t; // xpetra_mv problem type
  typedef Zoltan2::PartitioningProblem<xcrsGraph_t> xcrsGraph_problem_t; // xpetra_graph problem type
  typedef Zoltan2::PartitioningProblem<xcrsMatrix_t> xcrsMatrix_problem_t; // xpetra_matrix problem type
  typedef Zoltan2::PartitioningProblem<basic_vector_t> basicVector_problem_t; // vector problem type
  typedef Zoltan2::PartitioningProblem<pamgen_t> pamgen_problem_t; // pamgen mesh problem type

  int rank = comm->getRank();
  if(rank == 0)
    cout << "\nPeforming test: " << problem_parameters.get<string>("Name") << endl;
  
  
  ////////////////////////////////////////////////////////////
  // 1. get basic input adapter
  ////////////////////////////////////////////////////////////
  if(!problem_parameters.isParameter("InputAdapterParameters"))
    throw std::runtime_error("Input adapter parameters not provided");
  if(!problem_parameters.isParameter("Zoltan2Parameters"))
    throw std::runtime_error("Zoltan2 probnlem parameters not provided");
  
  const ParameterList &adapterPlist = problem_parameters.sublist("InputAdapterParameters");
  base_t * ia = AdapterForTests::getAdapterForInput(const_cast<UserInputForTests *>(&uinput), adapterPlist,comm); // a pointer to a basic type
  if(ia == nullptr)
  {
    if(rank == 0)
      cout << "Get adapter for input failed" << endl;
    return;
  }
  
  ////////////////////////////////////////////////////////////
  // 2. construct partitioning problem
  ////////////////////////////////////////////////////////////
  problem_t * problem;
  string adapter_name = adapterPlist.get<string>("inputAdapter"); // If we are here we have an input adapter, no need to check for one.
  // get Zoltan2 partion parameters
  ParameterList zoltan2_parameters = const_cast<ParameterList &>(problem_parameters.sublist("Zoltan2Parameters"));
  zoltan2_parameters.set("num_global_parts", comm->getSize());
  
  if(rank == 0){
    readPList(zoltan2_parameters, "Zoltan 2 Params:\n");
    cout <<"\n\n"<<endl;}
  
#ifdef HAVE_ZOLTAN2_MPI
  
  if(adapter_name == "BasicIdentifier"){
    problem = reinterpret_cast<problem_t * >(new basic_problem_t(reinterpret_cast<basic_id_t *>(ia),
                                                                 &zoltan2_parameters,
                                                                 MPI_COMM_WORLD));
  }else if(adapter_name == "XpetraMultiVector")
  {
    problem = reinterpret_cast<problem_t * >(new xpetra_mv_problem_t(reinterpret_cast<xpetra_mv_t *>(ia),
                                                                     &zoltan2_parameters,
                                                                     MPI_COMM_WORLD));
  }else if(adapter_name == "XpetraCrsGraph"){
    problem = reinterpret_cast<problem_t * >(new xcrsGraph_problem_t(reinterpret_cast<xcrsGraph_t *>(ia),
                                                                     &zoltan2_parameters,
                                                                     MPI_COMM_WORLD));
  }
  else if(adapter_name == "XpetraCrsMatrix")
  {
    problem = reinterpret_cast<problem_t * >(new xcrsMatrix_problem_t(reinterpret_cast<xcrsMatrix_t *>(ia),
                                                                      &zoltan2_parameters,
                                                                      MPI_COMM_WORLD));
  }else if(adapter_name == "BasicVector")
  {
    problem = reinterpret_cast<problem_t * >(new basicVector_problem_t(reinterpret_cast<basic_vector_t *>(ia),
                                                                       &zoltan2_parameters,
                                                                       MPI_COMM_WORLD));
  }else if(adapter_name == "PamgenMesh")
  {
    problem = reinterpret_cast<problem_t * >(new pamgen_problem_t(reinterpret_cast<pamgen_t *>(ia),
                                                                       &zoltan2_parameters,
                                                                       MPI_COMM_WORLD));
  }
  else
    throw std::runtime_error("Input adapter type not available, or misspelled.");
  
  
#else
  if(adapter_name == "BasicIdentifier"){
    problem = reinterpret_cast<problem_t * >(new basic_problem_t(reinterpret_cast<basic_id_t *>(ia),
                                                                 &zoltan2_parameters));
  }else if(adapter_name == "XpetraMultiVector")
  {
    problem = reinterpret_cast<problem_t * >(new xpetra_mv_problem_t(reinterpret_cast<xpetra_mv_t *>(ia),
                                                                     &zoltan2_parameters));
  }else if(adapter_name == "XpetraCrsGraph"){
    problem = reinterpret_cast<problem_t * >(new xcrsGraph_problem_t(reinterpret_cast<xcrsGraph_t *>(ia),
                                                                     &zoltan2_parameters));
  }
  else if(adapter_name == "XpetraCrsMatrix")
  {
    problem = reinterpret_cast<problem_t * >(new xcrsMatrix_problem_t(reinterpret_cast<xcrsMatrix_t *>(ia),
                                                                      &zoltan2_parameters));
  } else if(adapter_name == "BasicVector")
  {
    problem = reinterpret_cast<problem_t * >(new basicVector_problem_t(reinterpret_cast<basic_vector_t *>(ia),
                                                                       &zoltan2_parameters));
  }else if(adapter_name == "PamgenMesh")
  {
    problem = reinterpret_cast<problem_t * >(new pamgen_problem_t(reinterpret_cast<pamgen_t *>(ia),
                                                                  &zoltan2_parameters));
  }
  else
    throw std::runtime_error("Input adapter type not available, or misspelled.");
#endif
  
  ////////////////////////////////////////////////////////////
  // 3. Solve the problem
  ////////////////////////////////////////////////////////////
  reinterpret_cast<basic_problem_t *>(problem)->solve();
  if (rank == 0)
    cout << "Problem solved.\n" << endl;
  
  ////////////////////////////////////////////////////////////
  // 4. Print problem metrics
  ////////////////////////////////////////////////////////////
  
  if (comm->getRank() == 0)
  {
    // calculate pass fail based on imbalance
    if(rank == 0) cout << "analyzing metrics...\n" << endl;
    if(problem_parameters.isParameter("Metrics"))
    {
      reinterpret_cast<basic_problem_t *>(problem)->printMetrics(cout);
      ArrayRCP<const metric_t> metrics
      = reinterpret_cast<basic_problem_t *>(problem)->getMetrics();
      
      // get metric plist
      const ParameterList &metricsPlist = problem_parameters.sublist("Metrics");
      
      string test_name;
      bool all_tests_pass = true;
      for(int i = 0; i < metrics.size(); i++)
      {
        // print their names...
        ostringstream msg;
        test_name = metrics[i].getName();
        if(metricsPlist.isSublist(test_name))
        {
          if(!MinMaxTest(metrics[i], metricsPlist.sublist(test_name), msg))
            all_tests_pass = false;
          cout << msg.str() << endl;
          
        }
      }
      
      if(all_tests_pass) cout << "All tests PASSED." << endl;
      else cout << "Testing FAILED." << endl;
      
    }else{
      cout << "No test metrics provided." << endl;
      reinterpret_cast<basic_problem_t *>(problem)->printMetrics(cout);
    }
  }
  
  // 4b. timers
  if(zoltan2_parameters.isParameter("timer_output_stream"))
    reinterpret_cast<basic_problem_t *>(problem)->printTimers();
  
  ////////////////////////////////////////////////////////////
  // 5. Add solution to map for possible comparison testing
  ////////////////////////////////////////////////////////////
  ComparisonSource * comparison_source = new ComparisonSource;
  comparison_source->adapter = RCP<basic_id_t>(reinterpret_cast<basic_id_t *>(ia));
  comparison_source->problem = RCP<basic_problem_t>(reinterpret_cast<basic_problem_t *>(problem));
  comparison_source->problem_kind = problem_parameters.isParameter("kind") ? problem_parameters.get<string>("kind") : "?";
  comparison_source->adapter_kind = adapter_name;
  comparison_helper->AddSource(problem_parameters.name(), comparison_source);
  
  
  // write mesh solution
  auto sol = reinterpret_cast<basic_problem_t *>(problem)->getSolution();
  writePartionSolution(sol.getPartListView(), ia->getLocalNumIDs(), comm);

  ////////////////////////////////////////////////////////////
  // 6. Clean up
  ////////////////////////////////////////////////////////////
//  if(adapter_name == "XpetraCrsGraph")
//    delete reinterpret_cast<xcrsGraph_t *>(ia)->getCoordinateInput();
//  if(adapter_name == "XpetraCrsMatrix")
//    delete reinterpret_cast<xcrsMatrix_t *>(ia)->getCoordinateInput();
//  
//  delete ia;
//  delete reinterpret_cast<basic_problem_t *>(problem);
}


void readMesh(const UserInputForTests &uinput,const RCP<const Teuchos::Comm<int> > & comm)
{
  comm->barrier();
  PamgenMesh * mesh = const_cast<UserInputForTests *>(&uinput)->getPamGenMesh();
  printf("\n\nProc %d mesh report:\n", comm->getRank());
  
  int nodes, els;
  nodes = mesh->num_nodes;
  els = mesh->num_elem;
  
  printf("dimension: %d\n", mesh->num_dim);
  printf("local nodes: %d\n", nodes);
  printf("local elems: %d\n", els);
  
  int gnodes, gels;
  gnodes = mesh->num_nodes_global;
  gels = mesh->num_elems_global;
  printf("global nodes: %d\n", gnodes);
  printf("global elem: %d\n", gels);
  
  int blks = mesh->num_elem_blk;
  printf("num blocks: %d\n", blks);
  
  printf("\ncoordinates:\n");
  double * coord = mesh->coord;
  for (int i = 0; i < nodes; i++) {
    if(mesh->num_dim == 2)
    {
      printf("lid %d, gid %d: {%1.2f, %1.2f}\n",i,mesh->global_node_numbers[i],
             coord[i], coord[nodes+i]);
    }else
    {
      printf("lid %d, gid %d: {%1.2f, %1.2f, %1.2f}\n",i,mesh->global_node_numbers[i],
             coord[i], coord[nodes+i], coord[2*nodes+i]);
    }
  }
  
  
  int el_count = 0;
  printf("\nElements:\n");
  for(int i = 0; i < blks; i++)
  {
    int elb = mesh->elements[i];
    int npel = mesh->nodes_per_element[i];
    printf("blkid %d has %d els, with %d nodes/element\n", mesh->block_id[i], elb, npel);
    // nodes for el
    int * connect = mesh->elmt_node_linkage[i];
    
    for(int j = 0; j < elb; j++)
    {
      printf("element %d:{", el_count);
      for(int k = 0; k < npel; k++)
        printf("%d ", connect[j*npel + k]);
      printf("}\n");
      el_count++;
    }
    
  }
  
  el_count = 0;
  printf("\nElements Centroids:\n");
  for(int i = 0; i < blks; i++)
  {
    int elb = mesh->elements[i];
    for(int j = 0; j < elb; j++)
    {
      printf("element %d:{", mesh->global_element_numbers[el_count]);
      for(int k = 0; k < mesh->num_dim; k++)
        printf("%1.2f ", mesh->element_coord[el_count + k * mesh->num_elem]);
      printf("}\n");
      el_count++;
    }
  }
  
  comm->barrier();
}


void writeMesh(const UserInputForTests &uinput,const RCP<const Teuchos::Comm<int> > & comm)
{
  comm->barrier();
  std::ofstream file;
  char title[256];
  
  static const string path = "/Users/davidson/trilinosall/trilinosrepo/packages/zoltan2/test/driver/pamgen_mesh_data";
  
  PamgenMesh * mesh = const_cast<UserInputForTests *>(&uinput)->getPamGenMesh();
  
  int nodes, els;
  nodes = mesh->num_nodes;
  els = mesh->num_elem;
  
  int gnodes, gels;
  gnodes = mesh->num_nodes_global;
  gels = mesh->num_elems_global;
  
  int blks = mesh->num_elem_blk;
  
  sprintf(title, "/coordinates_%d", comm->getRank());
  file.open(path + string(title));
  double * coord = mesh->coord;
  for (int i = 0; i < nodes; i++) {
    file << coord[i] << "\t" << coord[nodes + i];
    if(mesh->num_dim == 3) file << "\t" << coord[2*nodes +i];
    file << "\n";
  }
  
  file.close();
  
  // write all elements
  sprintf(title, "/elements_%d", comm->getRank());
  file.open(path + string(title));
  for(int i = 0; i < blks; i++)
  {
    int elb = mesh->elements[i];
    int npel = mesh->nodes_per_element[i];
    // nodes for el
    int * connect = mesh->elmt_node_linkage[i];
    
    for(int j = 0; j < elb; j++)
    {
      for(int k = 0; k < npel-1; k++)
        file << connect[j*npel + k] << "\t";
      
      file << connect[j*npel + npel-1] << "\n";
    }
  }
  
  file.close();
  
  // write global element numbers
  sprintf(title, "/global_element_id_%d", comm->getRank());
  file.open(path + string(title));
  size_t el = 0;
  for(int i = 0; i < blks; i++)
  {
    int elb = mesh->elements[i];
    for(int j = 0; j < elb; j++)
    {
      file << mesh->global_element_numbers[el++] <<"\n";
    }
  }
  
  file.close();
  
  // write boundary faces
  if(mesh->num_dim == 3 && comm->getRank() == 0)
  {
    sprintf(title, "/element_map");
    file.open(path + string(title));
    file << 1 <<"\t"<<5<<"\t"<<8<<"\t"<<4<<"\n";
    file << 2 <<"\t"<<3<<"\t"<<7<<"\t"<<6<<"\n";
    file << 1 <<"\t"<<2<<"\t"<<6<<"\t"<<5<<"\n";
    file << 4 <<"\t"<<8<<"\t"<<7<<"\t"<<3<<"\n";
    file << 1 <<"\t"<<4<<"\t"<<3<<"\t"<<2<<"\n";
    file << 5 <<"\t"<<6<<"\t"<<7<<"\t"<<8<<"\n";
    file.close();
  }
  
  sprintf(title, "/element_centroids_%d",comm->getRank());
  file.open(path + string(title));
  el = 0;
  for(int i = 0; i < blks; i++)
  {
    int elb = mesh->elements[i];
    for(int j = 0; j < elb; j++)
    {
      file << "Element " << el+1 <<":\t";
      for(int k = 0; k < mesh->num_dim; k++)
        file << mesh->element_coord[el + k * mesh->num_elem] << "\t";
      file <<"\n";
      el++;
    }
  }
  file.close();
  // write info
  if(comm->getRank() == 0)
  {
    sprintf(title,"/mesh_info");
    file.open(path + string(title));
    file << comm->getSize() <<"\n" << nodes << "\n" << els << "\n" << blks;
    file.close();
  }
  comm->barrier();
}

void getConnectivityGraph(const UserInputForTests &uinput,const RCP<const Teuchos::Comm<int> > & comm)
{
  comm->barrier(); // wait for everyone
  int rank = comm->getRank();
  if(rank == 0) cout << "Making a graph from our pamgen mesh...." << endl;
  
  typedef Tpetra::Vector<>::local_ordinal_type scalar_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
  typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> crs_matrix_type;
  typedef Tpetra::CrsGraph<> crs_graph_type;
  
  // get mesh
  PamgenMesh * mesh = const_cast<UserInputForTests *>(&uinput)->getPamGenMesh();
  
  // get info for setting up map
  int local_nodes, local_els;
  local_nodes = mesh->num_nodes;
  local_els = mesh->num_elem;
  
  int global_nodes, global_els;
  global_nodes = mesh->num_nodes_global;
  global_els = mesh->num_elems_global;
  
  // make map with global elements assigned to this mesh
  const zgno_t idxBase = 0;
  Teuchos::ArrayView<int> g_el_ids(mesh->global_element_numbers,local_els);
  for(auto && v : g_el_ids) v--; // shif to idx base 0
  RCP<const map_type> range_map = rcp(new map_type(static_cast<Tpetra::global_size_t>(global_els),
                                                   g_el_ids,
                                                   idxBase,
                                                   comm));
  
  // make domain map
  Teuchos::ArrayView<int> g_node_ids(mesh->global_node_numbers,local_nodes);
  for(auto && v : g_node_ids) v--; // shif to idx base
  RCP<const map_type> domain_map = rcp(new map_type(static_cast<Tpetra::global_size_t>(global_nodes),
                                                    g_node_ids,
                                                    idxBase,
                                                    comm));
  
  
  // make a connectivity matrix
  Teuchos::RCP<crs_matrix_type> C = rcp(new crs_matrix_type(range_map,domain_map,0));
  // write all nodes per el to matrix
  int blks = mesh->num_elem_blk;
  // write all elements
  zlno_t el_no = 0;
  scalar_type one = static_cast<scalar_type>(1);
  for(int i = 0; i < blks; i++)
  {
    int el_per_block = mesh->elements[i];
    int nodes_per_el = mesh->nodes_per_element[i];
    int * connect = mesh->elmt_node_linkage[i];
    
    for(int j = 0; j < el_per_block; j++)
    {
      const global_ordinal_type gid = static_cast<global_ordinal_type>(g_el_ids[el_no]);
      //      const global_ordinal_type gid = domain_map->getGlobalElement(el_no);
      for(int k = 0; k < nodes_per_el; k++)
      {
        int g_node_i = g_node_ids[connect[j*nodes_per_el+k]-1];
        //        printf("inserting [%d, %d]\n", gid, g_node_i);
        C->insertGlobalValues(gid,
                              Teuchos::tuple<global_ordinal_type>(g_node_i),
                              Teuchos::tuple<zlno_t>(one));
      }
      el_no++;
    }
  }
  
  if(rank == 0) cout << "Call Fill complete..." << endl;
  C->fillComplete(domain_map, range_map);
  //  C->fillComplete();
  
  if(rank == 0) cout << "Calculating adjacency matrix of pamgen mesh: C*C' = A..." << endl;
  // Matrix multiply by Transpose to get El connectivity
  RCP<crs_matrix_type> A = rcp(new crs_matrix_type(range_map,0));
  Tpetra::MatrixMatrix::Multiply(*C, false, *C, true, *A);
  if(rank == 0) cout << "Completed Multiply" << endl;
//  if(rank == 0)
//  {
//    cout << "C: \n" << endl;
//    C->print(std::cout);
//    
//    cout <<"\nA:\n" << endl;
//    A->print(std::cout);
//  }
  // remove entris not adjacent
  // make graph
  RCP<crs_matrix_type> modA = rcp(new crs_matrix_type(range_map,0));
  
  if(rank == 0) cout << "Setting graph of connectivity..." << endl;
  for(zgno_t gid : range_map->getNodeElementList())
  {
    size_t numEntriesInRow = A->getNumEntriesInGlobalRow (gid);
    Array<crs_matrix_type::scalar_type> rowvals (numEntriesInRow);
    Array<global_ordinal_type> rowinds (numEntriesInRow);
    
    // modified
    Array<scalar_type> mod_rowvals;
    Array<global_ordinal_type> mod_rowinds;
    A->getGlobalRowCopy (gid, rowinds (), rowvals (), numEntriesInRow);
    for (size_t i = 0; i < numEntriesInRow; i++) {
//      if (rowvals[i] >= 2*(mesh->num_dim-1))
//      {
      if (rowvals[i] >= mesh->num_dim-1)
      {
        mod_rowvals.push_back(one);
        mod_rowinds.push_back(rowinds[i]);
      }
    }
    modA->insertGlobalValues(gid, mod_rowinds, mod_rowvals);
  }
  
  modA->fillComplete();
  
  // get graph
  RCP<const crs_graph_type> G = modA->getCrsGraph();
  
  if(rank == 0) cout << "Wrtiting graph to file..." << endl;
  comm->barrier();
  // Write G to file
  std::ofstream file;
  char title[256];
  static const string path = "/Users/davidson/trilinosall/trilinosrepo/packages/zoltan2/test/driver/pamgen_mesh_data";
  sprintf(title, "/neighbors_%d", rank);
  file.open(path + string(title));
  
  ArrayView<const global_ordinal_type> ggid = G->getRowMap()->getNodeElementList();
  el_no = 0;
  for(global_ordinal_type gid : ggid)
  {
    size_t numEntriesInRow = G->getNumEntriesInGlobalRow(gid);
    Array<scalar_type>         rowvals (numEntriesInRow);
    Array<global_ordinal_type> rowinds (numEntriesInRow);
    G->getGlobalRowCopy(gid,  rowinds, numEntriesInRow);
    file << g_el_ids[el_no++]+1;
    for (size_t i = 0; i < numEntriesInRow; i++) {
      file << "\t" << rowinds[i] + 1;
    }
    file <<"\n";
  }
  
  file.close();
}

int main(int argc, char *argv[])
{
  ////////////////////////////////////////////////////////////
  // (0) Set up MPI environment
  ////////////////////////////////////////////////////////////
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  
  int rank = comm->getRank(); // get rank
  
  ////////////////////////////////////////////////////////////
  // (1) Get and read the input file
  // the input file defines tests to be run
  ////////////////////////////////////////////////////////////
  string inputFileName("driver.xml"); // assumes a default input file exists
  if(argc > 1)
    inputFileName = argv[1]; // user has provided an input file
  
  ////////////////////////////////////////////////////////////
  // (2) Get All Input Parameter Lists
  ////////////////////////////////////////////////////////////
  queue<ParameterList> problems, comparisons;
  getParameterLists(inputFileName,problems, comparisons, comm);
  
  ////////////////////////////////////////////////////////////
  // (3) Get Input Data Parameters
  ////////////////////////////////////////////////////////////
  
  // assumes that first block will always be
  // the input block
  const ParameterList inputParameters = problems.front();
  if(inputParameters.name() != "InputParameters")
  {
    if(rank == 0)
      cout << "InputParameters not defined" << endl;
    return 1;
  }
  
  // get the user input for all tests
  UserInputForTests uinput(inputParameters,comm,true,true);
  problems.pop();
  comm->barrier();
  
  ////////////////////////////////////////////////////////////
  // (4) Perform all tests
  ////////////////////////////////////////////////////////////
  // pamgen debugging
    writeMesh(uinput,comm);
    getConnectivityGraph(uinput, comm);
  
  RCP<ComparisonHelper> comparison_manager = rcp(new ComparisonHelper);
  while (!problems.empty()) {
    run(uinput, problems.front(), comparison_manager, comm);
    problems.pop();
  }
  
  ////////////////////////////////////////////////////////////
  // (5) Compare solutions
  ////////////////////////////////////////////////////////////
  while (!comparisons.empty()) {
    
    comparison_manager->CompareSolutions(comparisons.front().get<string>("A"),
                                        comparisons.front().get<string>("B"),
                                        comm);
    
    comparisons.pop();
  }
  
  return 0;
}

// helper functions

void readXML(const XMLObject &xml, const string &title)
{
  cout << "\nReading XML object " << title << " ...." << endl;
  xml.print(cout , 5);
}

void readPList(const ParameterList &plist, const string &title, bool doc, bool unused)
{  
  cout << "\nReading parameter list: " << title << " ...." << endl;
  plist.print(cout, ParameterList::PrintOptions().showDoc(doc).indent(3).showTypes(true));
  
  if(unused)
  {
    cout << "\nUnused fields: " << title << " ...." << endl;
    plist.unused(cout);
  }
}
