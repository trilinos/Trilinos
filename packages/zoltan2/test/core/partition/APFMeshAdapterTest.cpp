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

/*! \file APFMeshAdapterTest.cpp
    \brief An example of partitioning a SCOREC APF mesh

    \author Created by G. Diamond, K. Devine.

*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <Zoltan2_APFMeshAdapter.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_HyperGraphModel.hpp>

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Hashtable.hpp"

// SCOREC includes
#ifdef HAVE_ZOLTAN2_PARMA
#include <parma.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <gmi_mesh.h>
#endif

#include <set>

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::ArrayView;

// Computes and prints ghost metrics (takes in a Hyper graph model)
template <typename Adapter>
void PrintGhostMetrics(Zoltan2::HyperGraphModel<Adapter>& mdl) {
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::scalar_t       scalar_t;
  typedef typename Adapter::offset_t       offset_t;
  typedef Zoltan2::StridedData<lno_t, scalar_t>  input_t;

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> wgts;
  mdl.getEdgeList(Ids,wgts);
  ArrayView<bool> isOwner;
  mdl.getOwnedList(isOwner);
  size_t numOwned = mdl.getLocalNumOwnedVertices();
  ArrayView<const gno_t> pins;
  ArrayView<const offset_t> offsets;
  mdl.getPinList(pins,offsets,wgts);

  std::set<gno_t> gids;
  for (size_t i=0;i<mdl.getLocalNumVertices();i++) {
    if (isOwner[i])
      gids.insert(Ids[i]);
  }
  std::set<gno_t> ghosts;
  gno_t num_ghosts=0;
  for (size_t i=0;i<mdl.getLocalNumPins();i++) {
    gno_t pin = pins[i];
    if (gids.find(pin)==gids.end()) {
      num_ghosts++;
      if (ghosts.find(pin)==ghosts.end())
        ghosts.insert(pin);
    }
  }
  std::cout<< "[METRIC] " << PCU_Comm_Self() << " Total number of ghosts in the hypergraph: " << num_ghosts << "\n"
           << "[METRIC] " << PCU_Comm_Self() << " Number of unique ghosts: " << ghosts.size() << "\n";
  gno_t unique_ghosts =ghosts.size();
  gno_t owned_and_ghosts =unique_ghosts+numOwned;
  gno_t max_o_and_g,min_o_and_g;
  gno_t max_ghosts,max_u_ghosts;
  gno_t min_ghosts,min_u_ghosts;
  max_ghosts = min_ghosts = num_ghosts;
  max_u_ghosts = min_u_ghosts = unique_ghosts;
  max_o_and_g = min_o_and_g = owned_and_ghosts;
  double avg_ghosts,avg_u_ghosts,avg_o_and_g;
  PCU_Add_Ints(&num_ghosts,1);
  PCU_Add_Ints(&unique_ghosts,1);
  PCU_Add_Ints(&owned_and_ghosts,1);
  PCU_Max_Ints(&max_ghosts,1);
  PCU_Max_Ints(&max_u_ghosts,1);
  PCU_Max_Ints(&max_o_and_g,1);
  PCU_Min_Ints(&min_ghosts,1);
  PCU_Min_Ints(&min_u_ghosts,1);
  PCU_Min_Ints(&min_o_and_g,1);
  avg_ghosts = num_ghosts*1.0/PCU_Comm_Peers();
  avg_u_ghosts = unique_ghosts*1.0/PCU_Comm_Peers();
  avg_o_and_g = owned_and_ghosts*1.0/PCU_Comm_Peers();
  if (!PCU_Comm_Self())
    std::cout<< "[METRIC] Global ghosts in the hypergraph (tot max min avg imb): "
             << num_ghosts<<" "<<max_ghosts<<" "<<min_ghosts<<" "<<avg_ghosts<<" "
             <<max_ghosts/avg_ghosts << "\n"
             << "[METRIC] Global unique ghosts (tot max min avg imb): "
             << unique_ghosts<<" "<<max_u_ghosts<<" "<<min_u_ghosts<<" "<<avg_u_ghosts<<" "
             <<max_u_ghosts/avg_u_ghosts << "\n"
             << "[METRIC] Global owned and ghosts  (tot max min avg imb): "
             << owned_and_ghosts<<" "<<max_o_and_g<<" "<<min_o_and_g<<" "<<avg_o_and_g<<" "
             <<max_o_and_g/avg_o_and_g << "\n";

}

/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > CommT = Tpetra::getDefaultComm();

  int me = CommT->getRank();
  //int numProcs = CommT->getSize();

  if (me == 0){
  std::cout 
    << "====================================================================\n" 
    << "|                                                                  |\n" 
    << "|                  Example: Partition APF Mesh                     |\n" 
    << "|                                                                  |\n"
    << "|  Questions? Contact  Karen Devine      (kddevin@sandia.gov),     |\n"
    << "|                      Erik Boman        (egboman@sandia.gov),     |\n"
    << "|                      Siva Rajamanickam (srajama@sandia.gov).     |\n"
    << "|                                                                  |\n"
    << "|  Zoltan2's website:  http://trilinos.sandia.gov/packages/zoltan2 |\n"
    << "|  Trilinos website:   http://trilinos.sandia.gov                  |\n"
    << "|                                                                  |\n"
    << "====================================================================\n";
  }


#ifdef HAVE_MPI
  if (me == 0) {
    std::cout << "PARALLEL executable \n";
  }
#else
  if (me == 0) {
    std::cout << "SERIAL executable \n";
  }
#endif

  /***************************************************************************/
  /******************************* GET INPUTS ********************************/
  /***************************************************************************/

  // default values for command-line arguments
  std::string meshFileName("4/");
  std::string modelFileName("torus.dmg");
  std::string action("parma");
  std::string parma_method("VtxElm");
  std::string output_loc("");
  int nParts = CommT->getSize();
  double imbalance = 1.1;
  int layers=2;
  int ghost_metric=false;
  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("meshfile", &meshFileName,
                 "Mesh file with APF specifications (.smb file(s))");
  cmdp.setOption("modelfile", &modelFileName,
                 "Model file with APF specifications (.dmg file)");
  cmdp.setOption("action", &action,
                 "Method to use:  mj, scotch, zoltan_rcb, parma or color");
  cmdp.setOption("parma_method", &parma_method,
                 "Method to use: Vertex, Element, VtxElm, VtxEdgeElm, Ghost, Shape, or Centroid ");
  cmdp.setOption("nparts", &nParts,
                 "Number of parts to create");
  cmdp.setOption("imbalance", &imbalance,
                 "Target imbalance for the partitioning method");
  cmdp.setOption("output", &output_loc,
                 "Location of new partitioned apf mesh. Ex: 4/torus.smb");
  cmdp.setOption("layers", &layers,
                 "Number of layers for ghosting");
  cmdp.setOption("ghost_metric", &ghost_metric,
                 "0 does not compute ghost metric otherwise compute both before and after");
  cmdp.parse(narg, arg);

  /***************************************************************************/
  /********************** GET CELL TOPOLOGY **********************************/
  /***************************************************************************/

  // Get dimensions
  //int dim = 3;

  /***************************************************************************/
  /***************************** GENERATE MESH *******************************/
  /***************************************************************************/

#ifdef HAVE_ZOLTAN2_PARMA

  if (me == 0) std::cout << "Generating mesh ... \n\n";

  //Setup for SCOREC
  PCU_Comm_Init();

  // Generate mesh with MDS
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(modelFileName.c_str(),meshFileName.c_str());
  apf::verify(m);

  //Data for APF MeshAdapter
  std::string primary="region";
  std::string adjacency="face";
  if (m->getDimension()==2) {
    primary="face";
    adjacency="edge";
  }
  bool needSecondAdj=false;

  // Set parameters for partitioning
  if (me == 0) std::cout << "Creating parameter list ... \n\n";

  Teuchos::ParameterList params("test params");
  params.set("timer_output_stream" , "std::cout");

  bool do_partitioning = false;
  if (action == "mj") {
    do_partitioning = true;
    params.set("debug_level", "basic_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("num_global_parts", nParts);
    params.set("algorithm", "multijagged");
    params.set("rectilinear", "yes");
  }
  else if (action == "scotch") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("objects_to_partition","mesh_elements");
    params.set("algorithm", "scotch");
    needSecondAdj=true;
  }
  else if (action == "zoltan_rcb") {
    do_partitioning = true;
    params.set("debug_level", "verbose_detailed_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "zoltan");
  }
  else if (action == "parma") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("algorithm", "parma");
    Teuchos::ParameterList &pparams = params.sublist("parma_parameters",false);
    pparams.set("parma_method",parma_method);
    pparams.set("step_size",1.1);
    if (parma_method=="Ghost") {
      pparams.set("ghost_layers",layers);
      pparams.set("ghost_bridge",m->getDimension()-1);
    }
    adjacency="vertex";
  }
  else if (action=="zoltan_hg") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("algorithm", "zoltan");
    params.set("num_global_parts", nParts);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","HYPERGRAPH");
    zparams.set("LB_APPROACH","PARTITION");
    //params.set("compute_metrics","yes");
    adjacency="vertex";
  }
  else if (action=="hg_ghost") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("algorithm", "zoltan");
    params.set("num_global_parts", nParts);
    params.set("hypergraph_model_type","ghosting");
    params.set("ghost_layers",layers);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","HYPERGRAPH");
    zparams.set("LB_APPROACH","PARTITION");
    zparams.set("PHG_EDGE_SIZE_THRESHOLD", "1.0");
    primary="vertex";
    adjacency="edge";
    needSecondAdj=true;
  }
  else if (action == "color") {
    params.set("debug_level", "verbose_detailed_status");
    params.set("debug_output_file", "kdd");
    params.set("debug_procs", "all");
  }
  Parma_PrintPtnStats(m,"before");

  // Creating mesh adapter
  if (me == 0) std::cout << "Creating mesh adapter ... \n\n";
  typedef Zoltan2::APFMeshAdapter<apf::Mesh2*> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
  typedef Zoltan2::MeshAdapter<apf::Mesh2*> baseMeshAdapter_t;

  double time_1=PCU_Time();
  inputAdapter_t *ia =
    new inputAdapter_t(*CommT, m,primary,adjacency,needSecondAdj);
  double time_2=PCU_Time();


  inputAdapter_t::scalar_t* arr =
    new inputAdapter_t::scalar_t[ia->getLocalNumOf(ia->getPrimaryEntityType())];
  for (size_t i=0;i<ia->getLocalNumOf(ia->getPrimaryEntityType());i++) {
    arr[i]=PCU_Comm_Self()+1;
  }

  const inputAdapter_t::scalar_t* weights=arr;
  ia->setWeights(ia->getPrimaryEntityType(),weights,1);


  if (ghost_metric) {
    const baseMeshAdapter_t *base_ia = dynamic_cast<const baseMeshAdapter_t*>(ia);
    Zoltan2::modelFlag_t graphFlags_;
    RCP<Zoltan2::Environment> env;
    try{
      env = rcp(new Zoltan2::Environment(params, Tpetra::getDefaultComm()));
    }
    Z2_FORWARD_EXCEPTIONS

    RCP<const Zoltan2::Environment> envConst = Teuchos::rcp_const_cast<const Zoltan2::Environment>(env);

    RCP<const baseMeshAdapter_t> baseInputAdapter_(base_ia,false);
    Zoltan2::HyperGraphModel<inputAdapter_t> model(baseInputAdapter_,envConst,CommT,
                                                 graphFlags_,Zoltan2::HYPEREDGE_CENTRIC);
    PrintGhostMetrics(model);
  }

  // create Partitioning problem
  double time_3 = PCU_Time();
  if (do_partitioning) {
    if (me == 0) std::cout << "Creating partitioning problem ... \n\n";

    Zoltan2::PartitioningProblem<inputAdapter_t> problem(ia, &params, CommT);

    // call the partitioner
    if (me == 0) std::cout << "Calling the partitioner ... \n\n";

    problem.solve();



    if (me==0) std::cout << "Applying Solution to Mesh\n\n";
    apf::Mesh2** new_mesh = &m;
    ia->applyPartitioningSolution(m,new_mesh,problem.getSolution());

    // create metric object
    RCP<quality_t> metricObject =
      rcp(new quality_t(ia, &params, CommT, &problem.getSolution()));

    if (!me) {
      metricObject->printMetrics(std::cout);
    }
  }
  else {
    if (me == 0) std::cout << "Creating coloring problem ... \n\n";

    Zoltan2::ColoringProblem<inputAdapter_t> problem(ia, &params);

    // call the partitioner
    if (me == 0) std::cout << "Calling the coloring algorithm ... \n\n";

    problem.solve();

    problem.printTimers();


  }

  double time_4=PCU_Time();

  //Destroy the adapter
  ia->destroy();
  delete [] arr;
  //Parma_PrintPtnStats(m,"after");

  if (ghost_metric) {
    inputAdapter_t ia2(*CommT, m,primary,adjacency,true);
    const baseMeshAdapter_t *base_ia = dynamic_cast<const baseMeshAdapter_t*>(&ia2);

    Zoltan2::modelFlag_t graphFlags_;
    RCP<Zoltan2::Environment> env;
    try{
      env = rcp(new Zoltan2::Environment(params, Tpetra::getDefaultComm()));
    }
    Z2_FORWARD_EXCEPTIONS
    RCP<const Zoltan2::Environment> envConst = Teuchos::rcp_const_cast<const Zoltan2::Environment>(env);
    RCP<const baseMeshAdapter_t> baseInputAdapter_(base_ia,false);
    Zoltan2::HyperGraphModel<inputAdapter_t> model(baseInputAdapter_, envConst, CommT,
                                                   graphFlags_,Zoltan2::HYPEREDGE_CENTRIC);

    PrintGhostMetrics(model);
    ia2.destroy();
  }

  if (output_loc!="") {
    m->writeNative(output_loc.c_str());
  }

  // delete mesh
  if (me == 0) std::cout << "Deleting the mesh ... \n\n";
  time_4-=time_3;
  time_2-=time_1;
  PCU_Max_Doubles(&time_2,1);
  PCU_Max_Doubles(&time_4,1);
  if (!me) {
    std::cout<<"\nConstruction time: "<<time_2<<"\n"
             <<"Problem time: " << time_4<<"\n\n";
  }
  //Delete the APF Mesh
  m->destroyNative();
  apf::destroyMesh(m);
  //End communications
  PCU_Comm_Free();
#endif

  if (me == 0)
    std::cout << "PASS" << std::endl;

  return 0;

}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
