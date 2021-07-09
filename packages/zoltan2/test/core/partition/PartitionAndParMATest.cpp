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

/*! \file PartitionAndParMA.cpp
    \brief Runs a partitioning algorithm followed by ParMA for cleanup on a SCOREC mesh

    \author Created by G. Diamond, K. Devine.

*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <Zoltan2_APFMeshAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

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

using Teuchos::ParameterList;
using Teuchos::RCP;

#ifdef HAVE_ZOLTAN2_PARMA

void runTest(RCP<const Teuchos::Comm<int> >& CommT, apf::Mesh2* m,std::string action,
             std::string parma_method,int nParts, double imbalance, std::string output_title );

#endif
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
    << "|  Pamgen's website:   http://trilinos.sandia.gov/packages/pamgen  |\n"
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
  std::string action("zoltan_hg");
  std::string parma_method("VtxElm");
  std::string output_loc("");
  int nParts = CommT->getSize();
  double imbalance=1.1;

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("meshfile", &meshFileName,
                 "Mesh file with APF specifications (.smb file(s))");
  cmdp.setOption("modelfile", &modelFileName,
                 "Model file with APF specifications (.dmg file)");
  cmdp.setOption("action", &action,
                 "Method to use:  mj, scotch, zoltan_rcb, parma or color");
  cmdp.setOption("parma_method", &parma_method,
                 "Method to use: Vertex, Edge, Element, VtxElm, VtxEdgeElm, ElmLtVtx, Ghost, or Shape ");
  cmdp.setOption("nparts", &nParts,
                 "Number of parts to create");
  cmdp.setOption("imbalance", &imbalance,
                 "Target Imbalance for first partitioner");
  cmdp.setOption("output", &output_loc,
                 "Location of new partitioned apf mesh. Ex: 4/torus.smb");
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

  runTest(CommT,m,action,parma_method,nParts,imbalance,"partition");

  runTest(CommT,m,"parma",parma_method,nParts,imbalance,"parma");




  if (output_loc!="") {
    m->writeNative(output_loc.c_str());
  }

  // delete mesh
  if (me == 0) std::cout << "Deleting the mesh ... \n\n";

  //Delete APF Mesh;
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

#ifdef HAVE_ZOLTAN2_PARMA

void runTest(RCP<const Teuchos::Comm<int> >& CommT, apf::Mesh2* m,std::string action,
             std::string parma_method,int nParts, double imbalance,std::string output_title) {
  //Get rank
  int me = CommT->getRank();

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

  if (action == "mj") {
    params.set("debug_level", "basic_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("num_global_parts", nParts);
    params.set("algorithm", "multijagged");
    params.set("rectilinear", "yes");
  }
  else if (action == "scotch") {
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "scotch");
    needSecondAdj=true;
  }
  else if (action == "zoltan_rcb") {
    params.set("debug_level", "verbose_detailed_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "zoltan");
  }
  else if (action == "parma") {
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("algorithm", "parma");
    Teuchos::ParameterList &pparams = params.sublist("parma_parameters",false);
    pparams.set("parma_method",parma_method);
    pparams.set("step_size",1.1);
    if (parma_method=="Ghost") {
      pparams.set("ghost_layers",3);
      pparams.set("ghost_bridge",m->getDimension()-1);
    }
    adjacency="vertex";
  }
  else if (action=="zoltan_hg") {
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("algorithm", "zoltan");
    params.set("num_global_parts", nParts);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","HYPERGRAPH");
    //params.set("compute_metrics","yes");
    adjacency="vertex";
  }

  //Print the stats of original mesh
  Parma_PrintPtnStats(m,output_title+"_before");

  // Creating mesh adapter
  if (me == 0) std::cout << "Creating mesh adapter ... \n\n";
  typedef Zoltan2::APFMeshAdapter<apf::Mesh2*> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
  typedef Zoltan2::MeshAdapter<apf::Mesh2*> baseMeshAdapter_t;

  double time_1 = PCU_Time();
  inputAdapter_t *ia =
    new inputAdapter_t(*CommT, m,primary,adjacency,needSecondAdj);
  double time_2 = PCU_Time();


  // create Partitioning problem
  if (me == 0) std::cout << "Creating partitioning problem ... \n\n";
  double time_3=PCU_Time();
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(ia, &params, CommT);

  // call the partitioner
  if (me == 0) std::cout << "Calling the partitioner ... \n\n";

  problem.solve();


  //apply the partitioning solution to the mesh
  if (me==0) std::cout << "Applying Solution to Mesh\n\n";
  apf::Mesh2** new_mesh = &m;
  ia->applyPartitioningSolution(m,new_mesh,problem.getSolution());
  new_mesh=NULL;
  double time_4=PCU_Time();

  //create metric object
  RCP<quality_t> metricObject =
    rcp(new quality_t(ia, &params, CommT, &problem.getSolution()));

  if (!me) {
    metricObject->printMetrics(std::cout);
  }

  //Print the stats after partitioning
  Parma_PrintPtnStats(m,output_title+"_after");
  ia->destroy();

  time_4-=time_3;
  time_2-=time_1;
  PCU_Max_Doubles(&time_2,1);
  PCU_Max_Doubles(&time_4,1);
  if (!me) {
    std::cout<<"\n"<<output_title<<"Construction time: "<<time_2<<"\n"
             <<output_title<<"Problem time: " << time_4<<"\n\n";
  }

}

#endif
