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

/*! \file RPIMeshAdapterTest.cpp
    \brief An example of partitioning a SCOREC mesh with RCB.

    \author Created by G. Diamond, K. Devine.

*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <Zoltan2_RPIMeshAdapter.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>

//Tpetra includes
#include "Tpetra_DefaultPlatform.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
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

using namespace std;
using Teuchos::ParameterList;
using Teuchos::RCP;

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
//Tpetra typedefs
typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;



/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  Teuchos::GlobalMPISession mpiSession(&narg, &arg,0);
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > CommT = platform.getComm();

  int me = CommT->getRank();
  //int numProcs = CommT->getSize();

  if (me == 0){
  cout 
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
    cout << "PARALLEL executable \n";
  }
#else
  if (me == 0) {
    cout << "SERIAL executable \n";
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

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("meshfile", &meshFileName,
                 "Mesh file with APF specifications (.smb file(s))");
  cmdp.setOption("modelfile", &modelFileName,
		 "Model file with APF specifications (.dmg file)");
  cmdp.setOption("action", &action,
                 "Method to use:  mj, scotch, zoltan_rcb, parma or color");
  cmdp.setOption("parma_method", &parma_method,
                 "Method to use: Vertex, Element, VtxElm, VtxEdgeElm, Ghost, or Shape ");
  cmdp.setOption("nparts", &nParts,
                 "Number of parts to create");
  cmdp.setOption("imbalance", &imbalance,
                 "Target imbalance for the partitioning method");
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

  if (me == 0) cout << "Generating mesh ... \n\n";

  //Setup for SCOREC
  PCU_Comm_Init();
  
  // Generate mesh with MDS
  double time_1=PCU_Time();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(modelFileName.c_str(),meshFileName.c_str());
  apf::verify(m);
  // Creating mesh adapter
  if (me == 0) cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::RPIMeshAdapter<apf::Mesh2*> inputAdapter_t;

  inputAdapter_t ia(*CommT, m);
  
  double time_2=PCU_Time();
  // Set parameters for partitioning
  if (me == 0) cout << "Creating parameter list ... \n\n";

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
      pparams.set("ghost_layers",3);
      pparams.set("ghost_bridge",m->getDimension()-1);
    }
    params.set("compute_metrics","yes");

  }
  else if (action=="zoltan_hg") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", imbalance);
    params.set("algorithm", "zoltan");
    params.set("num_global_parts", nParts);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","HYPERGRAPH");
    zparams.set("LB_APPROACH","REPARTITION");
    //params.set("compute_metrics","yes");

  }
  else if (action == "color") {
    params.set("debug_level", "verbose_detailed_status");
    params.set("debug_output_file", "kdd");
    params.set("debug_procs", "all");
  }
  Parma_PrintPtnStats(m,"before");
  // create Partitioning problem
  double time_3 = PCU_Time();
  if (do_partitioning) {
    if (me == 0) cout << "Creating partitioning problem ... \n\n";

    Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params, CommT);

    // call the partitioner
    if (me == 0) cout << "Calling the partitioner ... \n\n";

    problem.solve();



    if (me==0) cout << "Applying Solution to Mesh\n\n";
    apf::Mesh2** new_mesh = &m;
    ia.applyPartitioningSolution(m,new_mesh,problem.getSolution());
    
    

    if (!me) problem.printMetrics(cout);
  }
  else {
    if (me == 0) cout << "Creating coloring problem ... \n\n";

    Zoltan2::ColoringProblem<inputAdapter_t> problem(&ia, &params);

    // call the partitioner
    if (me == 0) cout << "Calling the coloring algorithm ... \n\n";

    problem.solve();

    problem.printTimers();


  }
  double time_4=PCU_Time();
  //if (!me)
  Parma_PrintPtnStats(m,"after");
  if (output_loc!="") {
    m->writeNative(output_loc.c_str());
  }

  // delete mesh
  if (me == 0) cout << "Deleting the mesh ... \n\n";
  time_4-=time_3;
  time_2-=time_1;
  PCU_Max_Doubles(&time_2,1);
  PCU_Max_Doubles(&time_4,1);
  if (!me) {
    std::cout<<"\nConstruction time: "<<time_2<<"\n"
	     <<"Problem time: " << time_4<<"\n\n";
  }
  //Delete_APF_Mesh();
  ia.destroy();
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
