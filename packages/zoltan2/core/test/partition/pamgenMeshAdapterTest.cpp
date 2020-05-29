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

/*! \file pamgenMeshAdapterTest.cpp
    \brief An example of partitioning pamgen coordinates with RCB.

    \author Created by V. Leung, K. Devine.

*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <Zoltan2_PamgenMeshAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Pamgen includes
#include "create_inline_mesh.h"

using Teuchos::ParameterList;
using Teuchos::ArrayRCP;

typedef Zoltan2::BasicUserTypes<> basic_user_t;

/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > CommT = Tpetra::getDefaultComm();

  int me = CommT->getRank();
  int numProcs = CommT->getSize();

  if (me == 0){
  std::cout
    << "====================================================================\n"
    << "|                                                                  |\n"
    << "|          Example: Partition Pamgen Hexahedral Mesh               |\n"
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
  /*************************** GET XML INPUTS ********************************/
  /***************************************************************************/

  // default values for command-line arguments
  std::string xmlMeshInFileName("Poisson.xml");
  std::string action("mj");
  int nParts = CommT->getSize();

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("xmlfile", &xmlMeshInFileName,
                 "XML file with PamGen specifications");
  cmdp.setOption("action", &action,
                 "Method to use:  mj, scotch, zoltan_rcb, zoltan_hg, hg_ghost, "
                 "parma or color");
  cmdp.setOption("nparts", &nParts,
                 "Number of parts to create");
  cmdp.parse(narg, arg);

  // Read xml file into parameter list
  ParameterList inputMeshList;

  if(xmlMeshInFileName.length()) {
    if (me == 0) {
      std::cout << "\nReading parameter list from the XML file \""
                <<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile(xmlMeshInFileName,
                                         Teuchos::inoutArg(inputMeshList));
    if (me == 0) {
      inputMeshList.print(std::cout,2,true,true);
      std::cout << "\n";
    }
  }
  else {
    std::cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
    return 5;
  }

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,
                                                             "meshInput");

  /***************************************************************************/
  /********************** GET CELL TOPOLOGY **********************************/
  /***************************************************************************/

  // Get dimensions
  int dim = 3;

  /***************************************************************************/
  /***************************** GENERATE MESH *******************************/
  /***************************************************************************/

  if (me == 0) std::cout << "Generating mesh ... \n\n";

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, me, numProcs, maxInt);

  // Creating mesh adapter
  if (me == 0) std::cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;

  inputAdapter_t *ia = new inputAdapter_t(*CommT, "region");
  ia->print(me);

  // Set parameters for partitioning
  if (me == 0) std::cout << "Creating parameter list ... \n\n";

  Teuchos::ParameterList params("test params");
  params.set("timer_output_stream" , "std::cout");

  bool do_partitioning = false;
  if (action == "mj") {
    do_partitioning = true;
    params.set("debug_level", "basic_status");
    params.set("imbalance_tolerance", 1.1);
    params.set("num_global_parts", nParts);
    params.set("algorithm", "multijagged");
    params.set("rectilinear", true); // bool parameter
  }
  else if (action == "scotch") {
    do_partitioning = true;
    params.set("debug_level", "verbose_detailed_status");
    params.set("imbalance_tolerance", 1.1);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "scotch");
  }
  else if (action == "zoltan_rcb") {
    do_partitioning = true;
    params.set("debug_level", "verbose_detailed_status");
    params.set("imbalance_tolerance", 1.1);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "zoltan");
  }
  else if (action == "zoltan_hg") {
    do_partitioning = true;
    params.set("debug_level", "verbose_detailed_status");
    params.set("imbalance_tolerance", 1.1);
    params.set("num_global_parts", nParts);
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "zoltan");
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","phg");
    zparams.set("FINAL_OUTPUT", "1");
  }
  else if (action=="hg_ghost") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", 1.1);
    params.set("algorithm", "zoltan");
    params.set("num_global_parts", nParts);
    params.set("hypergraph_model_type","ghosting");
    params.set("ghost_layers",2);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","HYPERGRAPH");
    zparams.set("LB_APPROACH","PARTITION");
    zparams.set("PHG_EDGE_SIZE_THRESHOLD", "1.0");
  }

  else if (action == "parma") {
    do_partitioning = true;
    params.set("debug_level", "basic_status");
    params.set("imbalance_tolerance", 1.05);
    params.set("algorithm", "parma");
    Teuchos::ParameterList &pparams = params.sublist("parma_parameters",false);
    pparams.set("parma_method","VtxElm");
  }
  else if (action=="zoltan_hg") {
    do_partitioning = true;
    params.set("debug_level", "no_status");
    params.set("imbalance_tolerance", 1.1);
    params.set("algorithm", "zoltan");
    params.set("num_global_parts", nParts);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("LB_METHOD","HYPERGRAPH");

  }

  else if (action == "color") {
    params.set("debug_level", "verbose_detailed_status");
    params.set("debug_output_file", "kdd");
    params.set("debug_procs", "all");
  }

  if(me == 0) std::cout << "Action: " << action << std::endl;
  // create Partitioning problem
  if (do_partitioning) {
    if (me == 0) std::cout << "Creating partitioning problem ... \n\n";

    Zoltan2::PartitioningProblem<inputAdapter_t> problem(ia, &params, CommT);

    // call the partitioner
    if (me == 0) std::cout << "Calling the partitioner ... \n\n";

    problem.solve();

    // create metric object

    Teuchos::RCP<quality_t> metricObject = 
      rcp(new quality_t(ia, &params, CommT, &problem.getSolution()));

    if (!me) {
      metricObject->printMetrics(std::cout);
    }
  }
  else {
    if (me == 0) std::cout << "Creating coloring problem ... \n\n";

    Zoltan2::ColoringProblem<inputAdapter_t> problem(ia, &params, CommT);

    // call the partitioner
    if (me == 0) std::cout << "Calling the coloring algorithm ... \n\n";

    problem.solve();

    problem.printTimers();
  }

  // delete mesh
  if (me == 0) std::cout << "Deleting the mesh ... \n\n";

  Delete_Pamgen_Mesh();

  if (me == 0)
    std::cout << "PASS" << std::endl;

  Kokkos::finalize(); // Test uses directory with kokkos

  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
