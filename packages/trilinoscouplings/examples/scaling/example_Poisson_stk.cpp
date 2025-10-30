// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   example_Poisson_stk.cpp
    \brief  Example solution of a Poisson equation on a hexahedral or
            tetrahedral mesh using nodal (Hgrad) elements.

            This example requires a hexahedral or tetrahedral mesh in Exodus
            format with a nodeset containing boundary nodes. STK is used to
            read the mesh and populate a mesh database, Intrepid is used to
            build the stiffness matrix and right-hand side, and ML is used
            to solve the resulting linear system.

    \verbatim

     Poisson system:

            div A grad u = f in Omega
                       u = g on  Gamma

       where
             A is a symmetric, positive definite material tensor
             f is a given source term


     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson C. Siefert.

    \remark Usage:
    \code   ./example_Poisson_stk --help  \endcode

    \remark Example requires a hexahedral or tetrahedral mesh in Exodus format
*/
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/


/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <unistd.h>
#include <vector>
#include <map>

#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_TpetraIntrepidPoissonExample.hpp"
#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"


// Teuchos includes
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_TimeMonitor.hpp>

// Intrepid2 includes
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_FECrsGraph.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Assembly_Helpers.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// MueLu includes
#include "MueLu.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

// Sacado includes
#include "Sacado.hpp"

// STK includes
#include "Ionit_Initializer.h"
#include "Ioss_SubSystem.h"

#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

#include "stk_util/parallel/Parallel.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Field.hpp"


/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
// Tpetra typedefs
typedef Tpetra::Map<> Map;

typedef Map::node_type::memory_space memory_space;
typedef Map::node_type::device_type device_type;

typedef shards::CellTopology ShardsCellTopology;
typedef Intrepid2::FunctionSpaceTools<device_type> Intrepid2FSTools;
typedef Intrepid2::CellTools<device_type> Intrepid2CTools;
typedef Intrepid2::ScalarView<double, memory_space> Intrepid2ScalarView;

using Teuchos::TimeMonitor;

// Number of dimensions
int spaceDim;

/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \return Value of the exact solution at (x,y,z)
 */
template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar&  x, const Scalar&  y, const Scalar&  z);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y,z)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point
 */
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z);


/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
           and material tensor.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y,z)
 */
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z);


/** \brief Computation of the material tensor at array of points in physical space.

    \param worksetMaterialValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
                                                with the values of the material tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        worksetMaterialValues,
                            const ArrayIn &   evaluationPoints);


/** \brief Computation of the source term at array of points in physical space.

    \param sourceTermValues           [out]     Rank-2 (C,P) array with the values of the source term
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints);

/** \brief Computation of the exact solution at array of points in physical space.

    \param exactSolutionValues        [out]     Rank-2 (C,P) array with the values of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints);


/** \brief Computation of the gradient of the exact solution at array of points in physical space.

    \param exactSolutionGradValues    [out]     Rank-3 (C,P,D) array with the values of the gradient of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints);


/**********************************************************************************/
/**************** FUNCTION DECLARATION FOR ML PRECONDITIONER *********************/
/**********************************************************************************/
using Tpetra_CrsMatrix = Tpetra::CrsMatrix<>;
using Tpetra_MultiVector = Tpetra::MultiVector<>;
using Tpetra_Vector = Tpetra::Vector<>;
int TestMultiLevelPreconditioner(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Teuchos::RCP<Tpetra_CrsMatrix>   & A,
                                 Teuchos::RCP<Tpetra_Vector> & xexact,
                                 Teuchos::RCP<Tpetra_MultiVector> & b,
                                 Teuchos::RCP<Tpetra_MultiVector> & uh,
				 Teuchos::RCP<Tpetra_MultiVector> & coords,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol);

/**********************************************************************************/
/************* FUNCTION DECLARATIONS FOR SIMPLE BASIS FACTORY *********************/
/**********************************************************************************/

/** \brief  Simple factory that chooses basis function based on cell topology.

    \param  cellTopology  [in]    Shards cell topology
    \param  order         [in]    basis function order, currently unused
    \param  basis         [out]   pointer to Intrepid basis

    \return Intrepid basis
 */

int getDimension(const ShardsCellTopology & cellTopology);

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

// Do the actual reading of the mesh database and
// creation and population of the MetaData and BulkData.
void mesh_read_write(const std::string &type,
                     const std::string &working_directory,
                     const std::string &filename,
                     stk::io::StkMeshIoBroker &broker,
                     int db_integer_size,
                     stk::io::HeartbeatType hb_type)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  std::string absoluteFileName = working_directory + "/" + filename;
  if (myrank == 0) {
    std::cout<<"Reading mesh: "<<absoluteFileName<<std::endl;
  }
  size_t input_index = broker.add_mesh_database(absoluteFileName, type, stk::io::READ_MESH);
  broker.set_active_mesh(input_index);
  // creates metadata
  broker.create_input_mesh();

  // commits the meta data
  broker.populate_bulk_data();

} //mesh_read_write


/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
template<class crs_matrix_type, class multivector_type1, class multivector_type2, class solution_type>
void Apply_Dirichlet_BCs(std::vector<int> &BCNodes, crs_matrix_type & A, multivector_type1 & x, multivector_type2 & b, solution_type & solution_values) {
  using SC = typename multivector_type1::scalar_type;
  using LO = typename multivector_type1::local_ordinal_type;
  int N=(int)BCNodes.size();
  Teuchos::ArrayRCP<SC> xdata = x.getDataNonConst(0);
  Teuchos::ArrayRCP<SC> bdata = b.getDataNonConst(0);
  if (b.getMap()->getComm()->getRank() == 0)
    std::cout<<"Apply in Dirichlet BCs to "<<BCNodes.size() << " nodes"<<std::endl;

  Tpetra::beginModify(A,b);

  for(int i=0; i<N; i++) {
    LO lrid = BCNodes[i];

    xdata[lrid]=bdata[lrid] = solution_values[i];

    size_t numEntriesInRow = A.getNumEntriesInLocalRow(lrid);
    typename crs_matrix_type::nonconst_local_inds_host_view_type cols("cols", numEntriesInRow);
    typename crs_matrix_type::nonconst_values_host_view_type vals("vals", numEntriesInRow);
    A.getLocalRowCopy(lrid, cols, vals, numEntriesInRow);

    for(int j=0; j<vals.extent_int(0); j++)
      vals(j) = (cols(j) == lrid) ? 1.0 : 0.0;

    A.replaceLocalValues(lrid, cols, vals);
  }

  Tpetra::endModify(A,b);

}



/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/


int main_(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using entity_type = stk::mesh::Entity;
  using Tpetra_Map = Tpetra::Map<>;
  using Tpetra_FEMultiVector = Tpetra::FEMultiVector<>;
  using Tpetra_FECrsMatrix = Tpetra::FECrsMatrix<>;
  using Tpetra_FECrsGraph = Tpetra::FECrsGraph<>;
  using SC = Tpetra_Vector::scalar_type;
  using LO = Tpetra_Vector::local_ordinal_type;
  using GO = Tpetra_Vector::global_ordinal_type;
  //using NO = Tpetra_Vector::node_type;


  const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;
  const stk::mesh::EntityRank ELEMENT_RANK = stk::topology::ELEMENT_RANK;

  RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
  RCP<const Teuchos::MpiComm<int> > MpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Comm);
  int numRanks = Comm->getSize();
  int MyPID = Comm->getRank();


  Teuchos::CommandLineProcessor clp(false);

  std::string optMeshFile = "unit_cube_10int_hex.exo";
  clp.setOption("mesh",  &optMeshFile, "Exodus hexahedra, tetrahedral, quadrilateral or triangle mesh file with nodeset define for boundary");
  std::string optXmlFile  = "";
  clp.setOption("xml",   &optXmlFile,  "xml file containing ML solver options");
  bool optPrintLocalStats = false; clp.setOption("localstats", "nolocalstats", &optPrintLocalStats, "print per-process statistics");
  // If matrixFilename is nonempty, dump the matrix to that file
  // in MatrixMarket format.
  std::string matrixFilename;
  clp.setOption ("matrixFilename", &matrixFilename, "If nonempty, dump the "
		  "generated matrix to that file in MatrixMarket format.");

  // If rhsFilename is nonempty, dump the rhs to that file
  // in MatrixMarket format.
  std::string rhsFilename;
  clp.setOption ("rhsFilename", &rhsFilename, "If nonempty, dump the "
		  "generated rhs to that file in MatrixMarket format.");

  // If rhsFilename is nonempty, dump the rhs to that file
  // in MatrixMarket format.
  std::string initialGuessFilename;
  clp.setOption ("initialGuessFilename", &initialGuessFilename, "If nonempty, dump the "
		  "generated initial guess to that file in MatrixMarket format.");

  // If coordsFilename is nonempty, dump the coords to that file
  // in MatrixMarket format.
  std::string coordsFilename;
  clp.setOption ("coordsFilename", &coordsFilename, "If nonempty, dump the "
		  "generated coordinates to that file in MatrixMarket format.");



  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      break;
  }

  if (MyPID == 0) {
    std::cout
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|              Example: Solve Poisson Equation                                |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid2's website: http://trilinos.sandia.gov/packages/intrepid2         |\n" \
    << "|  STK's website:       http://trilinos.github.io/stk.html                    |\n" \
    << "|  ML's website:        http://trilinos.sandia.gov/packages/ml                |\n" \
    << "|  Trilinos website:    http://trilinos.sandia.gov                            |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  }


  /**********************************************************************************/
  /********************************** GET XML INPUTS ********************************/
  /**********************************************************************************/

  // get xml file from command line if provided, otherwise use default
  std::string  xmlSolverInFileName(optXmlFile);

  // Read xml file into parameter list
  Teuchos::ParameterList inputSolverList;

  if(xmlSolverInFileName.length()) {
    if (MyPID == 0)
      std::cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n" << std::endl;
    Teuchos::updateParametersFromXmlFile (xmlSolverInFileName, Teuchos::ptr (&inputSolverList));
  }
  else if (MyPID == 0)
    std::cout << "Using default solver values ..." << std::endl;


  /**********************************************************************************/
  /*********************************** READ MESH ************************************/
  /**********************************************************************************/
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) Read Mesh")));

  stk::io::StkMeshIoBroker broker(*MpiComm->getRawMpiComm());
  broker.property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 180));
  broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "HSFC"));
  broker.property_add(Ioss::Property("COMPOSE_RESULTS", false));  //Note!  true results in an error in Seacas

  std::string type = "exodusii";
  char buf[1024];
  if (getcwd(buf,sizeof(buf)) == NULL)
    throw(std::runtime_error("Could not get current working directory"));
  std::string working_directory(buf);
  std::string filename(optMeshFile);
  int db_integer_size = 4;
  stk::io::HeartbeatType hb_type = stk::io::NONE;
  stk::initialize(&argc, &argv);
  mesh_read_write(type, working_directory, filename, broker, db_integer_size, hb_type);

  stk::mesh::BulkData &bulkData = broker.bulk_data();
  stk::mesh::MetaData &metaData = broker.meta_data();

  // Count number of local nodes, record GIDs for Tpetra map.
  Teuchos::Array<GO> ownedGIDs, ownedPlusSharedGIDs;
  stk::mesh::Selector locallyOwnedSelector = metaData.locally_owned_part();
  int numLocalNodes = stk::mesh::count_entities(bulkData, NODE_RANK, locallyOwnedSelector);

  stk::mesh::for_each_entity_run(bulkData, NODE_RANK, locallyOwnedSelector,
    [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity node)
    {
        ownedGIDs.push_back(mesh.identifier(node)-1);
    });
  ownedPlusSharedGIDs.assign(ownedGIDs.begin(), ownedGIDs.end());

  // Now record the shared-but-not-owned nodal GIDs
  {
    stk::mesh::Selector globallySharedSelector = metaData.globally_shared_part();
    globallySharedSelector &= !locallyOwnedSelector; //not owned
    stk::mesh::for_each_entity_run(bulkData, NODE_RANK, globallySharedSelector,
      [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity node)
      {
        ownedPlusSharedGIDs.push_back(mesh.identifier(node)-1);
      });
  }

  // Count # of local elements
  int numLocalElems = stk::mesh::count_entities(bulkData, ELEMENT_RANK, locallyOwnedSelector);

  if (optPrintLocalStats) {
    for (int i=0; i<numRanks; ++i) {
      if (MyPID == i) {
        std::cout << "(" << MyPID << ")    Number of local Elements: " << numLocalElems << std::endl
                  << "(" << MyPID << ")    Number of local Nodes   : " << numLocalNodes << std::endl << std::endl;
      }
      Comm->barrier();
    }
  }

  GO numGlobalNodes=0, numGlobalElements=0;
  Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, (GO)numLocalNodes, Teuchos::outArg(numGlobalNodes));
  Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, (GO)numLocalElems, Teuchos::outArg(numGlobalElements));
  if (MyPID == 0) {
    std::cout << "       Number of global Nodes   : " << numGlobalNodes << std::endl;
    std::cout << "       Number of global Elements: " << numGlobalElements << std::endl;
  }
  tm = Teuchos::null;


  /**********************************************************************************/
  /********************* BUILD MAPS/GRAPHS FOR GLOBAL SOLUTION **********************/
  /**********************************************************************************/
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("2) Build Maps/Graph")));

  RCP<Tpetra_Map> globalMapG = Teuchos::rcp(new Tpetra_Map( (Tpetra::global_size_t)numGlobalNodes,ownedGIDs(), (GO)0, Comm));
  RCP<Tpetra_Map> ownedPlusSharedMapG = Teuchos::rcp(new Tpetra_Map( (Tpetra::global_size_t)numGlobalNodes,ownedPlusSharedGIDs(), (GO)0, Comm));
  Kokkos::DualView<size_t*> nnzPerRowUpperBound("nnzbound",ownedPlusSharedGIDs.size());
  nnzPerRowUpperBound.template modify<typename Kokkos::DualView<size_t*>::host_mirror_space>();
  auto nnzPerRowUpperBound_h = nnzPerRowUpperBound.view_host();

  // Count the local elements and get node index upper bound
  stk::mesh::for_each_entity_run(bulkData, ELEMENT_RANK, locallyOwnedSelector,
    [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
    {
      stk::mesh::ConnectedEntities nodes = mesh.get_connected_entities(elem,NODE_RANK);

      // NOTE: This will substantially overcount the NNZ needed.  You should use hash table to be smarter
      for (unsigned inode = 0; inode < nodes.size(); ++inode) {
        GO GID = mesh.identifier(nodes[inode])-1;
        LO LID = ownedPlusSharedMapG->getLocalElement(GID);
        if (LID != Teuchos::OrdinalTraits<LO>::invalid())
          nnzPerRowUpperBound_h[LID]+=nodes.size();
      }
    });

  // Build the Graph
  RCP<Tpetra_FECrsGraph> StiffGraph = rcp(new Tpetra_FECrsGraph(globalMapG,ownedPlusSharedMapG,nnzPerRowUpperBound));
  Tpetra::beginAssembly(*StiffGraph);
  stk::mesh::for_each_entity_run(bulkData, ELEMENT_RANK, locallyOwnedSelector,
    [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
    {
      stk::mesh::ConnectedEntities nodes = mesh.get_connected_entities(elem,NODE_RANK);

      Teuchos::Array<GO> global_ids(nodes.size());
      bool foundOwnedNode = false;
      for (unsigned inode = 0; inode < nodes.size(); ++inode) {
        if (mesh.bucket(nodes[inode]).owned()) {
          foundOwnedNode = true;
        }
        GO GID = mesh.identifier(nodes[inode])-1;
        global_ids[inode]=GID;
      }

      if (!foundOwnedNode) {
        std::cout<<"Warning, element "<<mesh.identifier(elem)<<" on P"<<mesh.parallel_rank()<<" doesn't have any locally-owned nodes."<<std::endl;
      }

      for (unsigned inode = 0; inode < nodes.size(); ++inode) {
        StiffGraph->insertGlobalIndices(global_ids[inode],global_ids());
      }
    });

  Tpetra::endAssembly(*StiffGraph);
  tm = Teuchos::null;

  /**********************************************************************************/
  /******************** COMPUTE COORDINATES AND STATISTICS **************************/
  /**********************************************************************************/
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3) Compute Coordinates and Stats")));

  std::vector<entity_type> bcNodes;

  stk::mesh::ExodusTranslator exodusTranslator(bulkData);
  stk::mesh::PartVector nodeSets = exodusTranslator.get_node_set_parts();
  ShardsCellTopology cellType;
  for (const stk::mesh::Part* nodeSet : nodeSets) {

      stk::mesh::Selector bcNodeSelector = *nodeSet & locallyOwnedSelector;

      if(bcNodes.size() == 0)
        stk::mesh::get_entities(bulkData, NODE_RANK, bcNodeSelector, bcNodes);
      else {
        std::vector<entity_type> myBcNodes;
        stk::mesh::get_entities(bulkData, NODE_RANK, bcNodeSelector, myBcNodes);
        std::copy(myBcNodes.begin(),myBcNodes.end(),std::back_inserter(bcNodes));
      }

      if(MyPID==0) printf("Adding nodeset %s\n",nodeSet->name().c_str());
  }

  stk::mesh::PartVector elemBlocks = exodusTranslator.get_element_block_parts();
  for(const stk::mesh::Part* elemBlock : elemBlocks) {
      // Here the topology is defined from the mesh. Note that it is assumed
      // that all parts with elements (aka ELEMENT_RANK) have the same topology type
      auto myCell = stk::mesh::get_cell_topology(metaData.get_topology( *elemBlock ));
      if(myCell.getCellTopologyData()) {
        cellType =  myCell;
      }
  }

  int numNodesPerElem = cellType.getNodeCount();
  int numEdgesPerElem = cellType.getEdgeCount();

  if(MyPID==0) {
    std::cout<<"Cell Topology: "<<cellType.getName() << " ("<<cellType.getBaseName()<<")"<<std::endl;
  }

  //MachineLearningStatistics_Hex3D<SC,LO,GO,NO> MLStatistics(numGlobalElements);
  bool do_statistics = !strcmp(cellType.getName(),"Hexahedron_8") && (numRanks == 1);
  if(MyPID==0) std::cout<<"do_statistics = "<<do_statistics<<std::endl;

  // if no boundary node set was found give a warning
  int numLocalBCs = bcNodes.size();
  int numGlobalBCs = 0;
  Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, numLocalBCs, Teuchos::outArg(numGlobalBCs));

  if (MyPID == 0) {
    std::cout << "       Number of Dirichlet Nodes: " << numGlobalBCs << std::endl;
    if (numGlobalBCs == 0) {
      std::cout << std::endl
                << "     Warning! - No boundary node set found." << std::endl
                << "     Boundary conditions will not be applied correctly.\n"
                << std::endl;
    }
  }

  if (optPrintLocalStats) {
    for (int i=0; i<numRanks; ++i) {
      if (MyPID == i) {
        std::cout << "(" << MyPID << ")    Number of local b.c. nodes: " << bcNodes.size() << std::endl;
      }
      Comm->barrier();
    }
    if(MyPID==0) std::cout << std::endl;
  }

  tm = Teuchos::null;

  /**********************************************************************************/
  /********************************* GET CUBATURE ***********************************/
  /**********************************************************************************/
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("4) Assembly Prep")));
  // Define cubature of the specified degree for the cellType
  Intrepid2::DefaultCubatureFactory  cubFactory;
  int cubDegree = 2;

  RCP<Intrepid2::Cubature<device_type> > cellCubature = cubFactory.create<device_type>(cellType, cubDegree);

  int cubDim       = cellCubature -> getDimension();
  int numCubPoints = cellCubature -> getNumPoints();

  // Get numerical integration points and weights
  Intrepid2ScalarView cubPoints ("cubPoints", numCubPoints, cubDim);
  Intrepid2ScalarView cubWeights("cubWeights", numCubPoints);

  cellCubature -> getCubature(cubPoints, cubWeights);

  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  // Select basis from the cell topology
  int order = 1;
  spaceDim = getDimension(cellType);
  auto HGradBasis = Intrepid2::getBasis<Intrepid2::DerivedNodalBasisFamily<device_type> >(cellType, Intrepid2::FUNCTION_SPACE_HGRAD, order);


  int numFieldsG = HGradBasis->getCardinality();

  Intrepid2ScalarView basisValues("HGBValues", numFieldsG, numCubPoints);
  Intrepid2ScalarView basisGrads("HGBGrads", numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
  HGradBasis->getValues(basisValues, cubPoints, Intrepid2::OPERATOR_VALUE);
  HGradBasis->getValues(basisGrads, cubPoints, Intrepid2::OPERATOR_GRAD);


  /**********************************************************************************/
  /******************************** STASH COORDINATES *******************************/
  /**********************************************************************************/
  typedef stk::mesh::Field<double>  CoordFieldType;
  // get coordinates field
  CoordFieldType *coords = metaData.get_field<double>(NODE_RANK,"coordinates");

  // Put coordinates in multivector for output
  RCP<Tpetra_MultiVector> nCoord = rcp(new Tpetra_MultiVector(globalMapG,spaceDim));
  auto nCoord_h = nCoord->get2dViewNonConst();
  // Loop over elements

  const stk::mesh::BucketVector &localElementBuckets = bulkData.get_buckets(ELEMENT_RANK, locallyOwnedSelector);
  for (const stk::mesh::Bucket* elemBucketPtr : localElementBuckets) {
    for (stk::mesh::Entity elem : *elemBucketPtr) {
      stk::mesh::ConnectedEntities nodes = bulkData.get_connected_entities(elem,NODE_RANK);
      for (unsigned inode = 0; inode < nodes.size(); ++inode) {
        double *coord = stk::mesh::field_data(*coords, nodes[inode]);
        LO lid = globalMapG->getLocalElement((int)bulkData.identifier(nodes[inode]) -1);
        if(lid != Teuchos::OrdinalTraits<LO>::invalid()) {
          nCoord_h[0][lid] = coord[0];
          nCoord_h[1][lid] = coord[1];
          if(spaceDim==3)
            nCoord_h[2][lid] = coord[2];
        }
      }
    }
  } // end loop over elements


  /**********************************************************************************/
  /****************************** STATISTICS (Part I) *******************************/
  /**********************************************************************************/
  if(do_statistics) {
    Intrepid2::ScalarView<int, memory_space> elemToNode("elemToNode",numLocalElems,numNodesPerElem);
    Intrepid2::ScalarView<int, memory_space> elemToEdge("elemToEdge",numLocalElems,numEdgesPerElem);
    Intrepid2::ScalarView<double, memory_space> nodeCoord ("nodeCoord", numLocalNodes, spaceDim);
    Intrepid2::ScalarView<double, memory_space> sigmaVal("sigmaVal", numLocalElems);

    int elem_ct=0;
    std::map<std::pair<int,int>,int> local_edge_hash;
    std::vector<std::pair<int,int> > edge_vector;

    for (const stk::mesh::Bucket* elemBucketPtr : localElementBuckets) {
      for (stk::mesh::Entity elem : *elemBucketPtr) {
        stk::mesh::ConnectedEntities nodes = bulkData.get_connected_entities(elem,NODE_RANK);
        for (unsigned inode = 0; inode < nodes.size(); ++inode) {
          const double *coord = stk::mesh::field_data(*coords, nodes[inode]);
          LO lid = globalMapG->getLocalElement((int)bulkData.identifier(nodes[inode]) -1);
          elemToNode(elem_ct,inode) = lid;
          if(lid != -1) {
            nodeCoord(lid,0) = coord[0];
            nodeCoord(lid,1) = coord[1];
            if(spaceDim==3)
              nodeCoord(lid,2) = coord[2];
          }
        }//end node loop

        auto data = cellType.getCellTopologyData();
        for(unsigned iedge=0; iedge<cellType.getEdgeCount(); iedge++) {
          int n0 = data->edge[iedge].node[0];
          int n1 = data->edge[iedge].node[1];
          int lid0 = globalMapG->getLocalElement((int)bulkData.identifier(nodes[n0]) -1);
          int lid1 = globalMapG->getLocalElement((int)bulkData.identifier(nodes[n1]) -1);
          if(lid0 != -1 && lid1 != -1) {
            int lo = std::min(lid0,lid1);
            int hi = std::max(lid0,lid1);
            std::pair<int,int> key(lo,hi);
            if (local_edge_hash.find(key) == local_edge_hash.end()) {
              int new_edge_id = edge_vector.size();
              local_edge_hash[key] = new_edge_id;
              edge_vector.push_back(key);
              elemToEdge(elem_ct,iedge) = new_edge_id;
            }
            else {
              elemToEdge(elem_ct,iedge) = local_edge_hash[key];
            }
          }
        }//end edge loop

        sigmaVal(elem_ct) = 1;// Not doing sigma here

        elem_ct++;
      }//end element loop
    }//end bucket loop

    Intrepid2::ScalarView<int, memory_space> edgeToNode("edgeToNode", edge_vector.size(), 2);
    for(int i=0; i<(int)edge_vector.size(); i++) {
      edgeToNode(i,0) = edge_vector[i].first;
      edgeToNode(i,1) = edge_vector[i].second;
    }


    //MLStatistics.Phase1(elemToNode,elemToEdge,edgeToNode,nodeCoord,sigmaVal);

  }
  tm = Teuchos::null;

  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("5) Matrix/RHS Assembly")));
  RCP<Tpetra_FECrsMatrix> StiffMatrix = rcp(new Tpetra_FECrsMatrix(StiffGraph));

  // Define desired workset size and count how many worksets there are on this processor's mesh block
  //int desiredWorksetSize = numElems;                   // change to desired workset size!
  ////int desiredWorksetSize = 100;                      // change to desired workset size!
  //int numWorksets        = numElems/desiredWorksetSize;
  int desiredWorksetSize = numLocalElems;                   // change to desired workset size!
  //int desiredWorksetSize = 100;                      // change to desired workset size!
  int numWorksets        = numLocalElems/desiredWorksetSize;

  // When numLocalElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numLocalElems) numWorksets += 1;

  if (MyPID == 0) {
    std::cout << "       Desired workset size:                 " << desiredWorksetSize << std::endl;
    std::cout << "       Number of worksets (per processor):   " << numWorksets << std::endl << std::endl;
  }



  // Start assembly
  RCP<Tpetra_FEMultiVector> rhsVector = rcp (new Tpetra_FEMultiVector(globalMapG,StiffGraph->getImporter(),1));
  Tpetra::beginAssembly(*StiffMatrix,*rhsVector);

  // Right now, this loop only increments once:
  //   numWorkset = 1
  //   start      = 0
  //   end        = numLocalElems
  //   worksetSize = numLocalElems
  for(int workset = 0; workset < numWorksets; workset++) {

    // compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // when numLocalElems is not divisible by desiredWorksetSize, the last workset ends at numLocalElems
    worksetEnd   = (worksetEnd <= numLocalElems) ? worksetEnd : numLocalElems;

    // allocate the array for the cell nodes
    worksetSize  = worksetEnd - worksetBegin;
    Intrepid2ScalarView cellWorkset("cellWorkset", worksetSize, numNodesPerElem, spaceDim);

    // copy coordinates into cell workset
    int cellCounter = 0;
    for (const stk::mesh::Bucket* elemBucket : localElementBuckets) {
      for (stk::mesh::Entity elem : *elemBucket) {
        stk::mesh::ConnectedEntities nodes = bulkData.get_connected_entities(elem,NODE_RANK);
        for (unsigned inode = 0; inode < nodes.size(); ++inode) {
          double *coord = stk::mesh::field_data(*coords, nodes[inode]);
          cellWorkset(cellCounter, inode, 0) = coord[0];
          cellWorkset(cellCounter, inode, 1) = coord[1];
          if(spaceDim==3) cellWorkset(cellCounter, inode, 2) = coord[2];
        }
        cellCounter++;
      }
    }

    /**********************************************************************************/
    /*                                Allocate arrays                                 */
    /**********************************************************************************/

    // Containers for Jacobians, integration measure & cubature points in workset cells
    Intrepid2ScalarView worksetJacobian  ("worksetJacobian", worksetSize, numCubPoints, spaceDim, spaceDim);
    Intrepid2ScalarView worksetJacobInv  ("worksetJacobInv", worksetSize, numCubPoints, spaceDim, spaceDim);
    Intrepid2ScalarView worksetJacobDet  ("worksetJacobDet", worksetSize, numCubPoints);
    Intrepid2ScalarView worksetCubWeights("worksetCubWeights", worksetSize, numCubPoints);
    Intrepid2ScalarView worksetCubPoints ("worksetCubPoints", worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and them multiplied by cubature weights
    Intrepid2ScalarView worksetBasisValues        ("worksetBasisValues", worksetSize, numFieldsG, numCubPoints);
    Intrepid2ScalarView worksetBasisValuesWeighted("worksetBasisValuesWeighted", worksetSize, numFieldsG, numCubPoints);
    Intrepid2ScalarView worksetBasisGrads         ("worksetBasisGrads", worksetSize, numFieldsG, numCubPoints, spaceDim);
    Intrepid2ScalarView worksetBasisGradsWeighted ("worksetBasisGradsWeighted", worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
    Intrepid2ScalarView worksetDiffusiveFlux("worksetDiffusiveFlux", worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for material values and source term. Require user-defined functions
    Intrepid2ScalarView worksetMaterialVals ("worksetMaterialVals", worksetSize, numCubPoints, spaceDim, spaceDim);
    Intrepid2ScalarView worksetSourceTerm   ("worksetSourceTerm", worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization matrix and the right hand side
    Intrepid2ScalarView worksetStiffMatrix ("worksetStiffMatrix", worksetSize, numFieldsG, numFieldsG);
    Intrepid2ScalarView worksetRHS         ("worksetRHS", worksetSize, numFieldsG);


    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/

    Intrepid2CTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    Intrepid2CTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    Intrepid2CTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    // map cubature points to physical frame
    Intrepid2CTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

    // get A at cubature points
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // get source term at cubature points
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

    /**********************************************************************************/
    /*                         Compute Stiffness Matrix                               */
    /**********************************************************************************/

    // Transform basis gradients to physical frame:                        DF^{-T}(grad u)
    Intrepid2FSTools::HGRADtransformGRAD(worksetBasisGrads,
                                                worksetJacobInv,   basisGrads);

    // Compute integration measure for workset cells:                      Det(DF)*w = J*w
    Intrepid2FSTools::computeCellMeasure(worksetCubWeights,
                                                worksetJacobDet, cubWeights);


    // Multiply transformed (workset) gradients with weighted measure:     DF^{-T}(grad u)*J*w
    Intrepid2FSTools::multiplyMeasure(worksetBasisGradsWeighted,
                                             worksetCubWeights, worksetBasisGrads);


    // Compute material tensor applied to basis grads:                     A*(DF^{-T}(grad u)
    Intrepid2FSTools::tensorMultiplyDataField(worksetDiffusiveFlux,
                                                     worksetMaterialVals,
                                                     worksetBasisGrads);

    // Integrate to compute contribution to global stiffness matrix:      (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
    Intrepid2FSTools::integrate(worksetStiffMatrix,
                                       worksetBasisGradsWeighted,
                                       worksetDiffusiveFlux);

    /**********************************************************************************/
    /*                                   Compute RHS                                  */
    /**********************************************************************************/

    // Transform basis values to physical frame:                        clones basis values (u)
    Intrepid2FSTools::HGRADtransformVALUE(worksetBasisValues,
                                                 basisValues);

    // Multiply transformed (workset) values with weighted measure:     (u)*J*w
    Intrepid2FSTools::multiplyMeasure(worksetBasisValuesWeighted,
                                             worksetCubWeights, worksetBasisValues);

    // Integrate worksetSourceTerm against weighted basis function set:  f.(u)*J*w
    Intrepid2FSTools::integrate(worksetRHS,
                                       worksetSourceTerm,
                                       worksetBasisValuesWeighted);

    /**********************************************************************************/
    /***************************** STATISTICS (Part II) ******************************/
    /**********************************************************************************/
    /*
      if(do_statistics)
      MLStatistics.Phase2a(worksetJacobDet,worksetCubWeights);
    */

    /**********************************************************************************/
    /*                         Assemble into Global Matrix                            */
    /**********************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numLocalElems
    //JJH runs from 0 to (#local cells - 1)
    int worksetCellOrdinal = 0;
    for (const stk::mesh::Bucket* elemBucketPtr : localElementBuckets) {
      for (stk::mesh::Entity elem : *elemBucketPtr) {

        // Compute cell ordinal relative to the current workset

        stk::mesh::ConnectedEntities worksetNodes = bulkData.get_connected_entities(elem,NODE_RANK);

        // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
        for (int cellRow = 0; cellRow < numFieldsG; cellRow++) {

          int globalRow = bulkData.identifier(worksetNodes[cellRow]) - 1;
          double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);
          rhsVector->sumIntoGlobalValue(globalRow, 0, sourceTermContribution);

          // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
          for (int cellCol = 0; cellCol < numFieldsG; cellCol++){
            int globalCol = bulkData.identifier(worksetNodes[cellCol]) - 1;
            double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);

            Teuchos::Array<GO> col(1); col[0] = globalCol;
            Teuchos::Array<SC> val(1); val[0] = operatorMatrixContribution;
            StiffMatrix->sumIntoGlobalValues(globalRow, col(), val());

          }// end cell col loop

        }// end cell row loop

        worksetCellOrdinal++;
      }// end workset cell loop


    } //for (localElementBuckets

  }// end workset loop


  Tpetra::endAssembly(*StiffMatrix,*rhsVector);

/**********************************************************************************/
/***************************** STATISTICS (Part IIb) ******************************/
/**********************************************************************************/
  /*
  if(do_statistics){
        MLStatistics.Phase2b(Xpetra::toXpetra(StiffMatrix->getGraph()),Teuchos::rcp(&nCoord,false));
  }
  */
/**********************************************************************************/
/***************************** STATISTICS (Part III) ******************************/
/**********************************************************************************/
  /*
  if(do_statistics){
    MLStatistics.Phase3();
    Teuchos::ParameterList problemStatistics = MLStatistics.GetStatistics();
    if(MyPID==0) std::cout<<"*** Problem Statistics ***"<<std::endl<<problemStatistics<<std::endl;
  }
  */
  tm = Teuchos::null;

/**********************************************************************************/
/************************** DIRICHLET BC SETUP ************************************/
/**********************************************************************************/
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("6) Matrix/RHS Dirichlet BCs")));
  RCP<Tpetra_MultiVector> lhsVector = rcp(new Tpetra_MultiVector(globalMapG,1));

  std::vector<SC> solution_values; solution_values.reserve(bcNodes.size());
  std::vector<int> ownedBoundaryNodes; ownedBoundaryNodes.reserve(bcNodes.size());
  // Loop over boundary nodes
  for (unsigned i = 0; i < bcNodes.size(); i++) {
    int bcNodeId = bulkData.identifier(bcNodes[i]);
    int lid = globalMapG->getLocalElement((int) bcNodeId -1);
    if(lid != -1) {
      ownedBoundaryNodes.push_back(lid);

      // get coordinates for this node
      entity_type bcnode = bulkData.get_entity(NODE_RANK,bcNodeId);
      double * coord = stk::mesh::field_data(*coords, bcnode);

      // look up exact value of function on boundary
      double x  = coord[0];
      double y  = coord[1];
      double z  = (spaceDim==3) ? coord[2] : 0;
      solution_values.push_back(exactSolution(x, y, z));
    }
  } // end loop over boundary nodes

  // Apply the Dirichlet conditions to matrix, lhs and rhs
  Apply_Dirichlet_BCs(ownedBoundaryNodes,*StiffMatrix,*lhsVector,*rhsVector,solution_values);

  tm = Teuchos::null;



   // Optionally dump the matrix and/or its coords to files.
  {
    typedef Tpetra::MatrixMarket::Writer<Tpetra_CrsMatrix> writer_type;
    if (matrixFilename != "") {
      writer_type::writeSparseFile (matrixFilename, StiffMatrix);
    }
    if (rhsFilename != "") {
      writer_type::writeDenseFile (rhsFilename, rhsVector);
    }
    if (initialGuessFilename != "") {
      writer_type::writeDenseFile (initialGuessFilename, lhsVector);
    }
    if (coordsFilename != "") {
      writer_type::writeDenseFile (coordsFilename, nCoord);
    }
  }

  /**********************************************************************************/
  /*********************************** SOLVE ****************************************/
  /**********************************************************************************/
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("7) Analytic Solution Generation")));

  // Run the solver
  Teuchos::ParameterList MLList = inputSolverList;
  RCP<Tpetra_Vector> exactNodalVals = rcp(new Tpetra_Vector(globalMapG));
  auto exactNodalVals_h = exactNodalVals->getDataNonConst();
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  const stk::mesh::BucketVector &localNodeBuckets = bulkData.get_buckets(NODE_RANK, locallyOwnedSelector);
  for (const stk::mesh::Bucket* bucketPtr : localNodeBuckets) {
    for (stk::mesh::Entity node : *bucketPtr) {
      double * coord = stk::mesh::field_data(*coords, node);
      // look up exact value of function on boundary
      double x  = coord[0];
      double y  = coord[1];
      double z  = (spaceDim==3) ? coord[2] : 0;
      int gNodeId = bulkData.identifier(node) - 1;
      //exactNodalVals[0][gNodeId]=exactSolution(x, y, z);
      int lid = globalMapG->getLocalElement(gNodeId);
      exactNodalVals_h[lid]=exactSolution(x, y, z);
    }
  }

  char probType[10] = "laplace";

  Teuchos::RCP<Tpetra_CrsMatrix> Acrs = Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(StiffMatrix);
  Teuchos::RCP<Tpetra_MultiVector> rhsMV = Teuchos::rcp_dynamic_cast<Tpetra_MultiVector>(rhsVector);
  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("8) Linear solve")));
  TestMultiLevelPreconditioner(probType,             MLList,
                               Acrs,                 exactNodalVals,
                               rhsMV,                lhsVector, nCoord,
                               TotalErrorResidual,   TotalErrorExactSol);

   // Timer Output
   tm = Teuchos::null;
   const bool alwaysWriteLocal = false;
   const bool writeGlobalStats = true;
   const bool writeZeroTimers  = false;
   const bool ignoreZeroTimers = true;
   const std::string filter    = "";
   TimeMonitor::summarize(Comm.ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                          writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);

   return 0;

}
/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z) {

  // Patch test: tri-linear function is in the FE space and should be recovered
  return 1. + x + y + z + x*y + x*z + y*z + x*y*z;


  // Patch test - tet: function is in the FE space and should be recovered
  //     return 1. + x + y + z ;

  // Patch test - hex: tri-linear function is in the FE space and should be recovered
  // return 1. + x + y + z + x*y + x*z + y*z + x*y*z;

  // Analytic solution with homogeneous Dirichlet boundary data for [0 1]x[0 1]x[0 1]
  // return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  // return exp(x + y + z)/(1. + x*y + y*z + x*y*z);
}


template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar& x, const Scalar& y, const Scalar& z) {

  material[0][0] = 1.;
  material[0][1] = 0.;
  material[0][2] = 0.;
  //
  material[1][0] = 0.;
  material[1][1] = 1.;
  material[1][2] = 0.;
  //
  material[2][0] = 0.;
  material[2][1] = 0.;
  material[2][2] = 1.;
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/************ Grad of Exact Solution ****************/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z) {

  // To enable derivatives of the gradient (i.e., 2nd derivatives of the exact solution) need 2 levels of fad types
  Sacado::Fad::SFad<Scalar,3> fad_x = x;
  Sacado::Fad::SFad<Scalar,3> fad_y = y;
  Sacado::Fad::SFad<Scalar,3> fad_z = z;
  Sacado::Fad::SFad<Scalar,3> u;

  // Indicate the independent variables
  fad_x.diff(0,3);
  fad_y.diff(1,3);
  fad_z.diff(2,3);

  u = exactSolution(fad_x, fad_y, fad_z);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
  gradExact[2] = u.dx(2);
}

/************ Source Term (RHS) ****************/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z){

  Scalar u;
  Scalar grad_u[3];
  Scalar flux[3] = {0.0, 0.0, 0.0};
  Scalar material[3][3];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,3);
  y.diff(1,3);
  z.diff(2,3);

  // Get exact solution and its gradient
  u = exactSolution(x, y, z);
  exactSolutionGrad(grad_u, x, y, z);

  // Get material tensor
  materialTensor<Scalar>(material, x, y, z);

  // Compute total flux = (A.grad u)
  for(int i = 0; i < 3; i++){

    // Add diffusive flux
    for(int j = 0; j < 3; j++){
      flux[i] += material[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(A.grad u)
  f = -(flux[0].dx(0) + flux[1].dx(1) + flux[2].dx(2));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

/************ Material Tensor ****************/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        matTensorValues,
                             const ArrayIn &   evaluationPoints){

  int numWorksetCells  = evaluationPoints.extent_int(0);
  int numPoints        = evaluationPoints.extent_int(1);
  int spaceDim2         = evaluationPoints.extent_int(2);

  double material[3][3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = 0.0;
      if(spaceDim2==3) z = evaluationPoints(cell, pt, 2);

      materialTensor<double>(material, x, y, z);

      for(int row = 0; row < spaceDim2; row++){
        for(int col = 0; col < spaceDim2; col++){
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}

/************ Source Term (RHS) ****************/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.extent_int(0);
  int numPoints = evaluationPoints.extent_int(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      Sacado::Fad::SFad<double,3> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,3> y = evaluationPoints(cell, pt, 1);
      Sacado::Fad::SFad<double,3> z;
      if(spaceDim==3)
	z = evaluationPoints(cell, pt, 2);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,3> >(x, y, z).val();
    }
  }
}

/************ Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.extent_int(0);
  int numPoints = evaluationPoints.extent_int(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z=0.0;
      if(spaceDim==3)
	z = evaluationPoints(cell, pt, 2);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y, z);
    }
  }
}


/************ Grad of Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.extent_int(0);
  int numPoints = evaluationPoints.extent_int(1);
  int spaceDim2  = evaluationPoints.extent_int(2);

  double gradient[3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = 0.0;
      if(spaceDim2==3)
	z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<double>(gradient, x, y, z);

      for(int row = 0; row < spaceDim2; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}

/**********************************************************************************/
/******************************* TEST MueLu****************************************/
/**********************************************************************************/

// Test MueLu
int TestMultiLevelPreconditioner(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Teuchos::RCP<Tpetra_CrsMatrix>   & A,
                                 Teuchos::RCP<Tpetra_Vector> & xexact,
                                 Teuchos::RCP<Tpetra_MultiVector> & b,
                                 Teuchos::RCP<Tpetra_MultiVector> & uh,
				 Teuchos::RCP<Tpetra_MultiVector> & coords,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace TrilinosCouplings::TpetraIntrepidPoissonExample;
  using MT = Teuchos::ScalarTraits<Tpetra_MultiVector::scalar_type>::magnitudeType;


  // Solver params
  std::string solverName = MLList.get("solver","cg");
  int maxNumIters = MLList.get("Maximum Iterations",200);
  double tol = MLList.get("Convergence Tolerance",1e-10);
  int num_steps = MLList.get("Number of Time Steps",1);


  // Multigrid Hierarchy
  Teuchos::ParameterList mueluParams;
  if (MLList.isSublist("MueLu"))
    mueluParams = MLList.sublist("MueLu");
  // Xpetrify coordinates
  mueluParams.sublist("user data").set("Coordinates",Xpetra::toXpetra(coords));
  if(A->getMap()->getComm()->getRank()==0)
    std::cout<<"*** MueLu Params ***" <<std::endl<<mueluParams<<std::endl;

  RCP<Tpetra::Operator<> > Aop = Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(A);
  auto M = MueLu::CreateTpetraPreconditioner(Aop,mueluParams);

  // Do the linear solve(s).
  bool converged = false;
  int numItersPerformed = 0;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Solve", total_solve);
    solveWithBelos (converged, numItersPerformed, solverName, tol,
                    maxNumIters, num_steps, uh, A, b, Teuchos::null, M);
  }

  if (b->getMap()->getComm()->getRank() == 0) {
    std::cout<<"Total Iterations: "<<numItersPerformed<<std::endl;
  }

  // Compute ||X-X_exact||_2
  Teuchos::Array<MT> norm_x(1), norm_error(1);
  xexact->norm2 (norm_x());
  xexact->update (-1.0, *uh, 1.0);
  xexact->norm2 (norm_error());
  if (b->getMap()->getComm()->getRank() == 0) {
    std::cout << std::endl
              << "||X - X_exact||_2 / ||X_exact||_2 = " << norm_error[0] / norm_x[0]
              << std::endl;
  }

  TotalErrorExactSol += norm_error[0];
  TotalErrorResidual += norm_x[0];


  return( numItersPerformed );
}

/**********************************************************************************/
/**************************** SIMPLE BASIS FACTORY ********************************/
/**********************************************************************************/


int getDimension( const ShardsCellTopology & cellTopology) {
  switch (cellTopology.getKey()) {
  case shards::Tetrahedron<4>::key:
  case shards::Hexahedron<8>::key:
    return 3;
  case shards::Triangle<3>::key:
  case shards::Quadrilateral<4>::key:
    return 2;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(1,std::invalid_argument,
			       "Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4, Quadrilateral_4 or Triangle_3");
  }
}


int main(int argc, char *argv[]) {
  Kokkos::initialize(argc,argv);
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  main_(argc,argv);
  Kokkos::finalize();
}
