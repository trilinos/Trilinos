// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
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
                       u = g on Gamma

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

#include "TrilinosCouplings_config.h"

// Teuchos includes
#include <Teuchos_CommandLineProcessor.hpp>

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#if defined(HAVE_TRILINOSCOUPLINGS_MUELU) 
// MueLu includes
#include "MueLu.hpp"
#if defined(HAVE_MUELU_EPETRA)
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_EpetraOperator.hpp"
#endif
#endif

#ifdef HAVE_INTREPID_KOKKOSCORE
#include "Sacado.hpp"
#else
// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#endif

// STK includes
#include "Ionit_Initializer.h"
#include "Ioss_SubSystem.h"

#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

#include "stk_util/parallel/Parallel.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Comm.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/GetBuckets.hpp"
#include "stk_mesh/base/CreateAdjacentEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian3d, etc
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef shards::CellTopology             ShardsCellTopology;
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;
typedef Intrepid::FieldContainer<double> IntrepidFieldContainer;


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
int TestMultiLevelPreconditioner(char                        ProblemType[],
                                 Teuchos::ParameterList &    MLList,
                                 Epetra_CrsMatrix &          A,
                                 const Epetra_MultiVector &  xexact,
                                 Epetra_MultiVector &        b,
                                 Epetra_MultiVector &        uh,
				 Epetra_MultiVector & coords,
                                 double &                    TotalErrorResidual,
                                 double &                    TotalErrorExactSol);


/**********************************************************************************/
/************* FUNCTION DECLARATIONS FOR SIMPLE BASIS FACTORY *********************/
/**********************************************************************************/

/** \brief  Simple factory that chooses basis function based on cell topology.

    \param  cellTopology  [in]    Shards cell topology
    \param  order         [in]    basis function order, currently unused
    \param  basis         [out]   pointer to Intrepid basis

    \return Intrepid basis
 */

void getBasis(Teuchos::RCP<Intrepid::Basis<double,IntrepidFieldContainer > > &basis,
               const ShardsCellTopology & cellTopology,
               int order);

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
  size_t input_index = broker.add_mesh_database(absoluteFileName, type, stk::io::READ_MESH);
  broker.set_active_mesh(input_index);
  // creates metadata
  broker.create_input_mesh();

  // commits the meta data
  broker.populate_bulk_data();

} //mesh_read_write


/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using entity_type = stk::mesh::Entity;

  const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;
  const stk::mesh::EntityRank ELEMENT_RANK = stk::topology::ELEMENT_RANK;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int numRanks = Comm.NumProc();
  int MyPID = Comm.MyPID();

  Epetra_Time Time(Comm);
 
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
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  STK's website:      http://trilinos.sandia.gov/packages/stk                |\n" \
    << "|  ML's website:       http://trilinos.sandia.gov/packages/ml                 |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
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
  stk::io::StkMeshIoBroker broker(Comm.GetMpiComm());
  broker.property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 180));
  broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "rcb"));
  broker.property_add(Ioss::Property("COMPOSE_RESULTS", false));  //Note!  true results in an error in Seacas

  std::string type = "exodusii";
  char buf[1024];
  if (getcwd(buf,sizeof(buf)) == NULL)
    throw(std::runtime_error("Could not get current working directory"));
  std::string working_directory(buf);
  std::string filename(optMeshFile);
  int db_integer_size = 4;
  stk::io::HeartbeatType hb_type = stk::io::NONE;

  mesh_read_write(type, working_directory, filename, broker, db_integer_size, hb_type);

  stk::mesh::BulkData &bulkData = broker.bulk_data();
  stk::mesh::MetaData &metaData = broker.meta_data();

  // Count number of local nodes, record GIDs for Epetra map.
  std::vector<int> epetraGIDs;
  int numLocalNodes=0;
  stk::mesh::Selector locallyOwnedSelector = metaData.locally_owned_part(); //locally-owned
  const stk::mesh::BucketVector &localNodeBuckets = bulkData.get_buckets(NODE_RANK, locallyOwnedSelector);
  for (size_t bucketIndex = 0; bucketIndex < localNodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket &nodeBucket = *localNodeBuckets[bucketIndex];
    numLocalNodes += nodeBucket.size();
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      epetraGIDs.push_back(bulkData.identifier(node)-1);
    }
  }


  // Count number of local elements
  int numLocalElems=0;
  const stk::mesh::BucketVector &localElementBuckets = bulkData.get_buckets(ELEMENT_RANK, locallyOwnedSelector);
  for (size_t bucketIndex = 0; bucketIndex < localElementBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket &elementBucket = *localElementBuckets[bucketIndex];
    numLocalElems += elementBucket.size();
  }

  if (optPrintLocalStats) {
    for (int i=0; i<numRanks; ++i) {
      if (MyPID == i) {
        std::cout << "(" << MyPID << ")    Number of local Elements: " << numLocalElems << std::endl
                  << "(" << MyPID << ")       Number of local Nodes: " << numLocalNodes << std::endl << std::endl;
      }
      Comm.Barrier();
    }
  }

  int numGlobalNodes = 0;
  Comm.SumAll(&numLocalNodes,&numGlobalNodes,1);
  if (MyPID == 0)
    std::cout << "       Number of global Nodes: " << numGlobalNodes << std::endl;

  typedef stk::mesh::Field<double, stk::mesh::Cartesian>  CoordFieldType;
  // get coordinates field
  CoordFieldType *coords = metaData.get_field<CoordFieldType>(NODE_RANK,"coordinates");

  // get buckets containing entities of node rank
  stk::mesh::Selector nothingSelector;
  stk::mesh::Selector allSelector(!nothingSelector);
  stk::mesh::BucketVector const & nodeBuckets = bulkData.get_buckets( NODE_RANK,allSelector );
  std::vector<entity_type> bcNodes;

  // loop over all mesh parts
  const stk::mesh::PartVector & all_parts = metaData.get_parts();
  ShardsCellTopology cellType;
  for (stk::mesh::PartVector::const_iterator i  = all_parts.begin(); i != all_parts.end(); ++i) {

    stk::mesh::Part & part = **i ;

    // if part only contains nodes, then it is a node set
    //   ! this assumes that the only node set defined is the set
    //   ! of boundary nodes
    if (part.primary_entity_rank() == NODE_RANK) {
      stk::mesh::Selector partSelector(part);
      stk::mesh::Selector bcNodeSelector = partSelector & locallyOwnedSelector;
      stk::mesh::get_selected_entities(bcNodeSelector, nodeBuckets, bcNodes);
    }
    else if(part.primary_entity_rank() == ELEMENT_RANK) {
      // Here the topology is defined from the mesh. Note that it is assumed
      // that all parts with elements (aka ELEMENT_RANK) have the same topology type
      auto myCell = stk::mesh::get_cell_topology(metaData.get_topology( part ));
      if(myCell.getCellTopologyData()) {
	cellType =  myCell;
      }
    }
    
  } // end loop over mesh parts

  int numNodesPerElem = cellType.getNodeCount();  

  // if no boundary node set was found give a warning
  int numLocalBCs = bcNodes.size();
  int numGlobalBCs = 0;
  Comm.SumAll(&numLocalBCs,&numGlobalBCs,1);
  if (numGlobalBCs == 0) {
    if (MyPID == 0) {
      std::cout << std::endl
                << "     Warning! - No boundary node set found." << std::endl
                << "  Boundary conditions will not be applied correctly.\n"
                << std::endl;
    }
  }

  if (optPrintLocalStats) {
    if (MyPID == 0)
      std::cout << "       Number of global b.c. nodes: " << numGlobalBCs << std::endl;
    for (int i=0; i<numRanks; ++i) {
      if (MyPID == i) {
        std::cout << "(" << MyPID << ")    Number of local b.c. nodes: " << bcNodes.size() << std::endl;
      }
      Comm.Barrier();
    }
    if(MyPID==0) std::cout << std::endl;
  }

  if(MyPID==0) {
    std::cout << "Read mesh" << "                                   "
              << Time.ElapsedTime() << " seconds" << std::endl;
    Time.ResetStartTime();

  }

  /**********************************************************************************/
  /********************************* GET CUBATURE ***********************************/
  /**********************************************************************************/

  // Define cubature of the specified degree for the cellType
  Intrepid::DefaultCubatureFactory<double>  cubFactory;
  int cubDegree = 2;

  RCP<Intrepid::Cubature<double> > cellCubature = cubFactory.create(cellType, cubDegree);

  int cubDim       = cellCubature -> getDimension();
  int numCubPoints = cellCubature -> getNumPoints();

  // Get numerical integration points and weights
  IntrepidFieldContainer cubPoints (numCubPoints, cubDim);
  IntrepidFieldContainer cubWeights(numCubPoints);

  cellCubature -> getCubature(cubPoints, cubWeights);

  if(MyPID==0) {
    std::cout << "Getting cubature                            "
              << Time.ElapsedTime() << " seconds" << std::endl;
    Time.ResetStartTime();
  }

  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  // Select basis from the cell topology
  int order = 1;
  spaceDim = getDimension(cellType);
  RCP<Intrepid::Basis<double, IntrepidFieldContainer > >  HGradBasis;
  getBasis(HGradBasis, cellType, order);


  int numFieldsG = HGradBasis->getCardinality();
  IntrepidFieldContainer basisValues(numFieldsG, numCubPoints);
  IntrepidFieldContainer basisGrads(numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
  HGradBasis->getValues(basisValues, cubPoints, Intrepid::OPERATOR_VALUE);
  HGradBasis->getValues(basisGrads, cubPoints, Intrepid::OPERATOR_GRAD);

  if(MyPID==0) {
    std::cout << "Getting basis                               "
              << Time.ElapsedTime() << " seconds" << std::endl;
    Time.ResetStartTime();
  }


  /**********************************************************************************/
  /********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
  /**********************************************************************************/

  Epetra_Map   globalMapG(numGlobalNodes,numLocalNodes, &epetraGIDs[0], 0, Comm);
  Epetra_FECrsMatrix StiffMatrix(Copy, globalMapG, numFieldsG);
  Epetra_FEVector rhsVector(globalMapG);

  if(MyPID==0) {
    std::cout <<  "Build global maps                           "
              << Time.ElapsedTime() << " seconds" << std::endl;
    Time.ResetStartTime();
  }


  /**********************************************************************************/
  /******************************** STASH COORDINATES *******************************/
  /**********************************************************************************/
  // Put coordinates in multivector for output
  Epetra_MultiVector nCoord(globalMapG,spaceDim);    
  // Loop over elements
  
  for (size_t bucketIndex = 0; bucketIndex < localElementBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket &elemBucket = *localElementBuckets[bucketIndex];
    for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
      stk::mesh::Entity elem = elemBucket[elemIndex];
      //TODO (Optimization) It's assumed all elements are the same type, so this is constant.
      //TODO Therefore there's no need to do this everytime.
      unsigned numNodes = bulkData.num_nodes(elem);
      stk::mesh::Entity const* nodes = bulkData.begin_nodes(elem);
      for (unsigned inode = 0; inode < numNodes; ++inode) {
	double *coord = stk::mesh::field_data(*coords, nodes[inode]);
	int lid = globalMapG.LID((int)bulkData.identifier(nodes[inode]) -1);
	nCoord[0][lid] = coord[0];
	nCoord[1][lid] = coord[1];
	if(spaceDim==3)
	  nCoord[2][lid] = coord[2];
      }
    }      
  } // end loop over elements
  if(coordsFilename != "") {
    EpetraExt::MultiVectorToMatrixMarketFile(coordsFilename.c_str(),nCoord,0,0,false);
    if(MyPID==0) {Time.ResetStartTime();}
  }


  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

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
    std::cout << "\tDesired workset size:                 " << desiredWorksetSize << std::endl;
    std::cout << "\tNumber of worksets (per processor):   " << numWorksets << std::endl << std::endl;
    Time.ResetStartTime();
  }

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
    IntrepidFieldContainer cellWorkset(worksetSize, numNodesPerElem, spaceDim);

    // copy coordinates into cell workset
    int cellCounter = 0;
    for (size_t bucketIndex = 0; bucketIndex < localElementBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket &elemBucket = *localElementBuckets[bucketIndex];
      for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
        stk::mesh::Entity elem = elemBucket[elemIndex];
        //TODO (Optimization) It's assumed all elements are the same type, so this is constant.
        //TODO Therefore there's no need to do this everytime.
        unsigned numNodes = bulkData.num_nodes(elem);
        stk::mesh::Entity const* nodes = bulkData.begin_nodes(elem);
        for (unsigned inode = 0; inode < numNodes; ++inode) {
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
    IntrepidFieldContainer worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
    IntrepidFieldContainer worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
    IntrepidFieldContainer worksetJacobDet  (worksetSize, numCubPoints);
    IntrepidFieldContainer worksetCubWeights(worksetSize, numCubPoints);
    IntrepidFieldContainer worksetCubPoints (worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and them multiplied by cubature weights
    IntrepidFieldContainer worksetBasisValues        (worksetSize, numFieldsG, numCubPoints);
    IntrepidFieldContainer worksetBasisValuesWeighted(worksetSize, numFieldsG, numCubPoints);
    IntrepidFieldContainer worksetBasisGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
    IntrepidFieldContainer worksetBasisGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
    IntrepidFieldContainer worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for material values and source term. Require user-defined functions
    IntrepidFieldContainer worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);
    IntrepidFieldContainer worksetSourceTerm   (worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization matrix and the right hand side
    IntrepidFieldContainer worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
    IntrepidFieldContainer worksetRHS         (worksetSize, numFieldsG);

    if(MyPID==0) {
      std::cout << "Allocate arrays                             "
                << Time.ElapsedTime() << " seconds" << std::endl;
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/

    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    if(MyPID==0) {
      std::cout << "Calculate Jacobians                         "
                << Time.ElapsedTime() << " seconds" << std::endl;
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    // map cubature points to physical frame
    IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

    // get A at cubature points
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // get source term at cubature points
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

    if(MyPID==0) {
      std::cout << "Map to physical frame and get source term   "
                << Time.ElapsedTime() << " seconds" << std::endl;
      Time.ResetStartTime();
    }


    /**********************************************************************************/
    /*                         Compute Stiffness Matrix                               */
    /**********************************************************************************/

    // Transform basis gradients to physical frame:                        DF^{-T}(grad u)
    IntrepidFSTools::HGRADtransformGRAD<double>(worksetBasisGrads,
                                                worksetJacobInv,   basisGrads);

    // Compute integration measure for workset cells:                      Det(DF)*w = J*w
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeights,
                                                worksetJacobDet, cubWeights);


    // Multiply transformed (workset) gradients with weighted measure:     DF^{-T}(grad u)*J*w
    IntrepidFSTools::multiplyMeasure<double>(worksetBasisGradsWeighted,
                                             worksetCubWeights, worksetBasisGrads);


    // Compute material tensor applied to basis grads:                     A*(DF^{-T}(grad u)
    IntrepidFSTools::tensorMultiplyDataField<double>(worksetDiffusiveFlux,
                                                     worksetMaterialVals,
                                                     worksetBasisGrads);

    // Integrate to compute contribution to global stiffness matrix:      (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
    IntrepidFSTools::integrate<double>(worksetStiffMatrix,
                                       worksetBasisGradsWeighted,
                                       worksetDiffusiveFlux, Intrepid::COMP_BLAS);

    if(MyPID==0) {
      std::cout << "Compute stiffness matrix                    "
                << Time.ElapsedTime() << " seconds" << std::endl;
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*                                   Compute RHS                                  */
    /**********************************************************************************/

    // Transform basis values to physical frame:                        clones basis values (u)
    IntrepidFSTools::HGRADtransformVALUE<double>(worksetBasisValues,
                                                 basisValues);

    // Multiply transformed (workset) values with weighted measure:     (u)*J*w
    IntrepidFSTools::multiplyMeasure<double>(worksetBasisValuesWeighted,
                                             worksetCubWeights, worksetBasisValues);

    // Integrate worksetSourceTerm against weighted basis function set:  f.(u)*J*w
    IntrepidFSTools::integrate<double>(worksetRHS,
                                       worksetSourceTerm,
                                       worksetBasisValuesWeighted, Intrepid::COMP_BLAS);

    if(MyPID==0) {
      std::cout << "Compute right-hand side                     "
                << Time.ElapsedTime() << " seconds" << std::endl;
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*                         Assemble into Global Matrix                            */
    /**********************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numLocalElems
    //JJH runs from 0 to (#local cells - 1)
    int worksetCellOrdinal = 0;
    for (size_t bucketIndex = 0; bucketIndex < localElementBuckets.size(); ++bucketIndex) {

      stk::mesh::Bucket &elemBucket = *localElementBuckets[bucketIndex];
      for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {

        // Compute cell ordinal relative to the current workset

        // Get element entity from id of cell
        entity_type elem = elemBucket[elemIndex];
        entity_type const* worksetNodes = bulkData.begin_nodes(elem);

        // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
        for (int cellRow = 0; cellRow < numFieldsG; cellRow++) {

          int globalRow = bulkData.identifier(worksetNodes[cellRow]) - 1;
          double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);
          rhsVector.SumIntoGlobalValues(1, &globalRow, &sourceTermContribution);

          // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
          for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

            //int globalCol = worksetNodes[cellCol].entity().identifier() - 1;
            int globalCol = bulkData.identifier(worksetNodes[cellCol]) - 1;
            double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);
            StiffMatrix.InsertGlobalValues(1, &globalRow, 1, &globalCol, &operatorMatrixContribution);

          }// end cell col loop

        }// end cell row loop

      }// end workset cell loop
      worksetCellOrdinal++;

    } //for (size_t bucketIndex = 0; ...

  }// end workset loop

  StiffMatrix.GlobalAssemble();
  StiffMatrix.FillComplete();
  rhsVector.GlobalAssemble();

  if (MyPID==0) {
    std::cout << "Global assembly                             "
              << Time.ElapsedTime() << " seconds" << std::endl;
    Time.ResetStartTime();
  }

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/
  Epetra_Vector lhsVector(globalMapG,true);

  std::vector<int> ownedBoundaryNodes; ownedBoundaryNodes.reserve(bcNodes.size());
  // Loop over boundary nodes
  for (unsigned i = 0; i < bcNodes.size(); i++) {
    int bcNodeId = bulkData.identifier(bcNodes[i]);
    int lid = globalMapG.LID((int) bcNodeId -1);
    if(lid != -1) {
      ownedBoundaryNodes.push_back(lid);
      
      // get coordinates for this node
      entity_type bcnode = bulkData.get_entity(NODE_RANK,bcNodeId);
      double * coord = stk::mesh::field_data(*coords, bcnode);
      
      // look up exact value of function on boundary
      double x  = coord[0];
      double y  = coord[1];
      double z  = (spaceDim==3) ? coord[2] : 0;
      lhsVector[lid]=rhsVector[0][lid]=exactSolution(x, y, z);
    }
  } // end loop over boundary nodes


  // Zero out rows of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  ML_Epetra::Apply_BCsToMatrixRows(&(ownedBoundaryNodes[0]), ownedBoundaryNodes.size(), StiffMatrix);
  ML_Epetra::Remove_Zeroed_Rows(StiffMatrix,0.0);

  if(MyPID==0) {
    std::cout << "Adjust global matrix and rhs due to BCs     "
              << Time.ElapsedTime()
              << " seconds" << std::endl;
    Time.ResetStartTime();
  }

  if(matrixFilename != "") 
    EpetraExt::RowMatrixToMatlabFile(matrixFilename.c_str(),StiffMatrix);
  if(rhsFilename != "") 
    EpetraExt::MultiVectorToMatrixMarketFile(rhsFilename.c_str(),rhsVector,0,0,false);
  if(initialGuessFilename != "") 
    EpetraExt::MultiVectorToMatrixMarketFile(initialGuessFilename.c_str(),lhsVector,0,0,false);

  /**********************************************************************************/
  /*********************************** SOLVE ****************************************/
  /**********************************************************************************/

  // Run the solver
  Teuchos::ParameterList MLList = inputSolverList;
  Epetra_FEVector exactNodalVals(globalMapG);
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  for (size_t bucketIndex = 0; bucketIndex < localNodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket &nodeBucket = *localNodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * coord = stk::mesh::field_data(*coords, node);
      // look up exact value of function on boundary
      double x  = coord[0];
      double y  = coord[1];
      double z  = (spaceDim==3) ? coord[2] : 0;
      int gNodeId = bulkData.identifier(node) - 1;
      //exactNodalVals[0][gNodeId]=exactSolution(x, y, z);
      int lid = globalMapG.LID(gNodeId);
      exactNodalVals[0][lid]=exactSolution(x, y, z);
    }
  }

  exactNodalVals.GlobalAssemble();

  char probType[10] = "laplace";

  TestMultiLevelPreconditioner(probType,             MLList,
                               StiffMatrix,          exactNodalVals,
                               rhsVector,            lhsVector, nCoord,
                               TotalErrorResidual,   TotalErrorExactSol);

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

  // Patch test - tet: function is in the FE space and should be recovered
     return 1. + x + y + z ;

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

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double material[3][3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = 0.0;
      if(spaceDim==3) z = evaluationPoints(cell, pt, 2);

      materialTensor<double>(material, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        for(int col = 0; col < spaceDim; col++){
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

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

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

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

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

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double gradient[3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = 0.0;
      if(spaceDim==3)
	z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<double>(gradient, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}

/**********************************************************************************/
/******************************* TEST ML ******************************************/
/**********************************************************************************/

// Test ML
int TestMultiLevelPreconditioner(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,				 
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
				 Epetra_MultiVector & coords,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  Epetra_LinearProblem Problem(&A,&uh,&b);
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();

  Epetra_Time Time(A.Comm());

  // =================== //
  // call ML and AztecOO //
  // =================== //

  AztecOO solver(Problem);
  
  std::string precondType = "ML";
  if(MLList.isParameter("Preconditioner")) {
    precondType = MLList.get("Preconditioner",precondType);
    MLList.remove("Preconditioner");
  }

  RCP<Epetra_Operator> M;
  if(precondType == "ML") {
    ML_Epetra::SetDefaults("SA", MLList, 0, 0, false);
    M = rcp(new ML_Epetra::MultiLevelPreconditioner(A, MLList, true));
  }
#if defined(HAVE_TRILINOSCOUPLINGS_MUELU) && defined(HAVE_MUELU_EPETRA)
  else if(precondType == "MueLu") {
    // Multigrid Hierarchy
    Teuchos::ParameterList mueluParams;
    if (MLList.isSublist("MueLu"))
      mueluParams = MLList.sublist("MueLu");
    mueluParams.sublist("user data").set("Coordinates",Teuchos::rcpFromRef(coords));
    std::cout<<"*** MueLu Params ***" <<std::endl<<mueluParams<<std::endl;
    M = MueLu::CreateEpetraPreconditioner(Teuchos::rcpFromRef(A),mueluParams);
  }

#endif
  solver.SetPrecOperator(&*M);
    

  // tell AztecOO to use this preconditioner, then solve

  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(200, 1e-10);

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  double d = 0.0, d_tot = 0.0;
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);

  A.Comm().SumAll(&d,&d_tot,1);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A.Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);

  std::string msg = ProblemType;

  if (A.Comm().MyPID() == 0) {
    std::cout << msg << std::endl << "......Using " << A.Comm().NumProc() << " processes" << std::endl;
    std::cout << msg << "......||A x - b||_2 = " << Norm << std::endl;
    std::cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << std::endl;
    std::cout << msg << "......Total Time = " << Time.ElapsedTime() << std::endl << std::endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;

  return( solver.NumIters() );

}

/**********************************************************************************/
/**************************** SIMPLE BASIS FACTORY ********************************/
/**********************************************************************************/


void getBasis(Teuchos::RCP<Intrepid::Basis<double,IntrepidFieldContainer > > &basis,
               const ShardsCellTopology & cellTopology,
               int order)  {


 // select basis based on cell topology only for now, and assume first order basis
    switch (cellTopology.getKey()) {

       case shards::Tetrahedron<4>::key:
         basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_TET_C1_FEM<double, IntrepidFieldContainer > );
         break;

       case shards::Hexahedron<8>::key:
         basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, IntrepidFieldContainer > );
         break;

       case shards::Triangle<3>::key:
         basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, IntrepidFieldContainer > );
         break;

       case shards::Quadrilateral<4>::key:
         basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, IntrepidFieldContainer > );
         break;


       default:
         TEUCHOS_TEST_FOR_EXCEPTION(1,std::invalid_argument,
				    "Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4, Quadrilateral_4 or Triangle_3");

     }

}

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
