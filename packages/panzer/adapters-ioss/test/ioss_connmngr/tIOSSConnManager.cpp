// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <string>
#include <vector>
#include <map>

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

//#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_RCP.hpp"
#include <Teuchos_TestForException.hpp>

#include "PanzerCore_config.hpp"
#include "PanzerAdaptersIOSS_config.hpp"
#include "Panzer_IOSSConnManager.hpp"

#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IOSSDatabaseTypeManager.hpp"

#include "Shards_BasicTopologies.hpp"
#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "Ionit_Initializer.h"
#include "Ioss_DBUsage.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Property.h"
#include "Ioss_IOFactory.h"
#include "Ioss_DatabaseIO.h"

//#ifdef HAVE_MPI
//   #include "Epetra_MpiComm.h"
//#else
//   #include "Epetra_SerialComm.h"
//#endif

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_ioss {

namespace {

    struct Correct {

    	std::size_t numElementBlocks;
    	std::vector<shards::CellTopology> elementBlockTopologies;
    	std::vector<std::string> elementBlockIds;

        std::map<int,std::vector<double>> dofCoords;

        std::map<std::string,std::vector<int>> elementBlockMap;
        std::size_t ownedElementCount;
        std::vector<std::string> localElementBlockIds;
        std::vector<int> connectivitySize;
        std::map<int, std::vector<int>> connectivityLocalElementMap;

    };


    Correct correctQuad4NodalFP;
    Correct correctQuadHGradC2FP;
    Correct correctQuadHGradC3FP;
    Correct correctHexHGradC2FP;

    void setCorrectQuad4NodalFP(Teuchos::MpiComm<int> & comm);
    void setCorrectQuadHGradC2FP(Teuchos::MpiComm<int> & comm);
    void setCorrectQuadHGradC3FP(Teuchos::MpiComm<int> & comm);
    void setCorrectHexHGradC2FP(Teuchos::MpiComm<int> & comm);

}

TEUCHOS_UNIT_TEST(tIOSSConnManager, basic)
{

	// build global (or serial communicator)
	#ifdef HAVE_MPI
	  Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
	#else
	  // THIS_REALLY_DOES_NOT_WORK
	#endif

	const int MAX_PROCESSES = 3;

	int np = comm.getSize(); // number of processors
	TEUCHOS_TEST_FOR_EXCEPTION(np > MAX_PROCESSES, std::logic_error,
					"Error, Test was run with " << np << " processes,"
					<< "but is only valid for up to " << MAX_PROCESSES << " processes." << std::endl);

	typedef Kokkos::DynRankView<double, PHX::Device> FieldContainer;

    // For printing
	bool print = false;
    Teuchos::FancyOStream output(Teuchos::rcpFromRef(std::cout));
	output.setShowProcRank(true);
    if (print) output << std::endl;

    //Kokkos::initialize(argc, argv);

	// Create the correct versions of data structures representing the mesh
	// for comparison with the actual data structures extracted from the database.
    if (print) output << "Setting correct structures" << std::endl;
	setCorrectQuad4NodalFP(comm);
	setCorrectQuadHGradC2FP(comm);
	setCorrectQuadHGradC3FP(comm);
	setCorrectHexHGradC2FP(comm);

	// Initialize Ioss.
	if (print) output << "Initializing IOSS" << std::endl;
	Ioss::Init::Initializer();

	std::vector<std::string> filenames;
	std::vector<std::string> iossDatabaseTypes;
	std::vector<Correct *> corrects;
	std::vector<Teuchos::RCP<panzer::FieldPattern>> fieldPatterns; // use nullptr for the appropriate nodal field pattern.
	std::vector<Teuchos::RCP<panzer::Intrepid2FieldPattern>> intrepid2FieldPatterns;

	RCP<panzer::FieldPattern> fp;
	RCP<Intrepid2::Basis<PHX::Device,double,double> > basis;
	RCP<panzer::Intrepid2FieldPattern> i2fp;

	// Test 0
	filenames.push_back("rectangle.pg");
	iossDatabaseTypes.push_back("pamgen");
	corrects.push_back(&correctQuad4NodalFP);
	fp.reset(); // Set fp to null. Use the appropriate nodal field pattern.
	basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::Device,double,double>);
	i2fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	fieldPatterns.push_back(fp);
	intrepid2FieldPatterns.push_back(i2fp);

	// Test 1
	filenames.push_back("rectangle.pg");
	iossDatabaseTypes.push_back("pamgen");
	corrects.push_back(&correctQuadHGradC2FP);
	basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::Device,double,double>);
	fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	i2fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	fieldPatterns.push_back(fp);
	intrepid2FieldPatterns.push_back(i2fp);

	// Test 2
	filenames.push_back("rectangle.pg");
	iossDatabaseTypes.push_back("pamgen");
	corrects.push_back(&correctQuadHGradC3FP);
	basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(3, Intrepid2::POINTTYPE_EQUISPACED));
	fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	i2fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	fieldPatterns.push_back(fp);
	intrepid2FieldPatterns.push_back(i2fp);

	// Test 3
	filenames.push_back("brick.pg");
	iossDatabaseTypes.push_back("pamgen");
	corrects.push_back(&correctHexHGradC2FP);
	basis = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::Device,double,double>);
	fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	i2fp = rcp(new panzer::Intrepid2FieldPattern(basis));
	fieldPatterns.push_back(fp);
	intrepid2FieldPatterns.push_back(i2fp);

	size_t numTests = filenames.size();

	if (print) output << "Checking that containers defining tests all have the same sizes." << std::endl;
	TEUCHOS_TEST_FOR_EXCEPTION(filenames.size() != numTests, std::logic_error,
						"Error, filenames.size() =  " << filenames.size()
						<< "but there should be " << numTests << " tests defined." << std::endl);
	TEUCHOS_TEST_FOR_EXCEPTION(iossDatabaseTypes.size() != numTests, std::logic_error,
							"Error, iossDatabaseTypes.size() =  " << filenames.size()
							<< "but there should be " << numTests << " tests defined." << std::endl);
	TEUCHOS_TEST_FOR_EXCEPTION(corrects.size() != numTests, std::logic_error,
							"Error, corrects.size() =  " << filenames.size()
							<< "but there should be " << numTests << " tests defined." << std::endl);
	TEUCHOS_TEST_FOR_EXCEPTION(fieldPatterns.size() != numTests, std::logic_error,
							"Error, filenames.size() =  " << filenames.size()
							<< "but there should be " << numTests << " tests defined." << std::endl);
	TEUCHOS_TEST_FOR_EXCEPTION(intrepid2FieldPatterns.size() != numTests, std::logic_error,
							"Error, filenames.size() =  " << filenames.size()
							<< "but there should be " << numTests << " tests defined." << std::endl);

	for (size_t test = 0; test < numTests; ++test) {

	  if (print) output << std::endl;
	  if (print) output << "Starting test " << test << std::endl;
	  if (print) output << "---------------" << std::endl;

	  std::string filename = filenames[test];
	  std::string iossDatabaseType = iossDatabaseTypes[test];
	  Correct & correct = *(corrects[test]);
	  RCP<panzer::FieldPattern> fieldPattern = fieldPatterns[test];
	  RCP<panzer::Intrepid2FieldPattern> intrepid2FieldPattern = intrepid2FieldPatterns[test];

      // Check for a valid database type.
	  if (print) output << "Checking for a valid database type." << std::endl;
	  TEUCHOS_TEST_FOR_EXCEPTION(!IossDatabaseTypeManager::isValidType(iossDatabaseType), std::logic_error,
		  		"Error, " << iossDatabaseType
		  		<< " is not a valid IOSS database type." << std::endl
			  	<< "Valid types are " << IossDatabaseTypeManager::validTypeList());

	  // Create an empty property manager.
	  if (print) output << "Creating an empty property manager." << std::endl;
      Ioss::PropertyManager databaseProps;

      // Create the Ioss database.
      if (print) output << "Creating the Ioss database." << std::endl;
      Ioss::DatabaseIO * iossMeshDB = Ioss::IOFactory::create(iossDatabaseType, filename,
    		Ioss::READ_MODEL, MPI_COMM_WORLD, databaseProps);

      // Create the connectivity manager.
      if (print) output << "Creating the connectivity manager." << std::endl;
      IOSSConnManager connManager(iossMeshDB);

      // Create a clone
      bool cloneCreated = false;
      Teuchos::RCP<panzer::ConnManager> cmClone;
      if (IossDatabaseTypeManager::supportsMultipleOpenDatabases(iossDatabaseType)) {
        if (print) output << "Creating a clone." << std::endl;
        cmClone = connManager.noConnectivityClone(iossDatabaseType, databaseProps);
        cloneCreated = true;
      }

      // Check for correct number of element blocks.
      if (print) output << "Checking for correct number of element blocks." << std::endl;
      std::size_t numElementBlocks = connManager.numElementBlocks();
      TEST_EQUALITY(numElementBlocks, correct.numElementBlocks);
      if (cloneCreated) {
        std::size_t cloneNumElementBlocks = cmClone->numElementBlocks();
        TEST_EQUALITY(cloneNumElementBlocks, correct.numElementBlocks);
      }

      // Check for correct element block topologies.
      if (print) output << "Checking for correct element block topologies." << std::endl;
      std::vector<shards::CellTopology> elementBlockTopologies;
      connManager.getElementBlockTopologies(elementBlockTopologies);
      TEST_COMPARE_ARRAYS(elementBlockTopologies, correct.elementBlockTopologies);
      if (cloneCreated) {
        std::vector<shards::CellTopology> cloneElementBlockTopologies;
        cmClone->getElementBlockTopologies(cloneElementBlockTopologies);
        TEST_COMPARE_ARRAYS(cloneElementBlockTopologies, correct.elementBlockTopologies);
      }

      // Make field pattern
      RCP<panzer::FieldPattern> connFP;
      if (fieldPattern.is_null()) {
        if (print) output << "Creating a NodalFieldPattern." << std::endl;
        connFP = rcp(new panzer::NodalFieldPattern(elementBlockTopologies[0]));
      }
      else {
        if (print) output << "Using the specified FieldPattern." << std::endl;
        connFP = fieldPattern;
      }

      // Build the connectivity
      if (print) output << "Building the connectivity" << std::endl;
      connManager.buildConnectivity(*connFP);
      if (cloneCreated) {
        cmClone->buildConnectivity(*connFP);
      } // end if supportsMultipleOpenDatabases

      // Check for correct element block ids.
      if (print) output << "Checking for the correct element block ids." << std::endl;
      std::vector<std::string> elementBlockIds;
      connManager.getElementBlockIds(elementBlockIds);
      TEST_COMPARE_ARRAYS(elementBlockIds, correct.elementBlockIds)
      if (cloneCreated) {
        std::vector<std::string> cloneElementBlockIds;
        cmClone->getElementBlockIds(cloneElementBlockIds);
        TEST_COMPARE_ARRAYS(cloneElementBlockIds, correct.elementBlockIds)
      }

      // Check for correct local element ids for each block
      if (print) output << "Checking for correct local element ids for each block." << std::endl;
      std::vector<int> elementBlock;
      for (std::string elementBlockId : elementBlockIds) {
        elementBlock = connManager.getElementBlock(elementBlockId);
        TEST_COMPARE_ARRAYS(elementBlock, correct.elementBlockMap[elementBlockId])
        if (cloneCreated) {
          elementBlock = cmClone->getElementBlock(elementBlockId);
    	  TEST_COMPARE_ARRAYS(elementBlock, correct.elementBlockMap[elementBlockId])
        }
      }

      // Check for correct owned element count
      // Method does not exist for clone, which is of parent type.
      if (print) output << "Checking for correct owned element count." << std::endl;
      std::size_t ownedElementCount = connManager.getOwnedElementCount();
      TEST_EQUALITY(ownedElementCount, correct.ownedElementCount);


      // Check for correct block id for all local elements
      if (print) output << "Checking for correct block id for all local elements." << std::endl;
      std::vector<std::string> localElementBlockIds;
      for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
        localElementBlockIds.push_back(connManager.getBlockId(localElementId));
      }
      TEST_COMPARE_ARRAYS(localElementBlockIds, correct.localElementBlockIds);

      // Check for correct connectivity size for all local elements
      if (print) output << "Checking for correct connectivity size for all local elements." << std::endl;
      std::vector<int> connectivitySize;
      for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
        connectivitySize.push_back(connManager.getConnectivitySize(localElementId));
      }
      TEST_COMPARE_ARRAYS(connectivitySize, correct.connectivitySize);

      // Check for correct connectivity for all local elements
      if (print) output << "Checking for correct connectivity for all local elements." << std::endl;
      panzer::GlobalOrdinal * conn;
      for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
        conn = connManager.getConnectivity(localElementId);
        for (int i = 0; i < connectivitySize[localElementId]; ++i) {
          TEST_EQUALITY(conn[i], correct.connectivityLocalElementMap[localElementId][i]);
        }
      }

      // Check for correct dofCoords for all local elements
      if (print) output << "Checking for correct dofCoords for all local elements." << std::endl;
      std::vector<int> localCellIds;
      FieldContainer points;
      std::string elementBlockId;
      double tolerance = 1.0e-12;
      double coordinate;
      double correctCoordinate;
      int localCellId, subcellIndex, subcellCount, indexInConn;
      if (!intrepid2FieldPattern.is_null()) {
        for (size_t block = 0; block < numElementBlocks; ++block) {
          size_t dim = elementBlockTopologies[block].getDimension();
          elementBlockId = elementBlockIds[block];
          connManager.getDofCoords(elementBlockId, *intrepid2FieldPattern, localCellIds, points);
          for (size_t el = 0; el < localCellIds.size(); ++el) {
    	    localCellId = localCellIds[el];
            conn = connManager.getConnectivity(localCellId);
            indexInConn = 0;
            for (size_t subcellDim = 0; subcellDim <= dim; ++subcellDim) {
              subcellCount = intrepid2FieldPattern->getSubcellCount(subcellDim);
              for (int subcell = 0; subcell < subcellCount; ++subcell) {
                const std::vector<int> subcellIndices = intrepid2FieldPattern->getSubcellIndices(subcellDim, subcell);
                for (size_t id = 0; id < subcellIndices.size(); ++id, ++indexInConn) {
                  subcellIndex = subcellIndices[id];
                  for (size_t k = 0; k < dim; ++k) {
                    coordinate = points(el, subcellIndex, k);
                    correctCoordinate = correct.dofCoords[conn[indexInConn]][k];
                    TEST_FLOATING_EQUALITY(coordinate, correctCoordinate, tolerance);
                  }
                }
              }
            }
          }
        }
      }

      if (print) output << "Finished with test " << test << std::endl;

	}

    //Kokkos::finalize();
}

namespace {

  void setCorrectQuad4NodalFP(Teuchos::MpiComm<int> & comm) {

	Correct & correct = correctQuad4NodalFP;

    typedef std::pair<std::string, std::vector<int>> strVecIntPair;
    typedef std::pair<int, std::vector<int>> intVecIntPair;

	int np = comm.getSize(); // number of processors
	int rank = comm.getRank(); // processor rank

	// Correct number of element blocks
	correct.numElementBlocks = 2;

	// Correct element block topologies
	for (size_t b = 0; b < correct.numElementBlocks; ++b) {
	  correct.elementBlockTopologies.push_back(shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4>>()));
	}

	// Correct element block ids
    correct.elementBlockIds.push_back("block_1");
    correct.elementBlockIds.push_back("block_2");

    // Correct dof coordinates for all nodes
    int rows = 3;
    int columns = 11;
    int dim = 2;
    std::vector<double> dofCoordsThisNode(dim);
    int nd = 1;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < columns; ++j, ++nd) {
        dofCoordsThisNode[0] = 1.0 * j;
        dofCoordsThisNode[1] = 1.0 * i;
        correct.dofCoords[nd] = dofCoordsThisNode;
      }
    }

    // Correct owned element count, local element ids,
    // and block ids, connectivity sizes, and connectivities for local elements.

    const int NODES = 4;
    const int ELEMENTS = 20;
    int conn[ELEMENTS][NODES] =
                      {{1, 2, 13, 12}, {2, 3, 14, 13}, {3, 4, 15, 14}, {4, 5, 16, 15},
                       {5, 6, 17, 16}, {12, 13, 24, 23}, {13, 14, 25, 24}, {14, 15, 26, 25},
                       {15, 16, 27, 26}, {16, 17, 28, 27}, {6, 7, 18, 17}, {7, 8, 19, 18},
					   {8, 9, 20, 19}, {9, 10, 21, 20}, {10, 11, 22, 21}, {17, 18, 29, 28},
					   {18, 19, 30, 29}, {19, 20, 31, 30}, {20, 21, 32, 31}, {21, 22, 33, 32}};

    std::vector<int> connectivity;
    std::map<int, std::vector<int>> connectivityGlobalElementMap;
    int count = 1;
    for (int element = 0; element < ELEMENTS; ++element) {
      for (int index = 0; index < NODES; ++index) {
        connectivity.push_back(conn[element][index]);
      }
      connectivityGlobalElementMap.insert(intVecIntPair(count, connectivity));
      connectivity.clear();
      ++count;
    }

    std::vector<int> elementBlock;
    std::vector<int> elementLidToGid;

	if (np == 1) {
	  correct.ownedElementCount = 20;
	  correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, NODES);


	  // block_1
	  elementBlock.clear();
	  for (int i = 0; i <= 9; ++i) {
	    elementBlock.push_back(i);
	    correct.localElementBlockIds.push_back(correct.elementBlockIds[0]);
	  }
	  correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[0], elementBlock));
	  //block_2
	  elementBlock.clear();
	  for (int i = 10; i <= 19; ++i) {
	    elementBlock.push_back(i);
	    correct.localElementBlockIds.push_back(correct.elementBlockIds[1]);
	  }
	  correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[1], elementBlock));


	  int elLidToGid[20] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
	                        11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

	  for (size_t i = 0; i < correct.ownedElementCount; ++i) {
		elementLidToGid.push_back(elLidToGid[i]);
	    correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elementLidToGid[i]]));
	  }

	}
	else if (np == 2) {

	  if (rank == 0) {
		correct.ownedElementCount = 10;
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, NODES);

	    // block_1
	    elementBlock.clear();
	    for (int i = 0; i <= 4; ++i) {
	      elementBlock.push_back(i);
	      correct.localElementBlockIds.push_back(correct.elementBlockIds[0]);
	    }
	    correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[0], elementBlock));
	    //block_2
	    elementBlock.clear();
	    for (int i = 5; i <= 9; ++i) {
	      elementBlock.push_back(i);
	      correct.localElementBlockIds.push_back(correct.elementBlockIds[1]);
	    }
	    correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[1], elementBlock));

	    int elLidToGid[10] = {1, 2, 3, 4, 5, 11, 12, 13, 14, 15};

	    for (size_t i = 0; i < correct.ownedElementCount; ++i) {
	      elementLidToGid.push_back(elLidToGid[i]);
	      correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elementLidToGid[i]]));
	    }

	  }
	  else {
		correct.ownedElementCount = 10;
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, NODES);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <= 4; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[0]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[0], elementBlock));
		//block_2
		elementBlock.clear();
		for (int i = 5; i <= 9; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[1]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[1], elementBlock));

		int elLidToGid[10] = {6, 7, 8, 9, 10, 16, 17, 18, 19, 20};

		for (size_t i = 0; i < correct.ownedElementCount; ++i) {
	      elementLidToGid.push_back(elLidToGid[i]);
		  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elementLidToGid[i]]));
	    }

	  }
	}
	else {
	  if (rank == 0) {
		correct.ownedElementCount = 6;
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, NODES);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <=4; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[0]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[0], elementBlock));
		// block_2
		elementBlock.clear();
		for (int i = 5; i <=5; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[1]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[1], elementBlock));

		int elLidToGid[6] = {1, 2, 3, 4, 5, 11};

		for (size_t i = 0; i < correct.ownedElementCount; ++i) {
		  elementLidToGid.push_back(elLidToGid[i]);
	      correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elementLidToGid[i]]));
	    }

	  }
	  else if (rank == 1) {
		correct.ownedElementCount = 6;
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, NODES);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <= 1; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[0]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[0], elementBlock));
		//block_2
		elementBlock.clear();
		for (int i = 2; i <= 5; ++i) {
	      elementBlock.push_back(i);
	      correct.localElementBlockIds.push_back(correct.elementBlockIds[1]);
		}
	    correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[1], elementBlock));

	    int elLidToGid[6] = {6, 7, 12, 13, 14, 15};

	    for (size_t i = 0; i < correct.ownedElementCount; ++i) {
	      elementLidToGid.push_back(elLidToGid[i]);
	      correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elementLidToGid[i]]));
	    }

	  }
	  else {
		correct.ownedElementCount = 8;
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, NODES);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <= 2; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[0]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[0], elementBlock));
		//block_2
		elementBlock.clear();
		for (int i = 3; i <= 7; ++i) {
		  elementBlock.push_back(i);
		  correct.localElementBlockIds.push_back(correct.elementBlockIds[1]);
		}
		correct.elementBlockMap.insert(strVecIntPair(correct.elementBlockIds[1], elementBlock));

		int elLidToGid[8] = {8, 9, 10, 16, 17, 18, 19, 20};

		for (size_t i = 0; i < correct.ownedElementCount; ++i) {
	      elementLidToGid.push_back(elLidToGid[i]);
	      correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elementLidToGid[i]]));
	    }

	  }
	}

	return;

  } // end setCorrectQuad4NodalFP()


  void setCorrectQuadHGradC2FP(Teuchos::MpiComm<int> & comm) {

	typedef std::pair<int, std::vector<int>> intVecIntPair;

	int np = comm.getSize(); // number of processors
	int rank = comm.getRank(); // processor rank

	correctQuadHGradC2FP = correctQuad4NodalFP; // copy all members

  	Correct & correct = correctQuadHGradC2FP;

  	// Clear containers, otherwise pair insertions do not overwrite Quad4 Nodal values
  	correct.connectivitySize.clear();
  	correct.connectivityLocalElementMap.clear();
  	correct.dofCoords.clear();

  	const int NODES = 4;
  	const int CORNER_NODES = 4;

  	const int EDGES = 4;
  	const int IDS_PER_EDGE = 1;
  	const int EDGE_IDS = EDGES * IDS_PER_EDGE;

  	const int IDS_PER_CELL = 1;
  	const int CELL_IDS = IDS_PER_CELL;

  	const int CONNECTIVITY_SIZE = CORNER_NODES + EDGE_IDS + CELL_IDS;
  	const int ELEMENTS = 20;

  	const int EDGE_OFFSET = 34;
  	const int CELL_OFFSET = 86;

  	// Correct connectivity size
    correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, CONNECTIVITY_SIZE);

    // Correct connectivity
    int nodes[ELEMENTS][NODES] =
                      {{1, 2, 13, 12}, {2, 3, 14, 13}, {3, 4, 15, 14}, {4, 5, 16, 15},
                       {5, 6, 17, 16}, {12, 13, 24, 23}, {13, 14, 25, 24}, {14, 15, 26, 25},
                       {15, 16, 27, 26}, {16, 17, 28, 27}, {6, 7, 18, 17}, {7, 8, 19, 18},
					   {8, 9, 20, 19}, {9, 10, 21, 20}, {10, 11, 22, 21}, {17, 18, 29, 28},
					   {18, 19, 30, 29}, {19, 20, 31, 30}, {20, 21, 32, 31}, {21, 22, 33, 32}};

    int edges_np_1[ELEMENTS][EDGES] =
                      {{0, 3, 21, 1}, {2, 5, 23, 3}, {4, 7, 25, 5}, {6, 9, 27, 7},
                       {8, 11, 29, 9}, {21, 24, 42, 22}, {23, 26, 43, 24}, {25, 28, 44, 26},
                       {27, 30, 45, 28}, {29, 32, 46, 30}, {10, 13, 31, 11}, {12, 15, 33, 13},
                       {14, 17, 35, 15}, {16, 19, 37, 17}, {18, 20, 39, 19}, {31, 34, 47, 32},
                       {33, 36, 48, 34}, {35, 38, 49, 36}, {37, 40, 50, 38}, {39, 41, 51, 40}};

    int edges_np_2[ELEMENTS][EDGES] =
                      {{25, 27, 11, 26}, {0, 28, 12, 27}, {1, 2, 13, 28}, {29, 31, 35, 2},
					   {30, 3, 14, 31}, {11, 18, 44, 17}, {12, 19, 45, 18}, {13, 39, 46, 19},
					   {35, 40, 47, 39}, {14, 20, 48, 40}, {32, 5, 36, 3}, {4, 7, 37, 5},
					   {6, 33, 15, 7}, {8, 9, 16, 33}, {34, 10, 38, 9}, {36, 41, 49, 20},
					   {37, 21, 23, 41}, {15, 42, 50, 21}, {16, 43, 51, 42}, {38, 22, 24, 43}};

	int edges_np_3[ELEMENTS][EDGES] =
					  {{38, 18, 40, 17}, {0, 19, 5, 18}, {1, 3, 23, 19}, {2, 20, 6, 3},
		               {39, 4, 24, 20}, {40, 44, 46, 43}, {5, 45, 47, 44}, {23, 32, 35, 45},
                       {6, 33, 15, 32}, {24, 11, 49, 33}, {21, 22, 25, 4}, {7, 27, 30, 22},
			           {26, 8, 9, 27}, {28, 41, 10, 8}, {29, 42, 31, 41}, {25, 12, 50, 11},
					   {30, 13, 36, 12}, {9, 34, 16, 13}, {10, 14, 37, 34}, {31, 48, 51, 14}};


	// This is needed because the global ids output by Zoltan2::findUniqueGids() depend on the number of processes.
	// It would be nice if this method could be made to be independent of number of processes.
	int edges[ELEMENTS][EDGES];
	for (int element = 0; element < ELEMENTS; ++element) {
	  for (int edge = 0; edge < EDGES; ++edge) {
		switch (np) {
		  case 1:
	        edges[element][edge] = edges_np_1[element][edge];
	        break;
		  case 2:
		    edges[element][edge] = edges_np_2[element][edge];
		    break;
		  case 3:
			edges[element][edge] = edges_np_3[element][edge];
			break;
		  default:
			edges[element][edge] = -1;
		}
	  }
	}

    std::vector<int> connectivity;
    std::map<int, std::vector<int>> connectivityGlobalElementMap;
    int count = 1;

    for (int element = 0; element < ELEMENTS; ++element) {
      for (int corner_node = 0; corner_node < CORNER_NODES; ++corner_node) {
        connectivity.push_back(nodes[element][corner_node]);
      }
      for (int edge = 0; edge < EDGES; ++edge) {
        for (int id = 0; id < IDS_PER_EDGE; ++id) {
          connectivity.push_back(EDGE_OFFSET + edges[element][edge] * IDS_PER_EDGE + id);
        }
      }
      for (int id = 0; id < IDS_PER_CELL; ++id) {
        connectivity.push_back(CELL_OFFSET + (1+element) * IDS_PER_CELL + id);
      }
      connectivityGlobalElementMap.insert(intVecIntPair(count, connectivity));
      connectivity.clear();
      ++count;
    }

    std::vector<int> elementLidToGid;


    if (np == 1) {
      int elLidToGid[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    		                11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  	  for (size_t i = 0; i < correct.ownedElementCount; ++i) {
  		elementLidToGid.push_back(elLidToGid[i]);
  	    correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
  	  }
    }

    else if (np == 2) {
      if (rank == 0) {
        int elLidToGid[10] = {1, 2, 3, 4, 5, 11, 12, 13, 14, 15};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
          elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else {
        int elLidToGid[10] = {6, 7, 8, 9, 10, 16, 17, 18, 19, 20};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }

    }

    else {
      if (rank == 0) {
        int elLidToGid[6] = {1, 2, 3, 4, 5, 11};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else if (rank == 1) {
        int elLidToGid[6] = {6, 7, 12, 13, 14, 15};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
          elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else {
        int elLidToGid[8] = {8, 9, 10, 16, 17, 18, 19, 20};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
    }

    // Correct dof coordinates for all nodes
    int dim = 2;
    std::vector<double> dofCoordsThisLoc(dim);
    std::vector<double> dofCoordsThisCenter(dim);

    double eLen = 1.0;
    double nodeDeltas[ELEMENTS][2] = {{-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}};
    double edgeIdDeltas[EDGE_IDS][2] = {{0.0, -0.5}, {0.5, 0.0}, {0.0, 0.5}, {-0.5, 0.0}};
    double cellIdDeltas[IDS_PER_CELL][2] = {{0., 0.}};
    double elementCorners[ELEMENTS][2] = {{0., 0.}, {1., 0.}, {2., 0.}, {3., 0.}, {4., 0.},
    		                              {0., 1.}, {1., 1.}, {2., 1.}, {3., 1.}, {4., 1.},
                                          {5., 0.}, {6., 0.}, {7., 0.}, {8., 0.}, {9., 0.},
										  {5., 1.}, {6., 1.}, {7., 1.}, {8., 1.}, {9., 1.}};

    for (int element = 0; element < ELEMENTS; ++element) {
      for (int k = 0; k < dim; ++k) {
    	nodeDeltas[element][k] *= eLen;
    	edgeIdDeltas[element][k] *= eLen;
    	cellIdDeltas[element][k] *= eLen;
        elementCorners[element][k] *= eLen;
      }
    }

    for (int element = 0; element < ELEMENTS; ++element) {
      dofCoordsThisCenter[0] = elementCorners[element][0] + eLen/2.0;
      dofCoordsThisCenter[1] = elementCorners[element][1] + eLen/2.0;
      for (int cornerNode = 0; cornerNode < CORNER_NODES; ++cornerNode) {
    	for (int k = 0; k < dim; ++k) {
          dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + nodeDeltas[cornerNode][k];
    	}
        correct.dofCoords[nodes[element][cornerNode]] = dofCoordsThisLoc;
      }
      for (int edge = 0; edge < EDGES; ++edge) {
        for (int i = 0; i < IDS_PER_EDGE; ++i) {
    	  for (int k = 0; k < dim; ++k) {
            dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + edgeIdDeltas[edge*IDS_PER_EDGE+i][k];
    	  }
          correct.dofCoords[EDGE_OFFSET + IDS_PER_EDGE * edges[element][edge] + i] = dofCoordsThisLoc;
        }
      }
      for (int i = 0; i < IDS_PER_CELL; ++i) {
        for (int k = 0; k < dim; ++k) {
          dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + cellIdDeltas[i][k];
        }
        correct.dofCoords[CELL_OFFSET + IDS_PER_CELL * (element+1) + i] = dofCoordsThisLoc;
      }
    }

  }

  void setCorrectQuadHGradC3FP(Teuchos::MpiComm<int> & comm) {

	typedef std::pair<int, std::vector<int>> intVecIntPair;

	int np = comm.getSize(); // number of processors
	int rank = comm.getRank(); // processor rank

	correctQuadHGradC3FP = correctQuad4NodalFP; // copy all members

  	Correct & correct = correctQuadHGradC3FP;

  	// Clear containers, otherwise pair insertions do not overwrite Quad4 Nodal values
  	correct.connectivitySize.clear();
  	correct.connectivityLocalElementMap.clear();
  	correct.dofCoords.clear();

  	const int NODES = 4;
  	const int CORNER_NODES = 4;

  	const int EDGES = 4;
  	const int IDS_PER_EDGE = 2;
  	const int EDGE_IDS = EDGES * IDS_PER_EDGE;

  	const int IDS_PER_CELL = 4;
  	const int CELL_IDS = IDS_PER_CELL;

  	const int CONNECTIVITY_SIZE = CORNER_NODES + EDGE_IDS + CELL_IDS;
  	const int ELEMENTS = 20;

  	const int EDGE_OFFSET = 34;
  	const int CELL_OFFSET = 138;

    correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, CONNECTIVITY_SIZE);

    int nodes[ELEMENTS][NODES] =
                      {{1, 2, 13, 12}, {2, 3, 14, 13}, {3, 4, 15, 14}, {4, 5, 16, 15},
                       {5, 6, 17, 16}, {12, 13, 24, 23}, {13, 14, 25, 24}, {14, 15, 26, 25},
                       {15, 16, 27, 26}, {16, 17, 28, 27}, {6, 7, 18, 17}, {7, 8, 19, 18},
					   {8, 9, 20, 19}, {9, 10, 21, 20}, {10, 11, 22, 21}, {17, 18, 29, 28},
					   {18, 19, 30, 29}, {19, 20, 31, 30}, {20, 21, 32, 31}, {21, 22, 33, 32}};

    int edges_np_1[ELEMENTS][EDGES] =
                      {{0, 3, 21, 1}, {2, 5, 23, 3}, {4, 7, 25, 5}, {6, 9, 27, 7},
                       {8, 11, 29, 9}, {21, 24, 42, 22}, {23, 26, 43, 24}, {25, 28, 44, 26},
                       {27, 30, 45, 28}, {29, 32, 46, 30}, {10, 13, 31, 11}, {12, 15, 33, 13},
                       {14, 17, 35, 15}, {16, 19, 37, 17}, {18, 20, 39, 19}, {31, 34, 47, 32},
                       {33, 36, 48, 34}, {35, 38, 49, 36}, {37, 40, 50, 38}, {39, 41, 51, 40}};

    int edges_np_2[ELEMENTS][EDGES] =
                      {{25, 27, 11, 26}, {0, 28, 12, 27}, {1, 2, 13, 28}, {29, 31, 35, 2},
					   {30, 3, 14, 31}, {11, 18, 44, 17}, {12, 19, 45, 18}, {13, 39, 46, 19},
					   {35, 40, 47, 39}, {14, 20, 48, 40}, {32, 5, 36, 3}, {4, 7, 37, 5},
					   {6, 33, 15, 7}, {8, 9, 16, 33}, {34, 10, 38, 9}, {36, 41, 49, 20},
					   {37, 21, 23, 41}, {15, 42, 50, 21}, {16, 43, 51, 42}, {38, 22, 24, 43}};

	int edges_np_3[ELEMENTS][EDGES] =
					  {{38, 18, 40, 17}, {0, 19, 5, 18}, {1, 3, 23, 19}, {2, 20, 6, 3},
		               {39, 4, 24, 20}, {40, 44, 46, 43}, {5, 45, 47, 44}, {23, 32, 35, 45},
                       {6, 33, 15, 32}, {24, 11, 49, 33}, {21, 22, 25, 4}, {7, 27, 30, 22},
			           {26, 8, 9, 27}, {28, 41, 10, 8}, {29, 42, 31, 41}, {25, 12, 50, 11},
					   {30, 13, 36, 12}, {9, 34, 16, 13}, {10, 14, 37, 34}, {31, 48, 51, 14}};


	// This is needed because the global ids output by Zoltan2::findUniqueGids() depend on the number of processes.
	// It would be nice if this method could be made to be independent of number of processes.
	int edges[ELEMENTS][EDGES];
	for (int element = 0; element < ELEMENTS; ++element) {
	  for (int edge = 0; edge < EDGES; ++edge) {
		switch (np) {
		  case 1:
	        edges[element][edge] = edges_np_1[element][edge];
	        break;
		  case 2:
		    edges[element][edge] = edges_np_2[element][edge];
		    break;
		  case 3:
			edges[element][edge] = edges_np_3[element][edge];
			break;
		  default:
			edges[element][edge] = -1;
		}
	  }
	}

    std::vector<int> connectivity;
    std::map<int, std::vector<int>> connectivityGlobalElementMap;
    int count = 1;

    for (int element = 0; element < ELEMENTS; ++element) {
      for (int corner_node = 0; corner_node < CORNER_NODES; ++corner_node) {
        connectivity.push_back(nodes[element][corner_node]);
      }
      for (int edge = 0; edge < EDGES; ++edge) {
        for (int id = 0; id < IDS_PER_EDGE; ++id) {
          connectivity.push_back(EDGE_OFFSET + edges[element][edge] * IDS_PER_EDGE + id);
        }
      }
      for (int id = 0; id < IDS_PER_CELL; ++id) {
        connectivity.push_back(CELL_OFFSET + (1+element) * IDS_PER_CELL + id);
      }
      connectivityGlobalElementMap.insert(intVecIntPair(count, connectivity));
      connectivity.clear();
      ++count;
    }

    std::vector<int> elementLidToGid;


    if (np == 1) {
      int elLidToGid[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    		                11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  	  for (size_t i = 0; i < correct.ownedElementCount; ++i) {
  		elementLidToGid.push_back(elLidToGid[i]);
  	    correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
  	  }
    }

    else if (np == 2) {
      if (rank == 0) {
        int elLidToGid[10] = {1, 2, 3, 4, 5, 11, 12, 13, 14, 15};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
          elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else {
        int elLidToGid[10] = {6, 7, 8, 9, 10, 16, 17, 18, 19, 20};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }

    }

    else {
      if (rank == 0) {
        int elLidToGid[6] = {1, 2, 3, 4, 5, 11};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else if (rank == 1) {
        int elLidToGid[6] = {6, 7, 12, 13, 14, 15};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
          elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else {
        int elLidToGid[8] = {8, 9, 10, 16, 17, 18, 19, 20};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
    }

    // Correct dof coordinates for all nodes
    int dim = 2;
    std::vector<double> dofCoordsThisLoc(dim);
    std::vector<double> dofCoordsThisCenter(dim);

    double eLen = 1.0;
    double nodeDeltas[ELEMENTS][2] = {{-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}};
    double edgeIdDeltas[EDGE_IDS][2] = {{-1.0/6, -0.5}, {1.0/6, -0.5},
    		                            {0.5, -1.0/6}, {0.5, 1.0/6},
    										{-1.0/6, 0.5}, {1.0/6, 0.5},
										{-0.5, -1.0/6}, {-0.5, 1.0/6}};
    double cellIdDeltas[IDS_PER_CELL][2] = {{-1.0/6, -1.0/6}, {1.0/6, -1.0/6}, {-1.0/6, 1.0/6}, {1.0/6, 1.0/6}};

    double elementCorners[ELEMENTS][2] = {{0., 0.}, {1., 0.}, {2., 0.}, {3., 0.}, {4., 0.},
    		                              {0., 1.}, {1., 1.}, {2., 1.}, {3., 1.}, {4., 1.},
                                          {5., 0.}, {6., 0.}, {7., 0.}, {8., 0.}, {9., 0.},
										  {5., 1.}, {6., 1.}, {7., 1.}, {8., 1.}, {9., 1.}};

    for (int element = 0; element < ELEMENTS; ++element) {
      for (int k = 0; k < dim; ++k) {
    	nodeDeltas[element][k] *= eLen;
    	edgeIdDeltas[element][k] *= eLen;
    	cellIdDeltas[element][k] *= eLen;
        elementCorners[element][k] *= eLen;
      }
    }

    for (int element = 0; element < ELEMENTS; ++element) {
      dofCoordsThisCenter[0] = elementCorners[element][0] + eLen/2.0;
      dofCoordsThisCenter[1] = elementCorners[element][1] + eLen/2.0;
      for (int cornerNode = 0; cornerNode < CORNER_NODES; ++cornerNode) {
    	for (int k = 0; k < dim; ++k) {
          dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + nodeDeltas[cornerNode][k];
    	}
        correct.dofCoords[nodes[element][cornerNode]] = dofCoordsThisLoc;
      }
      for (int edge = 0; edge < EDGES; ++edge) {
        for (int i = 0; i < IDS_PER_EDGE; ++i) {
    	  for (int k = 0; k < dim; ++k) {
            dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + edgeIdDeltas[edge*IDS_PER_EDGE+i][k];
    	  }
          correct.dofCoords[EDGE_OFFSET + IDS_PER_EDGE * edges[element][edge] + i] = dofCoordsThisLoc;
        }
      }
      for (int i = 0; i < IDS_PER_CELL; ++i) {
        for (int k = 0; k < dim; ++k) {
          dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + cellIdDeltas[i][k];
        }
        correct.dofCoords[CELL_OFFSET + IDS_PER_CELL * (element+1) + i] = dofCoordsThisLoc;
      }
    }

  }


  void setCorrectHexHGradC2FP(Teuchos::MpiComm<int> & comm) {

	typedef std::pair<int, std::vector<int>> intVecIntPair;

	int np = comm.getSize(); // number of processors
	int rank = comm.getRank(); // processor rank

	correctHexHGradC2FP = correctQuadHGradC2FP; // copy all members

  	Correct & correct = correctHexHGradC2FP;

  	// Clear containers, otherwise pair insertions do not overwrite Quad4 Grad C2 values
  	correct.elementBlockTopologies.clear();
  	correct.connectivitySize.clear();
  	correct.connectivityLocalElementMap.clear();
  	correct.dofCoords.clear();

	// Correct element block topologies
	for (size_t b = 0; b < correct.numElementBlocks; ++b) {
	  correct.elementBlockTopologies.push_back(shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8>>()));
	}

  	const int NODES = 8;
  	const int CORNER_NODES = 8;

  	const int EDGES = 12;
  	const int IDS_PER_EDGE = 1;
  	const int EDGE_IDS = EDGES * IDS_PER_EDGE;

  	const int FACES = 6;
  	const int IDS_PER_FACE = 1;
  	const int FACE_IDS = FACES * IDS_PER_FACE;

  	const int IDS_PER_CELL = 1;
  	const int CELL_IDS = IDS_PER_CELL;

  	const int CONNECTIVITY_SIZE = CORNER_NODES + EDGE_IDS + FACE_IDS + CELL_IDS;
  	const int ELEMENTS = 20;

  	const int EDGE_OFFSET = 67;
  	const int FACE_OFFSET = 204;
  	const int CELL_OFFSET = 296;

  	// Correct connectivity size
    correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, CONNECTIVITY_SIZE);

    // Correct connectivity
    int nodes[ELEMENTS][NODES] = {{1, 2, 13, 12, 34, 35, 46, 45}, {2, 3, 14, 13, 35, 36, 47, 46},
    		                      {3, 4, 15, 14, 36, 37, 48, 47}, {4, 5, 16, 15, 37, 38, 49, 48},
								  {5, 6, 17, 16, 38, 39, 50, 49}, {12, 13, 24, 23, 45, 46, 57, 56},
								  {13, 14, 25, 24, 46, 47, 58, 57}, {14, 15, 26, 25, 47, 48, 59, 58},
								  {15, 16, 27, 26, 48, 49, 60, 59}, {16, 17, 28, 27, 49, 50, 61, 60},
								  {6, 7, 18, 17, 39, 40, 51, 50}, {7, 8, 19, 18, 40, 41, 52, 51},
								  {8, 9, 20, 19, 41, 42, 53, 52}, {9, 10, 21, 20, 42, 43, 54, 53},
								  {10, 11, 22, 21, 43, 44, 55, 54}, {17, 18, 29, 28, 50, 51, 62, 61},
								  {18, 19, 30, 29, 51, 52, 63, 62}, {19, 20, 31, 30, 52, 53, 64, 63},
								  {20, 21, 32, 31, 53, 54, 65, 64}, {21, 22, 33, 32, 54, 55, 66, 65}};

    int edges_np_1[ELEMENTS][EDGES] = {{0, 4, 32, 1, 85, 88, 106, 86, 2, 5, 37, 34}, {3, 7, 35, 4, 87, 90, 108, 88, 5, 8, 40, 37},
    		                           {6, 10, 38, 7, 89, 92, 110, 90, 8, 11, 43, 40}, {9, 13, 41, 10, 91, 94, 112, 92, 11, 14, 46, 43},
						        	   {12, 16, 44, 13, 93, 96, 114, 94, 14, 17, 49, 46}, {32, 36, 64, 33, 106, 109, 127, 107, 34, 37, 67, 65},
								       {35, 39, 66, 36, 108, 111, 128, 109, 37, 40, 69, 67}, {38, 42, 68, 39, 110, 113, 129, 111, 40, 43, 71, 69},
								       {41, 45, 70, 42, 112, 115, 130, 113, 43, 46, 73, 71}, {44, 48, 72, 45, 114, 117, 131, 115, 46, 49, 75, 73},
								       {15, 19, 47, 16, 95, 98, 116, 96, 17, 20, 52, 49}, {18, 22, 50, 19, 97, 100, 118, 98, 20, 23, 55, 52},
								       {21, 25, 53, 22, 99, 102, 120, 100, 23, 26, 58, 55}, {24, 28, 56, 25, 101, 104, 122, 102, 26, 29, 61, 58},
								       {27, 30, 59, 28, 103, 105, 124, 104, 29, 31, 63, 61}, {47, 51, 74, 48, 116, 119, 132, 117, 49, 52, 77, 75},
								       {50, 54, 76, 51, 118, 121, 133, 119, 52, 55, 79, 77}, {53, 57, 78, 54, 120, 123, 134, 121, 55, 58, 81, 79},
								       {56, 60, 80, 57, 122, 125, 135, 123, 58, 61, 83, 81}, {59, 62, 82, 60, 124, 126, 136, 125, 61, 63, 84, 83}};

    int edges_np_2[ELEMENTS][EDGES] = {{69, 72, 16, 70, 27, 96, 37, 28, 71, 73, 18, 85}, {0, 74, 17, 72, 95, 98, 38, 96, 73, 2, 86, 18},
    		                           {1, 3, 19, 74, 97, 100, 39, 98, 2, 4, 88, 86}, {75, 77, 87, 3, 99, 29, 106, 100, 4, 5, 21, 88},
						        	   {76, 6, 20, 77, 101, 103, 107, 29, 5, 79, 22, 21}, {16, 42, 117, 41, 37, 55, 62, 54, 85, 18, 119, 47},
								       {17, 43, 118, 42, 38, 131, 63, 55, 18, 86, 121, 119}, {19, 112, 120, 43, 39, 56, 134, 131, 86, 88, 48, 121},
								       {87, 113, 122, 112, 106, 57, 135, 56, 88, 21, 49, 48}, {20, 44, 123, 113, 107, 58, 64, 57, 21, 22, 125, 49},
								       {78, 8, 89, 6, 102, 31, 108, 103, 79, 9, 91, 22}, {7, 11, 90, 8, 30, 32, 109, 31, 9, 80, 92, 91},
								       {10, 81, 23, 11, 104, 34, 110, 32, 80, 13, 25, 92}, {12, 14, 24, 81, 33, 36, 111, 34, 13, 83, 26, 25},
								       {82, 15, 93, 14, 35, 105, 40, 36, 83, 84, 94, 26}, {89, 114, 124, 44, 108, 132, 65, 58, 22, 91, 126, 125},
								       {90, 45, 50, 114, 109, 59, 66, 132, 91, 92, 128, 126}, {23, 115, 127, 45, 110, 60, 67, 59, 92, 25, 51, 128},
								       {24, 116, 129, 115, 111, 61, 68, 60, 25, 26, 53, 51}, {93, 46, 52, 116, 40, 133, 136, 61, 26, 94, 130, 53}};

    int edges_np_3[ELEMENTS][EDGES] = {{94, 48, 100, 47, 102, 62, 16, 10, 0, 49, 56, 55}, {1, 50, 6, 48, 11, 103, 17, 62, 49, 95, 58, 56},
    		                           {2, 4, 57, 50, 12, 14, 18, 103, 95, 96, 8, 58}, {3, 51, 7, 4, 13, 105, 64, 14, 96, 52, 60, 8},
						        	   {97, 5, 59, 51, 104, 63, 65, 105, 52, 98, 101, 60}, {100, 110, 114, 109, 16, 32, 33, 31, 55, 56, 29, 28},
								       {6, 111, 115, 110, 17, 78, 34, 32, 56, 58, 74, 29}, {57, 82, 85, 111, 18, 129, 45, 78, 58, 8, 122, 74},
								       {7, 83, 39, 82, 64, 130, 133, 129, 8, 60, 40, 122}, {59, 35, 123, 83, 65, 88, 134, 130, 60, 101, 41, 40},
								       {53, 54, 61, 5, 106, 15, 19, 63, 98, 99, 9, 101}, {20, 67, 71, 54, 30, 75, 79, 15, 99, 68, 72, 9},
								       {66, 21, 25, 67, 116, 118, 80, 75, 68, 22, 27, 72}, {69, 107, 26, 21, 117, 77, 81, 118, 22, 23, 112, 27},
								       {70, 108, 73, 107, 76, 119, 120, 77, 23, 24, 113, 112}, {61, 36, 124, 35, 19, 89, 92, 88, 101, 9, 125, 41},
								       {71, 37, 86, 36, 79, 131, 135, 89, 9, 72, 126, 125}, {25, 84, 42, 37, 80, 90, 136, 131, 72, 27, 43, 126},
								       {26, 38, 87, 84, 81, 132, 93, 90, 27, 112, 128, 43}, {73, 121, 127, 38, 120, 91, 46, 132, 112, 113, 44, 128}};

    int faces_np_1[ELEMENTS][FACES] = {{1, 5, 32, 2, 0, 72}, {4, 8, 35, 5, 3, 73},
    		                           {7, 11, 38, 8, 6, 74}, {10, 14, 41, 11, 9, 75},
						        	   {13, 17, 44, 14, 12, 76}, {32, 36, 62, 33, 31, 82},
								       {35, 39, 63, 36, 34, 83}, {38, 42, 64, 39, 37, 84},
								       {41, 45, 65, 42, 40, 85}, {44, 48, 66, 45, 43, 86},
								       {16, 20, 47, 17, 15, 77}, {19, 23, 50, 20, 18, 78},
								       {22, 26, 53, 23, 21, 79}, {25, 29, 56, 26, 24, 80},
								       {28, 30, 59, 29, 27, 81}, {47, 51, 67, 48, 46, 87},
								       {50, 54, 68, 51, 49, 88}, {53, 57, 69, 54, 52, 89},
								       {56, 60, 70, 57, 55, 90}, {59, 61, 71, 60, 58, 91}};

    int faces_np_2[ELEMENTS][FACES] = {{44, 47, 15, 45, 43, 66}, {0, 2, 59, 47, 46, 67},
    		                           {1, 51, 60, 2, 48, 18}, {50, 5, 61, 51, 49, 19},
						        	   {4, 8, 62, 5, 3, 20}, {15, 72, 33, 70, 24, 87},
								       {59, 26, 34, 72, 71, 38}, {60, 28, 35, 26, 25, 39},
								       {61, 29, 82, 28, 27, 40}, {62, 30, 83, 29, 73, 88},
								       {7, 53, 16, 8, 6, 68}, {9, 11, 63, 53, 52, 21},
								       {10, 12, 17, 11, 54, 69}, {56, 14, 64, 12, 55, 22},
								       {57, 58, 65, 14, 13, 23}, {16, 31, 84, 30, 74, 41},
								       {63, 77, 85, 31, 75, 89}, {17, 32, 36, 77, 76, 42},
								       {64, 80, 86, 32, 78, 90}, {65, 81, 37, 80, 79, 91}};

    int faces_np_3[ELEMENTS][FACES] = {{57, 59, 34, 0, 56, 68}, {58, 61, 65, 59, 1, 5},
    		                           {28, 62, 66, 61, 60, 6}, {2, 31, 4, 62, 29, 36},
						        	   {3, 32, 35, 31, 30, 7}, {34, 40, 15, 12, 39, 42},
								       {65, 77, 41, 40, 76, 18}, {66, 44, 47, 77, 82, 25},
								       {4, 84, 88, 44, 43, 52}, {35, 19, 48, 84, 83, 53},
								       {64, 33, 67, 32, 63, 8}, {9, 70, 78, 33, 69, 80},
								       {10, 11, 13, 70, 37, 81}, {72, 74, 14, 11, 71, 16},
								       {73, 75, 79, 74, 38, 17}, {67, 86, 89, 19, 45, 26},
								       {78, 21, 49, 86, 85, 27}, {13, 46, 50, 21, 20, 91},
								       {14, 23, 90, 46, 22, 54}, {79, 24, 51, 23, 87, 55}};

	// This is needed because the global ids output by Zoltan2::findUniqueGids() depend on the number of processes.
	// It would be nice if this method could be made to be independent of number of processes.
	int edges[ELEMENTS][EDGES];
	int faces[ELEMENTS][FACES];
	for (int element = 0; element < ELEMENTS; ++element) {
	  for (int edge = 0; edge < EDGES; ++edge) {
		switch (np) {
		  case 1:
	        edges[element][edge] = edges_np_1[element][edge];
	        break;
		  case 2:
		    edges[element][edge] = edges_np_2[element][edge];
		    break;
		  case 3:
			edges[element][edge] = edges_np_3[element][edge];
			break;
		  default:
			edges[element][edge] = -1;
		}
	  }
		for (int face = 0; face < FACES; ++face) {
		  switch (np) {
		    case 1:
		      faces[element][face] = faces_np_1[element][face];
		      break;
			case 2:
			  faces[element][face] = faces_np_2[element][face];
			  break;
			case 3:
		      faces[element][face] = faces_np_3[element][face];
			  break;
			default:
		      faces[element][face] = -1;
		  }
	    }
	}

    std::vector<int> connectivity;
    std::map<int, std::vector<int>> connectivityGlobalElementMap;
    int count = 1;

    for (int element = 0; element < ELEMENTS; ++element) {
      for (int corner_node = 0; corner_node < CORNER_NODES; ++corner_node) {
        connectivity.push_back(nodes[element][corner_node]);
      }
      for (int edge = 0; edge < EDGES; ++edge) {
        for (int id = 0; id < IDS_PER_EDGE; ++id) {
          connectivity.push_back(EDGE_OFFSET + edges[element][edge] * IDS_PER_EDGE + id);
        }
      }
      for (int face = 0; face < FACES; ++face) {
        for (int id = 0; id < IDS_PER_FACE; ++id) {
          connectivity.push_back(FACE_OFFSET + faces[element][face] * IDS_PER_FACE + id);
        }
      }
      for (int id = 0; id < IDS_PER_CELL; ++id) {
        connectivity.push_back(CELL_OFFSET + (1+element) * IDS_PER_CELL + id);
      }
      connectivityGlobalElementMap.insert(intVecIntPair(count, connectivity));
      connectivity.clear();
      ++count;
    }

    std::vector<int> elementLidToGid;


    if (np == 1) {
      int elLidToGid[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    		                11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  	  for (size_t i = 0; i < correct.ownedElementCount; ++i) {
  		elementLidToGid.push_back(elLidToGid[i]);
  	    correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
  	  }
    }

    else if (np == 2) {
      if (rank == 0) {
        int elLidToGid[10] = {1, 2, 3, 4, 5, 11, 12, 13, 14, 15};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
          elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else {
        int elLidToGid[10] = {6, 7, 8, 9, 10, 16, 17, 18, 19, 20};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }

    }

    else {
      if (rank == 0) {
        int elLidToGid[6] = {1, 2, 3, 4, 5, 11};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else if (rank == 1) {
        int elLidToGid[6] = {6, 7, 12, 13, 14, 15};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
          elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
      else {
        int elLidToGid[8] = {8, 9, 10, 16, 17, 18, 19, 20};
    	for (size_t i = 0; i < correct.ownedElementCount; ++i) {
    	  elementLidToGid.push_back(elLidToGid[i]);
    	  correct.connectivityLocalElementMap.insert(intVecIntPair(i, connectivityGlobalElementMap[elLidToGid[i]]));
    	}
      }
    }

    // Correct dof coordinates for all nodes
    int dim = 3;
    std::vector<double> dofCoordsThisLoc(dim);
    std::vector<double> dofCoordsThisCenter(dim);

    double eLen = 1.0;
    double nodeDeltas[ELEMENTS][3] = {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5},
    		                          {-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}};

    double edgeIdDeltas[EDGE_IDS][3] = {{0.0, -0.5, -0.5}, {0.5, 0.0, -0.5}, {0.0, 0.5, -0.5}, {-0.5, 0.0, -0.5},
    		                            {0.0, -0.5, 0.5}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}, {-0.5, 0.0, 0.5},
                                        {-0.5, -0.5, 0.0}, {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0}, {-0.5, 0.5, 0.0}};

    double faceIdDeltas[FACE_IDS][3] = {{0.0, -0.5, 0.0}, {0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {-0.5, 0.0, 0.0},
									    {0.0, 0.0, -0.5}, {0.0, 0.0, 0.5},};

    double cellIdDeltas[IDS_PER_CELL][3] = {{0., 0.}};

    double elementCorners[ELEMENTS][3] = {{0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}, {3., 0., 0.}, {4., 0., 0.},
    		                              {0., 1., 0.}, {1., 1., 0.}, {2., 1., 0.}, {3., 1., 0.}, {4., 1., 0.},
                                          {5., 0., 0.}, {6., 0., 0.}, {7., 0., 0.}, {8., 0., 0.}, {9., 0., 0.},
										  {5., 1., 0.}, {6., 1., 0.}, {7., 1., 0.}, {8., 1., 0.}, {9., 1., 0.}};

    for (int element = 0; element < ELEMENTS; ++element) {
      for (int k = 0; k < dim; ++k) {
    	nodeDeltas[element][k] *= eLen;
    	edgeIdDeltas[element][k] *= eLen;
    	faceIdDeltas[element][k] *= eLen;
    	cellIdDeltas[element][k] *= eLen;
        elementCorners[element][k] *= eLen;
      }
    }

    for (int element = 0; element < ELEMENTS; ++element) {
      dofCoordsThisCenter[0] = elementCorners[element][0] + eLen/2.0;
      dofCoordsThisCenter[1] = elementCorners[element][1] + eLen/2.0;
      dofCoordsThisCenter[2] = elementCorners[element][2] + eLen/2.0;
      for (int cornerNode = 0; cornerNode < CORNER_NODES; ++cornerNode) {
    	for (int k = 0; k < dim; ++k) {
          dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + nodeDeltas[cornerNode][k];
    	}
        correct.dofCoords[nodes[element][cornerNode]] = dofCoordsThisLoc;
      }
      for (int edge = 0; edge < EDGES; ++edge) {
        for (int i = 0; i < IDS_PER_EDGE; ++i) {
    	  for (int k = 0; k < dim; ++k) {
            dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + edgeIdDeltas[edge*IDS_PER_EDGE+i][k];
    	  }
          correct.dofCoords[EDGE_OFFSET + IDS_PER_EDGE * edges[element][edge] + i] = dofCoordsThisLoc;
        }
      }
      for (int face = 0; face < FACES; ++face) {
        for (int i = 0; i < IDS_PER_FACE; ++i) {
    	  for (int k = 0; k < dim; ++k) {
            dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + faceIdDeltas[face*IDS_PER_FACE+i][k];
    	  }
          correct.dofCoords[FACE_OFFSET + IDS_PER_FACE * faces[element][face] + i] = dofCoordsThisLoc;
        }
      }
      for (int i = 0; i < IDS_PER_CELL; ++i) {
        for (int k = 0; k < dim; ++k) {
          dofCoordsThisLoc[k] = dofCoordsThisCenter[k] + cellIdDeltas[i][k];
        }
        correct.dofCoords[CELL_OFFSET + IDS_PER_CELL * (element+1) + i] = dofCoordsThisLoc;
      }
    }

  }
} // end empty namespace

} // end panzer_ioss namespace
