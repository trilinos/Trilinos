// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
	int rank = comm.getRank(); // processor rank

	TEUCHOS_TEST_FOR_EXCEPTION(np > MAX_PROCESSES, std::logic_error,
					"Error, Test was run with " << np << " processes,"
					<< "but is only valid for up to " << MAX_PROCESSES << " processes." << std::endl);

	typedef Kokkos::Experimental::DynRankView<double, PHX::Device> FieldContainer;
	typedef std::pair<std::string, std::vector<int>> strVecIntPair;
	typedef std::pair<int, std::vector<int>> intVecIntPair;

	std::string filename = "rectangle.pg";
	std::string iossDatabaseType = "pamgen";

    // For printing
    //Teuchos::FancyOStream output(Teuchos::rcpFromRef(std::cout));
	//output.setShowProcRank(true);
    //output << std::endl;

    //Kokkos::initialize(argc, argv);

	// Create the correct versions of data structures representing the mesh
	// for comparison with the actual data structures extracted from the database.

	// Correct number of element blocks
	std::size_t correctNumElementBlocks = 2;

	// Correct element block topologies
	std::vector<shards::CellTopology> correctElementBlockTopologies;
	for (size_t b = 0; b < correctNumElementBlocks; ++b) {
	  correctElementBlockTopologies.push_back(shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4>>()));
	}

	// Correct element block ids
    std::vector<std::string> correctElementBlockIds;
    correctElementBlockIds.push_back("block_1");
    correctElementBlockIds.push_back("block_2");

    // Corred dof coordinates for all nodes
    const int DIM = 2;
    const int ROWS = 3;
    const int COLUMNS = 11;
    double correctDOFCoords[ROWS*COLUMNS][DIM];
    int nd = 0;
    for (int i = 0; i < ROWS; ++i) {
      for (int j = 0; j < COLUMNS; ++j, ++nd) {
        correctDOFCoords[nd][0] = 1.0 * j;
        correctDOFCoords[nd][1] = 1.0 * i;
      }
    }

    // Correct owned element count, local element ids,
    // and block ids, connectivity sizes, and connectivities for local elements.
    std::size_t correctOwnedElementCount;
	std::vector<int> elementBlock;
	std::map<std::string,std::vector<int>> correctElementBlockMap;
    std::vector<std::string> correctLocalElementBlockIds;
    std::vector<int> correctConnectivitySize;
    const int QUAD_CONNECTIVITY_SIZE = 4;
    std::map<int, std::vector<int>> correctConnectivityGlobalElementMap;
    std::map<int, std::vector<int>> correctConnectivityLocalElementMap;
    std::vector<int> correctConnectivity;
    std::vector<int> correctElementLidToGid;

    correctConnectivity.push_back(1); correctConnectivity.push_back(2); correctConnectivity.push_back(13); correctConnectivity.push_back(12);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(1, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(2); correctConnectivity.push_back(3); correctConnectivity.push_back(14); correctConnectivity.push_back(13);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(2, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(3); correctConnectivity.push_back(4); correctConnectivity.push_back(15); correctConnectivity.push_back(14);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(3, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(4); correctConnectivity.push_back(5); correctConnectivity.push_back(16); correctConnectivity.push_back(15);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(4, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(5); correctConnectivity.push_back(6); correctConnectivity.push_back(17); correctConnectivity.push_back(16);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(5, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(12); correctConnectivity.push_back(13); correctConnectivity.push_back(24); correctConnectivity.push_back(23);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(6, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(13); correctConnectivity.push_back(14); correctConnectivity.push_back(25); correctConnectivity.push_back(24);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(7, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(14); correctConnectivity.push_back(15); correctConnectivity.push_back(26); correctConnectivity.push_back(25);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(8, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(15); correctConnectivity.push_back(16); correctConnectivity.push_back(27); correctConnectivity.push_back(26);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(9, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(16); correctConnectivity.push_back(17); correctConnectivity.push_back(28); correctConnectivity.push_back(27);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(10, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(6); correctConnectivity.push_back(7); correctConnectivity.push_back(18); correctConnectivity.push_back(17);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(11, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(7); correctConnectivity.push_back(8); correctConnectivity.push_back(19); correctConnectivity.push_back(18);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(12, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(8); correctConnectivity.push_back(9); correctConnectivity.push_back(20); correctConnectivity.push_back(19);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(13, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(9); correctConnectivity.push_back(10); correctConnectivity.push_back(21); correctConnectivity.push_back(20);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(14, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(10); correctConnectivity.push_back(11); correctConnectivity.push_back(22); correctConnectivity.push_back(21);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(15, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(17); correctConnectivity.push_back(18); correctConnectivity.push_back(29); correctConnectivity.push_back(28);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(16, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(18); correctConnectivity.push_back(19); correctConnectivity.push_back(30); correctConnectivity.push_back(29);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(17, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(19); correctConnectivity.push_back(20); correctConnectivity.push_back(31); correctConnectivity.push_back(30);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(18, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(20); correctConnectivity.push_back(21); correctConnectivity.push_back(32); correctConnectivity.push_back(31);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(19, correctConnectivity)); correctConnectivity.clear();

    correctConnectivity.push_back(21); correctConnectivity.push_back(22); correctConnectivity.push_back(33); correctConnectivity.push_back(32);
    correctConnectivityGlobalElementMap.insert(intVecIntPair(20, correctConnectivity)); correctConnectivity.clear();


	if (np == 1) {
	  correctOwnedElementCount = 20;
	  correctConnectivitySize.insert(correctConnectivitySize.end(), correctOwnedElementCount, QUAD_CONNECTIVITY_SIZE);


	  // block_1
	  elementBlock.clear();
	  for (int i = 0; i <= 9; ++i) {
	    elementBlock.push_back(i);
	    correctLocalElementBlockIds.push_back(correctElementBlockIds[0]);
	  }
	  correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[0], elementBlock));
	  //block_2
	  elementBlock.clear();
	  for (int i = 10; i <= 19; ++i) {
	    elementBlock.push_back(i);
	    correctLocalElementBlockIds.push_back(correctElementBlockIds[1]);
	  }
	  correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[1], elementBlock));

	  correctElementLidToGid.push_back(1);
	  correctElementLidToGid.push_back(2);
	  correctElementLidToGid.push_back(3);
	  correctElementLidToGid.push_back(4);
	  correctElementLidToGid.push_back(5);
	  correctElementLidToGid.push_back(6);
	  correctElementLidToGid.push_back(7);
	  correctElementLidToGid.push_back(8);
	  correctElementLidToGid.push_back(9);
	  correctElementLidToGid.push_back(10);
	  correctElementLidToGid.push_back(11);
	  correctElementLidToGid.push_back(12);
	  correctElementLidToGid.push_back(13);
	  correctElementLidToGid.push_back(14);
	  correctElementLidToGid.push_back(15);
	  correctElementLidToGid.push_back(16);
	  correctElementLidToGid.push_back(17);
	  correctElementLidToGid.push_back(18);
	  correctElementLidToGid.push_back(19);
	  correctElementLidToGid.push_back(20);

	  for (size_t i = 0; i < correctOwnedElementCount; ++i) {
	    correctConnectivityLocalElementMap.insert(intVecIntPair(i, correctConnectivityGlobalElementMap[correctElementLidToGid[i]]));
	  }

	}
	else if (np == 2) {

	  if (rank == 0) {
		correctOwnedElementCount = 10;
		correctConnectivitySize.insert(correctConnectivitySize.end(), correctOwnedElementCount, QUAD_CONNECTIVITY_SIZE);

	    // block_1
	    elementBlock.clear();
	    for (int i = 0; i <= 4; ++i) {
	      elementBlock.push_back(i);
	      correctLocalElementBlockIds.push_back(correctElementBlockIds[0]);
	    }
	    correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[0], elementBlock));
	    //block_2
	    elementBlock.clear();
	    for (int i = 5; i <= 9; ++i) {
	      elementBlock.push_back(i);
	      correctLocalElementBlockIds.push_back(correctElementBlockIds[1]);
	    }
	    correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[1], elementBlock));

	    correctElementLidToGid.push_back(1);
	    correctElementLidToGid.push_back(2);
	    correctElementLidToGid.push_back(3);
	    correctElementLidToGid.push_back(4);
	    correctElementLidToGid.push_back(5);
	    correctElementLidToGid.push_back(11);
	    correctElementLidToGid.push_back(12);
	    correctElementLidToGid.push_back(13);
	    correctElementLidToGid.push_back(14);
	    correctElementLidToGid.push_back(15);

	    for (size_t i = 0; i < correctOwnedElementCount; ++i) {
	      correctConnectivityLocalElementMap.insert(intVecIntPair(i, correctConnectivityGlobalElementMap[correctElementLidToGid[i]]));
	    }

	  }
	  else {
		correctOwnedElementCount = 10;
		correctConnectivitySize.insert(correctConnectivitySize.end(), correctOwnedElementCount, QUAD_CONNECTIVITY_SIZE);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <= 4; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[0]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[0], elementBlock));
		//block_2
		elementBlock.clear();
		for (int i = 5; i <= 9; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[1]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[1], elementBlock));

		correctElementLidToGid.push_back(6);
		correctElementLidToGid.push_back(7);
		correctElementLidToGid.push_back(8);
		correctElementLidToGid.push_back(9);
		correctElementLidToGid.push_back(10);
		correctElementLidToGid.push_back(16);
		correctElementLidToGid.push_back(17);
		correctElementLidToGid.push_back(18);
		correctElementLidToGid.push_back(19);
		correctElementLidToGid.push_back(20);

		for (size_t i = 0; i < correctOwnedElementCount; ++i) {
		  correctConnectivityLocalElementMap.insert(intVecIntPair(i, correctConnectivityGlobalElementMap[correctElementLidToGid[i]]));
	    }

	  }
	}
	else {
	  if (rank == 0) {
		correctOwnedElementCount = 6;
		correctConnectivitySize.insert(correctConnectivitySize.end(), correctOwnedElementCount, QUAD_CONNECTIVITY_SIZE);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <=4; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[0]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[0], elementBlock));
		// block_2
		elementBlock.clear();
		for (int i = 5; i <=5; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[1]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[1], elementBlock));

		correctElementLidToGid.push_back(1);
		correctElementLidToGid.push_back(2);
		correctElementLidToGid.push_back(3);
		correctElementLidToGid.push_back(4);
		correctElementLidToGid.push_back(5);
		correctElementLidToGid.push_back(11);

		for (size_t i = 0; i < correctOwnedElementCount; ++i) {
	      correctConnectivityLocalElementMap.insert(intVecIntPair(i, correctConnectivityGlobalElementMap[correctElementLidToGid[i]]));
	    }

	  }
	  else if (rank == 1) {
		correctOwnedElementCount = 6;
		correctConnectivitySize.insert(correctConnectivitySize.end(), correctOwnedElementCount, QUAD_CONNECTIVITY_SIZE);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <= 1; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[0]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[0], elementBlock));
		//block_2
		elementBlock.clear();
		for (int i = 2; i <= 5; ++i) {
	      elementBlock.push_back(i);
	      correctLocalElementBlockIds.push_back(correctElementBlockIds[1]);
		}
	    correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[1], elementBlock));

	    correctElementLidToGid.push_back(6);
	    correctElementLidToGid.push_back(7);
	    correctElementLidToGid.push_back(12);
	    correctElementLidToGid.push_back(13);
	    correctElementLidToGid.push_back(14);
	    correctElementLidToGid.push_back(15);

	    for (size_t i = 0; i < correctOwnedElementCount; ++i) {
	      correctConnectivityLocalElementMap.insert(intVecIntPair(i, correctConnectivityGlobalElementMap[correctElementLidToGid[i]]));
	    }

	  }
	  else {
		correctOwnedElementCount = 8;
		correctConnectivitySize.insert(correctConnectivitySize.end(), correctOwnedElementCount, QUAD_CONNECTIVITY_SIZE);

	    // block_1
		elementBlock.clear();
		for (int i = 0; i <= 2; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[0]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[0], elementBlock));
		//block_2
		elementBlock.clear();
		for (int i = 3; i <= 7; ++i) {
		  elementBlock.push_back(i);
		  correctLocalElementBlockIds.push_back(correctElementBlockIds[1]);
		}
		correctElementBlockMap.insert(strVecIntPair(correctElementBlockIds[1], elementBlock));

		correctElementLidToGid.push_back(8);
		correctElementLidToGid.push_back(9);
		correctElementLidToGid.push_back(10);
		correctElementLidToGid.push_back(16);
		correctElementLidToGid.push_back(17);
		correctElementLidToGid.push_back(18);
		correctElementLidToGid.push_back(19);
		correctElementLidToGid.push_back(20);

		for (size_t i = 0; i < correctOwnedElementCount; ++i) {
	      correctConnectivityLocalElementMap.insert(intVecIntPair(i, correctConnectivityGlobalElementMap[correctElementLidToGid[i]]));
	    }

	  }
	}

	// Initialize Ioss.
    Ioss::Init::Initializer();

    // Check for a valid database type.
	TEUCHOS_TEST_FOR_EXCEPTION(!IossDatabaseTypeManager::isValidType(iossDatabaseType), std::logic_error,
				"Error, " << iossDatabaseType
				<< " is not a valid IOSS database type." << std::endl
				<< "Valid types are " << IossDatabaseTypeManager::validTypeList());

	// Create an empty property manager.
    Ioss::PropertyManager databaseProps;

    // Create the Ioss database.
    Ioss::DatabaseIO * iossMeshDB = Ioss::IOFactory::create(iossDatabaseType, filename,
    		Ioss::READ_MODEL, MPI_COMM_WORLD, databaseProps);

    // Create the connectivity manager.
    IOSSConnManager<int> connManager(iossMeshDB);

    // Create a clone
    bool cloneCreated = false;
    Teuchos::RCP<panzer::ConnManagerBase<int>> cmClone;
    if (IossDatabaseTypeManager::supportsMultipleOpenDatabases(iossDatabaseType)) {
      cmClone = connManager.noConnectivityClone(iossDatabaseType, databaseProps);
      cloneCreated = true;
    }

    // Check for correct number of element blocks.
    std::size_t numElementBlocks = connManager.numElementBlocks();
    TEST_EQUALITY(numElementBlocks, correctNumElementBlocks);
    if (cloneCreated) {
      std::size_t cloneNumElementBlocks = cmClone->numElementBlocks();
      TEST_EQUALITY(cloneNumElementBlocks, correctNumElementBlocks);
    }

    // Check for correct element block topologies.
    std::vector<shards::CellTopology> elementBlockTopologies;
    connManager.getElementBlockTopologies(elementBlockTopologies);
    TEST_COMPARE_ARRAYS(elementBlockTopologies, correctElementBlockTopologies);
    if (cloneCreated) {
      std::vector<shards::CellTopology> cloneElementBlockTopologies;
      cmClone->getElementBlockTopologies(cloneElementBlockTopologies);
      TEST_COMPARE_ARRAYS(cloneElementBlockTopologies, correctElementBlockTopologies);
    }

    // Make nodal field pattern
    Teuchos::RCP<panzer::FieldPattern> iossNodalFP = rcp(new panzer::NodalFieldPattern(elementBlockTopologies[0]));

    // Build the connectivity
    connManager.buildConnectivity(*iossNodalFP);
    if (cloneCreated) {
      cmClone->buildConnectivity(*iossNodalFP);
    } // end if supportsMultipleOpenDatabases

    // Check for correct element block ids.
    std::vector<std::string> elementBlockIds;
    connManager.getElementBlockIds(elementBlockIds);
    TEST_COMPARE_ARRAYS(elementBlockIds, correctElementBlockIds)
    if (cloneCreated) {
      std::vector<std::string> cloneElementBlockIds;
      cmClone->getElementBlockIds(cloneElementBlockIds);
      TEST_COMPARE_ARRAYS(cloneElementBlockIds, correctElementBlockIds)
    }

    // Check for correct local element ids for each block
    for (std::string elementBlockId : elementBlockIds) {
      elementBlock = connManager.getElementBlock(elementBlockId);
      TEST_COMPARE_ARRAYS(elementBlock, correctElementBlockMap[elementBlockId])
      if (cloneCreated) {
        elementBlock = cmClone->getElementBlock(elementBlockId);
    	TEST_COMPARE_ARRAYS(elementBlock, correctElementBlockMap[elementBlockId])
      }
    }

    // Check for correct owned element count
    // Method does not exist for clone, which is of parent type.
    std::size_t ownedElementCount = connManager.getOwnedElementCount();
    TEST_EQUALITY(ownedElementCount, correctOwnedElementCount);


    // Check for correct block id for all local elements
    std::vector<std::string> localElementBlockIds;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      localElementBlockIds.push_back(connManager.getBlockId(localElementId));
    }
    TEST_COMPARE_ARRAYS(localElementBlockIds, correctLocalElementBlockIds);


    // Check for correct connectivity size for all local elements
    std::vector<int> connectivitySize;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      connectivitySize.push_back(connManager.getConnectivitySize(localElementId));
    }
    TEST_COMPARE_ARRAYS(connectivitySize, correctConnectivitySize);

    // Check for correct connectivity for all local elements
    int * conn;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      conn = connManager.getConnectivity(localElementId);
      for (int i = 0; i < connectivitySize[localElementId]; ++i) {
        TEST_EQUALITY(conn[i], correctConnectivityLocalElementMap[localElementId][i]);
      }
    }

    RCP<Intrepid2::Basis<PHX::Device,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::Device,double,double>);
    RCP<const panzer::Intrepid2FieldPattern> intrepid2fp = rcp(new panzer::Intrepid2FieldPattern(basis));
    std::vector<int> localCellIds;
    FieldContainer points;
    std::string elementBlockId;
    size_t size;
    double tolerance = 1.0e-12;
    double coordinate;
    double correctCoordinate;
    int localCellId;
    for (size_t block = 0; block < numElementBlocks; ++block) {
      elementBlockId = elementBlockIds[block];
      connManager.getDofCoords(elementBlockId, *intrepid2fp, localCellIds, points);
      for (size_t el = 0; el < localCellIds.size(); ++el) {
    	localCellId = localCellIds[el];
        size = connManager.getConnectivitySize(localCellId);
        conn = connManager.getConnectivity(localCellId);
        for (size_t node = 0; node < size; ++node) {
          for (size_t k = 0; k < elementBlockTopologies[block].getDimension(); ++k) {
            coordinate = points(el, node, k);
            correctCoordinate = correctDOFCoords[conn[node]-1][k];
            TEST_FLOATING_EQUALITY(coordinate, correctCoordinate, tolerance);
          }
        }
      }
    }

    //Kokkos::finalize();
}

}
