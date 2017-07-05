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

namespace {

    struct Correct {

        static const int DIM = 2;
        static const int ROWS = 3;
        static const int COLUMNS = 11;

    	std::size_t numElementBlocks;
    	std::vector<shards::CellTopology> elementBlockTopologies;
    	std::vector<std::string> elementBlockIds;

        double dofCoords[ROWS*COLUMNS][DIM];

        std::map<std::string,std::vector<int>> elementBlockMap;
        std::size_t ownedElementCount;
        std::vector<std::string> localElementBlockIds;
        std::vector<int> connectivitySize;
        std::map<int, std::vector<int>> connectivityLocalElementMap;

    };

    Correct correct;
    void setCorrect(Teuchos::MpiComm<int> & comm);

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

	typedef Kokkos::Experimental::DynRankView<double, PHX::Device> FieldContainer;

	std::string filename = "rectangle.pg";
	std::string iossDatabaseType = "pamgen";

    // For printing
    //Teuchos::FancyOStream output(Teuchos::rcpFromRef(std::cout));
	//output.setShowProcRank(true);
    //output << std::endl;

    //Kokkos::initialize(argc, argv);

	// Create the correct versions of data structures representing the mesh
	// for comparison with the actual data structures extracted from the database.
	setCorrect(comm);

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
    TEST_EQUALITY(numElementBlocks, correct.numElementBlocks);
    if (cloneCreated) {
      std::size_t cloneNumElementBlocks = cmClone->numElementBlocks();
      TEST_EQUALITY(cloneNumElementBlocks, correct.numElementBlocks);
    }

    // Check for correct element block topologies.
    std::vector<shards::CellTopology> elementBlockTopologies;
    connManager.getElementBlockTopologies(elementBlockTopologies);
    TEST_COMPARE_ARRAYS(elementBlockTopologies, correct.elementBlockTopologies);
    if (cloneCreated) {
      std::vector<shards::CellTopology> cloneElementBlockTopologies;
      cmClone->getElementBlockTopologies(cloneElementBlockTopologies);
      TEST_COMPARE_ARRAYS(cloneElementBlockTopologies, correct.elementBlockTopologies);
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
    TEST_COMPARE_ARRAYS(elementBlockIds, correct.elementBlockIds)
    if (cloneCreated) {
      std::vector<std::string> cloneElementBlockIds;
      cmClone->getElementBlockIds(cloneElementBlockIds);
      TEST_COMPARE_ARRAYS(cloneElementBlockIds, correct.elementBlockIds)
    }

    // Check for correct local element ids for each block
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
    std::size_t ownedElementCount = connManager.getOwnedElementCount();
    TEST_EQUALITY(ownedElementCount, correct.ownedElementCount);


    // Check for correct block id for all local elements
    std::vector<std::string> localElementBlockIds;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      localElementBlockIds.push_back(connManager.getBlockId(localElementId));
    }
    TEST_COMPARE_ARRAYS(localElementBlockIds, correct.localElementBlockIds);


    // Check for correct connectivity size for all local elements
    std::vector<int> connectivitySize;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      connectivitySize.push_back(connManager.getConnectivitySize(localElementId));
    }
    TEST_COMPARE_ARRAYS(connectivitySize, correct.connectivitySize);

    // Check for correct connectivity for all local elements
    int * conn;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      conn = connManager.getConnectivity(localElementId);
      for (int i = 0; i < connectivitySize[localElementId]; ++i) {
        TEST_EQUALITY(conn[i], correct.connectivityLocalElementMap[localElementId][i]);
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
            correctCoordinate = correct.dofCoords[conn[node]-1][k];
            TEST_FLOATING_EQUALITY(coordinate, correctCoordinate, tolerance);
          }
        }
      }
    }

    //Kokkos::finalize();
}

namespace {

  void setCorrect(Teuchos::MpiComm<int> & comm) {

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
    int nd = 0;
    for (int i = 0; i < Correct::ROWS; ++i) {
      for (int j = 0; j < Correct::COLUMNS; ++j, ++nd) {
        correct.dofCoords[nd][0] = 1.0 * j;
        correct.dofCoords[nd][1] = 1.0 * i;
      }
    }

    // Correct owned element count, local element ids,
    // and block ids, connectivity sizes, and connectivities for local elements.

    const int QUAD_CONNECTIVITY_SIZE = 4;
    const int ELEMENTS = 20;
    int conn[ELEMENTS][QUAD_CONNECTIVITY_SIZE] =
                      {{1, 2, 13, 12}, {2, 3, 14, 13}, {3, 4, 15, 14}, {4, 5, 16, 15},
                       {5, 6, 17, 16}, {12, 13, 24, 23}, {13, 14, 25, 24}, {14, 15, 26, 25},
                       {15, 16, 27, 26}, {16, 17, 28, 27}, {6, 7, 18, 17}, {7, 8, 19, 18},
					   {8, 9, 20, 19}, {9, 10, 21, 20}, {10, 11, 22, 21}, {17, 18, 29, 28},
					   {18, 19, 30, 29}, {19, 20, 31, 30}, {20, 21, 32, 31}, {21, 22, 33, 32}};

    std::vector<int> connectivity;
    std::map<int, std::vector<int>> connectivityGlobalElementMap;
    int count = 1;
    for (int element = 0; element < ELEMENTS; ++element) {
      for (int index = 0; index < QUAD_CONNECTIVITY_SIZE; ++index) {
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
	  correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, QUAD_CONNECTIVITY_SIZE);


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
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, QUAD_CONNECTIVITY_SIZE);

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
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, QUAD_CONNECTIVITY_SIZE);

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
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, QUAD_CONNECTIVITY_SIZE);

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
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, QUAD_CONNECTIVITY_SIZE);

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
		correct.connectivitySize.insert(correct.connectivitySize.end(), correct.ownedElementCount, QUAD_CONNECTIVITY_SIZE);

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

  } // end setCorrect()

} // end empty namespace

} // end panzer_ioss namespace
