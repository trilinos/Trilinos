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


#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

//#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_RCP.hpp"

#include "PanzerCore_config.hpp"
#include "PanzerAdaptersIOSS_config.hpp"
#include "Panzer_IOSSConnManager.hpp"

#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

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

	//int np = comm.getSize(); // number of processors
	//int rank = comm.getRank(); // processor rank

	typedef Kokkos::Experimental::DynRankView<double, PHX::Device> FieldContainer;

	std::string filename = "rectangle.pg";

    //Kokkos::initialize(argc, argv);

    Ioss::Init::Initializer();

    Ioss::PropertyManager databaseProps;

    Ioss::DatabaseIO * iossMeshDB = Ioss::IOFactory::create("pamgen", filename,
    		Ioss::READ_MODEL, MPI_COMM_WORLD, databaseProps);

    // Create the connectivity manager.
    IOSSConnManager<int> connManager(iossMeshDB);

    // Make nodal field pattern
    std::vector<shards::CellTopology> elementBlockTopologies;
    connManager.getElementBlockTopologies(elementBlockTopologies);
    Teuchos::RCP<panzer::FieldPattern> iossNodalFP = rcp(new panzer::NodalFieldPattern(elementBlockTopologies[0]));

    // Build the connectivity
    connManager.buildConnectivity(*iossNodalFP);

    // Build clone
    // Calling noConnectivityClone() causes a runtime error. Still need to fix this.
    std::string databaseType = "pamgen";
    Teuchos::RCP<panzer::ConnManagerBase<int>> clone = connManager.noConnectivityClone(databaseType, databaseProps);

    // Build connectivity of the clone
    clone->buildConnectivity(*iossNodalFP);

    // Print information about the connectivity manager.
    // Need to change from print statements to checks for exception.
    Teuchos::FancyOStream output(Teuchos::rcpFromRef(std::cout));
	output.setShowProcRank(true);

    std::vector<double> vec_double;
    std::vector<int> vec_int;
    int count = 0;

    std::size_t ownedElementCount = connManager.getOwnedElementCount();

    output << std::endl;

    std::vector<std::string> elementBlockIds;
    std::vector<int> elementBlock;
    std::vector<shards::CellTopology> blockTopologies;

    connManager.getElementBlockIds(elementBlockIds);
    connManager.getElementBlockTopologies(blockTopologies);

    output << "There are " << connManager.numElementBlocks() << " element blocks." << std::endl;

    count = 0;
    for (std::string elementBlockId : elementBlockIds) {
      output << "element block " << count
          << ", name = " << elementBlockId
		  << ", topology = " << blockTopologies[count].getName()
		  << std::endl;
      ++count;
      elementBlock = connManager.getElementBlock(elementBlockId);
      for (auto element : elementBlock) {
        output << element << ", ";
      }
      output << std::endl;
    }

    output << "Individual elements (" << ownedElementCount << "):" << std::endl;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      output << "local element = " << localElementId << ", "
          << "block = " << connManager.getBlockId(localElementId) << std::endl;
      }

    size_t size;
    int * conn;
    for (size_t localElementId = 0; localElementId < ownedElementCount; ++localElementId) {
      size = connManager.getConnectivitySize(localElementId);
      output << "local element " << localElementId << " conn ( "
          << size << ") = ";
      conn = connManager.getConnectivity(localElementId);
      for (size_t subcell = 0; subcell < size; ++subcell) {
        output << conn[subcell] << ", ";
      }
      output << std::endl;
    }

    TEST_EQUALITY(elementBlockTopologies[0],shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4>>()));
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::Device,double,double>);
    RCP<const panzer::Intrepid2FieldPattern> intrepid2fp = rcp(new panzer::Intrepid2FieldPattern(basis));
    std::vector<int> localCellIds;
    FieldContainer points;
    std::string elementBlockId;
    int localCellId;
    for (size_t block = 0; block < connManager.numElementBlocks(); ++block) {
      elementBlockId = elementBlockIds[block];
      output << "Coordinates for element block " << elementBlockId << ":" << std::endl;
      connManager.getDofCoords(elementBlockId, *intrepid2fp, localCellIds, points);
      for (size_t el = 0; el < localCellIds.size(); ++el) {
    	localCellId = localCellIds[el];
        output << "local element " << localCellId << ": " << std::endl;
        size = connManager.getConnectivitySize(localCellId);
        conn = connManager.getConnectivity(localCellId);
        for (size_t node = 0; node < size; ++node) {
          output << "   node " << conn[node] << ": (";
          for (size_t k = 0; k < blockTopologies[block].getDimension(); ++k) {
            output << points(el, node, k) << ", ";
          }
          output << ")" << std::endl;
        }
      }
    }

    //Kokkos::finalize();
}

}
