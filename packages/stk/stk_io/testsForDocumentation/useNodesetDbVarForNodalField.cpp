#include <gtest/gtest.h>                // for AssertHelper, ASSERT_EQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <string>                       // for allocator, operator+, etc
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeSet.h"               // for NodeSet
#include "Ioss_Region.h"                // for NodeSetContainer, Region
#include "Ioss_Utils.h"                 // for Utils
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
namespace Ioss { class DatabaseIO; }

namespace {

  TEST(StkMeshIoBrokerHowTo, useNodesetDbVarForNodalFields)
  {
    std::string resultsFilename = "nodeset_fields.results";
    std::string dbFieldName = "temp";
    std::string appFieldName = "temperature";
    
    MPI_Comm communicator = MPI_COMM_WORLD;

    size_t num_elems_per_edge = 9;  
    {
      //-BEGIN
      // ============================================================
      // INITIALIZATION
      std::string s_elems_per_edge = Ioss::Utils::to_string(num_elems_per_edge);

      //+ Create a generated mesh containg hexes and shells.
      std::string input_filename = s_elems_per_edge + "x" +
	s_elems_per_edge + "x" +
	s_elems_per_edge + "|shell:xyzXYZ";

      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(input_filename, "generated",
			      stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::MetaData &meta_data = stkIo.meta_data();
      stk::mesh::Field<double> &temperature = meta_data.
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,
						 appFieldName, 1);

      // ============================================================
      //+ Put the temperature field on the nodes of the shell parts.
      const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();
      stk::mesh::Selector shell_subset;
      for (size_t i=0; i < all_parts.size(); i++) {
	const stk::mesh::Part *part = all_parts[i];
	stk::topology topo = part->topology();
	if (topo == stk::topology::SHELL_QUAD_4) {
	  stk::mesh::put_field(temperature, *part);
	}
      }

      stkIo.populate_bulk_data();

      // Create the output...
      size_t fh = stkIo.create_output_mesh(resultsFilename, stk::io::WRITE_RESULTS);

      //+ The "temperature" field will be output on nodesets consisting
      //+ of the nodes of each part the field is defined on.
      stkIo.use_nodeset_for_part_nodes_fields(fh, true);
      stkIo.add_field(fh, temperature, dbFieldName);

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(),
			      stk::topology::NODE_RANK, nodes);
    
      // Add three steps to the database
      // For each step, the value of the field is the value 'time'
      for (size_t i=0; i < 3; i++) {
	double time = i;

	for(size_t inode=0; inode<nodes.size(); inode++) {
	  double *fieldDataForNode = stk::mesh::field_data(temperature, nodes[inode]);
	  if (fieldDataForNode)
	    *fieldDataForNode = time;
	}

	stkIo.begin_output_step(fh, time);
	stkIo.write_defined_output_fields(fh);
	stkIo.end_output_step(fh);
      }
      // Verification omitted...
      //-END
    }

    {
      // ==================================================
      // VERIFICATION
      // Verify that the output mesh has 6 nodesets (input has 0)
      // corresponding to the nodes of the 6 shell blocks.
      // Each nodeset should have a single variable named dbFieldName

      Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", resultsFilename,
							 Ioss::READ_MODEL, communicator);
      Ioss::Region ioRegion(iossDb);

      // Six nodesets should have been created -- 1 for each shell block in the mesh.
      const Ioss::NodeSetContainer &nsets = ioRegion.get_nodesets();
      size_t nset_count = 6;
      ASSERT_EQ(nsets.size(), nset_count);

      for (size_t i=0; i < nsets.size(); i++) {
	const Ioss::NodeSet *nset = nsets[i];
	// Each nodeset should have a field named 'temp'
	ASSERT_TRUE(nset->field_exists(dbFieldName));
      }
    }
    unlink(resultsFilename.c_str());
}

}
