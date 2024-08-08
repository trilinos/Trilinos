/*
 * Akri_OutputUtils.cpp
 *
 *  Created on: Nov 21, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_
#include <Akri_DiagWriter.hpp>
#include <Akri_OutputUtils.hpp>

#include <string>

#include <stk_io/DatabasePurpose.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Ioss_SubSystem.h>
#include <Ioss_PropertyManager.h>

namespace krino {

static void enable_io_parts(const stk::mesh::PartVector & parts)
{
  for (auto * part : parts)
    stk::io::put_io_part_attribute(*part);
}

void output_mesh_with_fields_and_properties(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, Ioss::PropertyManager &properties, stk::io::DatabasePurpose purpose)
{
  const stk::mesh::PartVector emptyIoParts = turn_off_output_for_empty_io_parts(mesh, outputSelector);

  stk::io::StkMeshIoBroker stkIo;
  stk::mesh::BulkData & workAroundNonConstMesh = const_cast<stk::mesh::BulkData &>(mesh);
  stkIo.set_bulk_data(workAroundNonConstMesh);

  const size_t outputFileIndex = stkIo.create_output_mesh(fileName, purpose, properties);

  const int filterDisconnectedNodes = true;
  if (filterDisconnectedNodes)
  {
    // Will filter out nodes that are themselves selected, but without any attached elements that are selected.
    // For example, if selector is BLOCK_1 | NODESET_1, a node will not be output if it is in NODESET_1 and not BLOCK_1.
    std::shared_ptr<Ioss::Region> ioRegion = stkIo.get_output_ioss_region(outputFileIndex);
    ioRegion->property_add(Ioss::Property(stk::io::s_ignoreDisconnectedNodes, filterDisconnectedNodes));
  }

  if (step > 0)
  {
    const stk::mesh::FieldVector fields = stkIo.bulk_data().mesh_meta_data().get_fields();
    for (stk::mesh::FieldBase * field : fields)
    {
      const Ioss::Field::RoleType * fieldRole = stk::io::get_field_role(*field);
      if (fieldRole == nullptr || *fieldRole == Ioss::Field::TRANSIENT)
        stkIo.add_field(outputFileIndex, *field);
    }
  }

  stkIo.set_active_selector(outputSelector);
  stkIo.set_subset_selector(outputFileIndex, outputSelector);

  stkIo.write_output_mesh(outputFileIndex);

  stkIo.begin_output_step(outputFileIndex, time);
  stkIo.write_defined_output_fields(outputFileIndex);
  stkIo.end_output_step(outputFileIndex);

  enable_io_parts(emptyIoParts);
}

void output_composed_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose)
{
  Ioss::PropertyManager properties;
  properties.add(Ioss::Property("COMPOSE_RESULTS", 1));
  output_mesh_with_fields_and_properties(mesh, outputSelector, fileName, step, time, properties, purpose);
}

void output_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose)
{
  Ioss::PropertyManager properties;
  output_mesh_with_fields_and_properties(mesh, outputSelector, fileName, step, time, properties, purpose);
}

std::string create_filename_from_base_filename(const std::string & baseFileName, const int numFileRevisions)
{
  size_t lastIndex = baseFileName.find_last_of(".");
  STK_ThrowRequire(lastIndex != std::string::npos);
  const std::string fileBaseName = baseFileName.substr(0, lastIndex);
  return create_file_name(fileBaseName, numFileRevisions);
}

std::string create_file_name(const std::string & fileBaseName, const int fileIndex)
{
  STK_ThrowAssert(fileIndex >= 0);
  if (fileIndex < 1)
    return fileBaseName + ".e";

  STK_ThrowAssert(fileIndex <= 9998);
  char counterChar[32] = {'\0'};
  std::snprintf(counterChar, 32, "%.4d", fileIndex+1);
  const std::string filename = fileBaseName + ".e-s" + std::string(counterChar);
  return filename;
}

template<class FACET>
void write_facets( const std::vector<FACET> & facets, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm)
{
  int numFacets = facets.size();
  const int nodesPerFacet = FACET::DIM;
  const int numNodes = numFacets * nodesPerFacet;
  const unsigned numProcs = stk::parallel_machine_size(comm);
  const int procId = stk::parallel_machine_rank(comm);

  // Type (ExodusII) is hard-wired at this time.
  Ioss::PropertyManager properties;
  properties.add(Ioss::Property("COMPOSE_RESULTS", 1));
  const std::string fileName = create_file_name(fileBaseName, fileIndex);
  properties.add(Ioss::Property("base_filename", create_file_name(fileBaseName, 0)));
  if (fileIndex > 0)
    properties.add(Ioss::Property("state_offset", fileIndex));
  Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", create_file_name(fileBaseName, fileIndex), Ioss::WRITE_RESULTS, comm, properties);
  STK_ThrowRequireMsg(db != nullptr && db->ok(), "ERROR: Could not open output database '" << fileName << "' of type 'exodus'\n");
  Ioss::Region io(db, "FacetRegion");

  // find offsets for consistent global numbering
  int elem_start = 0;
  if ( numProcs > 1 ) {
    std::vector< int > numFaces(numProcs);
    // Gather everyone's facet list sizes.  I don't think there is a stk call for this...
    MPI_Allgather( &numFacets, 1, MPI_INT,
                   &(numFaces[0]), 1, MPI_INT, comm );
    for (int i = 0; i<procId; ++i) {
      elem_start += numFaces[i];
    }
  }

  const std::string description = "level set interface facets";
  io.property_add(Ioss::Property("title", description));
  io.begin_mode(Ioss::STATE_DEFINE_MODEL);

  Ioss::NodeBlock *nb = new Ioss::NodeBlock(db, "nodeblock_1", numNodes, nodesPerFacet);
  io.add(nb);
  nb->property_add(Ioss::Property("locally_owned_count", numNodes));

  std::string el_type;
  if constexpr (FACET::DIM == 3)
    el_type = "trishell3";
  else
    el_type = "shellline2d2";

  Ioss::ElementBlock *eb = new Ioss::ElementBlock(db, "block_1", el_type, numFacets);
  io.add(eb);

  io.end_mode(Ioss::STATE_DEFINE_MODEL);
  io.begin_mode(Ioss::STATE_MODEL);

  // loop over elements and node to create maps
  {
    std::vector< int > nmap, emap;
    nmap.reserve(numNodes);
    emap.reserve(numFacets);

    for ( unsigned n=0, e=0; e<facets.size(); ++e ) {
      emap.push_back(elem_start + e + 1);
      for ( int j = 0; j < nodesPerFacet; ++j ) {
        nmap.push_back(nodesPerFacet*elem_start + n + 1);
        ++n;
      }
    }
    nb->put_field_data("ids", nmap);
    eb->put_field_data("ids", emap);
  }

  if (nb->get_database()->needs_shared_node_information())
  {
    std::vector<int> owningProcessor(numNodes, procId);
    nb->put_field_data("owning_processor", owningProcessor);
  }

  // generate coordinates
  {
    std::vector< double > xyz;
    xyz.reserve(FACET::DIM*numNodes);

    for ( auto&& facet : facets ) {
      for ( int j = 0; j < nodesPerFacet; ++j ) {
        const stk::math::Vector3d & vert = facet.facet_vertex(j);
        xyz.push_back(vert[0]);
        xyz.push_back(vert[1]);
        if (3 == FACET::DIM) xyz.push_back(vert[2]);
      }
    }
    nb->put_field_data("mesh_model_coordinates", xyz);
  }

  // generate connectivity
  {
    std::vector< int > conn;
    conn.reserve(numNodes);
    for ( int n = 0; n < numNodes; ++n ) {
      conn.push_back(nodesPerFacet*elem_start + n+1);
    }
    eb->put_field_data("connectivity", conn);
  }

  io.end_mode(Ioss::STATE_MODEL);

  // write fake transient
  io.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
  io.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

  io.begin_mode(Ioss::STATE_TRANSIENT);
  const int currentOutputStep = io.add_state(1.0*fileIndex);
  io.begin_state(currentOutputStep);
  io.end_state(currentOutputStep);
  io.end_mode(Ioss::STATE_TRANSIENT);
}

stk::mesh::PartVector turn_off_output_for_empty_io_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector)
{
  stk::mesh::PartVector emptyParts;
  for (auto * part : mesh.mesh_meta_data().get_parts())
  {
    if (stk::io::is_part_io_part(*part))
    {
      uint64_t numEntities = stk::mesh::count_selected_entities(*part & outputSelector, mesh.buckets(part->primary_entity_rank()));
      const uint64_t localNumEntities = numEntities;
      stk::all_reduce_sum(mesh.parallel(), &localNumEntities, &numEntities, 1);
      if(numEntities == 0)
      {
        krinolog << "Skipping output of empty part " << part->name() << stk::diag::dendl;
        emptyParts.push_back(part);
        stk::io::remove_io_part_attribute(*part);
      }
    }
  }
  return emptyParts;
}

// Explicit template instantiation

template void write_facets<Facet2d>( const std::vector<Facet2d> & facets, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm);
template void write_facets<Facet3d>( const std::vector<Facet3d> & facets, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm);

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_ */
