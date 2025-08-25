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
#include <Akri_config.hpp>

#include <string>

#include <stk_io/DatabasePurpose.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Ioss_SubSystem.h>
#include <Ioss_PropertyManager.h>

namespace krino {

bool is_parallel_io_enabled()
{
#ifdef KRINO_HAVE_PARALELL_IO
  return true;
#else
  return false;
#endif
}

static void enable_compose_results_if_supported(Ioss::PropertyManager & properties)
{
  if (is_parallel_io_enabled())
    properties.add(Ioss::Property("COMPOSE_RESULTS", 1));
  else
    krinolog << "Skipping parallel composed exodus output in krino because the necessary libraries are not enabled." << stk::diag::dendl;
}

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

static void fix_ownership_for_output_selector(stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector)
{
  fix_face_and_edge_ownership_to_assure_selected_owned_element(mesh, outputSelector);
  fix_node_ownership_to_assure_selected_owned_element(mesh, outputSelector);
}

void fix_ownership_and_output_composed_mesh_with_fields(stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose)
{
  fix_ownership_for_output_selector(mesh, outputSelector);
  output_composed_mesh_with_fields(mesh, outputSelector, fileName, step, time, purpose);
}

void output_composed_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose)
{
  Ioss::PropertyManager properties;
  enable_compose_results_if_supported(properties);
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
  enable_compose_results_if_supported(properties);
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

static std::ofstream open_stl_file_and_write_header(const std::string &filename, const unsigned long numTris)
{
  std::ofstream fstream;
  std::string header_info = "binary solid " + filename + "-output";
  char head[80];
  std::strncpy(head, header_info.c_str(), sizeof(head) - 1);

  fstream.open(filename.c_str(), std::ios::out | std::ios::binary);
  fstream.write(head, sizeof(head));

  fstream.write((char*) &numTris, 4);
  return fstream;
}

static void write_facet(std::ofstream & fstream, const Facet3d & facet)
{
  constexpr uint16_t volAttribute = 1;
  const std::array<stk::math::Float3d, 3> facetNodeLocs {facet.facet_vertex(0).data(), facet.facet_vertex(1).data(),
    facet.facet_vertex(2).data()};
  const stk::math::Float3d normal(facet.facet_normal().data());

  fstream.write((char*) &normal[0], 4);
  fstream.write((char*) &normal[1], 4);
  fstream.write((char*) &normal[2], 4);
  for(const stk::math::Float3d &nodeLoc : facetNodeLocs)
  {
    fstream.write((char*) &nodeLoc[0], 4);
    fstream.write((char*) &nodeLoc[1], 4);
    fstream.write((char*) &nodeLoc[2], 4);
  }
  fstream.write((char*) &volAttribute, 2);
}

std::string stl_parallel_file_name(const std::string &fileBaseName, const stk::ParallelMachine comm)
{
  const std::string filename = (stk::parallel_machine_size(comm)==1) ?
    (fileBaseName + ".stl") :
    (fileBaseName + "." + std::to_string(stk::parallel_machine_size(comm)) + "." + std::to_string(stk::parallel_machine_rank(comm)) + ".stl");
  return filename;
}

void write_stl(const std::string &filename, const std::vector<const Facet3d *> & facets)
{
  std::ofstream myfile = open_stl_file_and_write_header(filename, facets.size());

  for(const auto * facet : facets)
    write_facet(myfile, *facet);

  myfile.close();
}

void write_stl(const std::string &filename, const std::vector<Facet3d> & facets)
{
  std::ofstream myfile = open_stl_file_and_write_header(filename, facets.size());

  for(const auto & facet : facets)
    write_facet(myfile, facet);

  myfile.close();
}

void write_stl(const std::string &filename,
  const std::vector<stk::math::Vector3d> & vertices,
  const std::vector<std::array<unsigned,3>> & facetConnectivity)
{
  std::ofstream myfile = open_stl_file_and_write_header(filename, facetConnectivity.size());

  for(const auto & facetConn : facetConnectivity)
    write_facet(myfile, Facet3d(vertices[facetConn[0]], vertices[facetConn[1]], vertices[facetConn[2]]));

  myfile.close();
}

// Explicit template instantiation

template void write_facets<Facet2d>( const std::vector<Facet2d> & facets, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm);
template void write_facets<Facet3d>( const std::vector<Facet3d> & facets, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm);

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_ */
