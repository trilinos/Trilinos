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

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include <stk_io/DatabasePurpose.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetsFromSides.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Ioss_SubSystem.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_DynamicTopologyFileControl.h>
#include <stk_util/environment/FileUtils.hpp>

namespace krino {

bool is_parallel_io_enabled()
{
#ifdef KRINO_HAVE_PARALELL_IO
  return true;
#else
  return false;
#endif
}

bool does_file_exist(const std::string& filename)
{
  std::ifstream f(filename);
  return f.good();
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

std::string create_file_name(const std::string & fileBaseName, const int fileIndex)
{
  STK_ThrowAssert(fileIndex >= 0);
  std::string stepFileName = Ioss::DynamicTopologyFileControl::get_linear_database_filename(fileBaseName, fileIndex+1);
  stk::util::filename_substitution(stepFileName);
  return stepFileName;
}

void write_shells_for_surfaces(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const std::vector<stk::mesh::Part*> & surfaces, const stk::mesh::Part & activePart, const std::string & fileBaseName, const int fileIndex)
{
  std::vector<std::pair<stk::mesh::Selector, stk::mesh::Selector>>
      interfaceAndNegativeSideBlockSelectors;
  for (auto * surface : surfaces)
    interfaceAndNegativeSideBlockSelectors.emplace_back(*surface,
        stk::mesh::selectUnion(mesh.mesh_meta_data().get_blocks_touching_surface(surface)));
  write_facets(build_facet_mesh_from_interface_sides(
                   mesh, interfaceAndNegativeSideBlockSelectors, activePart, coordsField),
      fileBaseName,
      fileIndex);
}

// The results database and region for writing facets as step fileIndex of the files fileBaseName.
static std::unique_ptr<Ioss::Region> create_facets_region(const std::string & fileBaseName,
    const int fileIndex,
    const stk::ParallelMachine comm,
    Ioss::PropertyManager properties)
{
  // Type (ExodusII) is hard-wired at this time.
  enable_compose_results_if_supported(properties);
  const std::string fileName = create_file_name(fileBaseName, fileIndex);
  properties.add(Ioss::Property("base_filename", create_file_name(fileBaseName, 0)));
  if (fileIndex > 0)
    properties.add(Ioss::Property("state_offset", fileIndex));
  Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", fileName, Ioss::WRITE_RESULTS, comm, properties);
  STK_ThrowRequireMsg(db != nullptr && db->ok(), "ERROR: Could not open output database '" << fileName << "' of type 'exodus'\n");
  auto io = std::make_unique<Ioss::Region>(db, "FacetRegion");
  io->property_add(Ioss::Property("title", std::string("level set interface facets")));
  return io;
}

// The number of facets on the lower processors, for consistent global element numbering.
static int get_facet_offset(const int numFacets, const stk::ParallelMachine comm)
{
  const unsigned numProcs = stk::parallel_machine_size(comm);
  const int procId = stk::parallel_machine_rank(comm);
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
  return elem_start;
}

static std::string get_shell_element_type(const unsigned dim)
{
  return (3 == dim) ? "trishell3" : "shellline2d2";
}

static void write_fake_transient(Ioss::Region & io, const int fileIndex)
{
  io.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
  io.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

  io.begin_mode(Ioss::STATE_TRANSIENT);
  const int currentOutputStep = io.add_state(1.0 * fileIndex);
  io.begin_state(currentOutputStep);
  io.end_state(currentOutputStep);
  io.end_mode(Ioss::STATE_TRANSIENT);
}

void write_facets(const FacetMesh & facets, const std::string & fileBaseName, const int fileIndex)
{
  const stk::ParallelMachine comm = facets.comm();
  const unsigned dim = facets.spatial_dimension();
  const int procId = stk::parallel_machine_rank(comm);
  const int numFacets = facets.num_owned_facets();
  const size_t numNodes = facets.num_owned_facet_vertices();
  const bool sharedVertices = facets.has_shared_vertices();
  const std::vector<int> owners =
      sharedVertices ? facets.compute_vertex_owners() : std::vector<int>(numNodes, procId);
  const int64_t numOwnedNodes = std::count(owners.begin(), owners.end(), procId);

  // Vertex ids are 64 bit.
  Ioss::PropertyManager properties;
  properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
  properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
  std::unique_ptr<Ioss::Region> io =
      create_facets_region(fileBaseName, fileIndex, comm, properties);
  Ioss::DatabaseIO * db = io->get_database();

  const int elemStart = get_facet_offset(numFacets, comm);
  // Without shared vertices, each facet's vertices are numbered after those of the facets before
  // it, across the processors in order.
  auto node_id = [&](const unsigned v) -> int64_t
  { return sharedVertices ? int64_t(facets.vertex_id(v)) : int64_t(dim) * elemStart + v + 1; };

  io->begin_mode(Ioss::STATE_DEFINE_MODEL);

  Ioss::NodeBlock * nb = new Ioss::NodeBlock(db, "nodeblock_1", numNodes, dim);
  io->add(nb);
  nb->property_add(Ioss::Property("locally_owned_count", numOwnedNodes));

  Ioss::ElementBlock * eb =
      new Ioss::ElementBlock(db, "block_1", get_shell_element_type(dim), numFacets);
  io->add(eb);

  io->end_mode(Ioss::STATE_DEFINE_MODEL);
  io->begin_mode(Ioss::STATE_MODEL);

  std::vector<int64_t> nodeIds;
  std::vector<double> xyz;
  nodeIds.reserve(numNodes);
  xyz.reserve(dim * numNodes);
  for (size_t v = 0; v < numNodes; ++v)
  {
    nodeIds.push_back(node_id(v));
    const stk::math::Vector3d & coords = facets.vertex_coordinates(v);
    for (unsigned d = 0; d < dim; ++d)
      xyz.push_back(coords[d]);
  }
  nb->put_field_data("ids", nodeIds);
  if (db->needs_shared_node_information()) nb->put_field_data("owning_processor", owners);
  nb->put_field_data("mesh_model_coordinates", xyz);

  std::vector<int64_t> elemIds;
  std::vector<int64_t> conn;
  elemIds.reserve(numFacets);
  conn.reserve(static_cast<size_t>(dim) * static_cast<size_t>(numFacets));
  for (int f = 0; f < numFacets; ++f)
  {
    elemIds.push_back(elemStart + f + 1);
    for (unsigned v : facets.facet_vertices(f))
      conn.push_back(node_id(v));
  }
  eb->put_field_data("ids", elemIds);
  eb->put_field_data("connectivity", conn);

  io->end_mode(Ioss::STATE_MODEL);

  write_fake_transient(*io, fileIndex);
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


}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_ */
