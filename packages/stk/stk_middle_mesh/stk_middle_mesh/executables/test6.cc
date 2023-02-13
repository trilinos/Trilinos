#include "stk_io/DatabasePurpose.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_topology/topology.hpp"
#include <stk_io/StkMeshIoBroker.hpp>
#include <string>

void create_example_mesh(MPI_Comm comm, const std::string& fname)
{
  stk::io::StkMeshIoBroker stkIo(comm);
  auto index = stkIo.add_mesh_database("generated:3x3x3|sideset:X", stk::io::READ_MESH);
  stkIo.set_active_mesh(index);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  auto fh = stkIo.create_output_mesh(fname, stk::io::WRITE_RESULTS);
  stkIo.write_output_mesh(fh);
}

void print_node_coordinates(stk::mesh::MetaData& metaData, stk::mesh::BulkData& bulkData)
{
  const stk::mesh::FieldBase& coordField = *(metaData.coordinate_field());

  stk::mesh::Selector select             = *(metaData.get_part("surface_1"));
  const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, select);
  for (stk::mesh::Bucket* bucket : buckets)
  {
    std::cout << "number of values per node = " << stk::mesh::field_scalars_per_entity(coordField, *bucket)
              << std::endl;

    double* coordsB = static_cast<double*>(stk::mesh::field_data(coordField, *bucket));
    for (size_t i = 0; i < bucket->size(); i++)
      std::cout << "coordinates = " << coordsB[3 * i] << ", " << coordsB[3 * i + 1] << ", " << coordsB[3 * i + 2]
                << std::endl;
  }
}

void print_face_coordinates(stk::mesh::MetaData& metaData, stk::mesh::BulkData& bulkData)
{
  const stk::mesh::FieldBase& coordField = *(metaData.coordinate_field());

  stk::mesh::Selector select             = *(metaData.get_part("surface_1"));
  const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::FACE_RANK, select);

  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& face : *bucket)
    {
      std::cout << "\nface " << face << std::endl;
      for (auto it = bulkData.begin_nodes(face); it != bulkData.end_nodes(face); ++it)
      {
        double* coordsN = static_cast<double*>(stk::mesh::field_data(coordField, *it));
        std::cout << "node " << *it << " coords = " << coordsN[0] << ", " << coordsN[1] << ", " << coordsN[2]
                  << std::endl;
      }
    }
}

std::string get_unused_part_name(stk::mesh::MetaData& metaData, const std::string& baseName)
{
  const auto& parts = metaData.get_parts();
  std::string currName(baseName);
  int suffix         = 1;
  bool foundFreeName = false;
  while (!foundFreeName)
  {
    for (auto& part : parts)
      if (part->name() == currName)
      {
        currName = baseName + "_" + std::to_string(suffix++);
        continue;
      }

    foundFreeName = true;
  }

  return currName;
}

void print_part_names(stk::mesh::MetaData& metaData)
{
  const auto& parts = metaData.get_parts();
  std::cout << "Part names:" << std::endl;
  for (auto& part : parts)
    std::cout << part->name() << std::endl;
}

int get_idx(const int i, const int j, const int nx, const int ny)
{
  return i * ny + j;
}

void populate_new_surface(stk::mesh::MetaData& metaData, stk::mesh::BulkData& bulkData, stk::mesh::Part& part)
{
  const int nvertX = 3;
  const int nvertY = 3;
  const int nelX   = 2;
  const int nelY   = 2;
  const int nvert  = nvertX * nvertY;
  const int nelem  = nelX * nelY;

  bulkData.modification_begin();

  // generate elements
  std::vector<size_t> requests(metaData.entity_rank_count(), 0);
  requests[stk::topology::NODE_RANK]    = nvert;
  requests[stk::topology::ELEMENT_RANK] = nelem;
  std::vector<stk::mesh::Entity> entities;
  bulkData.generate_new_entities(requests, entities);

  // set topology and downward relations of faces
  std::vector<stk::mesh::Part*> partVec{&part};
  for (int i = 0; i < nelY; ++i)
    for (int j = 0; j < nelX; ++j)
    {
      std::cout << "element idx = " << get_idx(i, j, nelY, nelX) + nvert << std::endl;
      stk::mesh::Entity el = entities[get_idx(i, j, nelY, nelX) + nvert];
      bulkData.change_entity_parts(el, partVec);
      std::cout << "element " << i << ", " << j << " has verts " << get_idx(i, j, nvertY, nvertX) << ", "
                << get_idx(i, j + 1, nvertY, nvertX) << ", " << get_idx(i + 1, j + 1, nvertY, nvertX) << ", "
                << get_idx(i + 1, j, nvertY, nvertX) << std::endl;

      auto v1 = entities[get_idx(i, j, nvertY, nvertX)];
      auto v2 = entities[get_idx(i, j + 1, nvertY, nvertX)];
      auto v3 = entities[get_idx(i + 1, j + 1, nvertY, nvertX)];
      auto v4 = entities[get_idx(i + 1, j, nvertY, nvertX)];
      std::cout << "v1 = " << v1 << std::endl;
      std::cout << "v2 = " << v2 << std::endl;
      std::cout << "v3 = " << v3 << std::endl;
      std::cout << "v4 = " << v4 << std::endl;
      std::cout << "el = " << el << std::endl;

      bulkData.declare_relation(el, v1, 0);
      bulkData.declare_relation(el, v2, 1);
      bulkData.declare_relation(el, v3, 2);
      bulkData.declare_relation(el, v4, 3);
    }

  std::cout << "finished" << std::endl;

  bulkData.modification_end();
  std::cout << "after modification_end()" << std::endl;

  // set coordinates
  const stk::mesh::FieldBase& coordField = *(metaData.coordinate_field());
  const double dx                        = 1;
  const double dy                        = 1;
  const double z                         = 10;
  for (int i = 0; i < nvertY; ++i)
    for (int j = 0; j < nvertX; ++j)
    {
      double* coords = static_cast<double*>(stk::mesh::field_data(coordField, entities[get_idx(i, j, 3, 3)]));
      coords[0]      = j * dx;
      coords[1]      = i * dy;
      coords[2]      = z;
    }

  stk::io::put_io_part_attribute(part);
}

void populate_new_field(stk::mesh::MetaData& metaData, stk::mesh::BulkData& bulkData, stk::mesh::Part& part,
                        stk::mesh::Field<size_t>& field)
{
  stk::mesh::Selector partSel = part;

  const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEMENT_RANK, partSel);
  for (stk::mesh::Bucket* bucket : buckets)
  {
    std::cout << "number of values per node = " << stk::mesh::field_scalars_per_entity(field, *bucket) << std::endl;
    std::cout << "bucket size = " << bucket->size() << std::endl;

    size_t* data = stk::mesh::field_data(field, *bucket);
    for (size_t i = 0; i < bucket->size(); i++)
      data[i] = bulkData.identifier((*bucket)[i]);
  }
}

// TODO: need to test that all fields present on the input mesh also
//       get output
void write_mesh(const std::string& fname, stk::mesh::BulkData& bulkData, stk::mesh::FieldBase& field)
{
  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(bulkData);
  auto fh = stkIo.create_output_mesh(fname, stk::io::WRITE_RESULTS);
  stkIo.add_field(fh, field);
  stkIo.write_output_mesh(fh);

  stkIo.begin_output_step(fh, -1);
  stkIo.write_defined_output_fields(fh);
  stkIo.end_output_step(fh);
}

int main(int argc, char* argv[])
{
  using FieldType = stk::mesh::Field<size_t>;
  MPI_Init(&argc, &argv);
  if (argc > 2)
    std::cerr << "usage: " << argv[0] << " [fname]" << std::endl;

  MPI_Comm comm           = MPI_COMM_WORLD;
  const std::string fname = argc == 1 ? "example.exo" : argv[1];
  std::cout << "fname = " << fname << std::endl;
  if (argc == 1)
  {
    std::cout << "creating mesh" << std::endl;
    create_example_mesh(comm, fname);
  }

  const int spatialDim = 3;
  auto bulkDataPtr     = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(spatialDim).create();
  auto metaDataPtr     = bulkDataPtr->mesh_meta_data_ptr();

  auto& metaData = *metaDataPtr;
  auto& bulkData = *bulkDataPtr;
  stk::mesh::Part* partP;
  FieldType* gidFieldP;

  {
    stk::io::StkMeshIoBroker reader(comm);
    reader.set_bulk_data(bulkData);
    reader.add_mesh_database(fname, stk::io::READ_MESH);
    reader.create_input_mesh();

    std::cout << "\ninitially" << std::endl;
    print_part_names(metaData);
    auto partName = get_unused_part_name(metaData, "new_surface");
    partP         = &(metaData.declare_part_with_topology(partName, stk::topology::SHELL_QUAD_4));

    gidFieldP = &(metaData.declare_field<FieldType>(stk::topology::ELEM_RANK, "gid_field"));
    stk::mesh::put_field_on_mesh(*gidFieldP, *partP, 0);

    std::cout << "\nafter declaring new part" << std::endl;
    print_part_names(metaData);

    reader.populate_bulk_data();
    //    std::cout << "\nafter populating bulk data" << std::endl;
    //    print_part_names(meta_data);
  }

  // write_mesh("input.exo", bulk_data);

  stk::mesh::Selector selectAll = metaData.universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(selectAll, bulkData, entityCounts);
  std::cout << "number of elements = " << entityCounts[stk::topology::ELEMENT_RANK] << std::endl;
  for (int i = stk::topology::BEGIN_RANK; i < stk::topology::END_RANK; ++i)
    std::cout << "number of entities of rank " << i << " = " << entityCounts[i] << std::endl;

  // print_node_coordinates(meta_data, bulk_data);
  print_face_coordinates(metaData, bulkData);
  populate_new_surface(metaData, bulkData, *partP);
  populate_new_field(metaData, bulkData, *partP, *gidFieldP);
  write_mesh("output.exo", bulkData, *gidFieldP);
  return 0;
}
