/*
 * Akri_OutputUtils.cpp
 *
 *  Created on: Nov 21, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_
#include <Akri_OutputUtils.hpp>

#include <string>

#include <stk_io/DatabasePurpose.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Ioss_SubSystem.h>

namespace krino {

void output_field_with_mesh(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose)
{
  stk::io::StkMeshIoBroker stkIo;
  stk::mesh::BulkData & workAroundNonConstMesh = const_cast<stk::mesh::BulkData &>(mesh);
  stkIo.set_bulk_data(workAroundNonConstMesh);

  size_t outputFileIndex = stkIo.create_output_mesh(fileName, purpose);

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

  stkIo.set_active_selector(activePart);
  stkIo.set_subset_selector(outputFileIndex, activePart);

  stkIo.write_output_mesh(outputFileIndex);

  stkIo.begin_output_step(outputFileIndex, time);
  stkIo.write_defined_output_fields(outputFileIndex);
  stkIo.end_output_step(outputFileIndex);
}

std::string create_file_name(const std::string & fileBaseName, const int fileIndex)
{
  ThrowAssert(fileIndex >= 0);
  if (fileIndex < 1)
    return fileBaseName + ".e";

  ThrowAssert(fileIndex <= 9998);
  char counterChar[32] = {'\0'};
  std::snprintf(counterChar, 32, "%.4d", fileIndex+1);
  const std::string filename = fileBaseName + ".e-s" + std::string(counterChar);
  return filename;
}

void
write_facets( const int dim, const Faceted_Surface & facetedSurface, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm)
{
  int nfaces = facetedSurface.size();
  const int nodes_per_elem = dim;
  const int nnodes = nfaces * nodes_per_elem;
  const int nelems = nfaces * 1; //writing out faces as elements

  // Type (ExodusII) is hard-wired at this time.
  Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", create_file_name(fileBaseName, fileIndex), Ioss::WRITE_RESULTS);
  Ioss::Region io(db, "FacetRegion");

  // find offsets for consistent global numbering
  int elem_start = 0;
  if ( stk::parallel_machine_size(comm) > 1 ) {
    std::vector< int > num_faces;
    num_faces.resize( stk::parallel_machine_size(comm) );

    // Gather everyone's facet list sizes
    // I don't think there is a framework call for this...
    MPI_Allgather( &nfaces, 1, MPI_INT,
                   &(num_faces[0]), 1, MPI_INT, comm );
    for (int i = 0; i< stk::parallel_machine_rank(comm); ++i) {
      elem_start += num_faces[i];
    }
  }

  const std::string description = "level set interface facets";
  io.property_add(Ioss::Property("title", description));
  io.begin_mode(Ioss::STATE_DEFINE_MODEL);

  // if we have no elements bail now
  if ( 0 == nelems ) {
    io.end_mode(Ioss::STATE_DEFINE_MODEL);
    return;
  }

  Ioss::NodeBlock *nb = new Ioss::NodeBlock(db, "nodeblock_1", nnodes, dim);
  io.add(nb);

  std::string el_type;
  if (dim == 3)
    el_type = "trishell3";
  else
    el_type = "shellline2d2";

  Ioss::ElementBlock *eb = new Ioss::ElementBlock(db, "block_1", el_type, nelems);
  io.add(eb);

  io.end_mode(Ioss::STATE_DEFINE_MODEL);
  io.begin_mode(Ioss::STATE_MODEL);

  const FacetOwningVec & facets = facetedSurface.get_facets();

  // loop over elements and node to create maps
  {
    std::vector< int > nmap, emap;
    nmap.reserve(nnodes);
    emap.reserve(nelems);

    for ( unsigned n=0, e=0; e<facets.size(); ++e ) {
      emap.push_back(elem_start + e + 1);
      for ( int j = 0; j < nodes_per_elem; ++j ) {
        nmap.push_back(nodes_per_elem*elem_start + n + 1);
        ++n;
      }
    }
    nb->put_field_data("ids", nmap);
    eb->put_field_data("ids", emap);
  }

  // generate coordinates
  {
    std::vector< double > xyz;
    xyz.reserve(dim*nnodes);

    for ( auto&& facet : facets ) {
      for ( int j = 0; j < nodes_per_elem; ++j ) {
        const Vector3d & vert = facet->facet_vertex(j);
        xyz.push_back(vert[0]);
        xyz.push_back(vert[1]);
        if (3 == dim) xyz.push_back(vert[2]);
      }
    }
    nb->put_field_data("mesh_model_coordinates", xyz);
  }

  // generate connectivity
  {
    std::vector< int > conn;
    conn.reserve(nnodes);
    for ( int n = 0; n < nnodes; ++n ) {
      conn.push_back(nodes_per_elem*elem_start + n+1);
    }
    eb->put_field_data("connectivity", conn);
  }

  io.end_mode(Ioss::STATE_MODEL);
}

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_CPP_ */
