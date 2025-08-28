#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Ioss_IOFactory.h>
#include <Ioss_Utils.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <memory>
#include <string>

namespace krino {

MeshFromFile::MeshFromFile(const std::string & fileName, stk::ParallelMachine comm, const std::string & decompMethod, const bool useAllSidesForShells)
: myComm(comm)
{
  myIOBroker = std::make_unique<stk::io::StkMeshIoBroker>(comm);

  myIOBroker->property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 180));

  if (!decompMethod.empty())
    myIOBroker->property_add(Ioss::Property("DECOMPOSITION_METHOD", Ioss::Utils::uppercase(decompMethod)));

  myIOBroker->add_mesh_database(fileName, stk::io::READ_MESH);

  if (useAllSidesForShells)
  {
    std::shared_ptr<Ioss::Region> ioRegion = myIOBroker->get_input_ioss_region();
    if(ioRegion && 3 == ioRegion->get_property("spatial_dimension").get_int())
      myIOBroker->set_enable_all_face_sides_shell_topo(true);
  }

  myIOBroker->create_input_mesh();
  myIOBroker->add_all_mesh_fields_as_input_fields();
  myMeta = &myIOBroker->meta_data();

  AuxMetaData::create(*myMeta);
}

MeshFromFile::~MeshFromFile()
{
}

void MeshFromFile::populate_mesh(const stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  std::shared_ptr<stk::mesh::MetaData> sharedMetaWeWontDelete(myMeta,[](auto /*ptrWeWontDelete*/){});
  std::shared_ptr<stk::mesh::BulkData> sharedBulkThatWeGiveToIoBroker = stk::mesh::MeshBuilder(myComm).set_aura_option(autoAuraOption).create(sharedMetaWeWontDelete);
  myIOBroker->set_bulk_data( sharedBulkThatWeGiveToIoBroker );

  myIOBroker->populate_bulk_data();
  myMesh = & myIOBroker->bulk_data();

  int numSteps = myIOBroker->get_num_time_steps();
  if(numSteps>0)
    myIOBroker->read_defined_input_fields(numSteps);

  stk::mesh::create_exposed_block_boundary_sides(*myMesh, myMeta->universal_part(), {&AuxMetaData::get(*myMeta).exposed_boundary_part()});
  stk::mesh::create_interior_block_boundary_sides(*myMesh, myMeta->universal_part(), {&AuxMetaData::get(*myMeta).block_boundary_part()});
}

} // namespace krino
