#include <gtest/gtest.h>
#include <unistd.h>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include "stk_mesh/base/FieldBase.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_io/MeshField.hpp"

namespace
{

class QueryExoVars : public stk::unit_test_util::MeshFixture
{
protected:
  void read_meta(stk::io::StkMeshIoBroker &stkIo, const std::string &filename)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    stkIo.set_bulk_data(get_bulk());
    stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.create_input_mesh();
  }

  void expect_names(const std::vector<std::pair<std::string,std::string>> &goldNames,
                    const stk::io::FieldNameToPartVector &varNames)
  {
    ASSERT_EQ(goldNames.size(), varNames.size());
    for(size_t i=0; i<goldNames.size(); i++)
    {
      EXPECT_EQ(goldNames[i].first, varNames[i].first);
      EXPECT_EQ(goldNames[i].second, varNames[i].second->name());
    }
  }
};

TEST_F(QueryExoVars, nodeVars_getNames)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    stk::io::StkMeshIoBroker stkIo;
    read_meta(stkIo, "allTypesOfData.exo");
    expect_names({{"dispx","{UNIVERSAL}"},
                  {"dispy","{UNIVERSAL}"},
                  {"dispz","{UNIVERSAL}"},
                  {"rotx","{UNIVERSAL}"},
                  {"roty","{UNIVERSAL}"},
                  {"rotz","{UNIVERSAL}"}}, stkIo.get_nodal_var_names());
  }
}

TEST_F(QueryExoVars, elemVars_getNames)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    stk::io::StkMeshIoBroker stkIo;
    read_meta(stkIo, "elemData.exo");
    expect_names({{"vonmises","block_1"},
                  {"vonmises","block_11"}}, stkIo.get_elem_var_names());
  }
}

TEST_F(QueryExoVars, nodesetVars_getNames)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    stk::io::StkMeshIoBroker stkIo;
    read_meta(stkIo, "nodesetData.exo");
    expect_names({{"apressure","nodelist_2"},
                  {"dispx","nodelist_1"},
                  {"dispy","nodelist_1"},
                  {"dispz","nodelist_1"}}, stkIo.get_nodeset_var_names());
  }
}

TEST_F(QueryExoVars, sidesetVars_getNames)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    stk::io::StkMeshIoBroker stkIo;
    read_meta(stkIo, "allTypesOfData.exo");
    expect_names({{"appliedpressure_sideset_30","surface_hex8_quad4_30"},
                  {"appliedpressure_sideset_31","surface_hex8_quad4_31"},
                  {"appliedpressure_sideset_32","surface_hex8_quad4_32"}}, stkIo.get_sideset_var_names());
  }
}

}
