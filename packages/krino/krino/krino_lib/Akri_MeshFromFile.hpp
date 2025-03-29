#ifndef KRINO_KRINO_KRINO_LIB_AKRI_MESHFROMFILE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_MESHFROMFILE_HPP_

#include <Akri_MeshInterface.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <memory>

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace io { class StkMeshIoBroker; } }

namespace krino {

class MeshFromFile : public MeshInterface {
public:
  MeshFromFile(const std::string & fileName, stk::ParallelMachine comm, const std::string & decompMethod = "", const bool useAllSidesForShells = true);
  virtual ~MeshFromFile();

  virtual void populate_mesh(const stk::mesh::BulkData::AutomaticAuraOption auto_aura_option = stk::mesh::BulkData::AUTO_AURA) override;
  virtual stk::mesh::MetaData & meta_data() override { STK_ThrowAssert(myMeta); return *myMeta; }
  virtual const stk::mesh::MetaData & meta_data() const override { STK_ThrowAssert(myMeta); return *myMeta; }
  virtual stk::mesh::BulkData & bulk_data() override { STK_ThrowAssert(myMesh); return *myMesh; }
  virtual const stk::mesh::BulkData & bulk_data() const override { STK_ThrowAssert(myMesh); return *myMesh; }

private:
  stk::ParallelMachine myComm;
  std::unique_ptr<stk::io::StkMeshIoBroker> myIOBroker;
  stk::mesh::MetaData * myMeta{nullptr};
  stk::mesh::BulkData * myMesh{nullptr};
};

} // namespace krino

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_MESHFROMFILE_HPP_ */
