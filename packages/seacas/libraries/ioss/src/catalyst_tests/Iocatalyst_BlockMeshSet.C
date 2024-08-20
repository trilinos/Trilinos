// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_CopyDatabase.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_IOFactory.h>
#include <Ioss_MeshCopyOptions.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Utils.h>
#include <catalyst/Iocatalyst_DatabaseIO.h>
#include <catalyst/Iocatalyst_CatalystManager.h>
#include <catalyst_tests/Iocatalyst_BlockMeshSet.h>
#include <unordered_set>

namespace Iocatalyst {

  void BlockMeshSet::addBlockMesh(const BlockMesh &blockMesh) { bms.push_back(blockMesh); }

  int BlockMeshSet::getNumLocalPointsInMeshSet()
  {
    std::unordered_set<BlockMesh::ID> pids;
    for (auto bm : bms) {
      BlockMesh::IDList ids = bm.getPartitionPointIDs();
      for (auto id : ids) {
        auto gid = bm.getGlobalIDForPointID(id);
        if (pids.find(gid) == pids.end()) {
          pids.insert(gid);
        }
      }
    }
    return pids.size();
  }

  void BlockMeshSet::writeCatalystIOSSFile(IOSSparams &iop)
  {
    CatalystManager::getInstance().reset();
    iop.isCatalyst = true;

    //Cat Writes
    writeIOSSFile(iop);
    Ioss::PropertyManager cdbProps = Ioss::PropertyManager(iop.dbProps);
    cdbProps.add(Ioss::Property("CATALYST_CONDUIT_NODE", iop.getCatalystConduitNode()));

    //Cat Reads here
    Ioss::DatabaseIO *cdbi =
        Ioss::IOFactory::create(CATALYST_DATABASE_TYPE, CATALYST_DUMMY_DATABASE, Ioss::READ_RESTART,
                                Ioss::ParallelUtils::comm_world(), cdbProps);
    if (cdbi == nullptr || !cdbi->ok(true)) {
      return;
    }
    Ioss::PropertyManager properties = Ioss::PropertyManager(iop.dbProps);
    Ioss::DatabaseIO *cdbo = Ioss::IOFactory::create(iop.dbType, iop.fileName, Ioss::WRITE_RESULTS,
                                                     Ioss::ParallelUtils::comm_world(), properties);

    if (cdbo == nullptr || !cdbo->ok(true)) {
      std::ostringstream errmsg;
      errmsg << "Unable to open IOSS database " + iop.fileName + "\n";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region          cor(cdbo);
    Ioss::Region          cir(cdbi);
    Ioss::MeshCopyOptions options;
    options.data_storage_type = 1;
    Ioss::copy_database(cir, cor, options);
  }

  Ioss::DatabaseIO* BlockMeshSet::getCatalystDatabase(IOSSparams &iop)
  {
    CatalystManager::getInstance().reset();
    iop.isCatalyst = true;
    
    //Write to Cat database
    writeIOSSFile(iop);

    Ioss::PropertyManager cdbProps = Ioss::PropertyManager(iop.dbProps);

    //Get Conduit
    cdbProps.add(Ioss::Property("CATALYST_CONDUIT_NODE", iop.getCatalystConduitNode()));

    //Read to Cat Database
    Ioss::DatabaseIO *cdbi =
        Ioss::IOFactory::create(CATALYST_DATABASE_TYPE, CATALYST_DUMMY_DATABASE, Ioss::READ_RESTART,
                                Ioss::ParallelUtils::comm_world(), cdbProps);
    if (cdbi == nullptr || !cdbi->ok(true)) {
      return nullptr;
    }

    return cdbi;
  }

  void BlockMeshSet::writeIOSSFile(IOSSparams &iop)
  {
    openIOSSDatabase(iop);
    switchStateDefineModel(iop);
    switchStateModel(iop);
    switchStateDefineTransient(iop);
    switchStateTransient(iop);
    closeIOSSDatabase(iop);
  }

  void BlockMeshSet::openIOSSDatabase(IOSSparams &iop)
  {
    Ioss::PropertyManager properties = Ioss::PropertyManager(iop.dbProps);
    std::string           dbType = iop.dbType;
    if (iop.isCatalyst) {
      dbType = CATALYST_DATABASE_TYPE;
    }
    iop.databaseIO = Ioss::IOFactory::create(dbType, iop.fileName, Ioss::WRITE_RESULTS,
                                             Ioss::ParallelUtils::comm_world(), properties);

    if (iop.databaseIO == nullptr || !iop.databaseIO->ok(true)) {
      std::ostringstream errmsg;
      errmsg << "Unable to open IOSS database " + iop.fileName + "\n";
      IOSS_ERROR(errmsg);
    }
    iop.region = std::unique_ptr<Ioss::Region>(new Ioss::Region(iop.databaseIO));
  }

  void BlockMeshSet::closeIOSSDatabase(IOSSparams &iop)
  {
    iop.region.reset();
    iop.databaseIO = nullptr;
  }

  void BlockMeshSet::switchStateDefineModel(IOSSparams &iop)
  {
    iop.region->begin_mode(Ioss::STATE_DEFINE_MODEL);
    if (iop.isStructured()) {
      writeStructuredBlockDefinitions(iop);
    }
    else {
      writeUnstructuredBlockDefinitions(iop);
    }
    iop.region->end_mode(Ioss::STATE_DEFINE_MODEL);
  }

  void BlockMeshSet::switchStateModel(IOSSparams &iop)
  {
    iop.region->begin_mode(Ioss::STATE_MODEL);
    if (iop.isStructured()) {
      writeStructuredBlockBulkData(iop);
    }
    else {
      writeUnstructuredBlockBulkData(iop);
    }
    if (iop.isCatalyst) {
      saveConduitNode(iop);
    }
    iop.region->end_mode(Ioss::STATE_MODEL);
  }

  void BlockMeshSet::switchStateDefineTransient(IOSSparams &iop)
  {
    iop.region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
    if (iop.isStructured()) {
      writeStructuredTransientFieldDefinitions(iop);
    }
    else {
      writeUnstructuredTransientFieldDefinitions(iop);
    }
    iop.region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }
  void BlockMeshSet::switchStateTransient(IOSSparams &iop)
  {
    iop.region->begin_mode(Ioss::STATE_TRANSIENT);
    int tstep = iop.region->add_state(0.0);
    iop.region->begin_state(tstep);
    if (iop.isStructured()) {
      writeStructuredTransientBulkData(iop);
    }
    else {
      writeUnstructuredTransientBulkData(iop);
    }
    iop.region->end_state(tstep);
    if (iop.isCatalyst) {
      saveConduitNode(iop);
    }
    iop.region->end_mode(Ioss::STATE_TRANSIENT);
  }

  void BlockMeshSet::saveConduitNode(IOSSparams &iop)
  {
    auto c_node = reinterpret_cast<conduit_node *>(
        ((Iocatalyst::DatabaseIO *)iop.databaseIO)->get_catalyst_conduit_node());
    auto cpp_node = conduit_cpp::cpp_node(c_node);
    iop.conduitNode.set(cpp_node);
  }

  void BlockMeshSet::writeStructuredBlockDefinitions(IOSSparams &iop)
  {
    int spatialDims = 3;
    for (auto bm : bms) {
      Ioss::IJK_t parentOffsets = {{(int)bm.getPartitionStart().i, (int)bm.getPartitionStart().j,
                                    (int)bm.getPartitionStart().k}};
      Ioss::IJK_t globalSizes   = {
          {(int)bm.getExtents().i, (int)bm.getExtents().j, (int)bm.getExtents().k}};
      Ioss::IJK_t localSizes = {{(int)bm.getPartitionExtents().i, (int)bm.getPartitionExtents().j,
                                 (int)bm.getPartitionExtents().k}};
      Ioss::StructuredBlock *iossBlock =
          new Ioss::StructuredBlock(iop.databaseIO, getStructuredBlockName(bm.getID()), spatialDims,
                                    localSizes, parentOffsets, globalSizes);
      int node_count = (bm.getPartitionExtents().i + 1) * (bm.getPartitionExtents().j + 1) *
                       (bm.getPartitionExtents().k + 1);
      Ioss::NodeBlock *nodeBlock = new Ioss::NodeBlock(
          iop.databaseIO, getStructuredNodeBlockName(bm.getID()), node_count, spatialDims);
      iop.region->add(iossBlock);
      iop.region->add(nodeBlock);
    }
  }

  void BlockMeshSet::writeStructuredBlockBulkData(IOSSparams &iop)
  {
    std::vector<double> coordx;
    std::vector<double> coordy;
    std::vector<double> coordz;
    BlockMesh::Point    origin;
    origin.x = 0.0;
    origin.y = 0.0;
    origin.z = 0.0;
    for (auto bm : bms) {
      const int numI      = bm.getPartitionExtents().i + 1;
      const int numJ      = bm.getPartitionExtents().j + 1;
      const int numK      = bm.getPartitionExtents().k + 1;
      const int numPoints = numI * numJ * numK;

      coordx.resize(numPoints);
      coordy.resize(numPoints);
      coordz.resize(numPoints);

      for (int k = 0; k < numK; ++k) {
        const int kOffset = k * numI * numJ;
        for (int j = 0; j < numJ; ++j) {
          const int kjOffset = kOffset + j * numI;
          for (int i = 0; i < numI; ++i) {
            coordx[kjOffset + i] =
                i * bm.BLOCK_LENGTH + bm.getOrigin().i + bm.getPartitionStart().i;
            coordy[kjOffset + i] =
                j * bm.BLOCK_LENGTH + bm.getOrigin().j + bm.getPartitionStart().j;
            coordz[kjOffset + i] =
                k * bm.BLOCK_LENGTH + bm.getOrigin().k + bm.getPartitionStart().k;
          }
        }
      }

      auto iossBlock = iop.region->get_structured_block(getStructuredBlockName(bm.getID()));
      iossBlock->put_field_data("mesh_model_coordinates_x", coordx);
      iossBlock->put_field_data("mesh_model_coordinates_y", coordy);
      iossBlock->put_field_data("mesh_model_coordinates_z", coordz);
    }
  }

  void BlockMeshSet::writeStructuredTransientFieldDefinitions(IOSSparams &iop)
  {
    for (auto bm : bms) {
      //Modify this to access field dict in "bm" and populate ioss block with those fields (as well).
      auto iossBlock = iop.region->get_structured_block(getStructuredBlockName(bm.getID()));
      iossBlock->field_add(Ioss::Field(IOSS_CELL_FIELD, Ioss::Field::REAL, IOSS_SCALAR_STORAGE,
                                       Ioss::Field::TRANSIENT));
      iossBlock->get_node_block().field_add(Ioss::Field(
          IOSS_POINT_FIELD, Ioss::Field::REAL, IOSS_SCALAR_STORAGE, Ioss::Field::TRANSIENT));
    }
  }

  void BlockMeshSet::writeStructuredTransientBulkData(IOSSparams &iop)
  {

    std::vector<double> values;
    for (auto bm : bms) {
      values.clear();
      auto iossBlock = iop.region->get_structured_block(getStructuredBlockName(bm.getID()));
      for (int j = 0; j < iossBlock->get_field(IOSS_CELL_FIELD).raw_count(); j++) {
        values.push_back(bm.getPartition().id);
      }
      iossBlock->put_field_data(IOSS_CELL_FIELD, values);

      values.clear();
      for (int j = 0; j < iossBlock->get_node_block().get_field(IOSS_POINT_FIELD).raw_count();
           j++) {
        values.push_back(bm.getPartition().id);
      }
      iossBlock->get_node_block().put_field_data(IOSS_POINT_FIELD, values);
    }
  }

  void BlockMeshSet::writeUnstructuredBlockDefinitions(IOSSparams &iop)
  {
    int              spatialDims = 3;
    Ioss::NodeBlock *nodeBlock =
        new Ioss::NodeBlock(iop.databaseIO, "nodeblock", getNumLocalPointsInMeshSet(), spatialDims);
    iop.region->add(nodeBlock);
    for (auto bm : bms) {
      Ioss::ElementBlock *elemBlock = new Ioss::ElementBlock(
          iop.databaseIO, getUnstructuredBlockName(bm.getID()), "hex8", bm.getNumPartitionBlocks());
      iop.region->add(elemBlock);
    }
  }

  void BlockMeshSet::writeUnstructuredBlockBulkData(IOSSparams &iop)
  {
    Ioss::NodeBlock                  *nodeBlock = iop.region->get_node_block("nodeblock");
    std::vector<double>               coordx;
    std::vector<double>               coordy;
    std::vector<double>               coordz;
    BlockMesh::IDList                 globalPointIds;
    std::unordered_set<BlockMesh::ID> pids;
    for (auto bm : bms) {
      BlockMesh::IDList ids = bm.getPartitionPointIDs();
      for (auto id : ids) {
        auto gid = bm.getGlobalIDForPointID(id);
        if (pids.find(gid) == pids.end()) {
          BlockMesh::Point point = bm.getPointCoordsForPointID(id);
          coordx.push_back(point.x);
          coordy.push_back(point.y);
          coordz.push_back(point.z);
          globalPointIds.push_back(gid);
          pids.insert(gid);
        }
      }
    }
    nodeBlock->put_field_data("mesh_model_coordinates_x", coordx);
    nodeBlock->put_field_data("mesh_model_coordinates_y", coordy);
    nodeBlock->put_field_data("mesh_model_coordinates_z", coordz);
    nodeBlock->put_field_data("ids", globalPointIds);

    for (auto bm : bms) {
      Ioss::ElementBlock *elemBlock =
          iop.region->get_element_block(getUnstructuredBlockName(bm.getID()));

      std::vector<int>  connectivity(8 * bm.getNumPartitionBlocks());
      BlockMesh::IDList globalElemIds;
      BlockMesh::IDList ids = bm.getPartitionBlockIDs();

      for (int i = 0; i < ids.size(); i++) {
        BlockMesh::BlockConn conn = bm.getBlockConnectivityPointIDs(ids[i]);
        globalElemIds.push_back(bm.getGlobalIDForBlockID(ids[i]));
        for (int j = 0; j < conn.size(); j++) {
          connectivity[(i * 8) + j] = bm.getGlobalIDForPointID(conn[j]);
        }
      }
      elemBlock->put_field_data("connectivity", connectivity);
      elemBlock->put_field_data("ids", globalElemIds);
    }
  }

  void BlockMeshSet::writeUnstructuredTransientFieldDefinitions(IOSSparams &iop)
  {
    iop.region->field_add(Ioss::Field(IOSS_GLOBAL_FIELD, Ioss::Field::REAL, IOSS_SCALAR_STORAGE,
                                      Ioss::Field::TRANSIENT));
    for (auto bm : bms) {
      auto elemBlock = iop.region->get_element_block(getUnstructuredBlockName(bm.getID()));
      elemBlock->field_add(Ioss::Field(IOSS_CELL_FIELD, Ioss::Field::REAL, IOSS_SCALAR_STORAGE,
                                       Ioss::Field::TRANSIENT));
      auto nodeBlock = iop.region->get_node_block("nodeblock");
      nodeBlock->field_add(Ioss::Field(IOSS_POINT_FIELD, Ioss::Field::REAL, IOSS_SCALAR_STORAGE,
                                       Ioss::Field::TRANSIENT));
      
      writeUnstructuredAddedTransientFields(bm, iop);
    }
  }

  void BlockMeshSet::writeUnstructuredAddedTransientFields(BlockMesh bm, IOSSparams &iop)
  {
    writeUnstructuredAddedCellTransientFields(bm, iop);
    writeUnstructuredAddedPointTransientFields(bm, iop);
  }

  void BlockMeshSet::writeUnstructuredAddedCellTransientFields(BlockMesh bm, IOSSparams &iop)
  {
    auto cell_fields = bm.getTransientCellFieldMap();
    auto elemBlock = iop.region->get_element_block(getUnstructuredBlockName(bm.getID()));
    for (auto itr = cell_fields->begin(); itr != cell_fields->end(); ++itr) 
    { 
      elemBlock->field_add(Ioss::Field(itr->first, Ioss::Field::REAL, IOSS_SCALAR_STORAGE,
                                       Ioss::Field::TRANSIENT));
    } 
  }

  void BlockMeshSet::writeUnstructuredAddedPointTransientFields(BlockMesh bm, IOSSparams &iop)
  {
    auto point_fields = bm.getTransientPointFieldMap();
    auto nodeBlock = iop.region->get_node_block("nodeblock");
    for (auto itr = point_fields->begin(); itr != point_fields->end(); ++itr) 
    { 
      nodeBlock->field_add(Ioss::Field(itr->first, Ioss::Field::REAL, IOSS_SCALAR_STORAGE,
                                       Ioss::Field::TRANSIENT));
    } 
  }

  void BlockMeshSet::writeUnstructuredTransientBulkData(IOSSparams &iop)
  {
    std::vector<double> values;
    values.push_back(3.14);
    iop.region->put_field_data(IOSS_GLOBAL_FIELD, values);
    for (auto bm : bms) {
      values.clear();
      auto elemBlock = iop.region->get_element_block(getUnstructuredBlockName(bm.getID()));
      for (int j = 0; j < elemBlock->get_field(IOSS_CELL_FIELD).raw_count(); j++) {
        values.push_back(bm.getPartition().id);
      }
      elemBlock->put_field_data(IOSS_CELL_FIELD, values);

      auto nodeBlock = iop.region->get_node_block("nodeblock");
      values.clear();
      for (int j = 0; j < nodeBlock->get_field(IOSS_POINT_FIELD).raw_count(); j++) {
        values.push_back(bm.getPartition().id);
      }
      nodeBlock->put_field_data(IOSS_POINT_FIELD, values);

      writeUnstructuredAddedTransientFieldsBulkData(bm, iop);
    }
  }

  void BlockMeshSet::writeUnstructuredAddedTransientFieldsBulkData(BlockMesh bm, IOSSparams &iop)
  {
    writeUnstructuredAddedCellTransientFieldsBulkData(bm, iop);
    writeUnstructuredAddedPointTransientFieldsBulkData(bm, iop);
  }

  void BlockMeshSet::writeUnstructuredAddedCellTransientFieldsBulkData(BlockMesh bm, IOSSparams &iop)
  {
    auto cell_fields = bm.getTransientCellFieldMap();
    auto elemBlock = iop.region->get_element_block(getUnstructuredBlockName(bm.getID()));
    std::vector<double> values;
    for (auto itr = cell_fields->begin(); itr != cell_fields->end(); ++itr) 
    {
      int num_elements = elemBlock->get_field(itr->first).raw_count();
      for (int j = 0; j < num_elements; j++) {
        values.push_back(itr->second + j*0.1);
      }
      elemBlock->put_field_data(itr->first, values);
      values.clear();
    } 
  }

  void BlockMeshSet::writeUnstructuredAddedPointTransientFieldsBulkData(BlockMesh bm, IOSSparams &iop)
  {
    auto point_fields = bm.getTransientPointFieldMap();
    auto nodeBlock = iop.region->get_node_block("nodeblock");
    std::vector<double> values;
    for (auto itr = point_fields->begin(); itr != point_fields->end(); ++itr) 
    {
      int num_nodes = nodeBlock->get_field(itr->first).raw_count();
      for (int j = 0; j < num_nodes; j++) {
        values.push_back(itr->second + j*0.1);
      }
      nodeBlock->put_field_data(itr->first, values);
      values.clear();
    } 
  }

  std::string BlockMeshSet::getStructuredBlockName(int index)
  {
    return "StructuredBlock" + std::to_string(index);
  }

  std::string BlockMeshSet::getStructuredNodeBlockName(int index)
  {
    return "StructuredNodeBlock" + std::to_string(index);
  }

  std::string BlockMeshSet::getUnstructuredBlockName(int index)
  {
    return "UnstructuredBlock" + std::to_string(index);
  }

} // namespace Iocatalyst
