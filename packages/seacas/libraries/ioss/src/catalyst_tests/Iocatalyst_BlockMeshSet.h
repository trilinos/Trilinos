// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iocatalyst_export.h"

#include <Ioss_Region.h>
#include <catalyst.hpp>
#include <catalyst_tests/Iocatalyst_BlockMesh.h>
#include <cstddef>
#include <vector>

namespace Iocatalyst {

  class IOCATALYST_EXPORT BlockMeshSet
  {

  public:
    inline static const std::string CATALYST_DATABASE_TYPE  = "catalyst";
    inline static const std::string CATALYST_DUMMY_DATABASE = "dummy.db";

    class IOSSparams
    {
    public:
      IOSSparams(const std::string &fileName, const std::string &dbType, Ioss::PropertyManager dbProps = {})
          : fileName(fileName), dbType(dbType), databaseIO(nullptr), isCatalyst(false), dbProps(dbProps)
      {
      }
      bool              isStructured() { return dbType == CGNS_DATABASE_TYPE; }
      void              printCatalystConduitNode() { conduitNode.print_detailed(); }
      void             *getCatalystConduitNode() { return conduit_cpp::c_node(&conduitNode); }
      std::string       fileName;
      std::string       dbType;
      Ioss::DatabaseIO *databaseIO;
      bool              isCatalyst;
      std::unique_ptr<Ioss::Region> region;
      conduit_cpp::Node             conduitNode;
      Ioss::PropertyManager         dbProps;

    private:
      IOSSparams();
    };

    void addBlockMesh(const BlockMesh &blockMesh);
    void writeIOSSFile(IOSSparams &iop);
    void writeCatalystIOSSFile(IOSSparams &iop);
    Ioss::DatabaseIO* getCatalystDatabase(IOSSparams &iop);
    
    int  getNumLocalPointsInMeshSet();

    std::string getStructuredBlockName(int index);
    std::string getStructuredNodeBlockName(int index);

    std::string getUnstructuredBlockName(int index);

  private:
    std::vector<BlockMesh> bms;

    void openIOSSDatabase(IOSSparams &iop);
    void closeIOSSDatabase(IOSSparams &iop);

    void switchStateDefineModel(IOSSparams &iop);
    void switchStateModel(IOSSparams &iop);
    void switchStateDefineTransient(IOSSparams &iop);
    void switchStateTransient(IOSSparams &iop);

    void writeStructuredBlockDefinitions(IOSSparams &iop);
    void writeStructuredBlockBulkData(IOSSparams &iop);
    void writeStructuredTransientFieldDefinitions(IOSSparams &iop);
    void writeStructuredTransientBulkData(IOSSparams &iop);

    void writeUnstructuredBlockDefinitions(IOSSparams &iop);
    void writeUnstructuredBlockBulkData(IOSSparams &iop);
    void writeUnstructuredTransientFieldDefinitions(IOSSparams &iop);
    void writeUnstructuredAddedTransientFields(BlockMesh bm, IOSSparams &iop);
    void writeUnstructuredAddedCellTransientFields(BlockMesh bm, IOSSparams &iop);
    void writeUnstructuredAddedPointTransientFields(BlockMesh bm, IOSSparams &iop);
    void writeUnstructuredTransientBulkData(IOSSparams &iop);
    void writeUnstructuredAddedTransientFieldsBulkData(BlockMesh bm, IOSSparams &iop);
    void writeUnstructuredAddedCellTransientFieldsBulkData(BlockMesh bm, IOSSparams &iop);
    void writeUnstructuredAddedPointTransientFieldsBulkData(BlockMesh bm, IOSSparams &iop);

    void saveConduitNode(IOSSparams &iop);

    inline static const std::string CGNS_DATABASE_TYPE   = "cgns";
    inline static const std::string EXODUS_DATABASE_TYPE = "exodus";
    inline static const std::string IOSS_CELL_FIELD      = "cell";
    inline static const std::string IOSS_POINT_FIELD     = "point";
    inline static const std::string IOSS_GLOBAL_FIELD    = "global";
    inline static const std::string IOSS_SCALAR_STORAGE  = "scalar";
  };

} // namespace Iocatalyst
