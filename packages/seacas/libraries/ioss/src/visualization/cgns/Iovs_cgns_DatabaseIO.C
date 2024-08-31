// Copyright(C) 1999-2021, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <visualization/cgns/Iovs_cgns_DatabaseIO.h>
#include <visualization/utils/Iovs_Utils.h>

#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_State.h"
#include "Ioss_StructuredBlock.h"
#include "Ioss_Utils.h"

namespace Iovs_cgns {

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, Ioss_MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region,
                         Iovs::Utils::getInstance().getDatabaseOutputFilePath(filename, props),
                         db_usage, communicator, props)
  {

    Iovs::Utils::DatabaseInfo dbinfo;
    dbinfo.databaseFilename   = this->DBFilename;
    dbinfo.separatorCharacter = std::string(1, this->get_field_separator());
    dbinfo.parallelUtils      = &this->util();
    isIdOutputCreated         = false;

    Iovs::Utils::getInstance().checkDbUsage(db_usage);
    Iovs::Utils::getInstance().createDatabaseOutputFile(dbinfo);
    dbState = Ioss::STATE_UNKNOWN;
    Iovs::Utils::getInstance().writeToCatalystLogFile(dbinfo, props);
    this->catCGNSMesh = Iovs::Utils::getInstance().createCatalystCGNSMesh(dbinfo, props);
  }

  DatabaseIO::~DatabaseIO() { this->catCGNSMesh->Delete(); }

  bool DatabaseIO::begin_nl(Ioss::State /*state*/) { return true; }

  bool DatabaseIO::end_nl(Ioss::State state)
  {
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL: {
      write_meta_data();
      break;
    }
    default: break;
    }
    return true;
  }

  bool DatabaseIO::begin_state_nl(int state, double time)
  {
    this->catCGNSMesh->ReleaseMemory();
    if (!isIdOutputCreated) {
      createIdOutput();
      isIdOutputCreated = true;
    }
    this->catCGNSMesh->SetTimeData(time, state - 1);
    return true;
  }

  bool DatabaseIO::end_state_nl(int /*state*/, double /*time*/)
  {
    std::vector<int>         error_codes;
    std::vector<std::string> error_messages;

    this->catCGNSMesh->logMemoryUsageAndTakeTimerReading();
    this->catCGNSMesh->PerformCoProcessing(error_codes, error_messages);
    this->catCGNSMesh->logMemoryUsageAndTakeTimerReading();
    Iovs::Utils::getInstance().reportCatalystErrorMessages(error_codes, error_messages,
                                                           this->parallel_rank());
    return true;
  }

  void DatabaseIO::read_meta_data_nl() {}

  void DatabaseIO::write_meta_data()
  {
    const auto &structured_blocks = get_region()->get_structured_blocks();
    int         zone              = 1;
    for (const auto &sb : structured_blocks) {
      sb->property_update("zone", zone);
      zone++;
    }
  }

  void DatabaseIO::createIdOutput()
  {
    const auto &structured_blocks = get_region()->get_structured_blocks();
    for (const auto &sb : structured_blocks) {
      size_t               node_count = sb->get_property("node_count").get_int();
      std::vector<int64_t> ids;
      ids.resize(node_count);
      sb->get_cell_node_ids(Data(ids), true);
      bool isCellField = false;
      outputId("cell_node_ids", ids, isCellField, sb);

      size_t cell_count = sb->get_property("cell_count").get_int();
      ids.resize(cell_count);
      sb->get_cell_ids(Data(ids), true);
      isCellField = true;
      outputId("cell_ids", ids, isCellField, sb);

      auto                 zone_id = sb->get_property("zone").get_int();
      std::vector<int64_t> object_id(cell_count, zone_id);
      outputId("object_id", object_id, isCellField, sb);
    }
  }

  void DatabaseIO::outputId(const std::string idName, std::vector<int64_t> &ids, bool isCellField,
                            const Ioss::StructuredBlock *sb)
  {
    CatalystCGNSMeshBase::ZoneData zoneData;
    initZoneDataFromStructuredBlock(zoneData, sb);
    zoneData.data_name     = idName;
    zoneData.comp_count    = 1;
    zoneData.is_cell_field = isCellField;
    zoneData.size          = ids.size();
    zoneData.data_type     = zoneData.T_INT64;
    zoneData.data.p_int64  = Data(ids);

    this->catCGNSMesh->AddStructuredZoneData(zoneData);
  }

  void DatabaseIO::initZoneDataFromStructuredBlock(CatalystCGNSMeshBase::ZoneData &zoneData,
                                                   const Ioss::StructuredBlock    *sb) const
  {
    zoneData.zone_id   = sb->get_property("zone").get_int();
    zoneData.zone_name = sb->name();
    zoneData.ni        = sb->get_property("ni").get_int();
    zoneData.nj        = sb->get_property("nj").get_int();
    zoneData.nk        = sb->get_property("nk").get_int();
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t node_count = sb->get_property("node_count").get_int();
    size_t num_to_get = field.verify(data_size);
    auto   var_type   = field.transformed_storage();
    auto   ioss_type  = field.get_type();
    int    comp_count = var_type->component_count();
    void  *rdata      = num_to_get > 0 ? data : nullptr;

    bool is_cell_field = true;
    if (node_count == num_to_get) {
      is_cell_field = false;
    }

    CatalystCGNSMeshBase::ZoneData zoneData;
    initZoneDataFromStructuredBlock(zoneData, sb);
    zoneData.data_name     = field.get_name();
    zoneData.comp_count    = comp_count;
    zoneData.is_cell_field = is_cell_field;
    zoneData.size          = num_to_get;

    if (ioss_type == Ioss::Field::INTEGER) {
      zoneData.data_type  = zoneData.T_INT;
      zoneData.data.p_int = static_cast<int *>(data);
    }
    else if (ioss_type == Ioss::Field::INT64) {
      zoneData.data_type    = zoneData.T_INT64;
      zoneData.data.p_int64 = static_cast<int64_t *>(data);
    }
    else {
      zoneData.data_type     = zoneData.T_DOUBLE;
      zoneData.data.p_double = static_cast<double *>(data);
    }

    if (field.get_name() == "mesh_model_coordinates") {
      int phys_dimension = get_region()->get_property("spatial_dimension").get_int();

      std::vector<double> coord(num_to_get);

      // ========================================================================
      // Repetitive code for each coordinate direction; use a lambda to consolidate...
      auto coord_lambda = [=, &zoneData, &coord](const char *ordinate, int ordinal) {
        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'
        // Map to global coordinate position...
        for (size_t i = 0; i < num_to_get; i++) {
          coord[i] = static_cast<double *>(rdata)[phys_dimension * i + ordinal];
        }

        zoneData.data_name     = ordinate;
        zoneData.data.p_double = Data(coord);
        this->catCGNSMesh->AddStructuredZoneData(zoneData);
      };
      // ========================================================================

      coord_lambda("mesh_model_coordinates_x", 0);

      if (phys_dimension >= 2) {
        coord_lambda("mesh_model_coordinates_y", 1);
      }

      if (phys_dimension == 3) {
        coord_lambda("mesh_model_coordinates_z", 2);
      }
    }
    else {
      this->catCGNSMesh->AddStructuredZoneData(zoneData);
    }
    return num_to_get;
  }

} // namespace Iovs_cgns
