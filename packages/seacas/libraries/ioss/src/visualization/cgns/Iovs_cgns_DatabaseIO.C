// Copyright(C) 1999-2021 National Technology & Engineering Solutions
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
    dbinfo.parallelUtils = &this->util();

    Iovs::Utils::getInstance().checkDbUsage(db_usage);
    Iovs::Utils::getInstance().createDatabaseOutputFile(dbinfo);
    dbState           = Ioss::STATE_UNKNOWN;
    Iovs::Utils::getInstance().writeToCatalystLogFile(dbinfo, props);
    this->catCGNSMesh = Iovs::Utils::getInstance().createCatalystCGNSMesh(dbinfo, props);
  }

  DatabaseIO::~DatabaseIO() { this->catCGNSMesh->Delete(); }

  bool DatabaseIO::begin__(Ioss::State state) { return true; }

  bool DatabaseIO::end__(Ioss::State state)
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

  bool DatabaseIO::begin_state__(int state, double time)
  {
    this->catCGNSMesh->ReleaseMemory();
    this->catCGNSMesh->SetTimeData(time, state - 1);
    return true;
  }

  bool DatabaseIO::end_state__(int state, double time)
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

  void DatabaseIO::read_meta_data__() {}

  void DatabaseIO::write_meta_data()
  {
    this->catCGNSMesh->CreateBase(0, "Base");
    const auto &structured_blocks = this->get_region()->get_structured_blocks();
    int         base              = 0;
    int         zone              = 0;
    for (const auto &sb : structured_blocks) {
      sb->property_update("zone", zone);
      sb->property_update("base", base);
      zone++;
    }
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {

    Ioss::Field::RoleType role       = field.get_role();
    size_t                base       = sb->get_property("base").get_int();
    size_t                zone       = sb->get_property("zone").get_int();
    size_t                node_count = sb->get_property("node_count").get_int();
    size_t                num_to_get = field.verify(data_size);

    auto var_type               = field.transformed_storage();
    int  comp_count             = var_type->component_count();
    char field_suffix_separator = get_field_separator();

    bool is_cell_field = true;
    if (node_count == num_to_get) {
      is_cell_field = false;
    }

    double *rdata = num_to_get > 0 ? static_cast<double *>(data) : nullptr;

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x" ||
          field.get_name() == "mesh_model_coordinates_y" ||
          field.get_name() == "mesh_model_coordinates_z") {

        this->catCGNSMesh->AddStructuredZoneData(
            base, zone, sb->name(), field.get_name(), sb->get_property("ni").get_int(),
            sb->get_property("nj").get_int(), sb->get_property("nk").get_int(), comp_count,
            is_cell_field, field_suffix_separator, rdata, num_to_get);
      }
      else if (field.get_name() == "mesh_model_coordinates") {
        int phys_dimension = get_region()->get_property("spatial_dimension").get_int();

        std::vector<double> coord(num_to_get);

        // ========================================================================
        // Repetitive code for each coordinate direction; use a lambda to consolidate...
        auto coord_lambda = [=, &coord](const char *ordinate, int ordinal) {
          // Data required by upper classes store x0, y0, z0, ... xn,
          // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
          // ..., yn, z0, ..., zn so we have to allocate some scratch
          // memory to read in the data and then map into supplied
          // 'data'
          // Map to global coordinate position...
          for (size_t i = 0; i < num_to_get; i++) {
            coord[i] = rdata[phys_dimension * i + ordinal];
          }

          this->catCGNSMesh->AddStructuredZoneData(
              base, zone, sb->name(), ordinate, sb->get_property("ni").get_int(),
              sb->get_property("nj").get_int(), sb->get_property("nk").get_int(), comp_count,
              is_cell_field, field_suffix_separator, coord.data(), num_to_get);
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
    }
    else if (role == Ioss::Field::TRANSIENT) {
      this->catCGNSMesh->AddStructuredZoneData(
          base, zone, sb->name(), field.get_name(), sb->get_property("ni").get_int(),
          sb->get_property("nj").get_int(), sb->get_property("nk").get_int(), comp_count,
          is_cell_field, field_suffix_separator, rdata, num_to_get);
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {

    size_t num_to_get = field.verify(data_size);
    return num_to_get;
  }

} // namespace Iovs_cgns
