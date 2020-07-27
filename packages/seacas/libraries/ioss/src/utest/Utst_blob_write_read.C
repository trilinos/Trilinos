// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_Blob.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Property.h>
#include <Ioss_Region.h>
#include <Ioss_ScopeGuard.h>
#include <Ioss_Utils.h>

#include <Ionit_Initializer.h>

#include <algorithm>
#include <string>
#include <vector>

//--------------------------------------------------------------------
/*----------------------------------------------------------------------
 * IOSS Blob Example
 *
 */
void write_blob();
bool read_blob();

std::vector<double> generate_data(double time, size_t global_size, double offset, size_t local_size,
                                  size_t proc_offset)
{
  // Determine this ranks portion of the data
  std::vector<double> data;
  data.reserve(local_size);
  for (size_t i = 0; i < local_size; i++) {
    size_t ii = proc_offset + i;
    data.push_back(10.0 * ii + 100 * time + offset);
  }
  return data;
}

namespace {
  std::pair<size_t, size_t> get_blob_size(size_t global_size, size_t parallel_size, size_t my_rank)
  {
    size_t per_proc = global_size / parallel_size;
    size_t extra    = global_size % parallel_size;
    size_t count    = per_proc + (my_rank < extra ? 1 : 0);

    size_t offset = 0;
    if (my_rank < extra) {
      offset = (per_proc + 1) * my_rank;
    }
    else {
      offset = (per_proc + 1) * extra + per_proc * (my_rank - extra);
    }
    return std::make_pair(count, offset);
  }
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  Ioss::Init::Initializer io; // Initialize IOSS library.

  write_blob();
  bool diff = read_blob();
  if (diff) {
    std::cout << "Differences detected\n";
    return 1;
  }
  else {
    std::cout << "No Differences detected\n";
    return 0;
  }
}

void write_blob()
{
  std::cout << "***** Writing Blob Example File...\n";
  Ioss::PropertyManager properties;
  properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(
      "exodus", "ioss_blob_example.e", Ioss::WRITE_RESTART, (MPI_Comm)MPI_COMM_WORLD, properties);
  if (dbo == NULL || !dbo->ok(true)) {
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'region' owns 'dbo' pointer at this time
  Ioss::Region region(dbo, "example_region");
  region.begin_mode(Ioss::STATE_DEFINE_MODEL);

  const size_t b1_size = 100;
  const size_t b2_size = 200;
  const size_t b3_size = 57;

  // For example, we will just spread blob evenly across all ranks..
  size_t par_size = dbo->util().parallel_size();
  size_t my_rank  = dbo->util().parallel_rank();

  // Define a blob -- give name and size
  auto        count_offset = get_blob_size(b1_size, par_size, my_rank);
  Ioss::Blob *blob1        = new Ioss::Blob(dbo, "Tempus", count_offset.first);
  region.add(blob1);

  // NOTE: These properties are not needed for serial case, but don't cause problems
  blob1->property_add(Ioss::Property("processor_offset", (int64_t)count_offset.second));
  blob1->property_add(Ioss::Property("global_size", (int64_t)b1_size));

  count_offset      = get_blob_size(b2_size, par_size, my_rank);
  Ioss::Blob *blob2 = new Ioss::Blob(dbo, "Solver", count_offset.first);
  region.add(blob2);

  blob2->property_add(Ioss::Property("processor_offset", (int64_t)count_offset.second));
  blob2->property_add(Ioss::Property("global_size", (int64_t)b2_size));

  count_offset      = get_blob_size(b3_size, par_size, my_rank);
  Ioss::Blob *blob3 = new Ioss::Blob(dbo, "ABlob", count_offset.first);
  region.add(blob3);

  blob3->property_add(Ioss::Property("processor_offset", (int64_t)count_offset.second));
  blob3->property_add(Ioss::Property("global_size", (int64_t)b3_size));

  // These are "entity attributes" for blob1. Non-transient (constant) property
  // applied to the blob, not each entry in the blob.
  std::vector<double> offsets{1.0, 2.0, 0.0};
  std::vector<double> scales{10.5, 11.5, 17.5};
  blob1->property_add(Ioss::Property("Offset", offsets, Ioss::Property::ATTRIBUTE));
  blob1->property_add(Ioss::Property("Scale", scales, Ioss::Property::ATTRIBUTE));

  region.end_mode(Ioss::STATE_DEFINE_MODEL);

  region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

  // Aaadd some transient fields to the blobs.  There will be a value per entry in the blob.
  Ioss::Field x("X", Ioss::Field::BasicType::REAL, "scalar", Ioss::Field::RoleType::TRANSIENT);
  Ioss::Field dx("XDOT", Ioss::Field::BasicType::REAL, "scalar", Ioss::Field::RoleType::TRANSIENT);
  Ioss::Field ddx("XDDOT", Ioss::Field::BasicType::REAL, "scalar",
                  Ioss::Field::RoleType::TRANSIENT);

  blob1->field_add(x);
  blob1->field_add(dx);
  blob1->field_add(ddx);

  // Blobs can have different fields
  blob2->field_add(x);
  blob2->field_add(dx);

  blob3->field_add(x);

  // Reduction Fields -- Single value per blob per timestep
  Ioss::Field momentum("Momentum", Ioss::Field::BasicType::REAL, "vector_3d",
                       Ioss::Field::RoleType::REDUCTION);
  Ioss::Field ke("kinetic_energy", Ioss::Field::BasicType::REAL, "scalar",
                 Ioss::Field::RoleType::REDUCTION);

  blob1->field_add(ke);
  blob1->field_add(momentum);

  blob2->field_add(ke);
  blob2->field_add(momentum);

  blob3->field_add(ke);
  blob3->field_add(momentum);

  region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

  region.begin_mode(Ioss::STATE_TRANSIENT);
  const size_t num_ts = 10;
  for (size_t ts = 0; ts < num_ts; ts++) {
    double time = ts / 10.0;
    auto   step = region.add_state(time);
    region.begin_state(step);

    const auto &blobs = region.get_blobs();
    int         idx   = 0;
    for (const auto *blob : blobs) {
      // Dummy data for the fields.  All the same here...
      // Would be different in reality.
      const size_t size = blob->entity_count();

      // Get the global size and offset of this blob on this rank...
      size_t gl_size  = blob->get_property("global_size").get_int();
      size_t p_offset = blob->get_property("processor_offset").get_int();

      // Get the fields that are defined on this blob...
      Ioss::NameList fields;
      blob->field_describe(Ioss::Field::RoleType::TRANSIENT, &fields);
      std::sort(fields.begin(), fields.end()); // Just done for testing; not needed
      for (const auto &field : fields) {
        std::vector<double> data = generate_data(time, gl_size, idx++, size, p_offset);
        blob->put_field_data(field, data);
      }

      // Reduction fields...
      Ioss::NameList red_fields;
      blob->field_describe(Ioss::Field::RoleType::REDUCTION, &red_fields);
      for (const auto &field : red_fields) {
        std::vector<double> data = generate_data(time, 3, idx++, 3, 0);
        blob->put_field_data(field, data);
      }
    }
    region.end_state(step);
  }
  region.end_mode(Ioss::STATE_TRANSIENT);
  // File closed when `region` goes out of scope.
}

bool read_blob()
{
  std::cout << "\n***** Reading Blob Example File...\n";
  Ioss::PropertyManager properties;
  properties.add(Ioss::Property("DECOMPOSITION_METHOD", "linear"));
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(
      "exodus", "ioss_blob_example.e", Ioss::READ_RESTART, (MPI_Comm)MPI_COMM_WORLD, properties);
  if (dbi == NULL || !dbi->ok(true)) {
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'region' owns 'dbi' pointer at this time
  Ioss::Region region(dbi, "example_region");

  // Print a summary of the properties and fields on each blob...
  const auto &blobs = region.get_blobs();
  for (const auto *blob : blobs) {
    std::cout << "\nBlob " << blob->name() << " contains: " << blob->entity_count()
              << " item(s).\n";

    Ioss::Utils::info_property(blob, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ");
    Ioss::Utils::info_fields(blob, Ioss::Field::TRANSIENT, "\n\tTransient: ");
    Ioss::Utils::info_fields(blob, Ioss::Field::REDUCTION, "\n\tTransient (Reduction):  ", "\t");
  }

  size_t ts_count = 0;
  if (region.property_exists("state_count")) {
    ts_count = region.get_property("state_count").get_int();
  }
  std::cout << "\nFile contains " << ts_count << " timesteps.\n";

  std::vector<Ioss::NameList> all_fields;
  std::vector<Ioss::NameList> all_red_fields;

  for (const auto *blob : blobs) {
    // Get the names of the fields that are defined on this blob...
    Ioss::NameList fields;
    blob->field_describe(Ioss::Field::RoleType::TRANSIENT, &fields);
    std::sort(fields.begin(), fields.end()); // Just done for testing; not needed
    all_fields.push_back(fields);

    // Reduction fields...
    Ioss::NameList red_fields;
    blob->field_describe(Ioss::Field::RoleType::REDUCTION, &red_fields);
    all_red_fields.push_back(red_fields);
  }

  size_t par_size = dbi->util().parallel_size();

  bool diff = false;
  for (size_t step = 1; step <= ts_count; step++) {
    region.begin_state(step);
    double time = region.get_state_time(step); // Region steps are 1-based
    std::cout << "\n*** Step " << step << " is at time " << time << "\n";

    int                 idx = 0;
    std::vector<double> data;
    int                 offset = 0;
    for (const auto *blob : blobs) {
      // Get the global size and offset of this blob on this rank...
      // These are only needed for the comparison, not needed to just read data.
      size_t size     = blob->entity_count();
      size_t gl_size  = size;
      size_t p_offset = 0;
      if (par_size > 1) {
        gl_size  = blob->get_property("global_size").get_int();
        p_offset = blob->get_property("processor_offset").get_int();
      }

      const auto &fields = all_fields[idx];
      for (const auto &field : fields) {
        blob->get_field_data(field, data);

        // Compare with what was written to make sure read/write is ok.
        std::vector<double> gold = generate_data(time, gl_size, offset++, size, p_offset);
        if (data != gold) {
          std::cout << "Difference for field " << field << " on blob " << blob->name()
                    << " at step " << step << "\n";
          diff = true;
        }
      }

      const auto &red_fields = all_red_fields[idx++];
      for (const auto &field : red_fields) {
        blob->get_field_data(field, data);
        offset++; // Just for the testing part.
        std::cout << "\t\tReduction Field " << field << ", Value = ";

        size_t comp_count = blob->get_field(field).raw_storage()->component_count();
        for (size_t i = 0; i < comp_count; i++) {
          std::cout << data[i] << " ";
        }
        std::cout << "\n";
      }
    }
  }
  // File closed when `region` goes out of scope.
  return diff;
}
