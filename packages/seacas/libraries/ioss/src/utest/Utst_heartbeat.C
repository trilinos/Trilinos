#include "Ionit_Initializer.h"
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <fmt/core.h>
#include <math.h>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_Field.h"
#include "Ioss_IOFactory.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_State.h"

namespace {
  std::string tst_filename = "test.hb";

  Ioss::DatabaseIO *create_heartbeat(const std::string &filename)
  {
    Ioss::Init::Initializer init_db;

    Ioss::PropertyManager properties;
    properties.add(Ioss::Property("FULL_PRECISION", "yes"));
    properties.add(Ioss::Property("SHOW_LABELS", "yes"));
    properties.add(Ioss::Property("SHOW_LEGEND", "no"));
    properties.add(Ioss::Property("SHOW_TIME_STAMP", 1));
    properties.add(Ioss::Property("TIME_STAMP_FORMAT", "{%F %H:%M:%S}"));
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("heartbeat", filename, Ioss::WRITE_HEARTBEAT,
                                                    Ioss::ParallelUtils::comm_world(), properties);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }
    return dbo;
  }
  // BeginDocTest2
  TEST_CASE("Ioss::write_file")
  {
    auto *db = create_heartbeat(tst_filename);
    CHECK(db->ok());

    Ioss::Region region(db);

    region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
    region.field_add({"double", Ioss::Field::REAL, "scalar", Ioss::Field::TRANSIENT, 1});
    region.field_add({"integer", Ioss::Field::INTEGER, "scalar", Ioss::Field::TRANSIENT, 1});
    region.field_add({"vector_3d", Ioss::Field::REAL, "vector_3d", Ioss::Field::TRANSIENT, 1});
    region.field_add({"intvector", Ioss::Field::INTEGER, "vector_3d", Ioss::Field::TRANSIENT, 1});
    region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

    std::vector<double> field_data(3);
    std::vector<int>    int_data(3);
    region.begin_mode(Ioss::STATE_TRANSIENT);
    for (int i = 1; i < 10; i++) {
      region.add_state(i);
      region.begin_state(i);
      field_data[0] = (i % 2) ? i : -i;
      field_data[1] = i * i;
      field_data[2] = sqrt(i);
      int_data[0]   = (i % 2) ? -i : i;
      int_data[1]   = i / 2;
      int_data[2]   = i * i;
      region.put_field_data("double", field_data);
      region.put_field_data("integer", int_data);
      region.put_field_data("vector_3d", field_data);
      region.put_field_data("intvector", int_data);
      region.end_state(i);
    }
    region.end_mode(Ioss::STATE_TRANSIENT);
  }
} // namespace

int main(IOSS_MAYBE_UNUSED int argc, char **argv)
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif
  Catch::Session session; // There must be exactly one instance

  // Build a new parser on top of Catch2's
  using namespace Catch::Clara;
  auto cli = session.cli()                   // Get Catch2's command line parser
             | Opt(tst_filename, "filename") // bind variable to a new option, with a hint string
                   ["-f"]["--filename"]      // the option names it will respond to
             ("Filename to write heartbeat data to."); // description string for the help output

  // Now pass the new composite back to Catch2 so it uses that
  session.cli(cli);

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) // Indicates a command line error
    return returnCode;

  fmt::print("'{}'\n", tst_filename);
  return session.run();
}
