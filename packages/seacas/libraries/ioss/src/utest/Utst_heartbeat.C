#define DOCTEST_CONFIG_IMPLEMENT
#define DOCTEST_CONFIG_NO_SHORT_MACRO_NAMES
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include <doctest.h>
#include <string>

#include <Ionit_Initializer.h>
#include <Ioss_ScopeGuard.h>
#include <Ioss_SubSystem.h>

#include <fmt/format.h>

namespace {
  std::string tst_filename = "test.hb";

  Ioss::DatabaseIO *create_heartbeat(const std::string &filename)
  {
    Ioss::Init::Initializer init_db;

    Ioss::PropertyManager properties;
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
  DOCTEST_TEST_CASE("Ioss::write_file")
  {
    auto *db = create_heartbeat(tst_filename);
    DOCTEST_CHECK(db->ok());

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

  doctest::Context context;

  while (*++argv) {
    if (std::string(*argv) == "--filename") {
      tst_filename = *++argv;
      break;
    }
    fmt::print("'{}'\n", tst_filename);
  }
  return context.run();
}
