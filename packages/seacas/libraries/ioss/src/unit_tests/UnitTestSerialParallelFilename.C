// Copyright(C) 2026 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// Tests database filename decoding when the client specifies the
// "processor_count" and "my_processor" properties ("serial-parallel" mode,
// as used by the ParaView IOSS reader).  Most codes omit the ".N.r" suffix
// for a single-rank decomposition, but some (e.g. Sculpt) write
// "<name>.1.0"; both must be readable.

#include "gtest/gtest.h"

#include "Ionit_Initializer.h"
#include "Ioss_CopyDatabase.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_FileInfo.h"
#include "Ioss_IOFactory.h"
#include "Ioss_MeshCopyOptions.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {
  constexpr int64_t elem_count = 2; // "generated:1x1x2"

  Ioss_MPI_Comm self_comm()
  {
#ifdef SEACAS_HAVE_MPI
    return MPI_COMM_SELF;
#else
    return Ioss::ParallelUtils::comm_world();
#endif
  }

  Ioss::PropertyManager serial_parallel_props(int proc_count, int my_proc)
  {
    Ioss::PropertyManager props;
    props.add(Ioss::Property("processor_count", proc_count));
    props.add(Ioss::Property("my_processor", my_proc));
    return props;
  }

  // Write a small exodus mesh to `filename` (no filename decoding applied).
  void write_mesh(const std::string &filename)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO       *db_i = Ioss::IOFactory::create("generated", "1x1x2", Ioss::READ_MODEL,
                                                           self_comm(), Ioss::PropertyManager());
    ASSERT_NE(db_i, nullptr);
    Ioss::Region region_i(db_i, "input");

    Ioss::DatabaseIO *db_o = Ioss::IOFactory::create("exodus", filename, Ioss::WRITE_RESTART,
                                                     self_comm(), Ioss::PropertyManager());
    ASSERT_NE(db_o, nullptr);
    ASSERT_TRUE(db_o->ok(true));
    Ioss::Region region_o(db_o, "output");

    Ioss::MeshCopyOptions options{};
    options.data_storage_type = 1;
    Ioss::copy_database(region_i, region_o, options);
  }

  // Open `base` for reading with `props` and return the element count, or -1 if it can't be
  // opened.  `opened` receives the name of the file actually opened.
  int64_t read_elem_count(const std::string &base, const Ioss::PropertyManager &props,
                          std::string *opened = nullptr)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO       *db =
        Ioss::IOFactory::create("exodus", base, Ioss::READ_RESTART, self_comm(), props);
    if (db == nullptr || !db->ok(false)) {
      delete db;
      return -1;
    }
    if (opened != nullptr) {
      *opened = db->decoded_filename();
    }
    Ioss::Region region(db, "region");
    return region.get_property("element_count").get_int();
  }

  class SerialParallelFilename : public ::testing::Test
  {
  protected:
    void SetUp() override
    {
      if (Ioss::ParallelUtils(Ioss::ParallelUtils::comm_world()).parallel_size() != 1) {
        GTEST_SKIP() << "Serial-parallel filename tests run on a single rank";
      }
      const auto *info = ::testing::UnitTest::GetInstance()->current_test_info();
      base             = std::string("utst-spfn-") + info->name() + ".e";
    }

    void TearDown() override
    {
      for (const auto &file : created) {
        std::remove(file.c_str());
      }
    }

    void make(const std::string &filename)
    {
      write_mesh(filename);
      ASSERT_TRUE(Ioss::FileInfo(filename).exists());
      created.push_back(filename);
    }

    std::string              base;
    std::vector<std::string> created;
  };

  // A: only "<name>.1.0" exists -- the file written by Sculpt that ParaView fails to open.
  TEST_F(SerialParallelFilename, OpensSingleRankSuffixedFile)
  {
    make(base + ".1.0");
    std::string opened;
    EXPECT_EQ(read_elem_count(base, serial_parallel_props(1, 0), &opened), elem_count);
    EXPECT_EQ(opened, base + ".1.0");
  }

  // B: only the unsuffixed file exists; same properties.
  TEST_F(SerialParallelFilename, OpensSingleRankUnsuffixedFile)
  {
    make(base);
    std::string opened;
    EXPECT_EQ(read_elem_count(base, serial_parallel_props(1, 0), &opened), elem_count);
    EXPECT_EQ(opened, base);
  }

  // C: both exist; the unsuffixed file is preferred (existing behavior).
  TEST_F(SerialParallelFilename, PrefersUnsuffixedFileWhenBothExist)
  {
    make(base);
    make(base + ".1.0");
    std::string opened;
    EXPECT_EQ(read_elem_count(base, serial_parallel_props(1, 0), &opened), elem_count);
    EXPECT_EQ(opened, base);
  }

  // D: multi-rank decomposition read one rank at a time.
  TEST_F(SerialParallelFilename, OpensMultiRankFile)
  {
    make(base + ".4.2");
    std::string opened;
    EXPECT_EQ(read_elem_count(base, serial_parallel_props(4, 2), &opened), elem_count);
    EXPECT_EQ(opened, base + ".4.2");
  }

  // E: no processor properties -- plain serial read.
  TEST_F(SerialParallelFilename, OpensPlainFileWithoutProperties)
  {
    make(base);
    std::string opened;
    EXPECT_EQ(read_elem_count(base, Ioss::PropertyManager(), &opened), elem_count);
    EXPECT_EQ(opened, base);
  }

  // F: output with a processor count of 1 is still written without a suffix.
  TEST_F(SerialParallelFilename, OutputSingleRankIsUnsuffixed)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodus", base, Ioss::WRITE_RESTART, self_comm(),
                                                   serial_parallel_props(1, 0));
    ASSERT_NE(db, nullptr);
    created.push_back(base);
    created.push_back(base + ".1.0");
    EXPECT_EQ(db->decoded_filename(), base);
    delete db;
  }
} // namespace
