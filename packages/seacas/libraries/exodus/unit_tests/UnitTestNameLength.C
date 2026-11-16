// Copyright(C) 2026 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// Tests reading names whose length fills the `len_name` dimension of the
// database, i.e. names that are stored without a terminating NUL.  Older
// writers produce databases with `len_name = 32` holding 32-character
// names; these must be returned in full, without overrunning the caller's
// buffer.

#include <gtest/gtest.h>

#include <exodusII.h>
#include <netcdf.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {
  const std::string name32 = "LN012345678901234567890123456789";
  const std::string blk32  = "BLOCK_NAME_THAT_IS_32_CHARACTERS";
  const std::string name40 = "LONG_NAME_THAT_IS_EXACTLY_40_CHARACTERS_";

  constexpr char sentinel = '\x7f';

  // Create a minimal database with one element block and one global
  // variable whose names are written directly with netCDF so that the
  // `len_name` dimension can be smaller than exodus would choose.  If
  // `len_name` is 0, exodus picks its default dimension size.
  void create_file(const std::string &filename, int len_name, const std::string &var_name,
                   const std::string &blk_name)
  {
    int cpu_ws = 8;
    int io_ws  = 8;
    int exoid  = ex_create(filename.c_str(), EX_CLOBBER, &cpu_ws, &io_ws);
    ASSERT_GE(exoid, 0);

    if (len_name > 0) {
      int dimid = 0;
      nc_redef(exoid);
      ASSERT_EQ(nc_def_dim(exoid, "len_name", len_name, &dimid), NC_NOERR);
      nc_enddef(exoid);
    }

    ASSERT_EQ(ex_put_init(exoid, "name length test", 2, 4, 1, 1, 0, 0), EX_NOERR);
    std::vector<double> x{0.0, 1.0, 1.0, 0.0};
    std::vector<double> y{0.0, 0.0, 1.0, 1.0};
    ASSERT_EQ(ex_put_coord(exoid, x.data(), y.data(), nullptr), EX_NOERR);
    ASSERT_EQ(ex_put_block(exoid, EX_ELEM_BLOCK, 1, "QUAD4", 1, 4, 0, 0, 0), EX_NOERR);
    std::vector<int> conn{1, 2, 3, 4};
    ASSERT_EQ(ex_put_conn(exoid, EX_ELEM_BLOCK, 1, conn.data(), nullptr, nullptr), EX_NOERR);
    ASSERT_EQ(ex_put_variable_param(exoid, EX_GLOBAL, 1), EX_NOERR);

    // Write the names bypassing exodus so they can fill `len_name` with no NUL.
    int rootid = exoid; // No groups, so the file id is the root id.
    int varid  = 0;
    ASSERT_EQ(nc_inq_varid(rootid, "name_glo_var", &varid), NC_NOERR);
    size_t start[2] = {0, 0};
    size_t count[2] = {1, var_name.size()};
    ASSERT_EQ(nc_put_vara_text(rootid, varid, start, count, var_name.c_str()), NC_NOERR);

    ASSERT_EQ(nc_inq_varid(rootid, "eb_names", &varid), NC_NOERR);
    count[1] = blk_name.size();
    ASSERT_EQ(nc_put_vara_text(rootid, varid, start, count, blk_name.c_str()), NC_NOERR);

    int max_len = static_cast<int>(std::max(var_name.size(), blk_name.size()));
    ASSERT_EQ(nc_redef(rootid), NC_NOERR);
    ASSERT_EQ(nc_put_att_int(rootid, NC_GLOBAL, "maximum_name_length", NC_INT, 1, &max_len),
              NC_NOERR);
    ASSERT_EQ(nc_enddef(rootid), NC_NOERR);

    ASSERT_EQ(ex_close(exoid), EX_NOERR);
  }

  int open_file(const std::string &filename)
  {
    int   cpu_ws  = 8;
    int   io_ws   = 0;
    float version = 0.0;
    return ex_open(filename.c_str(), EX_READ, &cpu_ws, &io_ws, &version);
  }

  // Buffer of `max_name + 1` characters followed by sentinel bytes to detect overruns.
  struct NameBuffer
  {
    explicit NameBuffer(size_t max_name) : size(max_name + 1), data(size + 8, sentinel)
    {
      std::fill(data.begin(), data.begin() + size, '\0');
    }
    char       *get() { return data.data(); }
    std::string str() const { return std::string(data.data()); }
    bool        intact() const
    {
      return std::all_of(data.begin() + size, data.end(), [](char c) { return c == sentinel; });
    }
    size_t            size;
    std::vector<char> data;
  };

  class NameLength32 : public ::testing::Test
  {
  protected:
    static void SetUpTestSuite() { create_file(filename, 32, name32, blk32); }
    static void TearDownTestSuite() { std::remove(filename.c_str()); }

    void SetUp() override
    {
      exoid = open_file(filename);
      ASSERT_GE(exoid, 0);
      ASSERT_EQ(ex_set_max_name_length(exoid, 32), EX_NOERR);
    }
    void TearDown() override { ex_close(exoid); }

    static inline const std::string filename = "utst-len-name-32.exo";
    int                             exoid    = -1;
  };

  TEST_F(NameLength32, DimensionIs32)
  {
    EXPECT_EQ(ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH), 31);
  }

  TEST_F(NameLength32, GetVariableNamesReturnsFull32Characters)
  {
    NameBuffer buf(32);
    char      *names[] = {buf.get()};
    ASSERT_EQ(ex_get_variable_names(exoid, EX_GLOBAL, 1, names), EX_NOERR);
    EXPECT_EQ(buf.str(), name32);
    EXPECT_TRUE(buf.intact());
  }

  TEST_F(NameLength32, GetVariableNameReturnsFull32Characters)
  {
    NameBuffer buf(32);
    ASSERT_EQ(ex_get_variable_name(exoid, EX_GLOBAL, 1, buf.get()), EX_NOERR);
    EXPECT_EQ(buf.str(), name32);
    EXPECT_TRUE(buf.intact());
  }

  TEST_F(NameLength32, GetNameReturnsFull32Characters)
  {
    NameBuffer buf(32);
    ASSERT_EQ(ex_get_name(exoid, EX_ELEM_BLOCK, 1, buf.get()), EX_NOERR);
    EXPECT_EQ(buf.str(), blk32);
    EXPECT_TRUE(buf.intact());
  }

  TEST_F(NameLength32, GetNamesReturnsFull32Characters)
  {
    NameBuffer buf(32);
    char      *names[] = {buf.get()};
    ASSERT_EQ(ex_get_names(exoid, EX_ELEM_BLOCK, names), EX_NOERR);
    EXPECT_EQ(buf.str(), blk32);
    EXPECT_TRUE(buf.intact());
  }

  // A database with a large `len_name` holding names longer than the client
  // read length must be truncated to the read length without overrunning the
  // client buffer.
  class NameLongerThanReadLength : public ::testing::Test
  {
  protected:
    static void SetUpTestSuite() { create_file(filename, 0, name40, name40); }
    static void TearDownTestSuite() { std::remove(filename.c_str()); }

    void SetUp() override
    {
      exoid = open_file(filename);
      ASSERT_GE(exoid, 0);
    }
    void TearDown() override { ex_close(exoid); }

    static inline const std::string filename = "utst-len-name-long.exo";
    int                             exoid    = -1;
  };

  TEST_F(NameLongerThanReadLength, GetVariableNamesTruncatesToReadLength)
  {
    NameBuffer buf(32);
    char      *names[] = {buf.get()};
    ASSERT_GE(ex_get_variable_names(exoid, EX_GLOBAL, 1, names), EX_NOERR);
    EXPECT_EQ(buf.str(), name40.substr(0, 32));
    EXPECT_TRUE(buf.intact());
  }

  TEST_F(NameLongerThanReadLength, GetVariableNameTruncatesToReadLength)
  {
    NameBuffer buf(32);
    ASSERT_GE(ex_get_variable_name(exoid, EX_GLOBAL, 1, buf.get()), EX_NOERR);
    EXPECT_EQ(buf.str(), name40.substr(0, 32));
    EXPECT_TRUE(buf.intact());
  }

  TEST_F(NameLongerThanReadLength, GetNameTruncatesToReadLength)
  {
    NameBuffer buf(32);
    ASSERT_GE(ex_get_name(exoid, EX_ELEM_BLOCK, 1, buf.get()), EX_NOERR);
    EXPECT_EQ(buf.str(), name40.substr(0, 32));
    EXPECT_TRUE(buf.intact());
  }

  TEST_F(NameLongerThanReadLength, FullNameReturnedWhenReadLengthAllows)
  {
    ASSERT_EQ(ex_set_max_name_length(exoid, 40), EX_NOERR);
    NameBuffer buf(40);
    ASSERT_EQ(ex_get_variable_name(exoid, EX_GLOBAL, 1, buf.get()), EX_NOERR);
    EXPECT_EQ(buf.str(), name40);
    EXPECT_TRUE(buf.intact());
  }
} // namespace

int main(int argc, char **argv)
{
  ex_opts(EX_VERBOSE);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
