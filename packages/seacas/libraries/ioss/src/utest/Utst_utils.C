// Copyright(C) 1999-2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_ConcreteVariableType.h"
#include "Ioss_Field.h"
#include "Ioss_Utils.h"

TEST_CASE("number_width")
{
  SECTION("single digit")
  {
    for (int i = 0; i < 10; i++) {
      CHECK(1 == Ioss::Utils::number_width(i));
      CHECK(1 == Ioss::Utils::number_width(i, true));
    }
  }

  SECTION("double digit")
  {
    for (int i = 11; i < 100; i += 11) {
      CHECK(2 == Ioss::Utils::number_width(i));
      CHECK(2 == Ioss::Utils::number_width(i, true));
    }
  }

  SECTION("triple digit")
  {
    for (int i = 111; i < 1'000; i += 111) {
      CHECK(3 == Ioss::Utils::number_width(i));
      CHECK(3 == Ioss::Utils::number_width(i, true));
    }
  }

  SECTION("quad digit")
  {
    for (int i = 1111; i < 10'000; i += 1111) {
      CHECK(4 == Ioss::Utils::number_width(i));
      CHECK(5 == Ioss::Utils::number_width(i, true));
    }
  }

  SECTION("larger")
  {
    CHECK(6 == Ioss::Utils::number_width(999'999));
    CHECK(7 == Ioss::Utils::number_width(999'999, true));
    CHECK(7 == Ioss::Utils::number_width(1'000'000));
    CHECK(9 == Ioss::Utils::number_width(1'000'000, true));
    CHECK(10 == Ioss::Utils::number_width(1'111'111'111));
    CHECK(13 == Ioss::Utils::number_width(1'111'111'111, true));
#if !defined(__IOSS_WINDOWS__)
    CHECK(15 == Ioss::Utils::number_width(111'111'111'111, true));
#endif
  }
}

TEST_CASE("detect large int64")
{
  Ioss::StorageInitializer init_fields;
  static int64_t           max_double = 2LL << 53;

  Ioss::Field field{"test", Ioss::Field::INT64, "vector_2d", Ioss::Field::RoleType::TRANSIENT, 3};
  std::vector<int64_t> data{0, max_double, max_double + 1, max_double - 1};
  CHECK(true == Ioss::Utils::check_int_to_real_overflow(field, Data(data), data.size()));
}

#if !defined __NVCC__
TEST_CASE("str_equal")
{
  CHECK(Ioss::Utils::str_equal("", ""));
  CHECK(!Ioss::Utils::str_equal("", "a"));
  CHECK(!Ioss::Utils::str_equal("a", ""));
  CHECK(Ioss::Utils::str_equal("a", "a"));
  CHECK(Ioss::Utils::str_equal("A", "a"));
  CHECK(Ioss::Utils::str_equal("a", "A"));

  CHECK(Ioss::Utils::str_equal("longer_than_single_character", "longer_than_single_character"));
  CHECK(Ioss::Utils::str_equal("longer_than_single_character", "LONGER_THAN_SINGLE_CHARACTER"));
  CHECK(Ioss::Utils::str_equal("LONGER_THAN_SINGLE_CHARACTER", "longer_than_single_character"));
  CHECK(Ioss::Utils::str_equal("LONGER_THAN_SINGLE_CHARACTER", "LONGER_THAN_SINGLE_CHARACTER"));
  CHECK(Ioss::Utils::str_equal("LoNgEr_ThAn_SiNgLe_ChArAcTeR", "lOnGeR_tHaN_sInGlE_cHaRaCtEr"));

  CHECK(!Ioss::Utils::str_equal("Almost_The_Same", "almost_the_sam"));
}

TEST_CASE("substr_equal")
{
  CHECK(Ioss::Utils::substr_equal("", ""));
  CHECK(Ioss::Utils::substr_equal("", "a"));
  CHECK(!Ioss::Utils::substr_equal("a", ""));
  CHECK(Ioss::Utils::substr_equal("a", "a"));
  CHECK(Ioss::Utils::substr_equal("A", "a"));
  CHECK(Ioss::Utils::substr_equal("a", "A"));

  CHECK(!Ioss::Utils::substr_equal("prefix", "pref"));
  CHECK(Ioss::Utils::substr_equal("prefix", "PREFIX"));
  CHECK(Ioss::Utils::substr_equal("pre", "PREFIX"));
  CHECK(Ioss::Utils::substr_equal("pre", "prefix"));
  CHECK(Ioss::Utils::substr_equal("PRe", "prefix"));
  CHECK(Ioss::Utils::substr_equal("PRe", "PREFIX"));
}

TEST_CASE("format_id_list")
{
  SECTION("ids_empty")
  {
    std::string ret = Ioss::Utils::format_id_list({});
    CHECK(ret == std::string(""));
  }

  SECTION("ids_distinct")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 3, 5, 7, 9, 11});
    CHECK(ret == std::string("1, 3, 5, 7, 9, 11"));
  }

  SECTION("ids two element ranges")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 4, 5, 7, 8, 10, 11});
    CHECK(ret == std::string("1, 2, 4, 5, 7, 8, 10, 11"));
  }

  SECTION("ids one-range")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, "--");
    CHECK(ret == std::string("1--11"));
  }

  SECTION("ids one value")
  {
    std::string ret = Ioss::Utils::format_id_list({11}, "--");
    CHECK(ret == std::string("11"));
  }

  SECTION("ids two consecutive values")
  {
    std::string ret = Ioss::Utils::format_id_list({10, 11}, "--");
    CHECK(ret == std::string("10, 11"));
  }

  SECTION("ids two separated values")
  {
    std::string ret = Ioss::Utils::format_id_list({2, 11}, "--");
    CHECK(ret == std::string("2, 11"));
  }

  SECTION("ids two separated values pipe sep")
  {
    std::string ret = Ioss::Utils::format_id_list({2, 11}, "--");
    CHECK(ret == std::string("2, 11"));
  }

  SECTION("ids two ranges")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 3, 4, 6, 7, 8}, "--");
    CHECK(ret == std::string("1--4, 6--8"));
  }

  SECTION("ids two ranges pipe separated")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 3, 4, 6, 7, 8}, "--", "|");
    CHECK(ret == std::string("1--4|6--8"));
  }

  SECTION("ids large range")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    std::string ret = Ioss::Utils::format_id_list(range, "--");
    CHECK(ret == std::string("42--100041"));
  }

  SECTION("ids large range with singles")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range, "--");
    CHECK(ret == std::string("1, 43--100040, 100042"));
  }

  SECTION("ids large range with singles with to")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range);
    CHECK(ret == std::string("1, 43 to 100040, 100042"));
  }

  SECTION("ids large range with singles with colon")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range, ":");
    CHECK(ret == std::string("1, 43:100040, 100042"));
  }

  SECTION("ids large range with singles with words")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range, " to ", " and ");
    CHECK(ret == std::string("1 and 43 to 100040 and 100042"));
  }

#ifndef SEACAS_HAVE_MPI
  SECTION("detect unsorted two ids") { CHECK_THROWS(Ioss::Utils::format_id_list({2, 1})); }

  SECTION("detect unsorted ids") { CHECK_THROWS(Ioss::Utils::format_id_list({1, 2, 3, 4, 5, 1})); }

  SECTION("detect duplicate ids")
  {
    CHECK_THROWS(Ioss::Utils::format_id_list({1, 2, 3, 3, 4, 5, 6}));
  }
#endif
}
#endif
