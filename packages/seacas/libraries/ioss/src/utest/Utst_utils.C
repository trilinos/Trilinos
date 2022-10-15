// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_NO_SHORT_MACRO_NAMES
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include <doctest.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_ConcreteVariableType.h>
#include <Ioss_Utils.h>
#include <exception>
#include <numeric>
#include <vector>

DOCTEST_TEST_CASE("number_width")
{
  DOCTEST_SUBCASE("single digit")
  {
    for (int i = 0; i < 10; i++) {
      DOCTEST_REQUIRE(1 == Ioss::Utils::number_width(i));
      DOCTEST_REQUIRE(1 == Ioss::Utils::number_width(i, true));
    }
  }

  DOCTEST_SUBCASE("double digit")
  {
    for (int i = 11; i < 100; i += 11) {
      DOCTEST_REQUIRE(2 == Ioss::Utils::number_width(i));
      DOCTEST_REQUIRE(2 == Ioss::Utils::number_width(i, true));
    }
  }

  DOCTEST_SUBCASE("triple digit")
  {
    for (int i = 111; i < 1'000; i += 111) {
      DOCTEST_REQUIRE(3 == Ioss::Utils::number_width(i));
      DOCTEST_REQUIRE(3 == Ioss::Utils::number_width(i, true));
    }
  }

  DOCTEST_SUBCASE("quad digit")
  {
    for (int i = 1111; i < 10'000; i += 1111) {
      DOCTEST_REQUIRE(4 == Ioss::Utils::number_width(i));
      DOCTEST_REQUIRE(5 == Ioss::Utils::number_width(i, true));
    }
  }

  DOCTEST_SUBCASE("larger")
  {
    DOCTEST_REQUIRE(6 == Ioss::Utils::number_width(999'999));
    DOCTEST_REQUIRE(7 == Ioss::Utils::number_width(999'999, true));
    DOCTEST_REQUIRE(7 == Ioss::Utils::number_width(1'000'000));
    DOCTEST_REQUIRE(9 == Ioss::Utils::number_width(1'000'000, true));
    DOCTEST_REQUIRE(10 == Ioss::Utils::number_width(1'111'111'111));
    DOCTEST_REQUIRE(13 == Ioss::Utils::number_width(1'111'111'111, true));
#if !defined(__IOSS_WINDOWS__)
    DOCTEST_REQUIRE(15 == Ioss::Utils::number_width(111'111'111'111, true));
#endif
  }
}

DOCTEST_TEST_CASE("detect large int64")
{
  Ioss::StorageInitializer init_fields;
  static int64_t           max_double = 2LL << 53;

  Ioss::Field field{"test", Ioss::Field::INT64, "vector_2d", Ioss::Field::RoleType::TRANSIENT, 3};
  std::vector<int64_t> data{0, max_double, max_double + 1, max_double - 1};
  DOCTEST_REQUIRE(true == Ioss::Utils::check_int_to_real_overflow(field, data.data(), data.size()));
}

#if !defined __NVCC__
DOCTEST_TEST_CASE("str_equal")
{
  DOCTEST_REQUIRE(Ioss::Utils::str_equal("", ""));
  DOCTEST_REQUIRE(!Ioss::Utils::str_equal("", "a"));
  DOCTEST_REQUIRE(!Ioss::Utils::str_equal("a", ""));
  DOCTEST_REQUIRE(Ioss::Utils::str_equal("a", "a"));
  DOCTEST_REQUIRE(Ioss::Utils::str_equal("A", "a"));
  DOCTEST_REQUIRE(Ioss::Utils::str_equal("a", "A"));

  DOCTEST_REQUIRE(
      Ioss::Utils::str_equal("longer_than_single_character", "longer_than_single_character"));
  DOCTEST_REQUIRE(
      Ioss::Utils::str_equal("longer_than_single_character", "LONGER_THAN_SINGLE_CHARACTER"));
  DOCTEST_REQUIRE(
      Ioss::Utils::str_equal("LONGER_THAN_SINGLE_CHARACTER", "longer_than_single_character"));
  DOCTEST_REQUIRE(
      Ioss::Utils::str_equal("LONGER_THAN_SINGLE_CHARACTER", "LONGER_THAN_SINGLE_CHARACTER"));
  DOCTEST_REQUIRE(
      Ioss::Utils::str_equal("LoNgEr_ThAn_SiNgLe_ChArAcTeR", "lOnGeR_tHaN_sInGlE_cHaRaCtEr"));

  DOCTEST_REQUIRE(!Ioss::Utils::str_equal("Almost_The_Same", "almost_the_sam"));
}

DOCTEST_TEST_CASE("substr_equal")
{
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("", ""));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("", "a"));
  DOCTEST_REQUIRE(!Ioss::Utils::substr_equal("a", ""));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("a", "a"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("A", "a"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("a", "A"));

  DOCTEST_REQUIRE(!Ioss::Utils::substr_equal("prefix", "pref"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("prefix", "PREFIX"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("pre", "PREFIX"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("pre", "prefix"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("PRe", "prefix"));
  DOCTEST_REQUIRE(Ioss::Utils::substr_equal("PRe", "PREFIX"));
}

DOCTEST_TEST_CASE("format_id_list")
{
  DOCTEST_SUBCASE("ids_empty")
  {
    std::string ret = Ioss::Utils::format_id_list({});
    DOCTEST_REQUIRE(ret == std::string(""));
  }

  DOCTEST_SUBCASE("ids_distinct")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 3, 5, 7, 9, 11});
    DOCTEST_REQUIRE(ret == std::string("1, 3, 5, 7, 9, 11"));
  }

  DOCTEST_SUBCASE("ids two element ranges")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 4, 5, 7, 8, 10, 11});
    DOCTEST_REQUIRE(ret == std::string("1, 2, 4, 5, 7, 8, 10, 11"));
  }

  DOCTEST_SUBCASE("ids one-range")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, "--");
    DOCTEST_REQUIRE(ret == std::string("1--11"));
  }

  DOCTEST_SUBCASE("ids one value")
  {
    std::string ret = Ioss::Utils::format_id_list({11}, "--");
    DOCTEST_REQUIRE(ret == std::string("11"));
  }

  DOCTEST_SUBCASE("ids two consecutive values")
  {
    std::string ret = Ioss::Utils::format_id_list({10, 11}, "--");
    DOCTEST_REQUIRE(ret == std::string("10, 11"));
  }

  DOCTEST_SUBCASE("ids two separated values")
  {
    std::string ret = Ioss::Utils::format_id_list({2, 11}, "--");
    DOCTEST_REQUIRE(ret == std::string("2, 11"));
  }

  DOCTEST_SUBCASE("ids two separated values pipe sep")
  {
    std::string ret = Ioss::Utils::format_id_list({2, 11}, "--");
    DOCTEST_REQUIRE(ret == std::string("2, 11"));
  }

  DOCTEST_SUBCASE("ids two ranges")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 3, 4, 6, 7, 8}, "--");
    DOCTEST_REQUIRE(ret == std::string("1--4, 6--8"));
  }

  DOCTEST_SUBCASE("ids two ranges pipe separated")
  {
    std::string ret = Ioss::Utils::format_id_list({1, 2, 3, 4, 6, 7, 8}, "--", "|");
    DOCTEST_REQUIRE(ret == std::string("1--4|6--8"));
  }

  DOCTEST_SUBCASE("ids large range")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    std::string ret = Ioss::Utils::format_id_list(range, "--");
    DOCTEST_REQUIRE(ret == std::string("42--100041"));
  }

  DOCTEST_SUBCASE("ids large range with singles")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range, "--");
    DOCTEST_REQUIRE(ret == std::string("1, 43--100040, 100042"));
  }

  DOCTEST_SUBCASE("ids large range with singles with to")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range);
    DOCTEST_REQUIRE(ret == std::string("1, 43 to 100040, 100042"));
  }

  DOCTEST_SUBCASE("ids large range with singles with colon")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range, ":");
    DOCTEST_REQUIRE(ret == std::string("1, 43:100040, 100042"));
  }

  DOCTEST_SUBCASE("ids large range with singles with words")
  {
    std::vector<size_t> range(100'000);
    std::iota(range.begin(), range.end(), 42);
    range[0] = 1;
    range[range.size() - 1]++;
    std::string ret = Ioss::Utils::format_id_list(range, " to ", " and ");
    DOCTEST_REQUIRE(ret == std::string("1 and 43 to 100040 and 100042"));
  }

  DOCTEST_SUBCASE("detect unsorted two ids")
  {
    DOCTEST_CHECK_THROWS(Ioss::Utils::format_id_list({2, 1}));
  }

  DOCTEST_SUBCASE("detect unsorted ids")
  {
    DOCTEST_CHECK_THROWS(Ioss::Utils::format_id_list({1, 2, 3, 4, 5, 1}));
  }

  DOCTEST_SUBCASE("detect duplicate ids")
  {
    DOCTEST_CHECK_THROWS(Ioss::Utils::format_id_list({1, 2, 3, 3, 4, 5, 6}));
  }
}
#endif
