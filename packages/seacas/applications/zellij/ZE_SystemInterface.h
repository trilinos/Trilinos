// Copyright(C) 2021, 2022, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "Ioss_GetLongOpt.h" // for GetLongOption
#include <array>
#include <string> // for string

//! \file

using vector3d = std::array<double, 3>;

enum class Minimize { NONE = 0, UNIT = 1, OUTPUT = 2, ALL = 3 };

class SystemInterface
{
public:
  explicit SystemInterface(int my_rank = 0);

  bool parse_options(int argc, char **argv);

  int debug() const { return debugLevel_; }

  std::string lattice() const { return lattice_; }
  std::string decomp_method() const { return decompMethod_; }
  std::string sideset_surfaces() const { return sidesetSurfaces_; }
  std::string sideset_names() const { return sidesetNames_; }
  bool        ignore_internal_sidesets() const { return ignoreInternalSidesets_; }

  bool ints32bit() const { return ints32bit_; }
  bool use_netcdf4() const { return useNetcdf4_; }
  bool use_netcdf5() const { return useNetcdf5_; }

  int  compression_level() const { return compressionLevel_; }
  bool zlib() const { return zlib_; }
  bool szip() const { return szip_; }

  int  ranks() const { return ranks_; }
  int  start_rank() const { return startRank_; }
  int  rank_count() const { return rankCount_; }
  bool subcycle() const { return subcycle_; }

  bool     equivalence_nodes() const { return equivalenceNodes_; }
  Minimize minimize_open_files() const { return minimizeOpenFiles_; }

  double   scale_factor() const { return scaleFactor_; }
  vector3d offset() const { return offset_; }

  int skip() const { return skip_; }
  int repeat() const { return repeat_; }

  static void show_version();

  // Make this private eventually...
  std::string outputName_;

private:
  void enroll_options();

  Ioss::GetLongOption options_; //!< Options parsing

  std::string   lattice_{};
  std::string   decompMethod_{"HSFC"};
  std::string   sidesetSurfaces_{};
  std::string   sidesetNames_{};
  vector3d      offset_{0.0, 0.0, 0.0};
  double        scaleFactor_{1.0};
  int           myRank_{0};
  int           debugLevel_{0};
  int           compressionLevel_{0};
  int           ranks_{1};
  int           startRank_{0};
  int           rankCount_{0};
  int           skip_{0};
  int           repeat_{1};
  bool          ints32bit_{false};
  bool          useNetcdf4_{true};
  bool          useNetcdf5_{false};
  bool          zlib_{true};
  bool          szip_{false};
  bool          equivalenceNodes_{true};
  bool          subcycle_{false};
  bool          ignoreInternalSidesets_{false};
  enum Minimize minimizeOpenFiles_ { Minimize::NONE };
};
