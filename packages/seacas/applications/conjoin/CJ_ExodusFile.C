// Copyright(C) 1999-2021, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CJ_CodeTypes.h" // for StringIdVector, etc
#include "CJ_ExodusFile.h"
#include "CJ_SystemInterface.h"
#include "open_file_limit.h"
#include "smart_assert.h"

#include <climits>
#include <cstdlib>
#include <fmt/ostream.h>
#include <iomanip>
#include <iostream>
#include <string>

#include <exodusII.h>

std::vector<int>         Excn::ExodusFile::fileids_;
std::vector<std::string> Excn::ExodusFile::filenames_;
int                      Excn::ExodusFile::outputId_    = -1;
int                      Excn::ExodusFile::ioWordSize_  = 0;
int                      Excn::ExodusFile::cpuWordSize_ = 0;
std::string              Excn::ExodusFile::outputFilename_;
bool                     Excn::ExodusFile::keepOpen_          = false;
bool                     Excn::ExodusFile::usingChangeSets_   = false;
int                      Excn::ExodusFile::maximumNameLength_ = 32;
int                      Excn::ExodusFile::exodusMode_        = 0;

Excn::ExodusFile::ExodusFile(size_t which) : myLocation_(which)
{
  SMART_ASSERT(which < filenames_.size())(which)(filenames_.size());
  SMART_ASSERT(fileids_.size() == filenames_.size());
  if (!keepOpen_ && which != 0) {
    float version       = 0.0;
    int   cpu_word_size = cpuWordSize_;
    int   io_wrd_size   = ioWordSize_;
    SMART_ASSERT(fileids_[which] == -1)(which)(fileids_[which]);
    fileids_[which] = ex_open(filenames_[which].c_str(), EX_READ | exodusMode_, &cpu_word_size,
                              &io_wrd_size, &version);
    if (fileids_[which] < 0) {
      fmt::print(stderr, "ERROR: Cannot open file '{}' - exiting\n", filenames_[which]);
      exit(1);
    }
    ex_set_max_name_length(fileids_[which], maximumNameLength_);

    SMART_ASSERT(io_wrd_size == ioWordSize_);
    SMART_ASSERT(cpu_word_size == cpuWordSize_);
  }
}

bool Excn::ExodusFile::ints_64_bit() { return exodusMode_ & EX_ALL_INT64_API; }

int Excn::ExodusFile::output()
{
  SMART_ASSERT(outputId_ >= 0);
  return outputId_;
}

Excn::ExodusFile::operator int() const
{
  SMART_ASSERT(fileids_[myLocation_] >= 0);
  return fileids_[myLocation_];
}

Excn::ExodusFile::~ExodusFile()
{
  try {
    if (!keepOpen_ && myLocation_ != 0) {
      SMART_ASSERT(fileids_[myLocation_] > 0)(myLocation_)(fileids_[myLocation_]);
      ex_close(fileids_[myLocation_]);
      fileids_[myLocation_] = -1;
    }
  }
  catch (...) {
  }
}

void Excn::ExodusFile::close_all()
{
  if (usingChangeSets_) {
    ex_close(fileids_[0]);
  }
  else {
    for (auto &elem : fileids_) {
      if (elem > 0) {
        ex_close(elem);
      }
      elem = -1;
    }
  }
  ex_close(outputId_);
  outputId_ = -1;
}

bool Excn::ExodusFile::initialize(const SystemInterface &si)
{
  // See if we can keep files open
  size_t max_files = open_file_limit() - 1; // We also have an output exodus file...
  if (si.inputFiles_.size() <= max_files) {
    keepOpen_ = true;
    if ((si.debug() & 1) != 0) {
      fmt::print("Files kept open... (Max open = {})\n\n", max_files);
    }
  }
  else {
    keepOpen_ = false;
    fmt::print("Single file mode... (Max open = {})\n\n", max_files);
  }

  float   version                 = 0.0;
  int64_t overall_max_name_length = 32;

  if (si.inputFiles_.size() == 1) {
    // The file should contain multiple change sets which will be concatenated...
    int         cpu_word_size = sizeof(double);
    int         io_wrd_size   = sizeof(double);
    std::string name          = si.inputFiles_[0];
    int         exoid = ex_open(name.c_str(), EX_READ, &cpu_word_size, &io_wrd_size, &version);
    if (exoid < 0) {
      fmt::print(stderr, "ERROR: Cannot open file '{}'\n", name);
      return false;
    }

    if (auto num_change_sets = ex_inquire_int(exoid, EX_INQ_NUM_CHILD_GROUPS);
        num_change_sets > 1) {
      usingChangeSets_ = true;

      if (io_wrd_size < static_cast<int>(sizeof(float))) {
        io_wrd_size = sizeof(float);
      }

      ioWordSize_  = io_wrd_size;
      cpuWordSize_ = io_wrd_size;

      if (((ex_int64_status(exoid) & EX_ALL_INT64_DB) != 0) || si.ints_64_bit()) {
        exodusMode_ = EX_ALL_INT64_API;
      }

      overall_max_name_length =
          std::max(overall_max_name_length, ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH));

      // create exo names
      filenames_.resize(num_change_sets);
      fileids_.resize(num_change_sets, exoid);

      // Get names of change sets...
      auto group_name_length = ex_inquire_int(exoid, EX_INQ_GROUP_NAME_LEN);
      group_name_length      = std::max(int64_t(32), group_name_length);
      std::vector<char> group_name(group_name_length + 1, '\0');

      for (int i = 1; i <= num_change_sets; i++) {
        int   idum = 0;
        float rdum = 0.0;
        // Get name of this group...
        int ierr = ex_inquire(exoid + i, EX_INQ_GROUP_NAME, &idum, &rdum, group_name.data());
        if (ierr != EX_NOERR) {
          fmt::print(stderr, "ERROR: Could not get name for group {} in input file '{}'\n", i + 1,
                     name);
          return false;
        }
        filenames_[i - 1] = std::string(group_name.data());
        fileids_[i - 1]   = exoid + i;
        fmt::print("Part {}: '{}'\n", i, filenames_[i - 1]);
      }
    }
  }
  if (!usingChangeSets_) {
    // create exo names
    filenames_.resize(si.inputFiles_.size());
    fileids_.resize(si.inputFiles_.size(), -1);

    for (size_t p = 0; p < si.inputFiles_.size(); p++) {
      std::string name = si.inputFiles_[p];

      filenames_[p] = name;

      if (p == 0) {
        int cpu_word_size = sizeof(float);
        int io_wrd_size   = 0;
        int exoid = ex_open(filenames_[p].c_str(), EX_READ, &cpu_word_size, &io_wrd_size, &version);
        if (exoid < 0) {
          fmt::print(stderr, "ERROR: Cannot open file '{}'\n", filenames_[p]);
          return false;
        }

        overall_max_name_length = std::max(overall_max_name_length,
                                           ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH));

        if (((ex_int64_status(exoid) & EX_ALL_INT64_DB) != 0) || si.ints_64_bit()) {
          exodusMode_ = EX_ALL_INT64_API;
        }

        ex_close(exoid);

        if (io_wrd_size < static_cast<int>(sizeof(float))) {
          io_wrd_size = sizeof(float);
        }

        ioWordSize_  = io_wrd_size;
        cpuWordSize_ = io_wrd_size;
      }

      if (keepOpen_ || p == 0) {
        int io_wrd_size = 0;
        int mode        = EX_READ | exodusMode_;

        fileids_[p] = ex_open(filenames_[p].c_str(), mode, &cpuWordSize_, &io_wrd_size, &version);
        if (fileids_[p] < 0) {
          fmt::print(stderr, "ERROR: Cannot open file '{}'\n", filenames_[p]);
          return false;
        }
        if (auto num_change_sets = ex_inquire_int(fileids_[p], EX_INQ_NUM_CHILD_GROUPS);
            num_change_sets > 1) {
          fmt::print(stderr,
                     "ERROR: Cannot (yet) handle multiple input files containing change "
                     "sets.\n\tFile '{}' contains {} change sets.\n",
                     filenames_[p], num_change_sets);
          return false;
        }

        SMART_ASSERT(ioWordSize_ == io_wrd_size)(ioWordSize_)(io_wrd_size);
      }

      fmt::print("Part {}: '{}'\n", p + 1, name);
    }
  }

  maximumNameLength_ = overall_max_name_length;
  if (keepOpen_) {
    for (size_t p = 0; p < si.inputFiles_.size(); p++) {
      ex_set_max_name_length(fileids_[p], maximumNameLength_);
    }
  }
  else {
    ex_set_max_name_length(fileids_[0], maximumNameLength_);
  }
  return true;
}

bool Excn::ExodusFile::create_output(const SystemInterface &si)
// Create output file...
{
  outputFilename_ = si.outputName_;

  int mode = EX_CLOBBER | exodusMode_;
  if (si.ints_64_bit()) {
    mode |= EX_ALL_INT64_DB;
  }

  if (si.compress_data() > 0 || si.use_netcdf4() || si.szip()) {
    mode |= EX_NETCDF4;
  }

  fmt::print("Output:   '{}'\n", outputFilename_);
  outputId_ = ex_create(outputFilename_.c_str(), mode, &cpuWordSize_, &ioWordSize_);
  if (outputId_ < 0) {
    fmt::print(stderr, "ERROR: Cannot open file '{}'\n", outputFilename_);
    return false;
  }

  if (si.compress_data() > 0) {
    ex_set_option(outputId_, EX_OPT_COMPRESSION_LEVEL, si.compress_data());
    ex_set_option(outputId_, EX_OPT_COMPRESSION_SHUFFLE, 1);

    if (si.szip()) {
      ex_set_option(outputId_, EX_OPT_COMPRESSION_TYPE, EX_COMPRESS_SZIP);
    }
    else if (si.zlib()) {
      ex_set_option(outputId_, EX_OPT_COMPRESSION_TYPE, EX_COMPRESS_ZLIB);
    }
    else if (si.zstd()) {
      ex_set_option(outputId_, EX_OPT_COMPRESSION_TYPE, EX_COMPRESS_ZSTD);
    }
    else if (si.bz2()) {
      ex_set_option(outputId_, EX_OPT_COMPRESSION_TYPE, EX_COMPRESS_BZ2);
    }

    if (si.quantize()) {
      ex_set_option(outputId_, EX_OPT_QUANTIZE_NSD, si.quantize_nsd());
    }
  }

  fmt::print("IO Word size is {} bytes.\n", ioWordSize_);
  ex_set_max_name_length(outputId_, maximumNameLength_);
  return true;
}
