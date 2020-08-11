/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "EP_ExodusFile.h"
#include "EP_Internals.h"
#include "EP_ParallelDisks.h"
#include "EP_SystemInterface.h"
#include "fmt/ostream.h"
#include "smart_assert.h"
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if !defined(_WIN32)
#include <unistd.h>
#endif

#include <exodusII.h>

std::vector<int>         Excn::ExodusFile::fileids_;
std::vector<std::string> Excn::ExodusFile::filenames_;
int                      Excn::ExodusFile::processorCount_ = 0;
int                      Excn::ExodusFile::partCount_      = 0;
int                      Excn::ExodusFile::startPart_      = 0;
int                      Excn::ExodusFile::outputId_       = -1;
int                      Excn::ExodusFile::ioWordSize_     = 0;
int                      Excn::ExodusFile::cpuWordSize_    = 0;
int                      Excn::ExodusFile::mode64bit_      = 0;
std::string              Excn::ExodusFile::outputFilename_;
bool                     Excn::ExodusFile::keepOpen_          = false;
int                      Excn::ExodusFile::maximumNameLength_ = 32;

Excn::ExodusFile::ExodusFile(int processor) : myProcessor_(processor)
{
  SMART_ASSERT(processor < processorCount_)(processor)(processorCount_);
  SMART_ASSERT(fileids_.size() == (size_t)processorCount_);
  if (!keepOpen_ && processor != 0) {
    float version          = 0.0;
    int   cpu_word_size    = cpuWordSize_;
    int   io_word_size_var = ioWordSize_;
    int   mode             = EX_READ;
    mode |= mode64bit_;

    fileids_[processor] =
        ex_open(filenames_[processor].c_str(), mode, &cpu_word_size, &io_word_size_var, &version);
    if (fileids_[processor] < 0) {
      std::ostringstream errmsg;
      fmt::print(errmsg, "Cannot open file '{}' - exiting\n", filenames_[processor]);
      throw std::runtime_error(errmsg.str());
    }
    ex_set_max_name_length(fileids_[processor], maximumNameLength_);

    SMART_ASSERT(io_word_size_var == ioWordSize_);
    SMART_ASSERT(cpu_word_size == cpuWordSize_);
  }
}

int Excn::ExodusFile::output()
{
  SMART_ASSERT(outputId_ >= 0);
  return outputId_;
}

Excn::ExodusFile::operator int() const
{
  SMART_ASSERT(fileids_[myProcessor_] >= 0);
  return fileids_[myProcessor_];
}

Excn::ExodusFile::~ExodusFile()
{
  try {
    if (!keepOpen_ && myProcessor_ != 0) {
      ex_close(fileids_[myProcessor_]);
      fileids_[myProcessor_] = -1;
    }
  }
  catch (...) {
  }
}

void Excn::ExodusFile::close_all()
{
  for (int p = 0; p < partCount_; p++) {
    if (fileids_[p] >= 0) {
      ex_close(fileids_[p]);
      fileids_[p] = -1;
    }
  }
  if (outputId_ >= 0) {
    ex_close(outputId_);
    outputId_ = -1;
  }
}

void Excn::ExodusFile::unlink_temporary_files()
{
  for (int p = 0; p < partCount_; p++) {
    if (fileids_[p] >= 0) {
      ex_close(fileids_[p]);
      unlink(filenames_[p].c_str());
      fileids_[p] = -1;
    }
  }
  if (outputId_ >= 0) {
    ex_close(outputId_);
    outputId_ = -1;
  }
}

bool Excn::ExodusFile::initialize(const SystemInterface &si, int start_part, int part_count,
                                  int cycle, bool joining_subcycle)
{
  processorCount_ = si.processor_count(); // Total number processors
  partCount_      = part_count;           // Total parts we are processing
  startPart_      = start_part;           // Which one to start with
  SMART_ASSERT(partCount_ + startPart_ <= processorCount_)(partCount_)(startPart_)(processorCount_);

  // EPU always wants entity (block, set, map) ids as 64-bit quantities...
  mode64bit_ = EX_IDS_INT64_API;
  if (si.int64()) {
    mode64bit_ |= EX_ALL_INT64_API;

    // For output...
    mode64bit_ |= EX_ALL_INT64_DB;
  }

  // See if we can keep files open
  int max_files = get_free_descriptor_count();
  if (partCount_ <= max_files) {
    keepOpen_ = true;
    if (cycle == 0) {
      if ((si.debug() & 1) != 0) {
        fmt::print("Files kept open... (Max open = {})\n\n", max_files);
      }
    }
  }
  else {
    keepOpen_ = false;
    if (cycle == 0) {
      fmt::print("Single file mode... (Max open = {})\n"
                 "Consider using the -subcycle option for faster execution...\n\n",
                 max_files);
    }
  }

  fileids_.resize(processorCount_);
  filenames_.resize(processorCount_);

  std::string file_prefix   = si.basename();
  std::string exodus_suffix = si.exodus_suffix();

  std::string root_dir = si.root_dir();
  std::string sub_dir  = si.sub_dir();

  ParallelDisks disks;
  disks.Raid_Offset(si.raid_offset());
  disks.Number_of_Raids(si.raid_count());

  float version = 0.0;

  // create exo names
  for (int p = 0; p < partCount_; p++) {
    std::string name = file_prefix;
    if (!exodus_suffix.empty()) {
      name += "." + exodus_suffix;
    }

    int proc = p + startPart_;
    if (joining_subcycle) {
      name = si.output_filename();
      disks.rename_file_for_mp("", "", name, proc, processorCount_);
    }
    else {
      disks.rename_file_for_mp(root_dir, sub_dir, name, proc, processorCount_);
    }
    filenames_[p] = name;

    if (p == 0) {
      int cpu_word_size    = sizeof(float);
      int io_word_size_var = 0;
      int mode             = EX_READ;
      mode |= mode64bit_;
      int exoid = ex_open(filenames_[p].c_str(), mode, &cpu_word_size, &io_word_size_var, &version);
      if (exoid < 0) {
        fmt::print(stderr, "Cannot open file '{}'\n", filenames_[p]);
        return false;
      }

      int int64db = ex_int64_status(exoid) & EX_ALL_INT64_DB;
      if (int64db != 0) {
        // If anything stored on input db as 64-bit int, then output db will have
        // everything stored as 64-bit ints and all API functions will use 64-bit
        mode64bit_ |= EX_ALL_INT64_API;
        mode64bit_ |= EX_ALL_INT64_DB;
      }

      int max_name_length = static_cast<int>(ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH));
      if (max_name_length > maximumNameLength_) {
        maximumNameLength_ = max_name_length;
      }

      ex_close(exoid);

      if (io_word_size_var < static_cast<int>(sizeof(float))) {
        io_word_size_var = sizeof(float);
      }

      ioWordSize_  = io_word_size_var;
      cpuWordSize_ = io_word_size_var;
    }

    if (keepOpen_ || p == 0) {
      int io_word_size_var = 0;
      int mode             = EX_READ;
      // All entity ids (block, sets) are read/written as 64-bit...
      mode |= mode64bit_;
      fileids_[p] =
          ex_open(filenames_[p].c_str(), mode, &cpuWordSize_, &io_word_size_var, &version);
      if (fileids_[p] < 0) {
        fmt::print(stderr, "Cannot open file '{}'\n", filenames_[p]);
        return false;
      }
      ex_set_max_name_length(fileids_[p], maximumNameLength_);
      SMART_ASSERT(ioWordSize_ == io_word_size_var)(ioWordSize_)(io_word_size_var);
    }

    if (((si.debug() & 64) != 0) || p == 0 || p == partCount_ - 1) {
      fmt::print("[{}] Input({}): '{}'\n", cycle, p, name);
      if (((si.debug() & 64) == 0) && p == 0) {
        fmt::print("...\n");
      }
    }
  }

  if ((mode64bit_ & EX_ALL_INT64_DB) != 0) {
    if (cycle == 0) {
      fmt::print("Input files contain 8-byte integers.\n");
    }
    si.set_int64();
  }

  return true;
}

bool Excn::ExodusFile::create_output(const SystemInterface &si, int cycle)
// Create output file...
{
  std::string curdir        = si.cwd();
  std::string file_prefix   = si.basename();
  std::string output_suffix = si.output_suffix();

  outputFilename_ = file_prefix;
  if (!output_suffix.empty()) {
    outputFilename_ += "." + output_suffix;
  }

  if (!curdir.empty() && !Excn::is_path_absolute(outputFilename_)) {
    outputFilename_ = curdir + "/" + outputFilename_;
  }

  si.set_output_filename(outputFilename_);

  if (si.subcycle() > 1) {
    Excn::ParallelDisks::Create_IO_Filename(outputFilename_, cycle, si.subcycle());
  }

  // See if output file should be opened in netcdf4 format...
  // Did user specify it via -netcdf4 or -large_model argument...
  int mode = 0;

  if (si.compress_data() > 0) {
    // Force netcdf-4 if compression is specified...
    mode |= EX_NETCDF4;
  }
  else if (si.use_netcdf4()) {
    mode |= EX_NETCDF4;
  }
  else if (ex_large_model(fileids_[0]) == 1) {
    mode |= EX_LARGE_MODEL;
  }

  mode |= mode64bit_;

  if (si.int64()) {
    mode |= EX_ALL_INT64_DB;
    mode |= EX_ALL_INT64_API;
  }

  if (si.append()) {
    fmt::print("[{}] Output:   '{}' (appending)\n", cycle, outputFilename_);
    float version = 0.0;
    mode |= EX_WRITE;
    outputId_ = ex_open(outputFilename_.c_str(), mode, &cpuWordSize_, &ioWordSize_, &version);
  }
  else {
    mode |= EX_CLOBBER;
    fmt::print("[{}] Output:   '{}'\n", cycle, outputFilename_);
    outputId_ = ex_create(outputFilename_.c_str(), mode, &cpuWordSize_, &ioWordSize_);
  }
  if (outputId_ < 0) {
    fmt::print(stderr, "Cannot open file '{}'\n", outputFilename_);
    return false;
  }

  if (si.compress_data() > 0) {
    ex_set_option(outputId_, EX_OPT_COMPRESSION_LEVEL, si.compress_data());
    ex_set_option(outputId_, EX_OPT_COMPRESSION_SHUFFLE, 1);
  }

  // EPU Can add a name of "processor_id_epu" which is 16 characters long.
  // Make sure maximumNameLength_ is at least that long...

  if (maximumNameLength_ < 16) {
    maximumNameLength_ = 16;
  }
  ex_set_option(outputId_, EX_OPT_MAX_NAME_LENGTH, maximumNameLength_);

  int int_size = si.int64() ? 8 : 4;
  if (cycle == 0) {
    fmt::print("IO Word sizes: {} bytes floating point and {} bytes integer.\n", ioWordSize_,
               int_size);
  }
  return true;
}

size_t Excn::ExodusFile::get_free_descriptor_count()
{
// Returns maximum number of files that one process can have open
// at one time. (POSIX)
#if defined(_WIN32)
  int fdmax = _getmaxstdio();
#else
  int fdmax = sysconf(_SC_OPEN_MAX);
  if (fdmax == -1) {
    // POSIX indication that there is no limit on open files...
    fdmax = INT_MAX;
  }
#endif
  // File descriptors are assigned in order (0,1,2,3,...) on a per-process
  // basis.

  // Assume that we have stdin, stdout, stderr, and output exodus
  // file (4 total).

  return fdmax - 4;

  // Could iterate from 0..fdmax and check for the first
  // EBADF (bad file descriptor) error return from fcntl, but that takes
  // too long and may cause other problems.  There is a isastream(filedes)
  // call on Solaris that can be used for system-dependent code.
  //
  // Another possibility is to do an open and see which descriptor is
  // returned -- take that as 1 more than the current count of open files.
  //
}
