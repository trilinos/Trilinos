/*
 * Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include <ExodusFile.h>
#include <SystemInterface.h>
#include <ParallelDisks.h>
#include <smart_assert.h>
#include <stdlib.h>
#include <limits.h>
#include <to_string.h>

#include <iostream>
#include <iomanip>
#include <string>

#include <exodusII.h>

std::vector<int>         Excn::ExodusFile::fileids_;
std::vector<std::string> Excn::ExodusFile::filenames_;
int Excn::ExodusFile::processorCount_ = 0;
int Excn::ExodusFile::partCount_ = 0;
int Excn::ExodusFile::startPart_ = 0;
int Excn::ExodusFile::outputId_ = -1;
int Excn::ExodusFile::ioWordSize_ = 0;
int Excn::ExodusFile::cpuWordSize_ = 0;
std::string Excn::ExodusFile::outputFilename_;
bool Excn::ExodusFile::keepOpen_ = false;
int Excn::ExodusFile::maximumNameLength_ = 32;

namespace {
  int get_free_descriptor_count();
}

Excn::ExodusFile::ExodusFile(int processor)
  : myProcessor_(processor)
{
  SMART_ASSERT(processor < processorCount_)(processor)(processorCount_);
  SMART_ASSERT(fileids_.size() == (size_t)processorCount_);
  if (!keepOpen_ && processor != 0) {
    float version = 0.0;
    int cpu_word_size = cpuWordSize_;
    int io_word_size_var  = ioWordSize_;
    fileids_[processor] = ex_open(filenames_[processor].c_str(),
				  EX_READ, &cpu_word_size,
				  &io_word_size_var, &version);
    if (fileids_[processor] < 0) {
      std::cerr << "Cannot open file '" << filenames_[processor]
	   << "' - exiting" << std::endl;
      exit(1);
    }
    ex_set_max_name_length(fileids_[processor], maximumNameLength_);
    
    SMART_ASSERT(io_word_size_var  == ioWordSize_);
    SMART_ASSERT(cpu_word_size == cpuWordSize_);
  }
}

int Excn::ExodusFile::output()
{
  SMART_ASSERT(outputId_ >= 0);
  return outputId_;
}
  
Excn::ExodusFile::operator int () const
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
  } catch (...) {
  }
}

void Excn::ExodusFile::close_all()
{
  for(int p = 0; p < partCount_; p++) {
    ex_close(fileids_[p]);
    fileids_[p] = -1;
  }
  ex_close(outputId_);
  outputId_ = -1;
}

bool Excn::ExodusFile::initialize(const SystemInterface& si, int start_part, int part_count)
{
  processorCount_ = si.processor_count(); // Total number processors
  partCount_      = part_count;      // Total parts we are processing
  startPart_      = start_part;      // Which one to start with
  SMART_ASSERT(partCount_ + startPart_ <= processorCount_)(partCount_)(startPart_)(processorCount_);
  
  // See if we can keep files open 
  int max_files = get_free_descriptor_count();
  if (partCount_ <= max_files) {
    keepOpen_ = true;
    if (si.debug() & 1)
      std::cout << "Files kept open... (Max open = " << max_files << ")\n\n";
  } else {
    keepOpen_ = false;
    std::cout << "Single file mode... (Max open = " << max_files << ")\n"
	      << "Consider using the -subcycle option for faster execution...\n\n";
  }

  fileids_.resize(processorCount_);
  filenames_.resize(processorCount_);

  std::string curdir = si.cwd();
  std::string file_prefix = si.basename();
  std::string exodus_suffix = si.exodus_suffix();
  
  std::string root_dir = si.root_dir();
  std::string sub_dir  = si.sub_dir();

  ParallelDisks disks;
  disks.Raid_Offset(si.raid_offset());
  disks.Number_of_Raids(si.raid_count());

  float version = 0.0;

  // create exo names
  for(int p = 0; p < partCount_; p++) {
    std::string name = file_prefix + "." + exodus_suffix;

    int proc = p + startPart_;
    disks.rename_file_for_mp(root_dir, sub_dir, name, proc,
			     processorCount_); 
    filenames_[p] = name;

    if (p == 0) {
      int cpu_word_size  = sizeof(float);
      int io_word_size_var   = 0;
      int exoid = ex_open(filenames_[p].c_str(),
			  EX_READ, &cpu_word_size,
			  &io_word_size_var, &version);
      if (exoid < 0) {
	std::cerr << "Cannot open file '" << filenames_[p] << "'" << std::endl;
	return false;
      }

      int max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
      if (max_name_length > maximumNameLength_)
	maximumNameLength_ = max_name_length;

      ex_close(exoid);

      if (io_word_size_var < (int)sizeof(float))
	io_word_size_var = sizeof(float);

      ioWordSize_ = io_word_size_var;
      cpuWordSize_ = io_word_size_var;
    }

    if (keepOpen_ || p == 0) {
      int io_word_size_var   = 0;
      fileids_[p] = ex_open(filenames_[p].c_str(),
			    EX_READ, &cpuWordSize_,
			    &io_word_size_var, &version);
      if (fileids_[p] < 0) {
	std::cerr << "Cannot open file '" << filenames_[p] << "'" << std::endl;
	return false;
      }
      ex_set_max_name_length(fileids_[p], maximumNameLength_);
      SMART_ASSERT(ioWordSize_ == io_word_size_var)(ioWordSize_)(io_word_size_var);
    }
    
    if (si.debug()&2 || p==0 || p==partCount_-1) {
      std::cout << "Input(" << p << "): '" << name.c_str() << "'" << std::endl;
      if (!(si.debug()&2) && p==0)
	std::cout << "..." << std::endl;
    }
  }
  return true;
}

bool Excn::ExodusFile::create_output(const SystemInterface& si, int cycle)
  // Create output file...
{
  std::string curdir = si.cwd();
  std::string file_prefix = si.basename();
  std::string output_suffix = si.output_suffix();

  outputFilename_ = file_prefix + "." + output_suffix;
  if(curdir.length() && outputFilename_[0] != '/') {
    outputFilename_ = curdir + "/" + outputFilename_;
  }

  if (si.subcycle() > 1) {
    Excn::ParallelDisks::Create_IO_Filename(outputFilename_,
					    cycle, si.subcycle());
  }

  // See if output file should be opened in the large model format...
  // Did user specify it via -large_model argument...
  int mode = EX_CLOBBER;
  
  if (si.large_model() || ex_large_model(fileids_[0]) == 1) {
    mode |= EX_LARGE_MODEL;
  }
  
  if (si.append()) {
    std::cout << "Output:   '" << outputFilename_ << "' (appending)" << std::endl;
    float version = 0.0;
    mode = EX_WRITE;
    outputId_ = ex_open(outputFilename_.c_str(), mode, 
			&cpuWordSize_, &ioWordSize_, &version);
  } else {
    std::cout << "Output:   '" << outputFilename_ << "'" << std::endl;
    outputId_ = ex_create(outputFilename_.c_str(), mode,
			  &cpuWordSize_, &ioWordSize_);
  }
  if (outputId_ < 0) {
    std::cerr << "Cannot open file '" << outputFilename_ << "'" << std::endl;
    return false;
  }

  // EPU Can add a name of "processor_id_epu" which is 16 characters long.
  // Make sure maximumNameLength_ is at least that long...
  
  if (maximumNameLength_ < 16)
    maximumNameLength_ = 16;
  ex_set_max_name_length(outputId_, maximumNameLength_);

  std::cout << "IO Word size is " << ioWordSize_ << " bytes.\n";
  return true;
}

#if defined(__PUMAGON__)
#include <stdio.h>
#else
#include <unistd.h>
#endif

namespace {
  int get_free_descriptor_count()
  {
    // Returns maximum number of files that one process can have open
    // at one time. (POSIX)
#if defined(__PUMAGON__)
    int fdmax = FOPEN_MAX;
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
}
