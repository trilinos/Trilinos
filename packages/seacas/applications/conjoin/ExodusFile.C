// Copyright(C) 2009-2010 Sandia Corporation.
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include <ExodusFile.h>
#include <SystemInterface.h>
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
int Excn::ExodusFile::outputId_ = -1;
int Excn::ExodusFile::ioWordSize_ = 0;
int Excn::ExodusFile::cpuWordSize_ = 0;
std::string Excn::ExodusFile::outputFilename_;
bool Excn::ExodusFile::keepOpen_ = false;
int Excn::ExodusFile::maximumNameLength_ = 32;

namespace {
  int get_free_descriptor_count();
}

Excn::ExodusFile::ExodusFile(size_t which)
  : myLocation_(which)
{
  SMART_ASSERT(which < filenames_.size())(which)(filenames_.size());
  SMART_ASSERT(fileids_.size() == filenames_.size());
  if (!keepOpen_ && which != 0) {
    float version = 0.0;
    int cpu_word_size = cpuWordSize_;
    int io_wrd_size  = ioWordSize_;
    fileids_[which] = ex_open(filenames_[which].c_str(),
				  EX_READ, &cpu_word_size,
				  &io_wrd_size, &version);
    if (fileids_[which] < 0) {
      std::cerr << "Cannot open file '" << filenames_[which]
	   << "' - exiting" << std::endl;
      exit(1);
    }
    ex_set_max_name_length(fileids_[which], maximumNameLength_);

    SMART_ASSERT(io_wrd_size  == ioWordSize_);
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
  SMART_ASSERT(fileids_[myLocation_] >= 0);
  return fileids_[myLocation_];
}

Excn::ExodusFile::~ExodusFile()
{
  try {
    if (!keepOpen_ && myLocation_ != 0) {
      ex_close(fileids_[myLocation_]);
      fileids_[myLocation_] = -1;
    }
  } catch (...) {
  }
}

void Excn::ExodusFile::close_all()
{
  for(size_t p = 0; p < fileids_.size(); p++) {
    ex_close(fileids_[p]);
    fileids_[p] = -1;
  }
  ex_close(outputId_);
  outputId_ = -1;
}

bool Excn::ExodusFile::initialize(const SystemInterface& si)
{
  // See if we can keep files open 
  size_t max_files = get_free_descriptor_count();
  if (si.inputFiles_.size() <= max_files) {
    keepOpen_ = true;
    if (si.debug() & 1)
      std::cout << "Files kept open... (Max open = " << max_files << ")\n\n";
  } else {
    keepOpen_ = false;
    std::cout << "Single file mode... (Max open = " << max_files << ")\n"
	      << "Consider using the -subcycle option for faster execution...\n\n";
  }

  float version = 0.0;

  // create exo names
  filenames_.resize(si.inputFiles_.size());
  fileids_.resize(si.inputFiles_.size());
  
  int overall_max_name_length = 0;
  for(size_t p = 0; p < si.inputFiles_.size(); p++) {
    std::string name = si.inputFiles_[p];

    filenames_[p] = name;

    if (p == 0) {
      int cpu_word_size  = sizeof(float);
      int io_wrd_size   = 0;
      int exoid = ex_open(filenames_[p].c_str(),
			  EX_READ, &cpu_word_size,
			  &io_wrd_size, &version);
      if (exoid < 0) {
	std::cerr << "Cannot open file '" << filenames_[p] << "'" << std::endl;
	return false;
      }

      int max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
      if (max_name_length > overall_max_name_length)
	overall_max_name_length = max_name_length;

      ex_close(exoid);

      if (io_wrd_size < (int)sizeof(float))
	io_wrd_size = sizeof(float);

      ioWordSize_ = io_wrd_size;
      cpuWordSize_ = io_wrd_size;
    }

    if (keepOpen_ || p == 0) {
      int io_wrd_size   = 0;
      fileids_[p] = ex_open(filenames_[p].c_str(),
			    EX_READ, &cpuWordSize_,
			    &io_wrd_size, &version);
      if (fileids_[p] < 0) {
	std::cerr << "Cannot open file '" << filenames_[p] << "'" << std::endl;
	return false;
      }
      SMART_ASSERT(ioWordSize_ == io_wrd_size)(ioWordSize_)(io_wrd_size);
    }
    
    std::cout << "Part " << p+1 << ": '" << name.c_str() << "'" << std::endl;
  }

  maximumNameLength_ = overall_max_name_length;
  for(size_t p = 0; p < si.inputFiles_.size(); p++) {
    ex_set_max_name_length(fileids_[p], maximumNameLength_);
  }

  return true;
}

bool Excn::ExodusFile::create_output(const SystemInterface& si)
  // Create output file...
{
  outputFilename_ = si.outputName_;

  int mode = EX_CLOBBER;
  
  std::cout << "Output:   '" << outputFilename_ << "'" << std::endl;
  outputId_ = ex_create(outputFilename_.c_str(), mode,
			&cpuWordSize_, &ioWordSize_);
  if (outputId_ < 0) {
    std::cerr << "Cannot open file '" << outputFilename_ << "'" << std::endl;
    return false;
  }
  std::cout << "IO Word size is " << ioWordSize_ << " bytes.\n";
  ex_set_max_name_length(outputId_, maximumNameLength_);
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
