// Copyright(C) 2009-2010-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//     * Neither the name of NTESS nor the names of its
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
#ifndef SEACAS_ExodusFile_H
#define SEACAS_ExodusFile_H

#include "CJ_CodeTypes.h"
#include <string>
#include <vector>

namespace Excn {

  class SystemInterface;
  class ExodusFile
  {
  public:
    explicit ExodusFile(size_t which);
    ~ExodusFile();

    static bool initialize(const SystemInterface &si);
    static bool create_output(const SystemInterface &si);
    static void close_all();

    static int output();
    static int io_word_size() { return ioWordSize_; }
               operator int() const;
    static int max_name_length() { return maximumNameLength_; }

  private:
    size_t                          myLocation_;
    static std::vector<std::string> filenames_;
    static std::vector<int>         fileids_;
    static int                      outputId_;
    static int                      ioWordSize_;
    static int                      cpuWordSize_;
    static std::string              outputFilename_;
    static bool                     keepOpen_;
    static int                      maximumNameLength_;
    static int                      exodusMode_;

    // Disable copying and assignment...
    ExodusFile(const ExodusFile &);
    ExodusFile operator=(const ExodusFile &);
  };
} // namespace Excn
#endif /* SEACAS_ExodusFil_H */
