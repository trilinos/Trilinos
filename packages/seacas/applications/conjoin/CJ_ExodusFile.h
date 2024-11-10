// Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

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

    ExodusFile(const ExodusFile &)           = delete;
    ExodusFile operator=(const ExodusFile &) = delete;

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
  };
} // namespace Excn
