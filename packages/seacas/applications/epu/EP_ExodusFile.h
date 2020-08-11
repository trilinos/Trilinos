/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SEACAS_ExodusFile_H
#define SEACAS_ExodusFile_H

#include <string>
#include <vector>

namespace Excn {

  class SystemInterface;
  class ExodusFile
  {
  public:
    explicit ExodusFile(int processor);
    ~ExodusFile();
    ExodusFile(const ExodusFile &) = delete;
    ExodusFile operator=(const ExodusFile &) = delete;

    static bool initialize(const SystemInterface &si, int start_part, int part_count, int cycle,
                           bool joining_subcycle);
    static bool create_output(const SystemInterface &si, int cycle);
    static void close_all();

    static int    output();
    static int    io_word_size() { return ioWordSize_; }
                  operator int() const;
    static int    max_name_length() { return maximumNameLength_; }
    static size_t get_free_descriptor_count();
    static void   unlink_temporary_files();

  private:
    int                             myProcessor_;
    static std::vector<std::string> filenames_;
    static std::vector<int>         fileids_;
    static int                      processorCount_;
    static int                      partCount_;
    static int                      startPart_;
    static int                      outputId_;
    static int                      ioWordSize_;
    static int                      cpuWordSize_;
    static std::string              outputFilename_;
    static bool                     keepOpen_;
    static int                      maximumNameLength_;
    static int                      mode64bit_;
  };
} // namespace Excn
#endif /* SEACAS_ExodusFil_H */
