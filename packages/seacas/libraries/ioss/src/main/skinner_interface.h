/*
 * Copyright(C) 1999-2020, 2022, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <iosfwd> // for ostream
#include <limits>
#include <string> // for string
#include <vector>

#include "Ioss_GetLongOpt.h" // for GetLongOption

/** \brief A special namespace for the skinner demonstration program interFace.
 */
namespace Skinner {
  class Interface
  {
  public:
    Interface();

    bool parse_options(int argc, char **argv);

    bool ints_64_bit() const { return ints64Bit_; }

    bool no_output() const { return noOutput_; }
    bool output_transient() const { return outputTransient_; }

    std::string input_filename() const { return inputFile_; }
    std::string output_filename() const { return outputFile_; }
    std::string input_type() const { return inFiletype_; }
    std::string output_type() const { return outFiletype_; }

  private:
    void enroll_options();

    Ioss::GetLongOption options_;

    std::string inputFile_;
    std::string outputFile_;
    std::string inFiletype_{"unknown"};
    std::string outFiletype_{"unknown"};

  public:
    std::string decomp_method;
    std::string compose_output{"default"};
    double      maximum_time{std::numeric_limits<double>::max()};
    double      minimum_time{-std::numeric_limits<double>::max()};
    int         compression_level{0};
    bool        shuffle{false};
    bool        debug{false};
    bool        statistics{false};
    bool        ints64Bit_{false};
    bool        netcdf4_{false};
    bool        useFaceHashIds_{true};
    bool        noOutput_{false};
    bool        outputTransient_{false};
    bool        blocks_{false};

    //! If non-empty, then it is a list of times that should be transferred to the output file.
    std::vector<double> selected_times{};
  };
} // namespace Skinner
