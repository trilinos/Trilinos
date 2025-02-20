/*
 * Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <iosfwd> // for ostream
#include <string> // for string

#include "Ioss_GetLongOpt.h" // for GetLongOption
#include "io_info_lib_export.h"

/** \brief A special namespace for the io_info demonstration program interFace.
 */
namespace Info {
  class IO_INFO_LIB_EXPORT Interface
  {
  public:
    explicit Interface(std::string app_version);

    bool parse_options(int argc, char **argv);

    bool summary() const { return summary_; }
    bool check_node_status() const { return checkNodeStatus_; }
    bool compute_volume() const { return computeVolume_; }
    bool compute_bbox() const { return computeBBox_; }
    bool adjacencies() const { return adjacencies_; }
    bool ints_64_bit() const { return ints64Bit_; }
    bool list_change_sets() const { return listChangeSets_; }
    bool show_config() const { return showConfig_; }
    bool query_timesteps_only() const { return queryTimeOnly_; }
    bool field_details() const { return fieldDetails_; }

    int         surface_split_scheme() const { return surfaceSplitScheme_; }
    char        field_suffix_separator() const { return fieldSuffixSeparator_; }
    bool        use_generic_names() const { return useGenericNames_; }
    bool        disable_field_recognition() const { return disableFieldRecognition_; }
    std::string decomp_method() const { return decompMethod_; }
    std::string filename() const { return filename_; }
    std::string type() const { return filetype_; }
    std::string change_set_name() const { return changeSetName_; }
    std::string custom_field() const { return customField_; }

    std::string version{};

  private:
    void enroll_options();

    Ioss::GetLongOption options_;
    std::string         filetype_{"exodus"};
    std::string         filename_{};
    std::string         changeSetName_{};
    std::string         decompMethod_{};
    std::string         customField_{};

    bool checkNodeStatus_{false};
    bool computeVolume_{false};
    bool adjacencies_{false};
    bool ints64Bit_{false};
    bool computeBBox_{false};
    bool listChangeSets_{false};
    bool useGenericNames_{false};
    bool disableFieldRecognition_{false};
    bool showConfig_{false};
    bool summary_{false};
    bool queryTimeOnly_{false};
    bool fieldDetails_{false};

    char fieldSuffixSeparator_{'_'};

    int surfaceSplitScheme_{-1};
  };
} // namespace Info
