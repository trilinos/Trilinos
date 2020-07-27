/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef info_SystemInterface_h
#define info_SystemInterface_h

#include "Ioss_GetLongOpt.h" // for GetLongOption
#include <iosfwd>            // for ostream
#include <string>            // for string

/** \brief A special namespace for the io_info demonstration program interFace.
 */
namespace Info {
  class Interface
  {
  public:
    Interface();
    ~Interface();

    bool parse_options(int argc, char **argv);

    bool summary() const { return summary_; }
    bool check_node_status() const { return checkNodeStatus_; }
    bool compute_volume() const { return computeVolume_; }
    bool compute_bbox() const { return computeBBox_; }
    bool adjacencies() const { return adjacencies_; }
    bool ints_64_bit() const { return ints64Bit_; }
    bool list_groups() const { return listGroups_; }
    bool show_config() const { return showConfig_; }

    int         surface_split_scheme() const { return surfaceSplitScheme_; }
    char        field_suffix_separator() const { return fieldSuffixSeparator_; }
    bool        use_generic_names() const { return useGenericNames_; }
    bool        disable_field_recognition() const { return disableFieldRecognition_; }
    std::string decomp_method() const { return decompMethod_; }
    std::string filename() const { return filename_; }
    std::string type() const { return filetype_; }
    std::string groupname() const { return groupname_; }

    //! Dumps representation of data in this class to cerr

  private:
    void enroll_options();

    Ioss::GetLongOption options_;
    std::string         filetype_{"exodus"};
    std::string         filename_{};
    std::string         groupname_{};
    std::string         decompMethod_{};

    bool checkNodeStatus_{false};
    bool computeVolume_{false};
    bool adjacencies_{false};
    bool ints64Bit_{false};
    bool computeBBox_{false};
    bool listGroups_{false};
    bool useGenericNames_{false};
    bool disableFieldRecognition_{false};
    bool showConfig_{false};
    bool summary_{false};

    char fieldSuffixSeparator_{'_'};

    int surfaceSplitScheme_{1};
  };
} // namespace Info
#endif
