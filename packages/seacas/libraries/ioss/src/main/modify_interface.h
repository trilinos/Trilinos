/*
 * Copyright(C) 2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE file for details.
 */

#ifndef modify_SystemInterface_h
#define modify_SystemInterface_h

#include "Ioss_GetLongOpt.h" // for GetLongOption
#include <string>            // for string

/** \brief A special namespace for the io_info demonstration program interFace.
 */
namespace Modify {
  class Interface
  {
  public:
    Interface();
    ~Interface();

    bool parse_options(int argc, char **argv);

    std::string filename() const { return filename_; }
    std::string type() const { return filetype_; }
    bool        modify_existing_assembly() const { return allowModification_; }

  private:
    void enroll_options();

    Ioss::GetLongOption options_;
    std::string         filetype_{"exodus"};
    std::string         filename_{};
    bool                allowModification_{false};
  };
} // namespace Modify
#endif
