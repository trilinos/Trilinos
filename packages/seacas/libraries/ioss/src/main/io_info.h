/*
 * Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "io_info_lib_export.h"

namespace Info {
  class Interface;
} // namespace Info

namespace Ioss {
  class DatabaseIO;
  class Region;

  // internal to io_info
  IO_INFO_LIB_EXPORT void io_info_file_info(const Info::Interface &interFace);
  IO_INFO_LIB_EXPORT void io_info_change_set_info(Info::Interface &interFace);

  // for external calls
  IO_INFO_LIB_EXPORT void io_info_set_db_properties(const Info::Interface &interFace,
                                                    Ioss::DatabaseIO      *dbi);
  IO_INFO_LIB_EXPORT void io_info_file_info(const Info::Interface &interFace, Ioss::Region &region);
} // namespace Ioss
