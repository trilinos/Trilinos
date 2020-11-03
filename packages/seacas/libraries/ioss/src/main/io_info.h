/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef Ioss_io_info_h
#define Ioss_io_info_h

#include "info_interface.h"

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Utils.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#if defined(SEACAS_HAVE_EXODUS)
#include <exodusII.h>
#endif

#include "Ioss_Assembly.h"
#include "Ioss_Blob.h"
#include "Ioss_CommSet.h"
#include "Ioss_CoordinateFrame.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_StructuredBlock.h"
#include "Ioss_VariableType.h"

#include <cassert>

namespace Ioss {

  // internal to io_info
  void io_info_file_info(const Info::Interface &interFace);
  void io_info_group_info(Info::Interface &interFace);

  // for external calls
  void io_info_set_db_properties(const Info::Interface &interFace, Ioss::DatabaseIO *dbi);
  void io_info_file_info(const Info::Interface &interFace, Ioss::Region &region);
} // namespace Ioss

#endif
