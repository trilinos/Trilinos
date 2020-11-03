// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_DatabaseIO.h>
#include <kelpie/Kelpie.hh>

#include <vector>

#ifndef IOSS_IODW_UTILS_H
#define IOSS_IODW_UTILS_H

namespace Iodw {

  namespace Utils {

    struct Property
    {
      kelpie::Key KeyType;
    };

    struct Field
    {
      kelpie::Key KeyType;
    };

    struct ElemBlock
    {
    };

    struct RegionKeys
    {
      using Keys = std::vector<kelpie::Key>;

      Keys edge_blocks;
      Keys elem_blocks;
      Keys face_blocks;
      Keys node_blocks;

      Keys edge_sets;
      Keys elem_sets;
      Keys face_sets;
      Keys node_sets;

      // TODO: serialize
    };

    class IossToDW
    {
    public:
      void operator()(Ioss::DatabaseIO *dbi);
    };

  } // namespace Utils
} // namespace Iodw

#endif
