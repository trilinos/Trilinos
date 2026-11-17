// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "IossMesh.h"
#include "PartitionTypes.h"
#include <vector>

namespace utest_util {

  class Partition
  {
  public:
    Partition(IossMesh *mesh, const std::vector<EntityProc> &procAssign, const int nProc);

    virtual void write_nemesis_data(int exoid) const = 0;

    virtual size_t get_num_local_entities() const = 0;

    virtual int get_processor(size_t localId) const = 0;

    virtual size_t api_size() const = 0;

    virtual std::string type() const = 0;

    virtual ~Partition() = default;

  protected:
    IossMesh *m_mesh{nullptr};

    int                     m_numProcs;
    std::vector<EntityProc> m_entityPartition;

  protected:
    void verify_input_partition();
  };

  template <typename INT>
  void fill_node_element_connectivity(IossMesh                      *mesh,
                                      std::vector<std::vector<INT>> &surroundingElements,
                                      int                           &maxNumSurrroundingElements);

} // namespace utest_util
