/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
namespace Ioss {
  struct MeshCopyOptions
  {
    double minimum_time{0.0};
    double maximum_time{0.0};
    double delay{0.0};
    int    data_storage_type{0};
    bool   memory_statistics{false};
    bool   debug{false};
    bool   verbose{false};
    bool   ints_64_bit{false};
    bool   delete_timesteps{false};
    bool   reverse{false};     // Used for testing CGNS
    bool   add_proc_id{false}; // CGNS: Add proc_id field.
  };
} // namespace Ioss
