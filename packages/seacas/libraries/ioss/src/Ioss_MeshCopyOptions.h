namespace Ioss {
  struct MeshCopyOptions
  {
    bool   memory_statistics{false};
    bool   debug{false};
    bool   verbose{false};
    bool   ints_64_bit{false};
    bool   delete_timesteps{false};
    double minimum_time{0.0};
    double maximum_time{0.0};
    int    data_storage_type{0};
  };
}
